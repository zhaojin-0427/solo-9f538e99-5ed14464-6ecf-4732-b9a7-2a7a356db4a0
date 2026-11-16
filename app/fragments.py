#!/usr/bin/env python3
"""已配对 FASTQ 的片段长度分析（纯标准库，流式处理）。

输入是配对流水线发布的 ``paired_R1.fastq.gz`` / ``paired_R2.fastq.gz``。
逐对读取并校验：

  * 四行 FASTQ 记录（标题 @、第三行 +、序列非空、序列与质量等长）；
  * 两条记录同位 ID 相同（ID 取标题首个词，去掉末尾 /1、/2）；
  * 质量字符为 Phred+33（ASCII 33..126）；
  * 任一侧读长超过 300 bp 时整次分析拒绝。

片段长度估计：

  * 将 R2 反向互补，在 R1 后缀与 R2 前缀之间寻找长度 >= 12 bp 的重叠；
  * 重叠最多允许 2 个错配，且每个错配位置至少一端 Phred 质量 < 20；
  * 重叠范围内任一端为 N（或其他非 ACGT 字符）时该候选重叠不成立；
  * 只有恰好一个合法重叠长度时，片段长度 = len(R1) + len(R2) - 重叠长度；
    零个或多个合法重叠都记为未知。

内存约束：两个文件逐行读取，不保存读段；片段长度最大 600 bp，
用固定 601 个计数槽做直方图，内存占用与输入规模无关。
"""
import gzip
import hashlib
import json
import os
import zlib

# 读长硬上限；超过则整次分析拒绝
MAX_READ_LENGTH = 300
# 重叠长度下限与错配上限
MIN_OVERLAP = 12
MAX_MISMATCHES = 2
LOW_QUALITY_PHRED = 20
PHRED_OFFSET = 33
# Phred+33 字符小于 '5'（ASCII 53）即质量值 < 20
LOW_QUALITY_BYTE = PHRED_OFFSET + LOW_QUALITY_PHRED

# 单行保护上限（与配对流水线一致）
MAX_LINE_BYTES = 16 * 1024 * 1024
# 单侧解压后大小上限：256 MiB
DEFAULT_MAX_DECOMP_BYTES = 256 * 1024 * 1024
# 分布展示的分箱宽度（bp）
BIN_WIDTH = 10
_MAX_FRAGMENT = MAX_READ_LENGTH * 2

# 碱基编码：A=0 C=1 G=2 T=3（大小写等价），其余（N 等）=8；
# 质量低（<20）在 bit4 置位
_SEQ_ENC = bytearray([8]) * 256
for _i, _b in enumerate(b"ACGTacgt"):
    _SEQ_ENC[_b] = _i % 4
_SEQ_ENC = bytes(_SEQ_ENC)
_QLOW_ENC = bytes(0x10 if i < LOW_QUALITY_BYTE else 0 for i in range(256))
_COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

# 按读长预计算逐字节掩码（little-endian：第 0 个碱基在最低字节）
_ALL_BYTES = [((1 << (8 * n)) - 1) for n in range(MAX_READ_LENGTH + 1)]
_BASE_BITS = [int.from_bytes(b"\x03" * n, "little")
              for n in range(MAX_READ_LENGTH + 1)]
_LOWQ_BITS = [int.from_bytes(b"\x10" * n, "little")
              for n in range(MAX_READ_LENGTH + 1)]
_N_BITS = [int.from_bytes(b"\x08" * n, "little")
           for n in range(MAX_READ_LENGTH + 1)]
_ONE_BITS = [int.from_bytes(b"\x01" * n, "little")
             for n in range(MAX_READ_LENGTH + 1)]


class AnalysisError(Exception):
    """输入数据或资源限制导致的可预期失败；调用方不得发布新结论。"""


def file_signature(path):
    """返回文件身份指纹（用于绑定分析依据）；文件缺失返回 None。"""
    try:
        st = os.stat(path)
    except OSError:
        return None
    return _signature_from_stat(st)


def _signature_from_stat(st):
    return {
        "dev": st.st_dev,
        "ino": st.st_ino,
        "size": st.st_size,
        "mtime_ns": st.st_mtime_ns,
    }


def _signature_pair(path1, path2):
    return {"r1": file_signature(path1), "r2": file_signature(path2)}


def _sha256_file(path):
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def source_sha256(path1, path2):
    """重算源文件内容指纹（分块读取，内存与文件大小无关）。

    任一端缺失/不可读时返回 None（调用方据此判定过期/失败）。
    不做缓存：即使大小与 mtime 被保留，内容变化也必须立刻可见。
    """
    hashes = {}
    for key, path in (("r1", path1), ("r2", path2)):
        try:
            hashes[key] = _sha256_file(path)
        except OSError:
            return None
    return json.dumps(hashes, sort_keys=True)


class _HashingReader:
    """包裹原始文件对象，边读边计算压缩文件的 SHA-256。"""

    def __init__(self, raw, hasher):
        self.raw = raw
        self.hasher = hasher

    def read(self, n=-1):
        data = self.raw.read(n)
        if data:
            self.hasher.update(data)
        return data

    def readinto(self, buf):
        n = self.raw.readinto(buf)
        if n:
            self.hasher.update(bytes(memoryview(buf)[:n]))
        return n


def _readline(gz, side):
    try:
        line = gz.readline(MAX_LINE_BYTES + 1)
    except (OSError, EOFError, zlib.error) as e:
        raise AnalysisError(f"{side} 不是有效的 gzip 文件或数据已损坏: {e}") from e
    if not line:
        return None
    if len(line) > MAX_LINE_BYTES:
        raise AnalysisError(
            f"{side} 存在超过 {MAX_LINE_BYTES} 字节的超长行，记录非法")
    return line


def _extract_read_id(raw_header):
    """返回 (去 /1、/2 的 ID, 标题中声明的端号或 None)。"""
    body = raw_header[1:].strip()
    if not body:
        raise AnalysisError("读段标题为空")
    token = body.split(None, 1)[0]
    side = None
    if token.endswith(b"/1"):
        side, token = 1, token[:-2]
    elif token.endswith(b"/2"):
        side, token = 2, token[:-2]
    else:
        # 新式标题：@ID 1:N:0:XX / @ID 2:N:0:XX
        rest = body.split(None, 1)
        if len(rest) == 2 and rest[1][:1] in (b"1", b"2") and (
                len(rest[1]) == 1 or rest[1][1:2] in (b":", b"\t", b" ")):
            side = int(rest[1][:1])
    if not token:
        raise AnalysisError("读段 ID 为空")
    return token, side


def _read_record(gz, side, nrec, counter, max_bytes):
    """读取并校验一条四行记录；EOF 返回 None。"""
    lines = []
    for _ in range(4):
        line = _readline(gz, side)
        if line is None:
            break
        counter[0] += len(line)
        if counter[0] > max_bytes:
            raise AnalysisError(
                f"{side} 解压后大小超过上限 "
                f"({max_bytes // (1024 * 1024)} MiB)")
        lines.append(line)
    if not lines:
        return None
    if len(lines) < 4:
        raise AnalysisError(f"{side} 第 {nrec} 条记录不完整（文件末尾缺少行）")
    header, seq_line, plus_line, qual_line = lines
    if not header.startswith(b"@"):
        raise AnalysisError(f"{side} 第 {nrec} 条记录标题行未以 @ 开头")
    if not plus_line.startswith(b"+"):
        raise AnalysisError(f"{side} 第 {nrec} 条记录第三行未以 + 开头")
    seq = seq_line.rstrip(b"\r\n")
    qual = qual_line.rstrip(b"\r\n")
    if not seq:
        raise AnalysisError(f"{side} 第 {nrec} 条记录序列为空")
    if len(seq) != len(qual):
        raise AnalysisError(f"{side} 第 {nrec} 条记录序列与质量值长度不一致")
    if min(qual) < PHRED_OFFSET or max(qual) > 126:
        raise AnalysisError(
            f"{side} 第 {nrec} 条记录质量值不是合法的 Phred+33 字符")
    if len(seq) > MAX_READ_LENGTH:
        raise AnalysisError(
            f"{side} 第 {nrec} 条记录读长 {len(seq)} bp 超过 "
            f"{MAX_READ_LENGTH} bp 上限，整次分析拒绝")
    return _extract_read_id(header), seq, qual


def _encode(seq, qual):
    """序列与低质量标记按字节编码进一个大整数（第 i 个碱基占第 i 字节）。"""
    return (int.from_bytes(seq.translate(_SEQ_ENC), "little")
            | int.from_bytes(qual.translate(_QLOW_ENC), "little"))


def estimate_fragment_length(seq1, seq2, qual1, qual2):
    """估算一对读段的片段长度；无法唯一确定时返回 None。

    仅使用 R1 后缀 / revcomp(R2) 前缀的重叠；重叠须 >= 12 bp、
    最多 2 个错配、每个错配至少一端质量 < 20、重叠内不允许 N。
    """
    n1 = len(seq1)
    n2 = len(seq2)
    if n1 > MAX_READ_LENGTH or n2 > MAX_READ_LENGTH:
        raise AnalysisError(
            f"读长超过 {MAX_READ_LENGTH} bp 上限，整次分析拒绝")
    if n1 < MIN_OVERLAP or n2 < MIN_OVERLAP:
        return None

    r2rc = seq2.translate(_COMPLEMENT)[::-1]
    a = _encode(seq1, qual1)
    b = _encode(r2rc, qual2[::-1])

    # 重叠长度 L = n1 - d，d 为 R1 起点相对 revcomp(R2) 前端的位移。
    d_lo = max(0, n1 - n2)
    d_hi = n1 - MIN_OVERLAP
    found = None
    for d in range(d_lo, d_hi + 1):
        length = n1 - d
        all_mask = _ALL_BYTES[length]
        shifted = a >> (8 * d)
        ab = (shifted | b) & all_mask
        if ab & _N_BITS[length]:
            continue  # 重叠范围内含 N（或其他非 ACGT），不能支持重叠
        diff = (shifted ^ b) & all_mask
        lowq = ab & _LOWQ_BITS[length]  # 任一端质量 < 20 的字节标记（bit4）
        mismatches = ((diff | (diff >> 1)) & _ONE_BITS[length]).bit_count()
        if mismatches > MAX_MISMATCHES:
            continue
        # 每个错配至少一端质量 < 20：错配字节若两端皆高质量即不合法。
        # 低质量标记在 bit4，乘 3 把每个低质量字节的低两位都置为允许位。
        allowed = (lowq >> 4) * 3
        if (diff & _BASE_BITS[length]) & ~allowed:
            continue
        if found is not None:
            return None  # 多个合法重叠 -> 未知
        found = length
    if found is None:
        return None
    return n1 + n2 - found


def _bin_distribution(hist):
    bins = []
    for start in range(0, _MAX_FRAGMENT + 1, BIN_WIDTH):
        end = min(start + BIN_WIDTH - 1, _MAX_FRAGMENT)
        count = sum(hist[start:end + 1])
        if count:
            bins.append({"start": start, "end": end, "count": count})
    return bins


def analyze(r1_gz, r2_gz, threshold, max_bytes=DEFAULT_MAX_DECOMP_BYTES,
            expected_sig=None):
    """分析已发布的配对 R1/R2，返回结果字典；失败抛 AnalysisError。

    ``expected_sig`` 为请求时捕获的 ``{"r1": sig, "r2": sig}``；分析开始前
    或分析期间指纹变化都拒绝发布结论。
    """
    if not isinstance(threshold, int) or isinstance(threshold, bool):
        raise AnalysisError("阈值必须是整数")
    before = _signature_pair(r1_gz, r2_gz)
    if before["r1"] is None or before["r2"] is None:
        missing = "R1" if before["r1"] is None else "R2"
        raise AnalysisError(f"配对结果文件缺失（{missing}），未发布结论")
    if expected_sig is not None and before != expected_sig:
        raise AnalysisError("分析开始前源文件已变化，未发布结论")

    hash1 = hashlib.sha256()
    hash2 = hashlib.sha256()
    hist = [0] * (_MAX_FRAGMENT + 1)
    total = known = unknown = below = 0
    max_read = 0
    sides_seen = {1: set(), 2: set()}
    counter1, counter2 = [0], [0]
    try:
        with open(r1_gz, "rb") as raw1, open(r2_gz, "rb") as raw2:
            gz1 = gzip.GzipFile(filename=r1_gz, mode="rb",
                                fileobj=_HashingReader(raw1, hash1))
            gz2 = gzip.GzipFile(filename=r2_gz, mode="rb",
                                fileobj=_HashingReader(raw2, hash2))
            try:
                while True:
                    rec1 = _read_record(gz1, "R1", total + 1, counter1, max_bytes)
                    rec2 = _read_record(gz2, "R2", total + 1, counter2, max_bytes)
                    if rec1 is None and rec2 is None:
                        break
                    if rec1 is None or rec2 is None:
                        raise AnalysisError(
                            f"两端配对文件记录数不一致（第 {total + 1} 对）")
                    total += 1
                    (id1, side1), seq1, qual1 = rec1
                    (id2, side2), seq2, qual2 = rec2
                    if id1 != id2:
                        raise AnalysisError(
                            f"第 {total} 对读段 ID 不一致（{id1.decode('utf-8', 'replace')} "
                            f"vs {id2.decode('utf-8', 'replace')}），未发布结论")
                    if side1 is not None:
                        sides_seen[1].add(side1)
                    if side2 is not None:
                        sides_seen[2].add(side2)
                    if side1 == 2 or side2 == 1 or len(sides_seen[1]) > 1 \
                            or len(sides_seen[2]) > 1:
                        raise AnalysisError(
                            f"第 {total} 对读段端号异常（R1 文件端号 "
                            f"{sorted(sides_seen[1])}，R2 文件端号 "
                            f"{sorted(sides_seen[2])}），疑似两端文件错位，"
                            "未发布结论")
                    max_read = max(max_read, len(seq1), len(seq2))
                    try:
                        frag = estimate_fragment_length(seq1, seq2, qual1, qual2)
                    except AnalysisError as e:
                        raise AnalysisError(f"第 {total} 对读段：{e}") from e
                    if frag is None:
                        unknown += 1
                    else:
                        known += 1
                        hist[frag] += 1
                        if frag < threshold:
                            below += 1
            finally:
                gz1.close()
                gz2.close()
    except AnalysisError:
        raise
    except (OSError, EOFError, zlib.error) as e:
        raise AnalysisError(f"配对结果文件损坏: {e}") from e

    after = _signature_pair(r1_gz, r2_gz)
    if after != before:
        raise AnalysisError("分析期间源文件发生变化，未发布结论")

    def pct(n):
        return round(100.0 * n / total, 2) if total else 0.0

    known_lengths = [length for length in range(_MAX_FRAGMENT + 1)
                     if hist[length]]
    mean = (sum(length * hist[length] for length in known_lengths) / known
            if known else None)
    return {
        "threshold": threshold,
        "total_pairs": total,
        "known_pairs": known,
        "unknown_pairs": unknown,
        "below_threshold_pairs": below,
        "known_pct": pct(known),
        "unknown_pct": pct(unknown),
        "below_threshold_pct": pct(below),
        "min_length": known_lengths[0] if known_lengths else None,
        "max_length": known_lengths[-1] if known_lengths else None,
        "mean_length": round(mean, 2) if mean is not None else None,
        "max_read_length": max_read,
        "min_overlap": MIN_OVERLAP,
        "max_mismatches": MAX_MISMATCHES,
        "low_quality_phred": LOW_QUALITY_PHRED,
        "bin_width": BIN_WIDTH,
        "distribution": _bin_distribution(hist),
    }, {
        "source_sig": json.dumps(before, sort_keys=True),
        "source_sha256": json.dumps(
            {"r1": hash1.hexdigest(), "r2": hash2.hexdigest()}, sort_keys=True),
    }


if __name__ == "__main__":
    # 命令行自测: python3 -I fragments.py paired_R1.fastq.gz paired_R2.fastq.gz THRESHOLD
    import sys
    if len(sys.argv) != 4:
        sys.exit("用法: fragments.py paired_R1.fastq.gz paired_R2.fastq.gz THRESHOLD")
    try:
        result, _binding = analyze(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    except AnalysisError as e:
        sys.exit(f"分析失败: {e}")
    print(json.dumps(result, ensure_ascii=False, indent=2))
