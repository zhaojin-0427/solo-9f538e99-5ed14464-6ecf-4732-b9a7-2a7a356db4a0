# 双端 FASTQ 配对修复台

测序质控场景：收到**顺序错乱的双端读段**（R1/R2 两个 gzip FASTQ），本服务将其恢复为
可交给比对器的成对文件，并分离出缺端记录。Flask + SQLite，全流式处理，
容器内存限制 128 MiB 下可处理单端解压后 256 MiB 的输入。

## 功能

- 网页上传 R1、R2 两份 gzip FASTQ（每份解压后 ≤ 256 MiB，超限整批失败）；
- 读段 ID = 标题首个词，去掉末尾 `/1`、`/2`（兼容 `@ID/1` 与 `@ID 1:N:0:XX` 两种标题）；
- **整批校验**：同侧重复 ID、非法四行记录（标题非 `@`、第三行非 `+`、
  序列与质量长度不符、记录残缺、非 gzip 数据）一律整批失败——
  包括文件末尾才出现的错误（先全量校验建索引，后产出结果）；
- 按 ID 配对，两个成对输出均按 **R1 原顺序** 排列；
  缺端记录分别输出 `unpaired_R1` / `unpaired_R2`（各自保持原侧顺序）；
  原始标题、序列、质量值逐字节保留；
- 任务列表/详情页展示状态（排队中/处理中/完成/失败）、配对数与缺端数，
  完成后可下载四个 gzip 结果；
- **片段长度分析**（任务完成后在详情页发起）：逐对读取已发布的
  `paired_R1` / `paired_R2`，校验同位 ID、四行记录、Phred+33 质量字符与端号；
  R2 反向互补后寻找 R1 后缀与 R2 前缀 ≥ 12 bp 的重叠（≤ 2 个错配，
  每个错配至少一端质量 < 20，含 N 不成立）；仅有唯一合法重叠时按
  `len(R1)+len(R2)−重叠` 估算片段长度，零个或多个合法重叠均记未知；
  读长超 300 bp 整次拒绝。页面展示长度分布、低于所填阈值与未知的读段对占比；
- **可靠性**：结果先写入临时区，全部产出并 fsync、原子 rename、
  SQLite 登记成功后才开放下载；失败不保留任何部分结果；
  服务重启时自动清理中断任务的未完成产物，已完成结果仍可下载。
- **分析结论绑定依据**：每次分析记录源文件指纹（dev/inode/size/mtime）与
  SHA-256、阈值和发布时刻；文件缺失/损坏/两端错位/记录数或 ID 不一致、
  分析期间依据变化时均不发布新结论；源文件或阈值变化后旧结论在页面标记过期；
  并发/迟到的分析结果不能覆盖更新的分析。

## 快速开始（Docker Compose）

```bash
docker compose up --build
# 浏览器打开 http://localhost:8000
```

- 容器内存上限 128 MiB（`mem_limit`），数据持久化在命名卷 `fastq_data`
  （上传文件、结果文件、SQLite 登记库都在 `/data`）；
- 用 `samples/` 下的 `normal_r1.fastq.gz` / `normal_r2.fastq.gz` 即可体验
  乱序 + 缺端配对；`dup_*`、`badqual_*`、`truncated_*`、`notgzip_*`
  为必失败样例，可验证整批失败行为。

## 本地开发运行

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
cd app && DATA_DIR=../data PORT=8000 ../.venv/bin/python app.py
```

流水线也可脱离 Web 独立自测：

```bash
python3 -I app/pipeline.py R1.fastq.gz R2.fastq.gz WORK_DIR OUT_DIR
```

## 样例数据

`samples/make_samples.py` 可重新生成全部样例（确定性随机）：

| 样例 | 说明 | 预期 |
|---|---|---|
| `normal_r1/r2.fastq.gz` | 240 对乱序 + 16 条 R1 缺端 + 8 条 R2 缺端 | paired=240, r1_only=16, r2_only=8 |
| `dup_*` | R1 末尾出现重复 ID | 整批失败 |
| `badqual_*` | R1 末尾记录序列/质量长度不符 | 整批失败 |
| `truncated_*` | R1 末尾记录只有 3 行 | 整批失败 |
| `notgzip_*` | 内容不是 gzip | 整批失败 |
| `frag_r1/r2.fastq.gz` | 已配对、同位同序的片段分析样例：150 bp / 180 bp 各 100 对 + 4 对无重叠 | 阈值 180 时 known=200（below=100）、unknown=4 |

`frag_*` 可直接用片段分析 CLI 自测（不需要先跑配对）：

```bash
python3 -I app/fragments.py samples/frag_r1.fastq.gz samples/frag_r2.fastq.gz 180
```

## 片段长度分析规则

| 项目 | 规则 |
|---|---|
| 输入 | 仅读已发布的 `paired_R1.fastq.gz` / `paired_R2.fastq.gz` |
| 记录校验 | 四行结构、序列非空、序列/质量等长、质量字符 Phred+33（ASCII 33–126） |
| 配对校验 | 逐对同位 ID 相同；两端记录数一致；标题端号（`/1`、`/2` 或 `1:N:0`）与文件侧一致（错位即拒绝） |
| 读长上限 | 任意一条 > 300 bp 整次拒绝 |
| 重叠 | R2 反向互补后，R1 后缀 vs R2 前缀，长度 ≥ 12 bp |
| 错配 | ≤ 2 个；每个错配至少一端 Phred < 20 |
| N | 重叠范围内任一端含 N（或非 ACGT）该重叠不成立 |
| 片段长度 | 唯一合法重叠时 `len(R1)+len(R2)−重叠`；0 个或多个均记未知 |
| 统计 | 总对数、可估算数、低于阈值数、未知数及各占比；10 bp 分箱的分布 |
| 依据绑定 | 源文件指纹 + SHA-256 + 阈值 + 发布时刻写入 SQLite |
| 过期 | 源文件缺失/内容或元数据变化、任务不再是完成态、阈值变化后旧结论标记过期 |
| 迟到结果 | 以 request_token 校验，非当前分析的结果一律丢弃 |

## 输出文件

| 文件 | 内容 |
|---|---|
| `paired_R1.fastq.gz` | 有 mate 的 R1 记录，按 R1 原顺序 |
| `paired_R2.fastq.gz` | 对应的 R2 记录，与 paired_R1 同序 |
| `unpaired_R1.fastq.gz` | 缺端 R1（R2 中无配对），R1 原顺序 |
| `unpaired_R2.fastq.gz` | 缺端 R2（R1 中无配对），R2 原顺序 |

## 设计要点（为何能在 128 MiB 内存下工作）

1. **流式校验建索引**：逐行解压（`readline` 带 16 MiB 单行保护上限），
   原始记录字节顺序写入磁盘 records 文件，仅把 `(id, offset, length)`
   登记进临时 SQLite（唯一索引查重），全程不在内存保存读段集合；
2. **游标式配对**：按 R1 的 rowid（即原始顺序）分批游标遍历，
   对 R2 索引做点查，用 `os.pread` 按偏移取回记录写出四个 gzip；
3. **原子发布**：结果先写 `work/`，全部完成并 fsync 后 `rename` 进 `out/`，
   再提交 SQLite 登记；任一步失败即清理，下载接口只对 `done` 任务开放；
4. **启动恢复**：`queued`/`processing` 状态的任务在重启时标记失败并
   删除其 `work/`、`out/`；目录中无登记的任务目录一并清除。

实测：100 万对读段（每侧解压约 209 MiB）处理约 50 秒，
流水线峰值内存 17 MiB，含 Flask 的服务进程峰值约 39 MiB。

## 配置（环境变量）

| 变量 | 默认 | 说明 |
|---|---|---|
| `DATA_DIR` | `/data`（容器） | 数据目录（上传、结果、SQLite） |
| `PORT` | `8000` | 监听端口 |
| `MAX_DECOMP_BYTES` | `268435456` | 每份文件解压后上限（256 MiB） |
| `MAX_UPLOAD_BYTES` | 上值 + 2 MiB | 上传压缩包大小上限 |
| `SECRET_KEY` | 开发默认值 | Flask 会话密钥，生产应覆盖 |

## 目录结构

```
app/
  app.py            Flask 应用：页面、任务登记、后台处理、下载、启动恢复、
                    片段长度分析登记/发布/过期判定
  pipeline.py       流式配对流水线（校验/索引/配对/原子发布）
  fragments.py      片段长度分析（流式校验、重叠估算、直方图、依据指纹）
  templates/        页面模板
samples/            演示样例 + 生成脚本
tests/
  verify_pairing.py 独立校验脚本（比对输入与四个输出）
  make_big.py       大样本生成（压测用）
Dockerfile / docker-compose.yml / requirements.txt
```

## 已知限制

- 内置 Flask 开发服务器（threaded），适合演示与内部工具；
  生产部署建议换 gunicorn/waitress 并设置 `SECRET_KEY`；
- 仅支持标准四行 FASTQ（不支持序列折行的多行 FASTQ）；
  单行超过 16 MiB 视为非法记录（资源保护）；
- 末尾记录缺少换行符时输出会补一个换行（保证输出为合法 FASTQ），
  其余内容逐字节保留；
- 片段分析为演示级比对：只考虑“R1 后缀 / R2 前缀”这一种重叠方向，
  不处理末端修剪、adapter 或测序错误模型带来的复杂情形；
  重叠长度上限受读长约束（读长 ≤ 300 bp，片段长度 ≤ 600 bp）；
- 片段分析在后台线程执行（与配对共用单进程）；超大任务的页面轮询
  （详情页/状态接口）只读数据库，不会重复触发分析；
- 无任务删除/清理界面，历史任务与文件会一直保留在数据卷中。
