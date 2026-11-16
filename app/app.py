#!/usr/bin/env python3
"""双端 FASTQ 配对修复服务（Flask + SQLite）。

- 上传 R1/R2 两份 gzip FASTQ，后台流式配对；
- 全部结果文件发布且 SQLite 登记成功后才开放下载；
- 完成后可对已发布配对结果做片段长度分析（结论绑定源文件内容与阈值）；
- 失败不暴露任何部分结果；重启时清理中断任务的未完成产物。
"""
import json
import os
import shutil
import sqlite3
import threading
import time
import uuid

from flask import (Flask, abort, flash, g, jsonify, redirect, render_template,
                   request, send_file, url_for)
from werkzeug.exceptions import RequestEntityTooLarge
from werkzeug.utils import secure_filename

import fragments
import pipeline

DATA_DIR = os.path.abspath(os.environ.get(
    "DATA_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")))
JOBS_DIR = os.path.join(DATA_DIR, "jobs")
DB_PATH = os.path.join(DATA_DIR, "app.db")

# 解压后上限（默认 256 MiB，可用环境变量覆盖以便测试）
MAX_DECOMP_BYTES = int(os.environ.get(
    "MAX_DECOMP_BYTES", str(pipeline.DEFAULT_MAX_DECOMP_BYTES)))
# 上传压缩包上限：解压上限 + 2 MiB 余量（不可压缩数据 gzip 后略大于原始）
MAX_UPLOAD_BYTES = int(os.environ.get(
    "MAX_UPLOAD_BYTES", str(MAX_DECOMP_BYTES + 2 * 1024 * 1024)))

OUTPUTS = {
    "paired_r1": "配对 R1（按 R1 原顺序）",
    "paired_r2": "配对 R2（按 R1 原顺序）",
    "unpaired_r1": "缺端 R1（R2 中无配对）",
    "unpaired_r2": "缺端 R2（R1 中无配对）",
}

# 片段长度分析固定读取已发布的成对输出文件
FRAGMENT_INPUTS = {
    "r1": pipeline.OUTPUT_FILES["paired_r1"],
    "r2": pipeline.OUTPUT_FILES["paired_r2"],
}

# 阈值允许范围（bp）
MIN_THRESHOLD = 0
MAX_THRESHOLD = 1_000_000
# 任务详情页表单里的默认阈值
DEFAULT_FRAGMENT_THRESHOLD = 150

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES
app.secret_key = os.environ.get("SECRET_KEY", "fastq-pairer-dev-secret")


# ---------------- SQLite 登记库 ----------------

def connect_db():
    con = sqlite3.connect(DB_PATH, timeout=30)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA busy_timeout=30000")
    return con


def get_db():
    if "db" not in g:
        g.db = connect_db()
    return g.db


@app.teardown_appcontext
def close_db(_exc):
    db = g.pop("db", None)
    if db is not None:
        db.close()


def init_and_recover():
    """建表 + 启动恢复：中断任务标记失败并清理其未完成产物。"""
    os.makedirs(JOBS_DIR, exist_ok=True)
    con = connect_db()
    con.execute(
        """CREATE TABLE IF NOT EXISTS jobs(
               id TEXT PRIMARY KEY,
               created_at TEXT NOT NULL,
               status TEXT NOT NULL,          -- queued/processing/done/failed
               message TEXT NOT NULL DEFAULT '',
               r1_name TEXT NOT NULL DEFAULT '',
               r2_name TEXT NOT NULL DEFAULT '',
               paired INTEGER, r1_only INTEGER, r2_only INTEGER)""")
    con.commit()
    con.execute(
        """CREATE TABLE IF NOT EXISTS fragment_analyses(
               job_id TEXT PRIMARY KEY,
               status TEXT NOT NULL,          -- running/done/failed
               threshold INTEGER,             -- 本次尝试的阈值
               request_token TEXT NOT NULL DEFAULT '',
               attempt_sig TEXT NOT NULL DEFAULT '',
               message TEXT NOT NULL DEFAULT '',
               started_at TEXT NOT NULL DEFAULT '',
               finished_at TEXT NOT NULL DEFAULT '',
               published_json TEXT NOT NULL DEFAULT '',
               published_threshold INTEGER,
               published_sig TEXT NOT NULL DEFAULT '',
               published_sha256 TEXT NOT NULL DEFAULT '',
               published_at TEXT NOT NULL DEFAULT '')""")
    con.commit()
    # 兼容早于片段分析功能的数据库：补列（幂等）
    existing = {r["name"] for r in con.execute(
        "PRAGMA table_info(fragment_analyses)")}
    for column, ddl in (
            ("published_sha256", "TEXT NOT NULL DEFAULT ''"),
            ("published_threshold", "INTEGER"),
            ("published_sig", "TEXT NOT NULL DEFAULT ''"),
            ("published_at", "TEXT NOT NULL DEFAULT ''"),
            ("attempt_sig", "TEXT NOT NULL DEFAULT ''"),
            ("request_token", "TEXT NOT NULL DEFAULT ''"),
            ("started_at", "TEXT NOT NULL DEFAULT ''"),
            ("finished_at", "TEXT NOT NULL DEFAULT ''")):
        if column not in existing:
            con.execute(
                f"ALTER TABLE fragment_analyses ADD COLUMN {column} {ddl}")
    con.commit()
    interrupted = [r["id"] for r in con.execute(
        "SELECT id FROM jobs WHERE status IN ('queued','processing')")]
    for jid in interrupted:
        job_dir = os.path.join(JOBS_DIR, jid)
        shutil.rmtree(os.path.join(job_dir, "work"), ignore_errors=True)
        shutil.rmtree(os.path.join(job_dir, "out"), ignore_errors=True)
        con.execute(
            "UPDATE jobs SET status='failed', "
            "message='服务重启导致任务中断，未完成产物已清理' WHERE id=?", (jid,))
    # 分析线程随进程终止；未完成的分析标记失败（不发布任何结论）。
    # 若该任务在重启前已发布过旧结论，旧结论仍保留并按“过期”展示。
    con.execute(
        "UPDATE fragment_analyses SET status='failed', message=?, finished_at=? "
        "WHERE status='running'",
        ("服务重启导致分析中断，未发布新结论"
         + "（页面保留的旧结论按过期处理）", time.strftime("%Y-%m-%d %H:%M:%S")))
    con.commit()
    known = {r["id"] for r in con.execute("SELECT id FROM jobs")}
    for name in os.listdir(JOBS_DIR):
        if name not in known:
            shutil.rmtree(os.path.join(JOBS_DIR, name), ignore_errors=True)
    con.close()
    if interrupted:
        app.logger.info("已清理 %d 个中断任务", len(interrupted))


# ---------------- 后台处理 ----------------

def process_job(jid):
    """后台线程：跑流水线；全部结果发布 + SQLite 登记成功才算完成。"""
    job_dir = os.path.join(JOBS_DIR, jid)
    work = os.path.join(job_dir, "work")
    out = os.path.join(job_dir, "out")
    con = connect_db()
    try:
        con.execute("UPDATE jobs SET status='processing' WHERE id=?", (jid,))
        con.commit()
        os.makedirs(work, exist_ok=True)
        counts = pipeline.run(
            os.path.join(job_dir, "uploads", "r1.fastq.gz"),
            os.path.join(job_dir, "uploads", "r2.fastq.gz"),
            work, out, max_bytes=MAX_DECOMP_BYTES)
        # 结果文件已全部落盘，最后登记；若此处失败则走 except 删除结果
        con.execute(
            "UPDATE jobs SET status='done', paired=?, r1_only=?, r2_only=?, "
            "message='' WHERE id=?",
            (counts["paired"], counts["r1_only"], counts["r2_only"], jid))
        con.commit()
    except pipeline.PipelineError as e:
        con.rollback()
        shutil.rmtree(out, ignore_errors=True)  # 失败不暴露部分结果
        con.execute("UPDATE jobs SET status='failed', message=? WHERE id=?",
                    (str(e), jid))
        con.commit()
    except Exception as e:  # 意外错误（含登记失败）：同样清理结果
        con.rollback()
        shutil.rmtree(out, ignore_errors=True)
        con.execute("UPDATE jobs SET status='failed', message=? WHERE id=?",
                    (f"内部错误: {e!r}", jid))
        con.commit()
        app.logger.exception("任务 %s 处理异常", jid)
    finally:
        shutil.rmtree(work, ignore_errors=True)  # 中间产物（索引/记录文件）即清
        con.close()


# ---------------- 片段长度分析 ----------------

def _now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _source_paths(jid):
    job_dir = os.path.join(JOBS_DIR, jid, "out")
    return (os.path.join(job_dir, FRAGMENT_INPUTS["r1"]),
            os.path.join(job_dir, FRAGMENT_INPUTS["r2"]))


def _current_source_sig(jid):
    r1_path, r2_path = _source_paths(jid)
    return {"r1": fragments.file_signature(r1_path),
            "r2": fragments.file_signature(r2_path)}


def _sig_text(sig):
    return json.dumps(sig, sort_keys=True)


def _stale_reason(row, job):
    """判定已发布结论是否过期；返回原因文本，未过期返回空串。

    过期条件：源配对结果文件缺失/内容变化、任务不再是完成态、
    或当前请求的阈值与该结论绑定的阈值不同。
    """
    if not row or not row["published_json"]:
        return ""
    reasons = []
    if job is None or job["status"] != "done":
        reasons.append("任务状态已不再是完成")
    r1_path, r2_path = _source_paths(row["job_id"])
    sig = {"r1": fragments.file_signature(r1_path),
           "r2": fragments.file_signature(r2_path)}
    if sig["r1"] is None or sig["r2"] is None:
        reasons.append("源配对结果文件缺失")
    elif _sig_text(sig) != row["published_sig"]:
        reasons.append("源配对结果文件已变化")
    if not reasons:
        # 指纹一致时再以内容指纹复核（防止大小/时间戳被伪造或精度不足）
        sha = fragments.source_sha256(r1_path, r2_path)
        if sha is None or sha != row["published_sha256"]:
            reasons.append("源配对结果文件内容已变化")
    if (row["threshold"] is not None
            and row["published_threshold"] != row["threshold"]):
        reasons.append(
            f"阈值已变化（结论绑定 {row['published_threshold']} bp，"
            f"当前 {row['threshold']} bp）")
    return "；".join(reasons)


def _load_published(row, job):
    """返回 (结果 dict 或 None, 过期原因)。

    状态为 running/failed 时仍加载更早的已发布结论，但一律按过期展示
    （它已不是当前分析的结论）。
    """
    if not row or not row["published_json"]:
        return None, ""
    reason = _stale_reason(row, job)
    if not reason and row["status"] != "done":
        reason = "最近一次分析未成功，已发布结论不再代表当前依据"
    try:
        result = json.loads(row["published_json"])
    except ValueError:
        return None, reason or "已发布结论无法解析"
    return result, reason


def _record_analysis_failure(jid, token, message):
    con = connect_db()
    try:
        # 只更新仍然属于本次尝试的行；迟到失败不能覆盖新分析的状态字段
        con.execute(
            "UPDATE fragment_analyses SET status='failed', message=?, "
            "finished_at=? WHERE job_id=? AND request_token=?",
            (message, _now(), jid, token))
        con.commit()
    finally:
        con.close()


def process_fragment_analysis(jid, threshold, token, expected_sig, started_at):
    """后台线程：分析已发布配对结果；只有依据未变才发布结论。"""
    r1_path, r2_path = _source_paths(jid)
    try:
        result, binding = fragments.analyze(
            r1_path, r2_path, threshold, max_bytes=MAX_DECOMP_BYTES,
            expected_sig=expected_sig)
    except fragments.AnalysisError as e:
        _record_analysis_failure(jid, token, str(e))
        return
    except Exception as e:  # 意外错误：同样不发布结论
        app.logger.exception("任务 %s 片段分析异常", jid)
        _record_analysis_failure(jid, token, f"内部错误: {e!r}")
        return

    con = connect_db()
    try:
        con.execute("BEGIN IMMEDIATE")
        current = con.execute(
            "SELECT * FROM fragment_analyses WHERE job_id=?", (jid,)).fetchone()
        if current is None or current["request_token"] != token:
            con.rollback()  # 迟到结果：已被更新的分析取代，丢弃
            return
        job = con.execute(
            "SELECT status FROM jobs WHERE id=?", (jid,)).fetchone()
        if job is None or job["status"] != "done":
            con.rollback()
            con.execute(
                "UPDATE fragment_analyses SET status='failed', message=?, "
                "finished_at=? WHERE job_id=? AND request_token=?",
                ("任务不再是完成状态，未发布结论", _now(), jid, token))
            con.commit()
            return
        now_sig = _current_source_sig(jid)
        if now_sig["r1"] is None or now_sig["r2"] is None:
            con.rollback()
            con.execute(
                "UPDATE fragment_analyses SET status='failed', message=?, "
                "finished_at=? WHERE job_id=? AND request_token=?",
                ("源配对结果文件缺失，未发布结论", _now(), jid, token))
            con.commit()
            return
        if _sig_text(now_sig) != binding["source_sig"]:
            con.rollback()
            con.execute(
                "UPDATE fragment_analyses SET status='failed', message=?, "
                "finished_at=? WHERE job_id=? AND request_token=?",
                ("分析期间源文件发生变化，未发布结论", _now(), jid, token))
            con.commit()
            return
        # 元数据未变也不足以证明内容未变：复核内容指纹后再发布
        now_sha = fragments.source_sha256(r1_path, r2_path)
        if now_sha is None or now_sha != binding["source_sha256"]:
            con.rollback()
            con.execute(
                "UPDATE fragment_analyses SET status='failed', message=?, "
                "finished_at=? WHERE job_id=? AND request_token=?",
                ("分析期间源文件内容发生变化，未发布结论", _now(), jid, token))
            con.commit()
            return
        con.execute(
            "UPDATE fragment_analyses SET status='done', message='', "
            "finished_at=?, published_json=?, published_threshold=?, "
            "published_sig=?, published_sha256=?, published_at=? "
            "WHERE job_id=? AND request_token=?",
            (_now(), json.dumps(result, ensure_ascii=False), threshold,
             binding["source_sig"], binding["source_sha256"], _now(), jid, token))
        con.commit()
    except Exception:
        con.rollback()
        app.logger.exception("任务 %s 分析登记失败", jid)
        try:
            con.execute(
                "UPDATE fragment_analyses SET status='failed', message=?, "
                "finished_at=? WHERE job_id=? AND request_token=?",
                ("分析登记失败，未发布结论", _now(), jid, token))
            con.commit()
        except sqlite3.Error:
            pass
    finally:
        con.close()


def start_fragment_analysis(jid, threshold):
    """登记一次新分析并启动后台线程；返回 request_token。"""
    sig = _current_source_sig(jid)
    if sig["r1"] is None or sig["r2"] is None:
        raise fragments.AnalysisError("源配对结果文件缺失，无法分析")
    token = uuid.uuid4().hex
    db = get_db()
    now = _now()
    # 新分析只更新“当前尝试”字段；已发布的旧结论保留，供页面标记为过期。
    db.execute(
        """INSERT INTO fragment_analyses(
               job_id, status, threshold, request_token, attempt_sig,
               message, started_at, finished_at,
               published_json, published_threshold, published_sig,
               published_sha256, published_at)
           VALUES (?, 'running', ?, ?, ?, '', ?, '', '', NULL, '', '', '')
           ON CONFLICT(job_id) DO UPDATE SET
               status='running', threshold=excluded.threshold,
               request_token=excluded.request_token,
               attempt_sig=excluded.attempt_sig, message='',
               started_at=excluded.started_at, finished_at=''""",
        (jid, threshold, token, _sig_text(sig), now))
    db.commit()
    threading.Thread(
        target=process_fragment_analysis,
        args=(jid, threshold, token, sig, now), daemon=True).start()
    return token


# ---------------- 页面与接口 ----------------

@app.route("/")
def index():
    jobs = get_db().execute(
        "SELECT * FROM jobs ORDER BY created_at DESC, id DESC LIMIT 200").fetchall()
    active = any(j["status"] in ("queued", "processing") for j in jobs)
    return render_template("index.html", jobs=jobs, active=active,
                           max_mib=MAX_DECOMP_BYTES // (1024 * 1024))


@app.route("/jobs", methods=["POST"])
def create_job():
    f1 = request.files.get("r1")
    f2 = request.files.get("r2")
    if not f1 or not f1.filename or not f2 or not f2.filename:
        flash("请同时选择 R1 和 R2 两个 gzip FASTQ 文件", "error")
        return redirect(url_for("index"))
    jid = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
    uploads = os.path.join(JOBS_DIR, jid, "uploads")
    os.makedirs(uploads)
    f1.save(os.path.join(uploads, "r1.fastq.gz"))
    f2.save(os.path.join(uploads, "r2.fastq.gz"))
    db = get_db()
    db.execute(
        "INSERT INTO jobs(id, created_at, status, r1_name, r2_name) "
        "VALUES (?,?,?,?,?)",
        (jid, time.strftime("%Y-%m-%d %H:%M:%S"), "queued",
         secure_filename(f1.filename) or "r1.fastq.gz",
         secure_filename(f2.filename) or "r2.fastq.gz"))
    db.commit()
    threading.Thread(target=process_job, args=(jid,), daemon=True).start()
    return redirect(url_for("job_detail", jid=jid))


@app.route("/jobs/<jid>")
def job_detail(jid):
    db = get_db()
    job = db.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
    if job is None:
        abort(404)
    analysis_row = db.execute(
        "SELECT * FROM fragment_analyses WHERE job_id=?", (jid,)).fetchone()
    fragment_result, stale_reason = _load_published(analysis_row, job)
    analyzing = analysis_row is not None and analysis_row["status"] == "running"
    return render_template(
        "job.html", job=job, outputs=OUTPUTS,
        output_files=pipeline.OUTPUT_FILES,
        fragment_result=fragment_result, stale_reason=stale_reason,
        analysis_row=analysis_row, analyzing=analyzing,
        fragment_inputs=FRAGMENT_INPUTS,
        default_threshold=DEFAULT_FRAGMENT_THRESHOLD,
        min_threshold=MIN_THRESHOLD, max_threshold=MAX_THRESHOLD)


@app.route("/jobs/<jid>/fragments", methods=["POST"])
def analyze_fragments(jid):
    db = get_db()
    job = db.execute("SELECT status FROM jobs WHERE id=?", (jid,)).fetchone()
    if job is None:
        abort(404)
    if job["status"] != "done":
        flash("只有已完成配对的任务才能做片段长度分析", "error")
        return redirect(url_for("job_detail", jid=jid))
    raw = (request.form.get("threshold") or "").strip()
    try:
        threshold = int(raw, 10)
    except ValueError:
        flash("阈值必须是整数（bp）", "error")
        return redirect(url_for("job_detail", jid=jid))
    if not MIN_THRESHOLD <= threshold <= MAX_THRESHOLD:
        flash(f"阈值必须在 {MIN_THRESHOLD}..{MAX_THRESHOLD} bp 之间", "error")
        return redirect(url_for("job_detail", jid=jid))
    try:
        start_fragment_analysis(jid, threshold)
    except fragments.AnalysisError as e:
        flash(str(e), "error")
    return redirect(url_for("job_detail", jid=jid))


@app.route("/jobs/<jid>/fragments/status")
def fragment_status(jid):
    """供详情页轮询：返回当前分析状态（不触发新分析）。"""
    db = get_db()
    job = db.execute("SELECT * FROM jobs WHERE id=?", (jid,)).fetchone()
    if job is None:
        abort(404)
    row = db.execute(
        "SELECT * FROM fragment_analyses WHERE job_id=?", (jid,)).fetchone()
    result, stale_reason = _load_published(row, job)
    resp = jsonify({"status": "none"}) if row is None else jsonify({
        "status": row["status"],
        "threshold": row["threshold"],
        "message": row["message"],
        "stale_reason": stale_reason,
        "result": result,
    })
    resp.headers["Cache-Control"] = "no-store"
    return resp


@app.route("/jobs/<jid>/download/<key>")
def download(jid, key):
    if key not in OUTPUTS:
        abort(404)
    job = get_db().execute("SELECT status FROM jobs WHERE id=?", (jid,)).fetchone()
    if job is None:
        abort(404)
    if job["status"] != "done":
        # 未完成/失败任务一律不暴露产物
        abort(409)
    path = os.path.join(JOBS_DIR, jid, "out", pipeline.OUTPUT_FILES[key])
    if not os.path.isfile(path):
        abort(404)
    return send_file(path, mimetype="application/gzip", as_attachment=True,
                     download_name=f"{jid}_{pipeline.OUTPUT_FILES[key]}")


@app.route("/healthz")
def healthz():
    return "ok"


@app.errorhandler(RequestEntityTooLarge)
def too_large(_e):
    flash(f"上传文件超过压缩大小上限（{MAX_UPLOAD_BYTES // (1024 * 1024)} MiB）",
          "error")
    return redirect(url_for("index")), 413


if __name__ == "__main__":
    init_and_recover()
    port = int(os.environ.get("PORT", "8000"))
    app.run(host="0.0.0.0", port=port, threaded=True)
