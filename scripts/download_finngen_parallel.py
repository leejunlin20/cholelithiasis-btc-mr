#!/usr/bin/env python3
"""
FinnGen 汇总统计并行分块下载器 v2（修好 v1 的致命缺陷）。

v1 的问题（2026-09-16 实测）：
  - 24MB 大分块 + --max-time 600s：单连接仅 ~10-100 KB/s，24MB 必然超时 → 反复归零
  - curl --retry 会截断重来，没有真正的断点续传

v2 方案：
  - 4MB 小分块 + 64 并发（实测 ~11 MB/s；128 并发可达 25 MB/s）
  - 真续传：每次尝试用「剩余区间」下载到临时文件再追加，绝不截断已有字节
  - --speed-limit/--speed-time 快速掐掉僵死连接
  - 全部完成 → 按序拼接 → gzip -t 校验

用法:
  python3 download_finngen_parallel.py <输出目录> <endpoint1> [endpoint2 ...]
"""
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import threading

BASE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
CONCURRENCY = 64
CHUNK_MB = 4
MAX_ATTEMPTS = 12
ATTEMPT_MAX_TIME = 240      # 单次尝试上限（秒）
STALL_LIMIT = 5 * 1024      # 5 KB/s
STALL_TIME = 45             # 持续 45 秒低于该速度即断开重试

_lock = threading.Lock()
_progress = {"done_bytes": 0, "total": 0}


def head_size(url):
    args = ["curl", "--noproxy", "*", "-sIL", url]
    out = subprocess.run(args, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":")[1].strip())
    return None


def fetch_chunk(url, start, end, path):
    """把 [start,end] 写入 path，支持断点续传（多次尝试累加）。

    关键：curl 用 `-o -` 输出到 stdout，Python 以 append 模式直接写入 part 文件。
    这样**不需要任何临时文件**，也就不会触发沙箱的批量删除保护，
    且每次尝试都从「已有字节数」处继续，天然断点续传。
    """
    expect = end - start + 1
    for _ in range(MAX_ATTEMPTS):
        have = os.path.getsize(path) if os.path.exists(path) else 0
        if have >= expect:
            return expect
        if have > expect or not os.path.exists(path):   # 异常：重建
            have = 0
            with open(path, "wb"):
                pass
        args = ["curl", "--noproxy", "*", "-s",
                "--max-time", str(ATTEMPT_MAX_TIME),
                "--speed-limit", str(STALL_LIMIT), "--speed-time", str(STALL_TIME),
                "-r", f"{start + have}-{end}", "-o", "-", url]
        with open(path, "ab") as w:
            subprocess.run(args, stdout=w, check=False)
    return os.path.getsize(path) if os.path.exists(path) else 0


def download(endpoint, outdir):
    url = f"{BASE}/finngen_R12_{endpoint}.gz"
    final = os.path.join(outdir, f"finngen_R12_{endpoint}.gz")
    if os.path.exists(final):
        if subprocess.run(["gzip", "-t", final]).returncode == 0:
            print(f"[skip] {endpoint} 已存在且校验通过 ({os.path.getsize(final)/1e6:.0f} MB)", flush=True)
            return True
        os.remove(final)

    total = head_size(url)
    if not total:
        print(f"[error] {endpoint} 取不到文件大小", flush=True)
        return False

    chunk = CHUNK_MB * 1024 * 1024
    ranges = [(i, min(i + chunk - 1, total - 1)) for i in range(0, total, chunk)]
    pdir = os.path.join(outdir, f".parts_{endpoint}")
    os.makedirs(pdir, exist_ok=True)
    parts = [os.path.join(pdir, f"{i:05d}.part") for i in range(len(ranges))]
    print(f"[start] {endpoint}: {total/1e6:.0f} MB / {len(ranges)} 块 / {CONCURRENCY} 并发", flush=True)

    done = [0]
    lock = threading.Lock()

    def task(n):
        got = fetch_chunk(url, ranges[n][0], ranges[n][1], parts[n])
        with lock:
            done[0] += 1
            if done[0] % 25 == 0 or done[0] == len(parts):
                cur = sum(os.path.getsize(p) for p in parts if os.path.exists(p))
                print(f"  {endpoint}: {done[0]}/{len(parts)} 块已结束，{cur/1e6:.0f}/{total/1e6:.0f} MB", flush=True)
        return got

    with cf.ThreadPoolExecutor(CONCURRENCY) as ex:
        list(ex.map(task, range(len(ranges))))

    got = sum(os.path.getsize(p) for p in parts if os.path.exists(p))
    if got != total:
        print(f"[error] {endpoint} 不完整 {got}/{total} MB，保留分块供下次续传", flush=True)
        return False

    print(f"  拼接中 {endpoint} ...", flush=True)
    with open(final, "wb") as w:
        for p in parts:
            with open(p, "rb") as r:
                while True:
                    b = r.read(16 << 20)
                    if not b:
                        break
                    w.write(b)
    if subprocess.run(["gzip", "-t", final]).returncode == 0:
        print(f"[ok] {endpoint} 完成并通过 gzip 校验 ({os.path.getsize(final)/1e6:.0f} MB)", flush=True)
        # 注意：分块目录**不在此处删除**——沙箱的「每轮 50 个删除」保护会在
        # rmtree 196 个文件时直接终止进程，导致后续 endpoint 不再下载。
        # 完成后手动清理：rm -rf <outdir>/.parts_<ENDPOINT>
        print(f"[cleanup] 请稍后手动删除分块目录: {pdir}", flush=True)
        return True
    print(f"[error] {endpoint} gzip 校验失败", flush=True)
    return False


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    outdir, endpoints = sys.argv[1], sys.argv[2:]
    os.makedirs(outdir, exist_ok=True)
    ok = True
    for ep in endpoints:                        # 不短路：一个失败不影响下一个
        ok = download(ep, outdir) and ok
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
