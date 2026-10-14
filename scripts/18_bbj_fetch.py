#!/usr/bin/env python3
"""
BBJ（BioBank Japan, hum0197）汇总统计并发下载器
=================================================

背景
----
humandbs.dbcls.jp 单连接吞吐实测仅 ~7 KB/s（1.5 GB 需 60 小时），不可用。
与 FinnGen GCS 的情形相同：限速是**每 TCP 连接**的，并发 Range 可线性叠加。

做法
----
1. HEAD 拿 content-length；
2. 切成 1 MB 的 Range 请求，N 个线程并发；
3. 每个 chunk 落成独立文件（可断点续传：已存在且大小正确的 chunk 直接跳过）；
4. 全部完成后按序拼接，清理分片。

用法
----
  python scripts/18_bbj_fetch.py test          # 只测聚合吞吐
  python scripts/18_bbj_fetch.py get BBJ.CL    # 下载单个表型
  python scripts/18_bbj_fetch.py all           # 下载核心 4 个
"""
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://humandbs.dbcls.jp/files/hum0197"
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "bbj")
CHUNK = 1 << 20          # 1 MB
THREADS = 32
UA = "curl/8.4.0"

# 核心表型：东亚暴露/结局 + 欧洲与跨血统对照
PHENOS = ["BBJ.BtC", "BBJ.CL", "EUR.BtC", "EUR.CL", "META.BtC", "META.CL"]


def url_of(pheno: str) -> str:
    return f"{BASE}/hum0197.v3.{pheno}.v1.zip"


def _opener():
    # 不走代理：本机 HTTP(S)_PROXY 会拦截并返回 502
    return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    with _opener().open(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def fetch_range(url: str, beg: int, end: int, tries: int = 4) -> bytes:
    """取 [beg, end] 闭区间字节。"""
    last = None
    for k in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": UA, "Range": f"bytes={beg}-{end}"}
            )
            with _opener().open(req, timeout=180) as r:
                data = r.read()
            if len(data) == end - beg + 1:
                return data
            last = RuntimeError(f"短读 {len(data)} != {end - beg + 1}")
        except Exception as e:                                  # noqa: BLE001
            last = e
            time.sleep(1.5 * (k + 1))
    raise RuntimeError(f"range {beg}-{end} 失败: {last}")


def download(pheno: str, threads: int = THREADS) -> str:
    url = url_of(pheno)
    size = head_size(url)
    part_dir = os.path.join(OUTDIR, f".parts_{pheno}")
    os.makedirs(part_dir, exist_ok=True)
    out = os.path.join(OUTDIR, f"hum0197.v3.{pheno}.v1.zip")

    n_chunks = (size + CHUNK - 1) // CHUNK
    jobs = []
    for i in range(n_chunks):
        beg = i * CHUNK
        end = min(beg + CHUNK - 1, size - 1)
        fp = os.path.join(part_dir, f"{i:06d}")
        if os.path.exists(fp) and os.path.getsize(fp) == end - beg + 1:
            continue                                             # 断点续传
        jobs.append((i, beg, end, fp))

    print(f"[{pheno}] {size/2**30:.2f} GB, {n_chunks:,} chunks, 待下 {len(jobs):,}", flush=True)
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(fetch_range, url, b, e): (i, fp) for i, b, e, fp in jobs}
        for fu in as_completed(futs):
            i, fp = futs[fu]
            with open(fp, "wb") as fh:
                fh.write(fu.result())
            done += 1
            if done % 200 == 0 or done == len(jobs):
                el = time.time() - t0
                got = done * CHUNK / 2**20
                print(f"  {done:,}/{len(jobs):,}  {got:.0f} MB  {el:.0f}s  "
                      f"{got/max(el,1e-9):.1f} MB/s", flush=True)

    # 拼接
    with open(out, "wb") as fo:
        for i in range(n_chunks):
            with open(os.path.join(part_dir, f"{i:06d}"), "rb") as fi:
                fo.write(fi.read())
    assert os.path.getsize(out) == size, "拼接后大小不符"
    for f in os.listdir(part_dir):
        os.remove(os.path.join(part_dir, f))
    os.rmdir(part_dir)
    print(f"[{pheno}] 完成 → {out}  ({os.path.getsize(out)/2**30:.2f} GB, "
          f"{time.time()-t0:.0f}s)", flush=True)
    return out


def _stream_count(url: str, beg: int, end: int, deadline: float, acc: list) -> None:
    """从 beg 顺序读直到 deadline，把累计字节写入 acc（按时间窗测速用）。"""
    got = 0
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": UA, "Range": f"bytes={beg}-{end}"}
        )
        with _opener().open(req, timeout=120) as r:
            while time.time() < deadline:
                buf = r.read(1 << 18)          # 256 KB
                if not buf:
                    break
                got += len(buf)
    except Exception:                                             # noqa: BLE001
        pass
    acc.append(got)


def test(seconds: int = 20) -> None:
    """按固定时间窗测聚合吞吐：N 个线程各流式读取，统计窗口内总字节。"""
    url = url_of("EUR.CL")
    size = head_size(url)
    print(f"测试目标 {url}  ({size/2**30:.2f} GB), 时间窗 {seconds}s", flush=True)
    stride = 20 << 20                                       # 各线程错开 20 MB
    for n in (1, 8, 32, 64):
        deadline = time.time() + seconds
        acc: list[int] = []
        with ThreadPoolExecutor(max_workers=n) as ex:
            for i in range(n):
                b = (i * stride) % (size - (1 << 26))
                ex.submit(_stream_count, url, b, size - 1, deadline, acc)
        total = sum(acc)
        print(f"  {n:>2} 线程: {total/2**20:8.2f} MB / {seconds}s = "
              f"{total/2**20/seconds:6.2f} MB/s", flush=True)


if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    cmd = sys.argv[1] if len(sys.argv) > 1 else "test"
    if cmd == "test":
        test()
    elif cmd == "get":
        download(sys.argv[2])
    elif cmd == "all":
        for p in PHENOS[:2]:        # 东亚核心两份优先
            download(p)
    else:
        print(__doc__)
