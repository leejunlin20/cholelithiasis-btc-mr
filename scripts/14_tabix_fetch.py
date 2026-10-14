#!/usr/bin/env python3
"""
远端 tabix 高速抓取工具（针对 FinnGen R12 LD 矩阵这类超大 long-format BGZF 文件）。

背景 / 为什么需要它
------------------
GCS(Google Cloud Storage) 单 TCP 连接实测只有 ~115 KB/s，而 htslib(pysam.TabixFile)
走的是单连接流式读，导致读一个 200 kb 区域要 170 秒。但同一个桶用 48-64 并发
HTTP Range 可以跑到 ~10 MB/s。

做法
----
1. 下载 <file>.tbi 索引（几十 KB），按 tabix 二进制格式解析出每染色体的
   bin → chunks(virtual offsets) 与 linear index。
2. 用 reg2bins 找到覆盖目标区间的 bins，收集并合并 chunks，得到 .gz 的
   **压缩字节区间**（voffset 高 48 位即块起始偏移）。
3. 多线程 HTTP Range 并行下载这些字节段 → 在内存里拼接。
4. 本地 BGZF(=多成员 gzip) 解压，逐行回调。截断的尾块用 EOFError 容忍处理。

命令行
------
    python 14_tabix_fetch.py --url <gz_url> --chrom 22 --start .. --end .. --out x.tsv.gz

也可作为库 import: fetch_region(url, chrom, start, end) -> Iterator[str]
"""
from __future__ import annotations

import gzip
import io
import os
import struct
import sys
import time
import urllib.request
import urllib.error
import concurrent.futures as cf
from typing import Iterator, Optional

# ---------------------------------------------------------------- HTTP 基础

_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
_UA = {"User-Agent": "curl/8.4"}


def _http_get(url: str, headers: dict | None = None, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={**_UA, **(headers or {})})
    with _OPENER.open(req, timeout=timeout) as r:
        return r.read()


def _remote_size(url: str) -> int:
    """取远端对象大小（用 1 字节 Range 探测 Content-Range）。"""
    req = urllib.request.Request(url, headers={**_UA, "Range": "bytes=0-0"})
    with _OPENER.open(req, timeout=60) as r:
        cr = r.headers.get("Content-Range", "")
        if "/" in cr:
            return int(cr.rsplit("/", 1)[1])
        cl = r.headers.get("Content-Length")
        if cl:
            return int(cl)
    raise RuntimeError("无法获取远端文件大小: " + url)


# ---------------------------------------------------------------- tabix 索引

def parse_tbi(blob: bytes) -> tuple[dict, dict]:
    """解析已解压的 .tbi 字节流 → (header, {ref_name: {'bins':{}, 'intv':[]}})。"""
    if blob[:4] != b"TBI\x01":
        raise ValueError("不是 tabix 索引 (magic 不符)")
    off = 4
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack_from("<8i", blob, off)
    off += 32
    names = blob[off:off + l_nm].split(b"\x00")
    off += l_nm
    refs: dict = {}
    for r in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", blob, off)
        off += 4
        bins: dict[int, list[tuple[int, int]]] = {}
        for _ in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<Ii", blob, off)
            off += 8
            chunks = list(struct.unpack_from("<%dQ" % (2 * n_chunk), blob, off))
            off += 16 * n_chunk
            bins[bin_id] = [(chunks[i], chunks[i + 1]) for i in range(0, len(chunks), 2)]
        (n_intv,) = struct.unpack_from("<i", blob, off)
        off += 4
        intv = list(struct.unpack_from("<%dQ" % n_intv, blob, off))
        off += 8 * n_intv
        refs[names[r].decode()] = {"bins": bins, "intv": intv}
    header = {"format": fmt, "col_seq": col_seq, "col_beg": col_beg,
              "col_end": col_end, "meta": meta, "skip": skip}
    return header, refs


def _reg2bins(beg: int, end: int) -> list[int]:
    """tabix 标准 reg2bins（0-based half-open）→ 覆盖该区间的所有 bin 号。"""
    if beg < 0:
        beg = 0
    if end <= beg:
        end = beg + 1
    end -= 1
    bins = [0]
    for start, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        lo = start + (beg >> shift)
        hi = start + (end >> shift)
        bins.extend(range(lo, hi + 1))
    return bins


def region_chunks(ref: dict, start: int, end: int) -> list[list[int]]:
    """返回覆盖 [start,end) 的、合并后的 [(byte_start, byte_end), ...]（字节为单位）。"""
    chunks: list[tuple[int, int]] = []
    for b in _reg2bins(start, end):
        chunks.extend(ref["bins"].get(b, []))
    if not chunks:
        return []
    chunks.sort()
    merged: list[list[int]] = []
    for beg, fin in chunks:
        if merged and beg <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], fin)
        else:
            merged.append([beg, fin])
    out = []
    for beg, fin in merged:
        b0 = beg >> 16                      # 块起始（BGZF 成员边界）
        b1 = (fin >> 16) + 1                # 结束块之后的第一个字节
        out.append([b0, b1])
    return out


# ---------------------------------------------------------------- BGZF 下载

def _download_range(url: str, b0: int, b1: int, threads: int = 48,
                    chunk: int = 1024 * 1024) -> bytes:
    """并发分块下载 [b0, b1) 字节并拼接。"""
    total = b1 - b0
    if total <= 0:
        return b""
    n = (total + chunk - 1) // chunk
    parts: list[Optional[bytes]] = [None] * n

    def one(i: int) -> tuple[int, bytes]:
        s = b0 + i * chunk
        e = min(s + chunk, b1) - 1
        for attempt in range(4):
            try:
                return i, _http_get(url, {"Range": f"bytes={s}-{e}"}, timeout=90)
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(1.5 * (attempt + 1))
        raise RuntimeError("unreachable")

    with cf.ThreadPoolExecutor(min(threads, n)) as ex:
        for i, data in ex.map(one, range(n)):
            parts[i] = data
    return b"".join(p for p in parts if p)


def bgzf_lines(buf: bytes) -> Iterator[str]:
    """按 BGZF(=多成员 gzip) 解压 buf，逐行 yield。容忍尾部截断。"""
    stream = io.BytesIO(buf)
    dec = gzip.GzipFile(fileobj=stream)
    pending = b""
    while True:
        try:
            blk = dec.read(8 * 1024 * 1024)
        except (EOFError, OSError, gzip.BadGzipFile):
            blk = b""
        if not blk:
            break
        pending += blk
        *lines, pending = pending.split(b"\n")
        for ln in lines:
            yield ln.decode("utf-8", "replace")
    if pending:
        yield pending.decode("utf-8", "replace")


class RemoteTabix:
    """远端 .gz + .tbi 的轻量查询器（按基因组区间取记录，并发下载）。"""

    def __init__(self, gz_url: str, tbi_url: str | None = None, verbose: bool = True):
        self.gz_url = gz_url
        self.tbi_url = tbi_url or (gz_url + ".tbi")
        t0 = time.time()
        raw = _http_get(self.tbi_url, timeout=120)
        self.header, self.refs = parse_tbi(gzip.decompress(raw))
        self.size = _remote_size(gz_url)
        if verbose:
            print(f"  [tabix] 索引就绪 {len(self.refs)} 条染色体, "
                  f"目标 {self.size/1e6:.0f} MB, 耗时 {time.time()-t0:.1f}s", flush=True)

    def chunks_for(self, chrom: str, start: int, end: int) -> list[list[int]]:
        ref = self.refs.get(str(chrom))
        if ref is None:
            return []
        return region_chunks(ref, start, end)

    def fetch(self, chrom: str, start: int, end: int, threads: int = 48,
              merge_gap: int = 2 * 1024 * 1024) -> Iterator[str]:
        """yield 原始行（未做 pos 过滤），覆盖 [start, end) 区间。"""
        cs = self.chunks_for(chrom, start, end)
        if not cs:
            return
        # 合并靠得很近的块，减少 Range 请求次数
        merged: list[list[int]] = []
        for b0, b1 in cs:
            if merged and b0 - merged[-1][1] <= merge_gap:
                merged[-1][1] = max(merged[-1][1], b1)
            else:
                merged.append([b0, b1])
        for b0, b1 in merged:
            b0 = max(0, b0)
            b1 = min(self.size, b1)
            if b1 - b0 > 400 * 1024 * 1024:
                raise RuntimeError(f"请求字节段过大 ({ (b1-b0)/1e6:.0f} MB)，区域可能过宽")
            buf = _download_range(self.gz_url, b0, b1, threads=threads)
            for line in bgzf_lines(buf):
                yield line


# ---------------------------------------------------------------- CLI

def _cli() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="远端 tabix 区间抓取")
    ap.add_argument("--url", required=True)
    ap.add_argument("--tbi", default=None)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    tb = RemoteTabix(a.url, a.tbi)
    t0 = time.time()
    n = 0
    fh = gzip.open(a.out, "wt") if a.out else sys.stdout
    cs = tb.chunks_for(a.chrom, a.start, a.end)
    print(f"  [chunks] {len(cs)} 段, 合计 {sum(b-a for a,b in cs)/1e6:.1f} MB",
          file=sys.stderr, flush=True)
    for line in tb.fetch(a.chrom, a.start, a.end, threads=a.threads):
        fh.write(line + "\n")
        n += 1
    if a.out:
        fh.close()
    print(f"  写出 {n:,} 行, 耗时 {time.time()-t0:.1f}s", file=sys.stderr, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
