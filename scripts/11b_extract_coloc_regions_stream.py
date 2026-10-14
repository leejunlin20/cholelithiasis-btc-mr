#!/usr/bin/env python3
"""
11b_extract_coloc_regions_stream.py
-----------------------------------
绕过本机磁盘慢读，直接从 GCS 并行流式下载 FinnGen R12 summary stats，
在内存中解压并抽取 COLOC 区域内的变异。

背景（2026-09-17 实测）：
  - 本机 Python `os.read` 全文件：冷读≈1 MB/s（661MB 需 >10 min），缓存命中则 1388 MB/s，
    且缓存不保留 → 反复全量扫描不可行
  - `gzcat`（macOS BSD gzip）更慢：936s 墙钟 / 仅 9.8s CPU
  - 而从 GCS 并行分块下载实测 ≈11 MB/s（4MB 块 × 64 并发）
  故：网络获取 → 内存解压 → 本地只写结果（约 40MB/trait）

用法:
  python 11b_extract_coloc_regions_stream.py all
  python 11b_extract_coloc_regions_stream.py CHOLE
"""
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys
import threading
import time
import zlib

BASE = os.environ.get("MR_BASE", ".")
REG = os.path.join(BASE, "results/coloc/regions.tsv")
TMP = "/tmp/fg_coloc"
OUTDIR = os.path.join(BASE, "data/coloc")

STORAGE = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release"
TRAITS = {
    "CHOLE": "K11_CHOLELITH",
    "BTC":   "C3_BILIARY_GALLBLADDER_EXALLC",
    "BMI":   "BMI_IRN",
}

CHUNK = 4 * 1024 * 1024
CONC = 64
WINDOW = 128
MAX_ATTEMPTS = 12
HEADER = (b"region_id\tchrom\tpos\tref\talt\trsids\tnearest_genes\t"
          b"pval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls\n")


def get_size(url):
    out = subprocess.run(["curl", "--noproxy", "*", "-sIL", url],
                         capture_output=True, text=True).stdout
    for line in out.splitlines():
        if line.lower().startswith("content-length:"):
            return int(line.split(":")[1].strip())
    return None


def fetch(url, start, end):
    """下载 [start,end] 字节区间，返回 bytes；失败时抛异常。"""
    want = end - start + 1
    for _ in range(MAX_ATTEMPTS):
        args = ["curl", "--noproxy", "*", "-s", "--max-time", "180",
                "--speed-limit", "5120", "--speed-time", "45",
                "-r", f"{start}-{end}", "-o", "-", url]
        r = subprocess.run(args, capture_output=True)
        if r.returncode == 0 and len(r.stdout) == want:
            return r.stdout
    raise RuntimeError(f"chunk {start}-{end} 下载失败")


def load_regions(path):
    regs = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            regs.setdefault(p[1].encode(), []).append((int(p[2]), int(p[3]), p[0].encode()))
    for c in regs:
        regs[c].sort()
    return regs


class LineFilter:
    """按行切分解压输出，只保留落在 COLOC 区域内的行。"""

    def __init__(self, regs, fout):
        self.regs = regs
        self.fout = fout
        self.tail = b""
        self.first = True
        self.n_total = 0
        self.n_kept = 0
        self.obuf = []
        self.cur_chrom = None
        self.arr = ()

    def emit(self, data):
        buf = self.tail + data
        start = 0
        while True:
            nl = buf.find(b"\n", start)
            if nl < 0:
                break
            self._line(buf[start:nl + 1])
            start = nl + 1
        self.tail = buf[start:]
        if len(self.obuf) >= 20000:
            self.flush()

    def _line(self, line):
        if self.first:
            self.first = False
            if line.startswith(b"#"):
                self.fout.write(b"region_id\t" + line)
                return
            self.fout.write(HEADER)
        self.n_total += 1
        i = line.find(b"\t")
        if i <= 0:
            return
        chrom = line[:i]
        j = line.find(b"\t", i + 1)
        if j <= i:
            return
        try:
            pos = int(line[i + 1:j])
        except ValueError:
            return
        if chrom != self.cur_chrom:
            self.cur_chrom = chrom
            self.arr = self.regs.get(chrom, ())
        if self.arr:
            for s, e, rid in self.arr:
                if s <= pos <= e:
                    self.obuf.append(rid + b"\t" + line)
                    self.n_kept += 1
                    return

    def flush(self):
        if self.obuf:
            self.fout.write(b"".join(self.obuf))
            self.obuf.clear()

    def close(self):
        if self.tail:
            self._line(self.tail)
        self.flush()


def run(trait):
    endpoint = TRAITS[trait]
    url = f"{STORAGE}/finngen_R12_{endpoint}.gz"
    tmp_out = os.path.join(TMP, f"{trait}_regions.tsv")
    final_out = os.path.join(OUTDIR, f"{trait}_regions.tsv")
    os.makedirs(TMP, exist_ok=True)
    os.makedirs(OUTDIR, exist_ok=True)

    total = get_size(url)
    nchunks = (total + CHUNK - 1) // CHUNK
    regs = load_regions(REG)

    t0 = time.time()
    print(f"[{trait}] {total/1e6:.0f} MB / {nchunks} 块 / {CONC} 并发 → 开始", flush=True)

    d = zlib.decompressobj(31)
    finished = [0]
    lock = threading.Lock()

    with open(tmp_out, "wb") as fout:
        lf = LineFilter(regs, fout)
        with cf.ThreadPoolExecutor(CONC) as ex:
            futs = {}
            nxt_sub = 0
            nxt_con = 0
            while nxt_con < nchunks:
                while nxt_sub < nchunks and len(futs) < WINDOW:
                    s = nxt_sub * CHUNK
                    e = min(s + CHUNK, total) - 1
                    futs[nxt_sub] = ex.submit(fetch, url, s, e)
                    nxt_sub += 1
                data = futs.pop(nxt_con).result()
                nxt_con += 1
                with lock:
                    finished[0] += 1
                    if finished[0] % 25 == 0 or finished[0] == nchunks:
                        print(f"  [{trait}] {finished[0]}/{nchunks} 块, 已保留 {lf.n_kept:,} 变异, "
                              f"{time.time()-t0:.0f}s", flush=True)
                while data:
                    outb = d.decompress(data)
                    if outb:
                        lf.emit(outb)
                    rest = d.unused_data
                    if d.eof:
                        d = zlib.decompressobj(31)
                        data = rest
                    else:
                        data = b""
        lf.close()

    shutil.move(tmp_out, final_out)
    print(f"[{trait}] 完成: 扫描 {lf.n_total:,} 行 → 保留 {lf.n_kept:,} 变异, "
          f"耗时 {time.time()-t0:.0f}s → {final_out}", flush=True)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for t in (TRAITS.keys() if which == "all" else [which]):
        run(t)
