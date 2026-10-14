#!/usr/bin/env python3
"""
17_count_independent_signals.py —— 用 LD clumping 判定区域内独立信号数。

为什么还要用 clumping 与 susie 并行
-----------------------------------
【历史背景，已被修正】初版曾观察到 susie 在 R69（PNPLA3）2002 个 SNP 上分出
9 个"独立"信号、BTC 侧 ELBO 在 L=2 跳 +5225 / L=3 跌 −4919 震荡不收敛，一度
据此认为"稀疏 LD 让 susie 不可采信"。**该判断的根因已查明并修复**：FinnGen LD
矩阵(PLINK2 --r-unphased)的 r 符号以 **minor allele** 为参考，而 summary stats
的 beta 以 **alt** 为参考；af_alt>0.5 的位点两者反号（R69 占 16.9%）。修正符号后
susie 收敛到 1 个 credible set、PP.H4=0.994，与 coloc.abf 的 0.997 一致。
详见 scripts/16_coloc_susie.R 与 finngen-opengwas-data skill §2c-2。

因此本脚本**不是** susie 的替代品，而是**独立于 susie 的交叉验证**：clumping 只
依赖"两两 LD + p 值排序"，不做似然搜索，对 LD 稀疏与轻度非正定稳健，且是判断
"区域内是否存在多个独立信号"的标准做法。两者一致时结论最可靠。

⚠️ 注意：clumping 用 |r|² 判冗余，对符号不敏感；但**条件 z 必须用修正符号后的 r**，
否则 (z − r·z_lead) 会被反号污染。本脚本已加入等位基因方向修正。

产出 results/coloc_susie/independent_signals.tsv
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

BASE = os.environ.get("MR_BASE", ".")
OUT = os.path.join(BASE, "results/coloc_susie")
REGIONS = {
    "R69": ("PNPLA3", 22),
    "R51": ("SERPINA1", 14),
    "R67": ("HNF4A", 20),
    "R07": ("ABCG5/ABCG8 (negative control)", 2),
    "R54": ("FTO (positive control)", 16),
}
TRAITS = {1: "Cholelithiasis", 2: "Biliary tract cancer"}
P_GWS = 5e-8
R2_THRESHOLDS = [0.10, 0.01]     # 主用 0.10（|r|>0.316），敏感性 0.01（|r|>0.1）


def load_ld(rid: str, n: int, flip: np.ndarray | None = None):
    """读取稀疏 LD 并还原为 n×n 方阵。

    flip: 长度 n 的 ±1 向量，用于修正等位基因方向（见文件头说明）。
          等位基因翻转是相似变换 R → D·R·D，只改符号不改 |r|，
          故 clumping（用 |r|²）不受影响，但条件 z 必须用它。
    """
    p = f"/tmp/ld/{rid}_ld_long.tsv"
    if not os.path.exists(p):
        p = os.path.join(BASE, "data/ld", f"{rid}_ld_long.tsv.gz")
    d = pd.read_csv(p, sep="\t", header=None, names=["i1", "i2", "r", "r2"],
                    compression="gzip" if p.endswith(".gz") else None)
    R = np.eye(n)
    i1 = d.i1.values.astype(int)
    i2 = d.i2.values.astype(int)
    R[i1, i2] = d.r.values
    R[i2, i1] = d.r.values
    R = (R + R.T) / 2
    np.fill_diagonal(R, 1.0)
    if flip is not None:
        R = R * np.outer(flip, flip)
        np.fill_diagonal(R, 1.0)
    return R, len(d)


def clump(R, p, order_idx, r2_th):
    """贪心 clumping：按 p 升序取 lead，剔除与它 r²>r2_th 的位点。返回 [(lead, n_members)]。"""
    assigned = np.zeros(len(p), bool)
    groups = []
    for i in order_idx:
        if assigned[i]:
            continue
        members = np.where(np.abs(R[i, :]) ** 2 > r2_th)[0]
        members = members[~assigned[members]]
        groups.append((int(i), len(members)))
        assigned[np.abs(R[i, :]) ** 2 > r2_th] = True
    return groups


def main():
    rows = []
    for rid, (gene, chrom) in REGIONS.items():
        sf = os.path.join(BASE, "data/ld", f"{rid}_snps.tsv.gz")
        if not os.path.exists(sf):
            print(f"跳过 {rid}: 无 SNP 文件")
            continue
        s = pd.read_csv(sf, sep="\t")
        n = len(s)
        # ★ 等位基因方向修正：LD 的 r 以 minor allele 为参考，beta 以 alt 为参考
        flip = np.where(s["af"].values < 0.5, 1.0, -1.0)
        R, npairs = load_ld(rid, n, flip=flip)
        lo, hi = int(s.pos.min()), int(s.pos.max())
        print(f"\n=== {rid} {gene}  (chr{chrom}:{lo}-{hi}, {n} SNP, LD 配对 {npairs:,}) ===")
        print(f"    等位基因方向修正: {int((flip < 0).sum())}/{n} 个位点翻转 (af_alt>=0.5)")

        for which, tname in TRAITS.items():
            p = s[f"p_{which}"].values
            z = (s[f"beta_{which}"] / s[f"se_{which}"]).values
            sig = np.where(p < P_GWS)[0]
            if len(sig) == 0:
                print(f"  {tname:<22}: 无全基因组显著 SNP (min p = {p.min():.2e})")
                for r2t in R2_THRESHOLDS:
                    rows.append(dict(region=rid, gene=gene, trait=tname, r2_threshold=r2t,
                                     n_gws=0, n_independent=0, lead_snp="", lead_pos="",
                                     lead_p=float(p.min()), group_sizes="", lead_zs=""))
                continue
            order = sig[np.argsort(p[sig])]
            for r2t in R2_THRESHOLDS:
                gs = clump(R, p, order, r2t)
                lead_pos = [int(s.pos[i]) for i, _ in gs]
                # 条件 z：控制首个 lead 后，排除 |r|>0.99 的位点以避开 1/sqrt(1-r²) 数值爆炸
                k = order[0]
                rl = R[k, :]
                ok = np.abs(rl) < 0.99
                zc = np.full(n, np.nan)
                zc[ok] = (z[ok] - rl[ok] * z[k]) / np.sqrt(1 - rl[ok] ** 2)
                max_zc = float(np.nanmax(np.abs(zc))) if np.any(ok) else float("nan")
                if r2t == R2_THRESHOLDS[0]:
                    print(f"  {tname:<22}: GW 显著 {len(sig):>4} 个 SNP → 独立信号 "
                          f"{len(gs)} 个 (r²>{r2t})")
                    for j, (i, m) in enumerate(gs[:6], 1):
                        print(f"      {j}. {s.snp[i]:<28} pos={int(s.pos[i]):>10} "
                              f"p={p[i]:.2e} z={z[i]:+6.2f} 组内 {m} 个 SNP")
                    if len(gs) > 6:
                        print(f"      ... 另有 {len(gs)-6} 个")
                    print(f"      控制首要 lead 后: max|z_cond| = {max_zc:.2f} "
                          f"({'残余信号' if max_zc > 3.9 else '单信号可解释'})")
                rows.append(dict(region=rid, gene=gene, trait=tname, r2_threshold=r2t,
                                 n_gws=len(sig), n_independent=len(gs),
                                 lead_snp=",".join(s.snp[i] for i, _ in gs),
                                 lead_pos=",".join(str(int(s.pos[i])) for i, _ in gs),
                                 lead_p=float(p[order[0]]),
                                 group_sizes=",".join(str(m) for _, m in gs),
                                 lead_zs=",".join(f"{z[i]:+.2f}" for i, _ in gs),
                                 max_abs_z_cond=round(float(max_zc), 2)))

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "independent_signals.tsv"), sep="\t", index=False)
    print(f"\n已写出 {os.path.join(OUT,'independent_signals.tsv')}")


if __name__ == "__main__":
    main()
