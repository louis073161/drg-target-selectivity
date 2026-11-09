#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 FinnGen 汇总统计聚合成基因层面指标（GRCh38，基因体 ±10 kb）。

指标
----
n_snp      窗口内 SNP 数
p_min      窗口内最小 P
mlogp_max  对应 -log10(P)
n_gws      达到全基因组显著（P<5e-8）的 SNP 数
z_max      最大 |beta/sebeta|
mean_chi2  平均 z²（受 LD 影响，仅作横向比较，同时输出按全基因组均值归一后的值）

用法：python gene_level.py finngen_R12_G6_MIGRAINE.gz [窗口kb]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
COORD = ROOT / "gene_coords_grch38.tsv"
WIN = 10_000
CHUNK = 3_000_000


def load_gwas(path):
    """分块读取，只留必要列，按染色体分组并排序。"""
    cols = ["#chrom", "pos", "pval", "beta", "sebeta"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=cols, chunksize=CHUNK,
                          dtype={"#chrom": str, "pos": np.int64,
                                 "pval": np.float64, "beta": np.float64,
                                 "sebeta": np.float64}):
        ch = ch.dropna(subset=["pval"])
        ch = ch[(ch.pval > 0) & (ch.pval <= 1)]
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True)
    df["chrom"] = df["#chrom"].str.replace("chr", "", regex=False)
    df = df[df.chrom.isin([str(i) for i in range(1, 23)] + ["X"])]
    return df


def main():
    if len(sys.argv) < 2:
        sys.exit("需要 GWAS 文件名")
    src = ROOT / sys.argv[1]
    global WIN
    if len(sys.argv) > 2:
        WIN = int(float(sys.argv[2]) * 1000)
    tag = src.name.replace("finngen_R12_", "").replace(".gz", "")

    coords = pd.read_csv(COORD, sep="\t")
    coords = coords[coords.chrom.isin([str(i) for i in range(1, 23)] + ["X"])]

    print("[读入] %s ..." % src.name, flush=True)
    g = load_gwas(src)
    print("  SNP %d，染色体 %d" % (len(g), g.chrom.nunique()), flush=True)

    out = []
    for chrom, sub in g.groupby("chrom", sort=False):
        sub = sub.sort_values("pos")
        pos = sub.pos.values
        pv = sub.pval.values
        beta = sub.beta.values
        se = sub.sebeta.values
        z = np.abs(beta / np.where(se > 0, se, np.nan))
        cz = z ** 2

        genes = coords[coords.chrom == chrom]
        if genes.empty:
            continue
        lo = np.searchsorted(pos, genes.start.values - WIN, side="left")
        hi = np.searchsorted(pos, genes.end.values + WIN, side="right")
        for (gene, a, b, gs, ge) in zip(genes.gene, lo, hi, genes.start, genes.end):
            n = b - a
            if n <= 0:
                out.append((gene, chrom, gs, ge, 0, np.nan, np.nan, 0, np.nan, np.nan))
                continue
            pp = pv[a:b]
            k = int(np.argmin(pp))
            out.append((gene, chrom, gs, ge, int(n), float(pp[k]),
                        float(-np.log10(pp[k])), int((pp < 5e-8).sum()),
                        float(np.nanmax(z[a:b])), float(np.nanmean(cz[a:b]))))

    res = pd.DataFrame(out, columns=["gene", "chrom", "start", "end", "n_snp",
                                     "p_min", "mlogp_max", "n_gws", "z_max", "mean_chi2"])
    # 全基因组平均 chi2 归一（近似 MAGMA 的 LD 校正思路）
    gmean = float(np.nanmean(res.mean_chi2))
    res["chi2_norm"] = res.mean_chi2 / gmean
    res.attrs = {}
    dst = ROOT / ("genelvl_%s.csv" % tag)
    res.to_csv(dst, index=False)
    print("  全基因组平均 χ² = %.3f" % gmean)
    print("  基因 %d，有 SNP 的 %d，达显著(任一 SNP P<5e-8) 的基因 %d"
          % (len(res), int((res.n_snp > 0).sum()), int((res.n_gws > 0).sum())))
    print("  -> %s" % dst)


if __name__ == "__main__":
    main()
