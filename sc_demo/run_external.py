#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""在 GSE310000（两例人 DRG snRNA-seq）上跑与发现队列同一套 ΔAUC。

与 replicate_external.py 的区别：本数据是整块 DRG 的 snRNA-seq，簇里混有
卫星胶质（SGC）、施万、成纤维等非神经元成分。为了与发现队列（13 类神经元
+ 1 类巨噬细胞）可比，这里先按标志基因把神经元簇挑出来，再分别输出：
  A) 全部簇（与发表口径同构：簇内全基因组中位数作零分布）
  B) 仅神经元簇（更贴近"感觉神经元亚型"这一命题）
两份结果都落盘，便于在论文里如实交代两种口径的差别。
"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "sc_demo"))
from replicate_external import TARGETS, DISCOVERY, delta_auc  # noqa: E402

ROOT = Path(__file__).resolve().parent
NEURON = ["NEFL", "NEFH", "PRPH", "STMN2", "ELAVL2", "POU4F1", "SCN9A", "SCN10A", "TUBB3"]
GLIA = ["KCNJ10", "SLC1A3", "GFAP", "MPZ", "SOX10", "CD74", "COL1A1", "DCN"]


def score(ad, genes):
    g = ad[:, [x for x in genes if x in ad.var_names]].X
    g = g.toarray() if hasattr(g, "toarray") else np.asarray(g)
    z = (g - g.mean(0)) / (g.std(0) + 1e-9)
    return z.mean(1)


def main():
    ad = sc.read_h5ad(ROOT / "GSE310000.h5ad")
    print("[读入] %d 细胞 × %d 基因" % ad.shape)
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=10)
    print("[质控] %d 细胞 × %d 基因" % ad.shape)

    raw = ad.X
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
    if np.allclose(raw[:20], np.round(raw[:20])):
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        print("[预处理] 原始计数 → CP10K + log1p")

    sub = ad.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
    sub = sub[:, sub.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(sub, resolution=1.0, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False)
    ad.obs["cl"] = sub.obs["cl"].astype(str).values
    ad.obs["neu_score"] = score(ad, NEURON)
    ad.obs["glia_score"] = score(ad, GLIA)

    prof = ad.obs.groupby("cl", observed=True).agg(
        n=("cl", "size"), neu=("neu_score", "mean"), glia=("glia_score", "mean"))
    prof["神经元"] = prof.neu > prof.glia
    print("\n[簇概况]\n%s" % prof.round(2).to_string())

    tg = [(g, n) for g, n in TARGETS if g in ad.var_names]

    out_all, sizes_all, p50_all = delta_auc(ad, tg)
    out_all.round(4).to_csv(ROOT / "ext_all_clusters.csv", index=False, encoding="utf-8-sig")
    print("\n[A] 全部 %d 个簇：全基因组 maxΔAUC 中位数 %.3f" % (len(sizes_all), p50_all))

    neu_cl = prof.index[prof["神经元"]].tolist()
    adn = ad[ad.obs.cl.isin(neu_cl)].copy()
    print("[B] 神经元簇 %d 个，细胞 %d 个" % (len(neu_cl), adn.n_obs))
    out_neu, sizes_neu, p50_neu = delta_auc(adn, tg)
    out_neu.round(4).to_csv(ROOT / "ext_neuron_clusters.csv", index=False, encoding="utf-8-sig")
    print("[B] 全基因组 maxΔAUC 中位数 %.3f" % p50_neu)

    for tag, out in (("全部簇", out_all), ("神经元簇", out_neu)):
        rows = []
        for g, name in tg:
            r = out[out.gene == g]
            if r.empty:
                continue
            r = r.iloc[0]
            dd, dp = DISCOVERY.get(g, (np.nan, np.nan))
            rows.append({"基因": g, "外部maxΔAUC": r.max_dAUC, "外部百分位": r.percentile,
                         "外部最高簇": r.best_cluster, "外部检出%": r.detect_pct,
                         "发现百分位": dp})
        df = pd.DataFrame(rows).dropna(subset=["发现百分位"]).sort_values("发现百分位", ascending=False)
        rho, p = spearmanr(df["外部百分位"], df["发现百分位"])
        print("\n===== %s 与发现队列对照（ρ=%.3f, p=%.4f）=====" % (tag, rho, p))
        print(df.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
