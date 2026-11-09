#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""外部队列补充分析：簇标志基因、两样本拆分复算、图 5。"""
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "sc_demo"))
from replicate_external import TARGETS, DISCOVERY, delta_auc  # noqa: E402

ROOT = Path(__file__).resolve().parent
NEURON = ["NEFL", "NEFH", "PRPH", "STMN2", "ELAVL2", "POU4F1", "SCN9A", "SCN10A", "TUBB3"]
GLIA = ["KCNJ10", "SLC1A3", "GFAP", "MPZ", "SOX10", "CD74", "COL1A1", "DCN"]
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def score(ad, genes):
    g = ad[:, [x for x in genes if x in ad.var_names]].X
    g = g.toarray() if hasattr(g, "toarray") else np.asarray(g)
    z = (g - g.mean(0)) / (g.std(0) + 1e-9)
    return z.mean(1)


def build():
    ad = sc.read_h5ad(ROOT / "GSE310000.h5ad")
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=10)
    X = ad.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    if np.allclose(X[:20], np.round(X[:20])):
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
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
    return ad


def main():
    ad = build()
    prof = ad.obs.groupby("cl", observed=True).agg(
        n=("cl", "size"), neu=("neu_score", "mean"), glia=("glia_score", "mean"))
    prof["神经元"] = prof.neu > prof.glia

    sc.tl.rank_genes_groups(ad, "cl", method="wilcoxon", n_genes=8)
    mk = pd.DataFrame(ad.uns["rank_genes_groups"]["names"])
    prof["top_genes"] = ["/".join(mk[c].head(6)) for c in prof.index]
    prof = prof.sort_values("n", ascending=False)
    prof.round(3).to_csv(ROOT / "ext_cluster_markers.csv", encoding="utf-8-sig")
    print("\n[簇注释]\n%s" % prof.sort_index().to_string())

    nei = prof.index[prof["神经元"]].tolist()
    adn = ad[ad.obs.cl.isin(nei)].copy()
    tg = [(g, n) for g, n in TARGETS if g in adn.var_names]

    # ---- 两样本拆分 ----
    rows = {}
    for s in adn.obs["sample"].cat.categories:
        sub = adn[adn.obs["sample"] == s].copy()
        out, _, p50 = delta_auc(sub, tg)
        rows[s] = out.set_index("gene")["percentile"]
        print("[样本 %s] %d 细胞，全基因组 maxΔAUC 中位数 %.3f" % (s, sub.n_obs, p50))
    split = pd.DataFrame(rows)
    split["极差"] = split.max(1) - split.min(1)
    split.round(1).to_csv(ROOT / "ext_donor_split.csv", encoding="utf-8-sig")
    print("\n[两样本百分位]\n%s" % split.round(1).to_string())

    # ---- 图 5 ----
    out, _, p50 = delta_auc(adn, tg)
    rows2 = []
    for g, name in tg:
        r = out[out.gene == g].iloc[0]
        rows2.append({"基因": g, "外部百分位": r.percentile, "外部maxΔAUC": r.max_dAUC,
                      "发现百分位": DISCOVERY.get(g, (np.nan, np.nan))[1],
                      "外部检出": r.detect_pct, "簇": r.best_cluster})
    df = pd.DataFrame(rows2).dropna(subset=["发现百分位"]).sort_values("发现百分位")
    df.round(2).to_csv(ROOT / "ext_compare_neuron.csv", index=False, encoding="utf-8-sig")
    rho, p = spearmanr(df["外部百分位"], df["发现百分位"])

    fig, ax = plt.subplots(1, 2, figsize=(7.8, 3.72),
                           gridspec_kw={"width_ratios": [1.06, 1.0]})
    ax[0].scatter(df["发现百分位"], df["外部百分位"], s=64, c="#2b6cb0",
                  edgecolor="white", linewidth=0.9, zorder=3)
    crowd = df[(df["发现百分位"] >= 88) & (df["外部百分位"] >= 85)].sort_values(
        "外部百分位", ascending=False)
    for _, r in df.iterrows():
        lo = r["外部检出"]
        col = "#c53030" if lo < 20 else "#2b6cb0"
        ax[0].scatter(r["发现百分位"], r["外部百分位"], s=64, c=col,
                      edgecolor="white", linewidth=0.9, zorder=4)
    # 右上角拥挤区：标签移到绘图区右侧排成一列，用细引线连回数据点
    y0 = 100.0
    for k, (_, r) in enumerate(crowd.iterrows()):
        yl = y0 - k * 9.0
        ax[0].annotate(r["基因"], (r["发现百分位"], r["外部百分位"]),
                       xytext=(105, yl), textcoords="data", fontsize=7.6,
                       color="#2d3748", va="center", ha="left",
                       arrowprops=dict(arrowstyle="-", lw=0.42, color="#cbd5e0",
                                       shrinkA=2, shrinkB=4))
    for _, r in df.drop(crowd.index).iterrows():
        ax[0].annotate(r["基因"], (r["发现百分位"], r["外部百分位"]), fontsize=7.6,
                       xytext=(4.5, -3), textcoords="offset points", color="#2d3748")
    ax[0].plot([30, 100], [30, 100], ls=":", lw=1, c="#a0aec0")
    ax[0].set_xlim(30, 120); ax[0].set_ylim(0, 106)
    ax[0].set_xticks([40, 60, 80, 100])
    ax[0].set_xlabel("发现队列百分位（GSE249746）", fontsize=9.2)
    ax[0].set_ylabel("外部队列百分位（GSE310000）", fontsize=9.2)
    ax[0].set_title("A　两队列靶点定位的一致性　ρ＝%.2f" % rho, fontsize=10, loc="left")
    ax[0].tick_params(labelsize=8.2)

    ax[1].hist(out.max_dAUC, bins=70, color="#cbd5e0", edgecolor="none")
    ax[1].axvline(p50, ls="--", lw=1, c="#718096")
    ax[1].text(p50 + 0.009, ax[1].get_ylim()[1] * 0.935, "中位 %.3f" % p50, fontsize=8.2,
               color="#2d3748", bbox=dict(fc="white", ec="none", alpha=0.92, pad=1.6))
    for g in ["SCN11A", "SCN10A", "TRPV1", "P2RX3", "OPRM1", "ADORA1", "CNR1"]:
        v = out.loc[out.gene == g, "max_dAUC"]
        if len(v):
            ax[1].axvline(v.iloc[0], lw=1.1, c="#c53030", alpha=0.85)
    ax[1].set_xlabel("maxΔAUC（外部队列全部基因）", fontsize=9.2)
    ax[1].set_ylabel("基因数", fontsize=9.2)
    ax[1].set_title("B　外部队列的全基因组基准（红线示靶点）", fontsize=10, loc="left")
    ax[1].tick_params(labelsize=8.2)
    fig.tight_layout()
    fig.savefig(ROOT / "fig5_external_replication.png", dpi=400, facecolor="white", bbox_inches="tight")
    print("\n[输出] fig5_external_replication.png；ρ=%.3f p=%.4f" % (rho, p))
    print(df.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
