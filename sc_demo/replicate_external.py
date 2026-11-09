#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""外部独立重复：把与发现队列完全相同的 ΔAUC 流程，套用到第二套人类 DRG 数据上。

设计要点
--------
外部队列不需要与本队列注释同名。ΔAUC 的关键量是「一个基因的最大 ΔAUC 在
该数据集自身的全基因组分布中处于什么位置」，这一步只依赖簇的划分粒度，
不依赖簇的命名。因此脚本默认在外部数据上自行聚类，再输出：
  1) 每个靶点的 maxΔAUC、全基因组百分位、最高富集簇及其标志基因；
  2) 该数据集全基因组 maxΔAUC 分布的中位数（与发现队列的 0.112 对照）；
  3) 与发现队列百分位向量的一致性（Spearman ρ）。
若外部数据自带注释（--label-col），则直接使用，省去聚类。

用法
----
  python replicate_external.py --h5ad GSE201586_human.h5ad
  python replicate_external.py --h5ad ext.h5ad --label-col cell_type --figure
  python replicate_external.py --10x-h5 raw_feature_bc_matrix.h5 --resolution 1.2

输入要求：行为细胞、列为基因的表达矩阵（.h5ad / 10x .h5 / .mtx 目录）。
若矩阵为原始 UMI，脚本会自动执行 normalize_total(1e4) + log1p。
"""
import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

TARGETS = [
    ("SCN11A", "Nav1.9"), ("SCN10A", "Nav1.8"), ("TRPV1", "TRPV1"), ("P2RX3", "P2X3"),
    ("OPRM1", "μ 阿片受体"), ("OPRD1", "δ 阿片受体"), ("OPRK1", "κ 阿片受体"),
    ("ADORA1", "腺苷 A1 受体"), ("CNR1", "大麻素 CB1"), ("OPRL1", "痛敏肽受体"),
    ("CNR2", "大麻素 CB2"), ("GPR55", "GPR55"),
]

# 发现队列（GSE249746, 1 129 个细胞, 16 969 个基因）已发表值，用于并列比较
DISCOVERY = {
    "SCN11A": (0.495, 99.95), "TRPV1": (0.424, 98.8), "P2RX3": (0.406, 98.5),
    "SCN10A": (0.397, 98.2), "OPRM1": (0.379, 97.6), "OPRD1": (0.335, 95.8),
    "CNR1": (0.238, 87.3), "OPRL1": (0.119, 53.5), "ADORA1": (0.117, 52.4),
    "OPRK1": (0.096, 41.3), "CNR2": (np.nan, np.nan), "GPR55": (np.nan, np.nan),
}
DISCOVERY_P50 = 0.112


def load(path: Path, label_col=None, resolution=1.0, min_cells=200, min_genes=200):
    import scanpy as sc

    if path.suffix == ".h5ad":
        ad = sc.read_h5ad(path)
    elif path.suffix in (".h5", ".hdf5"):
        ad = sc.read_10x_h5(str(path))
        ad.var_names_make_unique()
    elif path.is_dir():
        ad = sc.read_10x_mtx(str(path), var_names="gene_symbols", make_unique=True)
    else:
        sys.exit("无法识别的输入格式：%s（支持 .h5ad / 10x .h5 / 10x mtx 目录）" % path)

    # 基本过滤
    sc.pp.filter_cells(ad, min_genes=min_genes)
    sc.pp.filter_genes(ad, min_cells=3)
    print("[输入] %d 个细胞 × %d 个基因（过滤后）" % ad.shape)

    X = ad.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    if not np.issubdtype(X.dtype, np.floating):
        X = X.astype(np.float32)
    # 未归一化的原始计数：整数且最大值远大于 30
    if np.allclose(X[:50], np.round(X[:50])) and X.max() > 30:
        print("[预处理] 判定为原始计数，执行 CP10K + log1p")
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
    else:
        print("[预处理] 判定为已归一化，直接使用")

    if label_col and label_col in ad.obs:
        print("[聚类] 使用自带注释列：%s（%d 类）" % (label_col, ad.obs[label_col].nunique()))
        ad.obs["cl"] = ad.obs[label_col].astype(str)
    else:
        # 聚类在 HVG 子集上做，标签回填到保留全部基因的原矩阵（ΔAUC 需要全基因组参照）
        sub = ad.copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
        sub = sub[:, sub.var.highly_variable].copy()
        sc.pp.scale(sub, max_value=10)
        sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
        sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
        sc.tl.leiden(sub, resolution=resolution, key_added="cl", flavor="igraph",
                     n_iterations=2, directed=False)
        ad.obs["cl"] = sub.obs["cl"].astype(str).values
        print("[聚类] leiden(res=%.2f) 得 %d 个簇" % (resolution, ad.obs.cl.nunique()))
    return ad


def delta_auc(ad, targets, min_detect=0.05):
    """与发现队列逐行一致的实现：全基因秩变换 → AUC → 各簇经验零分布 → ΔAUC。"""
    X = ad.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X, dtype=np.float32)
    genes = ad.var_names.to_numpy()
    det = (X > 0).mean(0)
    keep = np.where(det >= min_detect)[0]
    pan = [ad.var_names.get_loc(g) for g, _ in targets if g in ad.var_names]
    keep = np.union1d(keep, np.array(pan, dtype=int)) if len(pan) else keep
    E = X[:, keep].astype(np.float32)
    EG = genes[keep]

    R = rankdata(E, axis=0)                    # 每个基因在全部细胞中的秩
    cl = ad.obs["cl"].to_numpy()
    cats = sorted(pd.unique(cl), key=lambda c: -int((cl == c).sum()))
    n_all = len(cl)

    def auc_cols(mask):
        n1 = int(mask.sum())
        n0 = n_all - n1
        return (R[mask].mean(0) - (n1 + 1) / 2) / n0

    base, auc_by = {}, {}
    for c in cats:
        a = auc_cols(cl == c)
        auc_by[c] = a
        base[c] = np.median(a)                 # 该簇的全基因 AUC 中位数

    dmax = np.full(len(EG), -np.inf)
    best = np.empty(len(EG), dtype=object)
    for c in cats:
        d = auc_by[c] - base[c]
        upd = d > dmax
        dmax[upd] = d[upd]
        best[upd] = c

    df = pd.DataFrame({"gene": EG, "max_dAUC": dmax, "best_cluster": best,
                       "detect_pct": det[keep] * 100})
    df["percentile"] = df.max_dAUC.rank(pct=True) * 100
    return df, {c: int((cl == c).sum()) for c in cats}, float(np.median(dmax))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, help=".h5ad / 10x .h5 / 10x mtx 目录")
    ap.add_argument("--10x-h5", dest="h5", type=Path, help="10x 输出目录或 h5")
    ap.add_argument("--label-col", default=None, help="外部数据自带的注释列名")
    ap.add_argument("--resolution", type=float, default=1.0)
    ap.add_argument("--out", type=Path, default=Path("external_replication.csv"))
    ap.add_argument("--figure", action="store_true")
    a = ap.parse_args()
    path = a.h5ad or a.h5
    if path is None:
        sys.exit("请用 --h5ad 或 --10x-h5 指定外部数据集")

    ad = load(path, a.label_col, a.resolution)
    tg = [(g, n) for g, n in TARGETS if g in ad.var_names]
    missing = [g for g, _ in TARGETS if g not in ad.var_names]
    if missing:
        print("[警告] 数据中缺少：%s（按缺失处理）" % ", ".join(missing))

    df, sizes, p50 = delta_auc(ad, tg)
    print("\n[外部队列] 簇大小：%s" % sizes)
    print("[外部队列] 全基因组 maxΔAUC 中位数 = %.3f（发现队列 %.3f）" % (p50, DISCOVERY_P50))

    rows = []
    for g, name in tg:
        r = df[df.gene == g]
        if r.empty:
            continue
        r = r.iloc[0]
        dd, dp = DISCOVERY.get(g, (np.nan, np.nan))
        rows.append({"基因": g, "靶点": name,
                     "外部maxΔAUC": r.max_dAUC, "外部百分位": r.percentile,
                     "外部最高簇": r.best_cluster, "外部检出率%": r.detect_pct,
                     "发现maxΔAUC": dd, "发现百分位": dp,
                     "百分位差": r.percentile - dp})
    out = pd.DataFrame(rows).sort_values("发现百分位", ascending=False)

    ok = out.dropna(subset=["发现百分位"])
    rho, p = spearmanr(ok["外部百分位"], ok["发现百分位"])
    print("\n" + out.round(2).to_string(index=False))
    print("\n百分位一致性 Spearman ρ = %.3f (p = %.4f, n = %d)" % (rho, p, len(ok)))
    hi = ok[ok.发现百分位 >= 95]
    lo = ok[ok.发现百分位 <= 60]
    print("原判为高选择性（发现百分位≥95）的 %d 个靶点，外部百分位落在 %s"
          % (len(hi), ", ".join("%.1f" % v for v in hi.外部百分位)))
    print("原判为无选择性（发现百分位≤60）的 %d 个靶点，外部百分位落在 %s"
          % (len(lo), ", ".join("%.1f" % v for v in lo.外部百分位)))

    out.round(4).to_csv(a.out, index=False, encoding="utf-8-sig")
    print("\n已写出：%s" % a.out)

    if a.figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        matplotlib.rcParams["axes.unicode_minus"] = False
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), dpi=300)
        ax[0].scatter(ok.发现百分位, ok.外部百分位, s=70, c="#2b6cb0",
                      edgecolor="white", linewidth=1, zorder=3)
        for _, r in ok.iterrows():
            ax[0].annotate(r.基因, (r.发现百分位, r.外部百分位), fontsize=7.5,
                           xytext=(4, 4), textcoords="offset points")
        ax[0].plot([0, 100], [0, 100], ls=":", c="#888", lw=1)
        ax[0].set_xlabel("发现队列百分位"); ax[0].set_ylabel("外部队列百分位")
        ax[0].set_title("A  百分位一致性  ρ = %.2f" % rho, fontsize=10, loc="left")
        ax[0].set_xlim(30, 102); ax[0].set_ylim(30, 102)
        ax[1].hist(df.max_dAUC, bins=80, color="#cbd5e0", edgecolor="none")
        for g, name in tg:
            v = df.loc[df.gene == g, "max_dAUC"]
            if len(v):
                ax[1].axvline(v.iloc[0], lw=1.2, color="#c53030")
        ax[1].set_xlabel("maxΔAUC（外部队列全部基因）"); ax[1].set_ylabel("基因数")
        ax[1].set_title("B  外部队列的全基因组基准", fontsize=10, loc="left")
        fig.tight_layout()
        fig.savefig("fig5_external_replication.png", bbox_inches="tight")
        print("已写出：fig5_external_replication.png")


if __name__ == "__main__":
    main()
