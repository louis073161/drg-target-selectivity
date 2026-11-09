#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
人类 DRG 图谱侧独立分析 —— 步骤 1：重建、聚类、标志基因
固定随机种子，全流程可复现。
"""
import warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse

t0 = time.time()
SEED = 0
np.random.seed(SEED)

print("[1] 读入原始计数矩阵 ...", flush=True)
df = pd.read_csv("GSE249746_Expression_matrix_raw_counts.csv.gz", index_col=0)
df.columns = [c.strip('"') for c in df.columns]
print("    原始维度:", df.shape, flush=True)

adata = ad.AnnData(sparse.csr_matrix(df.values.astype("float32").T))
adata.obs_names = df.columns.tolist()
adata.var_names = df.index.astype(str).tolist()
adata.var_names_make_unique()
del df

# 条形码结构 N2-RL5-1：第1段=供体，第2段=取材节段
adata.obs["donor"] = [b.split("-")[0] if "-" in b else "NA" for b in adata.obs_names]
adata.obs["level"] = [b.split("-")[1] if b.count("-") == 2 else "NA" for b in adata.obs_names]
adata.obs["donor"] = adata.obs["donor"].astype("category")
adata.obs["level"] = adata.obs["level"].astype("category")
print("    供体构成:", adata.obs["donor"].value_counts().to_dict(), flush=True)
print("    节段构成:", adata.obs["level"].value_counts().to_dict(), flush=True)

adata.layers["counts"] = adata.X.copy()
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)

print("[2] 质控 ...", flush=True)
print("    QC 前:", adata.shape, flush=True)
print("    pct_mt 分位:", np.percentile(adata.obs.pct_counts_mt, [0, 5, 50, 95, 100]).round(1), flush=True)
min_genes, min_counts, max_mt = 1000, 5000, 50
m = (adata.obs.n_genes_by_counts > min_genes) & (adata.obs.total_counts > min_counts) & (adata.obs.pct_counts_mt < max_mt)
adata = adata[m].copy()
print("    QC 后:", adata.shape, "| 阈值: n_genes>%d, counts>%d, pct_mt<%d" % (min_genes, min_counts, max_mt), flush=True)

print("[3] 标准化 / HVG / PCA / 邻居图 ...", flush=True)
adata.layers["counts"] = adata.layers["counts"]  # 保留原始计数
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
hv = adata[:, adata.var.highly_variable].copy()
n_pcs = min(30, hv.n_obs - 1, hv.n_vars - 1)
sc.pp.scale(hv, max_value=10)
sc.tl.pca(hv, n_comps=n_pcs, svd_solver="arpack", random_state=SEED)
adata.obsm["X_pca"] = hv.obsm["X_pca"]
print("    PC 数:", n_pcs, "| 方差解释前5:", hv.uns["pca"]["variance_ratio"][:5].round(3), flush=True)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden", random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
print("    Leiden 簇数:", adata.obs.leiden.nunique(), flush=True)

print("[4] 各簇标志基因 (wilcoxon) ...", flush=True)
adata.obs["leiden"] = adata.obs["leiden"].astype(str)
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=12, use_raw=True)
res = adata.uns["rank_genes_groups"]
sz = adata.obs.leiden.value_counts()
for c in sorted(res["names"].dtype.names, key=lambda x: int(x)):
    genes = [res["names"][c][i] for i in range(12)]
    lfc = [res["logfoldchanges"][c][i] for i in range(12)]
    top = ", ".join("%s(%+.1f)" % (g, l) for g, l in zip(genes[:8], lfc[:8]))
    print("  簇%-2s n=%-4d %s" % (c, sz[c], top), flush=True)

adata.write_h5ad("atlas_s1_clustered.h5ad")
print("[完成] %.1fs -> atlas_s1_clustered.h5ad" % (time.time() - t0), flush=True)
