# -*- coding: utf-8 -*-
"""小鼠侧参照系 · 第一步：与人类完全同参重聚类。

人类侧参数（atlas_s1_build.py）：HVG 2000 (seurat)、scale max 10、PCA 30、
neighbors k=15、Leiden res=1.0、seed 0。小鼠侧原先用 res=1.2、k=15，
此处改为同参，并附 0.8 / 1.2 两档作敏感性。
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import sparse

SEED = 0
np.random.seed(SEED)
sc.settings.verbosity = 1

a = sc.read_h5ad("mdrg_neuron_59739.h5ad")
X = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X, dtype=np.float32)
a.X = X.astype(np.float32)
a.var_names = [v for v in a.var_names]
print("输入:", a.shape, "| 细胞:", a.n_obs, flush=True)

# X 已是 log1p(FPKM)
hv = a.copy()
sc.pp.highly_variable_genes(hv, n_top_genes=2000, flavor="seurat")
hv = hv[:, hv.var.highly_variable].copy()
sc.pp.scale(hv, max_value=10)
n_pcs = min(30, hv.n_obs - 1, hv.n_vars - 1)
sc.tl.pca(hv, n_comps=n_pcs, svd_solver="arpack", random_state=SEED)
a.obsm["X_pca"] = hv.obsm["X_pca"]
sc.pp.neighbors(a, n_neighbors=15, n_pcs=n_pcs, random_state=SEED)
sc.tl.umap(a, random_state=SEED)

for res in (0.8, 1.0, 1.2):
    key = "leiden_%.1f" % res
    sc.tl.leiden(a, resolution=res, key_added=key, random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    a.obs[key] = a.obs[key].astype(str)
    print("\n===== res=%.1f  %d 簇 =====" % (res, a.obs[key].nunique()), flush=True)
    print(a.obs[key].value_counts().sort_index().to_dict(), flush=True)

# 主分析取 res=1.0
a.obs["leiden"] = a.obs["leiden_1.0"].values

PANEL = ["Calca", "Tac1", "Nmb", "Sst", "Mrgprd", "Mrgpra3", "Mrgprb4", "Nppb", "Zcchc12",
         "Th", "Fam19a4", "Ntrk1", "Ntrk2", "Ntrk3", "Gfra1", "Gfra2", "Gfra3", "Ret",
         "Piezo2", "Kcnv1", "Nefh", "Nefl", "Prph", "Scn10a", "Scn11a", "Trpv1", "Trpm8",
         "Pvalb", "Pou4f1", "Pou4f2", "Pou4f3", "Slc17a8", "Casq2", "Smr2", "Prokr2",
         "Baiap2l1", "Wnt7a", "Lpar3", "Il31ra", "Mrgprx1", "P2rx3", "Oprm1", "Oprd1",
         "Oprk1", "Adora1", "Cnr1", "Sox10", "Mpz", "Fabp7", "Snap25"]

for key in ("leiden_0.8", "leiden_1.0", "leiden_1.2"):
    cl = a.obs[key].astype(str).values
    cats = sorted(pd.unique(cl), key=int)
    have = [g for g in PANEL if g in a.var_names]
    D = pd.DataFrame(index=cats)
    for g in have:
        D[g] = pd.Series(a[:, g].X.ravel()).groupby(cl).mean()
    print("\n===== 诊断基因均值 (%s) =====" % key, flush=True)
    pd.set_option("display.width", 400, "display.max_columns", 100)
    print(D.round(2).to_string(), flush=True)
    D.round(3).to_csv("mouse_ref_diag_%s.csv" % key, encoding="utf-8-sig")

    # 标志基因
    b = a.copy()
    b.obs["cl"] = cl
    sc.tl.rank_genes_groups(b, "cl", method="wilcoxon", n_genes=10, use_raw=False)
    r = pd.DataFrame(b.uns["rank_genes_groups"]["names"])
    sz = pd.Series(cl).value_counts()
    print("\n-- %s top10 --" % key, flush=True)
    for c in cats:
        print("  簇%-2s n=%-4d %s" % (c, sz[c], ", ".join(r[c][:10])), flush=True)

a.write_h5ad("mdrg_reclustered_59739.h5ad")
print("\n已存 mdrg_reclustered_59739.h5ad", flush=True)
