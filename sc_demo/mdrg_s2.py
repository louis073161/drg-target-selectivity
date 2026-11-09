# -*- coding: utf-8 -*-
"""小鼠 DRG：神经元纯化后重聚类（step2）"""
import warnings, sys
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse

a = sc.read_h5ad("mdrg_clustered_59739.h5ad")
X = a.X
X = np.asarray(X.todense()) if sparse.issparse(X) else np.asarray(X, dtype=np.float32)
a.X = X

def gmean(g):
    return float(a[:, g].X.mean()) if g in a.var_names else 0.0

markers = {"Snap25": gmean("Snap25"), "Ptprc": gmean("Ptprc"), "Col1a1": gmean("Col1a1"),
           "Acta2": gmean("Acta2"), "Hbb-b2": gmean("Hbb-b2"), "Fabp7": gmean("Fabp7"),
           "Mpz": gmean("Mpz"), "Sox10": gmean("Sox10")}
print("全库均值:", {k: round(v, 2) for k, v in markers.items()}, flush=True)

# 神经元判定：神经标志高 且 非免疫/非间质/非髓鞘/非红系
neu = np.ones(a.n_obs, bool)
for g, thr in (("Ptprc", 1.0), ("Col1a1", 2.0), ("Acta2", 2.0), ("Hbb-b2", 2.0), ("Mpz", 2.0)):
    if g in a.var_names:
        neu &= ~(a[:, g].X.ravel() > thr)
if "Snap25" in a.var_names:
    neu &= a[:, "Snap25"].X.ravel() > 0.5
print("神经元 %d / %d" % (neu.sum(), a.n_obs), flush=True)
a = a[neu].copy()

sc.pp.filter_cells(a, min_genes=1500)
sc.pp.filter_genes(a, min_cells=10)
print("质控后:", a.shape, flush=True)

sc.pp.highly_variable_genes(a, n_top_genes=2000)
b = a[:, a.var.highly_variable].copy()
sc.pp.scale(b, max_value=10)
npcs = min(30, b.n_obs - 1, b.n_vars - 1)
sc.tl.pca(b, n_comps=npcs, svd_solver="arpack")
sc.pp.neighbors(b, n_neighbors=15, n_pcs=npcs)
for res in (1.2,):
    sc.tl.leiden(b, resolution=res, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(b)
a.obs["leiden"] = b.obs.leiden.values
a.obsm["X_umap"] = b.obsm["X_umap"]
print("簇数:", a.obs.leiden.nunique(), a.obs.leiden.value_counts().sort_index().to_dict(), flush=True)

sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon", n_genes=10, use_raw=False)
r = pd.DataFrame(a.uns["rank_genes_groups"]["names"])
sz = a.obs.leiden.value_counts()
print("\n=== 各簇 top10 ===", flush=True)
for c in sorted(r.columns, key=lambda x: int(x)):
    print("簇%2s (n=%3d): %s" % (c, sz[c], ", ".join(r[c][:10])), flush=True)

PANEL = ["Calca", "Tac1", "Nmb", "Sst", "Mrgprd", "Mrgpra3", "Mrgprb4", "Nppb", "Zcchc12",
         "Th", "Fam19a4", "Ntrk1", "Ntrk2", "Ntrk3", "Gfra1", "Gfra2", "Gfra3", "Ret",
         "Piezo2", "Kcnv1", "Nefh", "Nefl", "Scn10a", "Scn11a", "Trpv1", "Trpm8",
         "Pvalb", "Pou4f1", "Pou4f2", "Pou4f3", "Slc17a8", "Casq2", "Smr2", "Prokr2",
         "Baiap2l1", "Wnt7a", "Lpar3", "Il31ra", "Mrgprx1", "P2rx3", "Oprm1", "Oprd1",
         "Oprk1", "Adora1", "Cnr1", "Cnr2", "Gpr55", "Cd74"]
have = [g for g in PANEL if g in a.var_names]
D = pd.DataFrame(index=sorted(a.obs.leiden.unique(), key=int))
A = a.X
for g in have:
    D[g] = pd.Series(A[:, a.var_names.get_loc(g)]).groupby(a.obs.leiden.values).mean()
pd.set_option("display.width", 320, "display.max_columns", 90)
print("\n=== 诊断基因 ===", flush=True)
print(D.round(2).to_string(), flush=True)
D.round(3).to_csv("mdrg_diag_neuron.csv", encoding="utf-8-sig")
a.write_h5ad("mdrg_neuron_59739.h5ad")
print("\n保存 mdrg_neuron_59739.h5ad", flush=True)
