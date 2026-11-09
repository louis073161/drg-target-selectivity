# -*- coding: utf-8 -*-
"""小鼠 DRG：与人类侧同流程的聚类与标志基因反查（step1）"""
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy import sparse

DS = sys.argv[1] if len(sys.argv) > 1 else "155622"
np.random.seed(0)

PANEL = ["Calca", "Tac1", "Nmb", "Sst", "Mrgprd", "Mrgpra3", "Mrgprb4", "Nppb", "Zcchc12",
         "Th", "Fam19a4", "Ntrk1", "Ntrk2", "Ntrk3", "Gfra1", "Gfra2", "Gfra3", "Ret",
         "Piezo2", "Kcnv1", "Nefh", "Nefl", "Prph", "Scn10a", "Scn11a", "Trpv1", "Trpm8",
         "Pvalb", "Pou4f1", "Pou4f2", "Pou4f3", "Slc17a8", "Casq2", "Smr2", "Prokr2",
         "Baiap2l1", "Cldn9", "Wnt7a", "Lpar3", "Il31ra", "S100b", "Fabp7", "Gfap",
         "Cd74", "C1qc", "Aif1", "Ptprc", "Sox10", "Mpz", "P2rx3", "Oprm1", "Oprd1",
         "Oprk1", "Adora1", "Cnr1", "Cnr2", "Gpr55"]

if DS == "155622":
    df = pd.read_csv("GSE155622_Smart-seq2_FPKM.txt.gz", sep="\t", index_col=0)
    meta = pd.read_csv("GSE155622_Smart-seq2_FPKM_metadata.txt.gz", sep="\t", index_col=0)
    X = df.T
    X = X.loc[[c for c in meta.index if c in X.index]]
    meta = meta.loc[X.index]
    obs = pd.DataFrame({"condition": meta.Condition.values,
                        "author_cluster": meta.clusters.values,
                        "ib4": meta.IB4.values}, index=X.index)
else:
    raw = pd.read_csv("GSE59739_DataTable.txt.gz", sep="\t", index_col=0).T
    hdr = raw.columns[:4].tolist()           # 元数据列
    obs = raw[hdr].copy()
    X = raw.drop(columns=hdr)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0)
    obs.columns = ["sess", "temp", "sex", "content"]

X = X.astype(np.float32)
print("表达矩阵(细胞×基因):", X.shape, flush=True)
print("元数据列:", obs.columns.tolist(), flush=True)

adata = ad.AnnData(sparse.csr_matrix(X.values), obs=obs, var=pd.DataFrame(index=X.columns))
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.layers["fpkm"] = adata.X.copy()

# 过滤：检出基因数过低或为空的细胞
adata.var["det"] = np.asarray((adata.X > 0).mean(0)).ravel()
adata = adata[:, adata.var.det >= 0.01].copy()
sc.pp.filter_cells(adata, min_genes=500)
print("过滤后:", adata.shape, flush=True)
if "content" in adata.obs:
    print("Content 分布:", adata.obs.content.value_counts().to_dict(), flush=True)

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars - 1))
a2 = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(a2, max_value=10)
npcs = min(30, a2.n_obs - 1, a2.n_vars - 1)
sc.tl.pca(a2, n_comps=npcs, svd_solver="arpack")
sc.pp.neighbors(a2, n_neighbors=12, n_pcs=npcs)
sc.tl.leiden(a2, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(a2)
adata.obs["leiden"] = a2.obs.leiden.values
adata.obsm["X_umap"] = a2.obsm["X_umap"]
print("簇数:", adata.obs.leiden.nunique(), "| 各簇细胞数:", adata.obs.leiden.value_counts().sort_index().to_dict(), flush=True)

# 各簇标志基因
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=12, use_raw=False)
r = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])
sz = adata.obs.leiden.value_counts()
print("\n=== 各簇 top12 标志基因 ===", flush=True)
for c in sorted(r.columns, key=lambda x: int(x)):
    print("簇%2s (n=%3d): %s" % (c, sz[c], ", ".join(r[c][:12])), flush=True)

# 诊断 panel 均值
have = [g for g in PANEL if g in adata.var_names]
D = pd.DataFrame(index=sorted(adata.obs.leiden.unique(), key=int))
A = adata.X
A = np.asarray(A.todense()) if sparse.issparse(A) else np.asarray(A, dtype=np.float32)
for g in have:
    j = adata.var_names.get_loc(g)
    D[g] = pd.Series(A[:, j]).groupby(adata.obs.leiden.values).mean()
pd.set_option("display.width", 300, "display.max_columns", 80)
print("\n=== 诊断基因均值 ===", flush=True)
print(D.round(2).to_string(), flush=True)
D.round(3).to_csv("mdrg_diag_panel_%s.csv" % DS, encoding="utf-8-sig")
adata.write_h5ad("mdrg_clustered_%s.h5ad" % DS)
print("\n已保存 mdrg_clustered_%s.h5ad" % DS, flush=True)
