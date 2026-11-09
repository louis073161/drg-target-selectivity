# -*- coding: utf-8 -*-
"""小鼠侧 ΔAUC：与人类侧 atlas_s3_stats.py 逐行同口径"""
import sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import rankdata, mannwhitneyu
from scipy import sparse

H5 = sys.argv[1]
MAPFILE = sys.argv[2]          # csv: leiden,cell_type,cell_class
OUT = sys.argv[3]

TARGETS = ["OPRM1", "OPRD1", "OPRK1", "ADORA1", "ADORA2A", "CNR1", "CNR2", "GPR55",
           "SCN10A", "SCN11A", "TRPV1", "TRPM8", "CALCA", "NMB", "MRGPRX1", "P2RX3"]

a = sc.read_h5ad(H5)
ann = pd.read_csv(MAPFILE)
mp = dict(zip(ann.leiden.astype(str), ann.cell_type))
cl = dict(zip(ann.leiden.astype(str), ann.cell_class))
a.obs["leiden"] = a.obs["leiden"].astype(str)
a.obs["cell_type"] = a.obs.leiden.map(mp).astype("category")
a.obs["cell_class"] = a.obs.leiden.map(cl).astype("category")
a = a[a.obs.cell_type.notna()].copy()
print("细胞数", a.n_obs, "| 类型数", a.obs.cell_type.nunique(), flush=True)
print(a.obs.cell_type.value_counts().to_dict(), flush=True)

# 统一为人类基因符号（小鼠仅首字母大写）
a.var["sym"] = [v.upper() for v in a.var_names]
sym2idx = {}
for i, s in enumerate(a.var.sym):
    sym2idx.setdefault(s, i)

Xall = a.X
Xall = np.asarray(Xall.todense()) if sparse.issparse(Xall) else np.asarray(Xall, dtype=np.float32)
det = (Xall > 0).mean(0)
keep_idx = np.where(det >= 0.05)[0]
pan_idx = [sym2idx[g] for g in TARGETS if g in sym2idx]
keep_idx = np.union1d(keep_idx, np.array(pan_idx, dtype=int))
E = Xall[:, keep_idx].astype(np.float32)
EG = [a.var.sym[i] for i in keep_idx]
print("[秩变换] 基因 %d / %d" % (len(EG), a.n_vars), flush=True)

R = rankdata(E, axis=0)
ct = a.obs.cell_type.values
cats = list(a.obs.cell_type.cat.categories)
n_all = len(ct)


def auc_matrix(Rm, mask):
    n1 = mask.sum(); n0 = n_all - n1
    return (Rm[mask].mean(0) - (n1 + 1) / 2) / n0


auc_all, base = {}, {}
for c in cats:
    m = ct == c
    auc_all[c] = auc_matrix(R, m)
    base[c] = np.median(auc_all[c])
print("[零分布] AUC 中位数:", {k: round(v, 3) for k, v in base.items()}, flush=True)

gidx = {g: i for i, g in enumerate(EG)}
rows = []
for g in TARGETS:
    if g not in gidx:
        print("  未检出基因:", g, flush=True); continue
    j = gidx[g]
    v = E[:, j]; glob = v.mean()
    pct_rank = (E.mean(0) < glob).mean() * 100
    for c in cats:
        m = ct == c
        if m.sum() < 5:
            continue
        rows.append(dict(gene=g, cell_type=c, n_cells=int(m.sum()), mean=v[m].mean(),
                         pct_expr=float((v[m] > 0).mean() * 100), global_mean=glob,
                         auroc=auc_all[c][j], delta_auc=auc_all[c][j] - base[c],
                         p=mannwhitneyu(v[m], v[~m], alternative="greater").pvalue,
                         abun_pct=pct_rank))
per = pd.DataFrame(rows)
per["q"] = np.nan
for g, idx in per.groupby("gene").groups.items():
    p = per.loc[idx, "p"].values; n = len(p); order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    per.loc[idx, "q"] = np.clip(q, 0, 1)

per.to_csv(OUT, index=False, encoding="utf-8-sig")
print("\n=== 最强富集类型 ===", flush=True)
for g in sorted(per.gene.unique()):
    r = per[per.gene == g].sort_values("delta_auc", ascending=False)
    x = r.iloc[0]
    print("  %-8s %-14s ΔAUC=%.3f 均值=%.2f 丰度P%.0f q=%.2g | 前二: %s"
          % (g, x.cell_type, x.delta_auc, x["mean"], x.abun_pct, x.q,
             "/".join(r.cell_type.head(2))), flush=True)
a.write_h5ad(H5.replace(".h5ad", "_named.h5ad"))
