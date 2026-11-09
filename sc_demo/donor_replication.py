# -*- coding: utf-8 -*-
"""供体内部复算：在每名供体内独立重算 ΔAUC 与全基因组百分位
用于回答"这篇的结论会不会只是某一两个供体带来的"
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
from scipy.stats import rankdata, spearmanr

A = sc.read_h5ad("atlas_s3_final.h5ad")
A = A[A.obs.cell_class != "非神经元"].copy()
X = A.X.toarray() if sparse.issparse(A.X) else np.asarray(A.X, dtype=np.float32)
X = X.astype(np.float64)
genes = np.array(A.var_names)
det_all = (X > 0).mean(0)
keep = np.where(det_all >= 0.05)[0]
print("神经元 %d，基因 %d" % (A.n_obs, len(keep)))

TARGETS = ["SCN10A", "SCN11A", "P2RX3", "TRPV1", "OPRM1", "OPRD1",
           "CNR1", "ADORA1", "OPRK1", "CALCA"]
gidx = {g: i for i, g in enumerate(genes)}
donors = sorted(A.obs.donor.astype(str).unique())
print("供体:", {d: int((A.obs.donor.astype(str) == d).sum()) for d in donors})


def run(mask_cells, label):
    """在给定细胞子集上计算全基因组 maxΔAUC 与靶点百分位"""
    Xs = X[mask_cells][:, keep]
    ct = A.obs.cell_type.astype(str).values[mask_cells]
    types = [t for t in pd.unique(ct) if (ct == t).sum() >= 8]   # 至少 8 个细胞
    if len(types) < 5:
        return None
    Xt = Xs.T
    R = rankdata(Xt, axis=1)
    n_all = Xt.shape[1]
    AUC = np.zeros((len(keep), len(types)))
    for j, t in enumerate(types):
        m1 = (ct == t); n1 = int(m1.sum()); n2 = n_all - n1
        AUC[:, j] = (R[:, m1].sum(1) - n1 * (n1 + 1) / 2) / (n1 * n2)
    D = AUC - np.median(AUC, axis=0, keepdims=True)
    mx = D.max(1)
    rows = []
    for g in TARGETS:
        if g not in gidx:
            continue
        p = gidx[g] - np.searchsorted(keep, gidx[g]) if False else None
        pos = np.where(keep == gidx[g])[0]
        if len(pos) == 0:
            continue
        k = pos[0]
        rows.append({"set": label, "gene": g, "max_dAUC": round(float(mx[k]), 3),
                     "pct": round(float((mx < mx[k]).mean() * 100), 1),
                     "top_type": types[int(D[k].argmax())]})
    return pd.DataFrame(rows)


out = [run(np.ones(A.n_obs, bool), "全部")]
for d in donors:
    out.append(run((A.obs.donor.astype(str) == d).values, d))
res = pd.concat([r for r in out if r is not None], ignore_index=True)
res.to_csv("donor_replication.csv", index=False, encoding="utf-8-sig")

pd.set_option("display.width", 200)
piv = res.pivot_table(index="gene", columns="set", values="pct")
piv = piv.reindex([g for g in TARGETS if g in piv.index])
print("\n=== 各靶点在全基因组分布中的百分位（按供体分别计算）===")
print(piv.to_string())

print("\n=== 最高富集类型 ===")
print(res.pivot_table(index="gene", columns="set", values="top_type", aggfunc="first")
      .reindex([g for g in TARGETS if g in piv.index]).to_string())

# 供体间一致性：靶点百分位向量的 Spearman 相关
dl = [c for c in piv.columns if c != "全部"]
print("\n=== 供体间百分位向量的一致性（Spearman）===")
for i in range(len(dl)):
    for j in range(i + 1, len(dl)):
        r = spearmanr(piv[dl[i]], piv[dl[j]]).correlation
        print("  %s vs %s : rho = %.2f" % (dl[i], dl[j], r))

print("\n=== 关键判断 ===")
for g in ["ADORA1", "SCN10A", "SCN11A", "TRPV1"]:
    if g in piv.index:
        v = piv.loc[g, dl].astype(float)
        print("  %-7s 各供体百分位 %s | 极差 %.1f" % (g, list(v.round(1)), v.max() - v.min()))
print("\n已写出 donor_replication.csv")
