# -*- coding: utf-8 -*-
"""重建跨物种对照表：人类取自 atlas（与图谱一致），小鼠取自原始 FPKM 矩阵
（避免 h5ad 的基因过滤删掉 Oprd1/Oprk1/Cnr2/Gpr55）。
检出率与 ΔAUC 必须来自同一矩阵。
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import pearsonr
import scanpy as sc
from scipy import sparse

CLS = ["NP", "PEP", "C-LTMR", "A-LTMR", "Aδ", "Propr"]

# ---------- 人类 ----------
h = pd.read_csv("cross_species_delta_auc.csv")
h = h[h.species == "human"].copy()
H = h.pivot_table(index="gene", columns="cls", values="dAUC")
print("人类：", H.shape, "| 类型", [c for c in H.columns])

# 人类检出率：在图谱神经元上重算
a = sc.read_h5ad("../sc_demo/atlas_s3_final.h5ad")
a = a[a.obs.cell_class != "非神经元"].copy()
X = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X)
sym = {str(v).upper(): i for i, v in enumerate(a.var_names)}
hdet = {g: float((X[:, sym[g]] > 0).mean() * 100) for g in H.index if g in sym}

# ---------- 小鼠 ----------
m = pd.read_csv("mouse_delta_auc_raw.csv")
M = m.pivot_table(index="gene", columns="cls", values="dAUC")
md = m.groupby("gene").pct.max()
print("小鼠：", M.shape, "| 类型", [c for c in M.columns])

rows = []
for g in sorted(set(H.index) & set(M.index)):
    shared = [c for c in CLS if c in H.columns and c in M.columns]
    hv = H.loc[g, shared].astype(float).values
    mv = M.loc[g, shared].astype(float).values
    ok = ~(np.isnan(hv) | np.isnan(mv))
    r = pearsonr(hv[ok], mv[ok])[0] if ok.sum() >= 3 else np.nan
    agree = int((np.sign(hv[ok]) == np.sign(mv[ok])).sum())
    hd, md_ = hdet.get(g, np.nan), float(md.get(g, np.nan))
    lo = np.nanmin([hd, md_])
    grade = "可定量比较" if lo >= 20 else ("仅比方向" if lo >= 5 else "不可评价（检出<5%）")
    rows.append({
        "gene": g, "human_detect_pct": round(hd, 1), "mouse_detect_pct": round(md_, 1),
        **{"H_" + c: (round(float(H.loc[g, c]), 3) if (c in H.columns and not pd.isna(H.loc[g, c])) else np.nan) for c in CLS},
        **{"M_" + c: (round(float(M.loc[g, c]), 3) if (c in M.columns and not pd.isna(M.loc[g, c])) else np.nan) for c in CLS},
        "H_top": H.loc[g, shared].astype(float).idxmax(),
        "H_topd": round(float(np.nanmax(hv)), 3),
        "M_top": M.loc[g, shared].astype(float).idxmax(),
        "M_topd": round(float(np.nanmax(mv)), 3),
        "sign_concord": "%d/%d" % (agree, int(ok.sum())),
        "pearson_r": round(float(r), 2) if not pd.isna(r) else np.nan,
        "grade": grade,
    })
D = pd.DataFrame(rows)
D.to_csv("cross_species_final_v2.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.width", 260)
print()
print(D[["gene", "human_detect_pct", "mouse_detect_pct", "H_top", "H_topd",
         "M_top", "M_topd", "sign_concord", "pearson_r", "grade"]].to_string(index=False))

allr = [x for x in D.pearson_r if not pd.isna(x)]
print("\n可计算 r 的基因数 %d，中位 r = %.2f" % (len(allr), np.median(allr)))
print("已写出 cross_species_final_v2.csv")
