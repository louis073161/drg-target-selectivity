# -*- coding: utf-8 -*-
"""小鼠侧参照系 · 第二步（对齐表 4 口径）：全基因组 maxΔAUC 经验分布。

口径对齐
--------
表 4 的 ΔAUC 是这样算的（cross_compare.py + rebuild_cross.py）：
  · 细胞：全部（小鼠 433 个，含 TH-C；人类 1 074 个，含 Aβ-LTMR 与损伤应激）
  · 零分布：每个类内、全部基因的 AUC 中位数
  · 取最大值时只在这 5 个两类共有的类上取：NP / PEP / C-LTMR / Aδ / Propr
本脚本沿用同一口径：ΔAUC 在全部细胞、全部类上算，max 只在这 5 类上取。
这样表 4 中的每个数值都能在自己物种的分布里定位，百分位才可比。
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
from scipy.stats import rankdata

TARG = ["SCN10A", "SCN11A", "P2RX3", "CALCA", "OPRM1", "TRPV1", "ADORA1", "CNR1",
        "MRGPRX1", "OPRD1", "OPRK1", "TRPM8"]
SHARED = ["NP", "PEP", "C-LTMR", "Aδ", "Propr"]

H2C = {"hNP1": "NP", "hPEP.TAC1+": "PEP", "hPEP.NMB+": "PEP", "hPEP.KIT": "PEP",
       "hPEP.CHRNA7": "PEP", "hSST": "PEP", "hC.LTMR": "C-LTMR", "hAδ.LTMR": "Aδ",
       "hAβ.LTMR": "A-LTMR", "hAβ.LTMR.SPP1+": "A-LTMR", "hPropr": "Propr",
       "hATF3": "损伤应激", "hTRPM8": "C-LTMR"}
M2C = {"0": "TH-C", "1": "Propr", "2": "PEP", "3": "NP",
       "4": "Aδ", "5": "PEP", "6": "C-LTMR", "7": "C-LTMR"}


def matrix(path, cls_map, src, drop_non_neuron):
    a = sc.read_h5ad(path)
    if drop_non_neuron:
        a = a[a.obs["cell_class"] != "非神经元"].copy()
    a.obs["cls"] = a.obs[src].astype(str).map(cls_map)
    a = a[a.obs.cls.notna()].copy()
    X = a.X.toarray() if sparse.issparse(a.X) else np.asarray(a.X, dtype=np.float32)
    sym = np.array([str(v).upper() for v in a.var_names])
    _, first = np.unique(sym, return_index=True)
    return a, X[:, np.sort(first)].astype(np.float64), sym[np.sort(first)]


def distribution(X, sym, cls, det_min=0.05):
    det = (X > 0).mean(0)
    keep = det >= det_min
    Xk, G = X[:, keep], sym[keep]
    R = rankdata(Xk.T, method="average", axis=1)
    n = Xk.shape[0]
    cats = sorted(pd.unique(cls))
    dauc = np.empty((len(G), len(cats)))
    for j, c in enumerate(cats):
        sub = cls == c
        n1 = int(sub.sum()); n2 = n - n1
        auc = (R[:, sub].sum(1) - n1 * (n1 + 1) / 2) / (n1 * n2)
        dauc[:, j] = auc - np.median(auc)
    D = pd.DataFrame(dauc, index=G, columns=cats)
    mx = D[[c for c in SHARED if c in D.columns]].max(1)
    st = dict(n_gene=int(len(G)), n_cell=int(n), n_class=len(cats),
              P50=float(np.percentile(mx, 50)), P75=float(np.percentile(mx, 75)),
              P90=float(np.percentile(mx, 90)), P95=float(np.percentile(mx, 95)),
              P99=float(np.percentile(mx, 99)))
    return mx, st, D, pd.Series(det[keep], index=G)


def table4_percentile(mx, D, detk, label):
    rows = []
    for t in TARG:
        if t not in mx.index:
            rows.append(dict(数据集=label, 基因=t, 最高类="", maxΔAUC=np.nan,
                             百分位=np.nan, 检出率=np.nan))
            continue
        v = float(mx[t])
        best = D.loc[t, [c for c in SHARED if c in D.columns]].astype(float).idxmax()
        rows.append(dict(数据集=label, 基因=t, 最高类=best, maxΔAUC=round(v, 3),
                         百分位=round(float((mx < v).mean() * 100), 1),
                         检出率=round(float(detk[t] * 100), 1)))
    return pd.DataFrame(rows)


out = []
md = sc.read_h5ad("mdrg_reclustered_59739.h5ad")

for res in (0.8, 1.0, 1.2):
    key = "leiden_%.1f" % res
    a, X, sym = matrix("mdrg_reclustered_59739.h5ad", M2C, key, False)
    cls = a.obs.cls.values
    mx, st, D, detk = distribution(X, sym, cls)
    print("\n===== 小鼠 res=%.1f | %d 细胞 %d 类 %d 基因 =====" %
          (res, st["n_cell"], st["n_class"], st["n_gene"]))
    print("   P50=%.3f P75=%.3f P90=%.3f P95=%.3f P99=%.3f" %
          (st["P50"], st["P75"], st["P90"], st["P95"], st["P99"]))
    T = table4_percentile(mx, D, detk, "小鼠 res=%.1f" % res)
    print(T.to_string(index=False))
    T.to_csv("mouse_ref_targets_res%.1f.csv" % res, index=False, encoding="utf-8-sig")
    pd.DataFrame({"max_dAUC": mx}).round(5).to_csv(
        "mouse_ref_genomewide_res%.1f.csv" % res, encoding="utf-8-sig")
    out.append(dict(物种="小鼠", 方案="res=%.1f" % res,
                    **{k: v for k, v in st.items() if k.startswith("P") or k == "n_gene"}))

# ---------- 人类同口径 ----------
ah, hx, hs = matrix("../sc_demo/atlas_s3_final.h5ad", H2C, "cell_type", True)
hcls = ah.obs.cls.values
mx, st, D, detk = distribution(hx, hs, hcls)
print("\n===== 人类（同口径）| %d 细胞 %d 类 %d 基因 =====" % (st["n_cell"], st["n_class"], st["n_gene"]))
print("   P50=%.3f P75=%.3f P90=%.3f P95=%.3f P99=%.3f" %
      (st["P50"], st["P75"], st["P90"], st["P95"], st["P99"]))
Th = table4_percentile(mx, D, detk, "人类")
print(Th.to_string(index=False))
Th.to_csv("human_collapsed_targets.csv", index=False, encoding="utf-8-sig")
pd.DataFrame({"max_dAUC": mx}).round(5).to_csv("human_collapsed_genomewide.csv", encoding="utf-8-sig")
out.append(dict(物种="人类", 方案="合并大类",
                **{k: v for k, v in st.items() if k.startswith("P") or k == "n_gene"}))

S = pd.DataFrame(out)
print("\n===== 汇总 =====")
print(S.round(3).to_string(index=False))
S.to_csv("mouse_ref_distribution_summary.csv", index=False, encoding="utf-8-sig")
