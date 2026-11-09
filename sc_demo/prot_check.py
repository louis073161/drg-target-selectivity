# -*- coding: utf-8 -*-
"""表 10 与 3.13 节的蛋白层核对：从公开人 DRG 蛋白组复现全部数字。

输入  ext/prot/ST11-ExpressionMatrix.csv
      PRIDE PXD061129（Barry et al., Pain 2025;166:2386-2402）的定量矩阵，
      DIA-PASEF，timsTOF HT，8 名供体的 ganglia 与 nerve root 成对样本。

输入  data/01_人类发现队列_GSE249746/genomewide_delta_auc.csv
      本文主分布（16 969 个基因 × 13 个神经元类型）。

输出  table10_proteome.csv      表 10 每一行的检出供体数、CV、蛋白丰度百分位
      adora1_by_class.csv        3.13 节 ADORA1 逐类检出率与类间相关

用法  python prot_check.py    （工作目录为 sc_demo/）
"""
import numpy as np
import pandas as pd
from scipy import stats

PROT = "../ext/prot/ST11-ExpressionMatrix.csv"
GW = "../data/01_人类发现队列_GSE249746/genomewide_delta_auc.csv"
BYCLASS = "../data/01_人类发现队列_GSE249746/atlas_target_by_class.csv"

# 表 3 的 12 个靶点，顺序与正文一致
TARGETS = ["TRPM8", "SCN11A", "SST", "TRPV1", "CALCA", "P2RX3",
           "SCN10A", "OPRM1", "OPRD1", "CNR1", "ADORA1", "OPRK1"]

# ---------- 1. 蛋白层 ----------
d = pd.read_csv(PROT)
gang = [c for c in d.columns if c.endswith(".ganglia")]
nerve = [c for c in d.columns if c.endswith(".nerve.root")]
d = d.dropna(subset=["genes"])
p = d.groupby("genes")[gang].mean()          # 同一基因多蛋白条目取均值
prot = p.mean(axis=1)                        # 8 名供体的平均 log2 强度

rows = []
for g in TARGETS:
    if g in p.index:
        v = p.loc[g]
        n = int(v.notna().sum())
        cv = float(v.std() / v.mean() * 100) if n > 1 else np.nan
        ab = float(prot[g])
        pct = float((prot.dropna() < ab).mean() * 100)
    else:
        n, cv, ab, pct = 0, np.nan, np.nan, np.nan
    rows.append([g, n, len(gang), round(cv, 1) if n > 1 else "—",
                 round(ab, 2) if n else "—", round(pct, 1) if n else "—"])
tab = pd.DataFrame(rows, columns=["gene", "n_donor_detected", "n_donor_total",
                                  "CV_pct", "log2_abundance", "abundance_pctile"])
tab.to_csv("table10_proteome.csv", index=False)

print("== 表 10 蛋白层 ==")
print(tab.to_string(index=False))
print("全蛋白组检出基因数 %d；8 供体 log2 强度的 CV 中位 %.1f%%"
      % (len(prot.dropna()), 100 * (p.std(axis=1) / p.mean(axis=1)).median()))

# 正文引用但未列入表 10 者（SCN9A）
extra = []
for gname in ["SCN9A"]:
    if gname in p.index:
        v = p.loc[gname]
        extra.append([gname, int(v.notna().sum()), len(gang),
                      round(float(v.std() / v.mean() * 100), 1)])
print("\n== 正文另引 ==")
print(pd.DataFrame(extra, columns=["gene", "detected", "total", "CV_pct"]).to_string(index=False))

# 全基因组尺度：蛋白丰度与 maxΔAUC 的关系
g = pd.read_csv(GW, index_col=0)
mx = g.max(axis=1)
common = [x for x in prot.index if x in mx.index and pd.notna(prot[x])]
a, b = prot[common].values, mx[common].values
ok = ~(np.isnan(a) | np.isnan(b))
r = stats.spearmanr(a[ok], b[ok])
print("蛋白丰度 ~ maxΔAUC：n=%d，Spearman ρ=%.3f，P=%.1e" % (ok.sum(), r.statistic, r.pvalue))

# ---------- 2. 原位层：ADORA1 逐类检出率 ----------
c = pd.read_csv(BYCLASS)
c = c[c.cell_class != "非神经元"]             # 与正文一致，只留 9 个神经元类
piv = c.pivot_table(index="cell_class", columns="gene", values="pct_expr")
out = piv[["ADORA1", "TRPV1"]].sort_values("ADORA1", ascending=False)
out["ADORA1"] = out["ADORA1"].round(1)
out["TRPV1"] = out["TRPV1"].round(1)
out.to_csv("adora1_by_class.csv")

print("\n== ADORA1 逐类检出率（9 个神经元类）==")
print(out.to_string())
v = out["ADORA1"]
print("极差 %.1f 个百分点（%.1f ~ %.1f）" % (v.max() - v.min(), v.min(), v.max()))
rr = stats.spearmanr(out["ADORA1"], out["TRPV1"])
print("ADORA1 与 TRPV1 的类间相关：ρ=%+.2f，P=%.2f（n=9，只作方向性对照）"
      % (rr.statistic, rr.pvalue))
