# -*- coding: utf-8 -*-
"""检出率（及基因长度）匹配的条件零分布：
把每个靶点放回"同检出率、同基因长度"的基因中定位，
用以排除"低表达基因的 ΔAUC 被系统性压缩"这一混杂。"""
import numpy as np, pandas as pd
from scipy import stats

r = pd.read_csv("gene_selectivity_ranked.csv", index_col=0).dropna(
    subset=["detect_pct", "max_dAUC"])

# 基因长度（GRCh38，end - start）
g = pd.read_csv("../ext/gwas/gene_coords_grch38.tsv", sep="\t")
g["len"] = (g.end - g.start).clip(lower=1)
gl = g.drop_duplicates("gene").set_index("gene")["len"]
r["len"] = gl.reindex(r.index)

d = r.detect_pct.values
a = r.max_dAUC.values
genes = np.array(r.index)

TARGETS = ["TRPM8", "SCN11A", "SST", "TRPV1", "CALCA", "P2RX3", "SCN10A",
           "OPRM1", "OPRD1", "CNR1", "ADORA1", "OPRK1", "SCN9A", "MRGPRX1"]

K = 800
rows = []
for name in TARGETS:
    if name not in set(genes):
        continue
    i = int(np.where(genes == name)[0][0])
    dd, aa = d[i], a[i]

    # --- 1D：检出率最近 800 个 ---
    near = [j for j in np.argsort(np.abs(d - dd)) if j != i][:K]
    pct_1d = float((a[near] < aa).mean() * 100)

    # --- 2D：检出率 + log 基因长度 ---
    ok = r["len"].notna().values
    idx = np.where(ok)[0]
    z1 = (d[ok] - d[ok].mean()) / d[ok].std()
    z2 = (np.log10(r["len"].values[ok]) - np.log10(r["len"].values[ok]).mean()) / \
         np.log10(r["len"].values[ok]).std()
    dist = np.hypot(z1 - z1[list(idx).index(i)], z2 - z2[list(idx).index(i)])
    order = [idx[k] for k in np.argsort(dist) if idx[k] != i][:K]
    pct_2d = float((a[order] < aa).mean() * 100)

    m = np.abs(d - dd) <= 5.0
    m[i] = False
    pct_win = float((a[m] < aa).mean() * 100) if m.sum() >= 50 else np.nan

    rows.append([name, round(aa, 3), round(dd, 1),
                 round(float((a < aa).mean() * 100), 1),
                 round(pct_1d, 1), round(pct_2d, 1), int(m.sum()), round(pct_win, 1)])

out = pd.DataFrame(rows, columns=["基因", "maxΔAUC", "检出率%", "全基因组百分位",
                                  "匹配百分位", "匹配百分位_含长度", "±5pp内基因数", "窗口百分位"])
out.to_csv("matched_null.csv", index=False, encoding="utf-8-sig")
print(out.to_string(index=False))

print("\n检出率 ~ maxΔAUC  Spearman ρ = %.3f" % stats.spearmanr(d, a).statistic)
print("表达量 ~ maxΔAUC  Spearman ρ = %.3f"
      % stats.spearmanr(r.mean_expr.values, a).statistic)
print("两种匹配口径的秩相关 ρ = %.3f"
      % stats.spearmanr(out["匹配百分位"], out["匹配百分位_含长度"]).statistic)
print("匹配 vs 全基因组 秩相关 ρ = %.3f"
      % stats.spearmanr(out["全基因组百分位"], out["匹配百分位_含长度"]).statistic)

# ---- TOST 等价界敏感性 ----
eq = pd.read_csv("equivalence_tests.csv")


def tost(delta, shift, se):
    return max(stats.norm.sf((delta - shift) / se),
               stats.norm.cdf((-delta - shift) / se))


print("\nTOST 在不同等价界下的 P：")
for delta in (0.05, 0.066, 0.08, 0.10):
    row = eq[eq["基因"] == "ADORA1"].iloc[0]
    k = eq[eq["基因"] == "OPRK1"].iloc[0]
    p1 = tost(delta, row["偏移"], row["SE"])
    pk = tost(delta, k["偏移"], k["SE"])
    print("  δ=%.3f  A1R P=%.4f %s   κ P=%.4f %s"
          % (delta, p1, "等效" if p1 < 0.05 else "—", pk, "等效" if pk < 0.05 else "—"))
