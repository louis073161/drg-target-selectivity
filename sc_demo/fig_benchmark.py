# -*- coding: utf-8 -*-
"""图2：全基因组经验基准与靶点定位（A 分布 / B 定位 / C 检出率匹配）"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib
from matplotlib import font_manager
for _f in ["C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/simhei.ttf"]:
    try:
        font_manager.fontManager.addfont(_f); break
    except Exception:
        pass
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

s = pd.read_csv("gene_selectivity_ranked.csv", index_col=0)
v = s.max_dAUC.values
P = {q: float(np.percentile(v, q)) for q in (50, 75, 90, 95, 99)}

TARGETS = [
    ("TRPM8",  "亚型标志",       "#8C8C8C"),
    ("SST",    "亚型标志",       "#8C8C8C"),
    ("CALCA",  "亚型标志",       "#8C8C8C"),
    ("SCN11A", "高选择性靶点",   "#0F6E56"),
    ("TRPV1",  "高选择性靶点",   "#0F6E56"),
    ("P2RX3",  "高选择性靶点",   "#0F6E56"),
    ("SCN10A", "高选择性靶点",   "#0F6E56"),
    ("SCN9A",  "高选择性靶点",   "#0F6E56"),
    ("OPRM1",  "高选择性靶点",   "#0F6E56"),
    ("OPRD1",  "高选择性靶点",   "#0F6E56"),
    ("CNR1",   "中等选择性",     "#BA7517"),
    ("ADORA1", "无亚型选择性",   "#A32D2D"),
    ("OPRK1",  "无亚型选择性",   "#A32D2D"),
]
rows = []
for g, grp, col in TARGETS:
    if g in s.index:
        r = s.loc[g]
        rows.append((g, grp, col, float(r.max_dAUC),
                     float((v < r.max_dAUC).mean() * 100), r.best_type,
                     float(r.detect_pct)))
D = pd.DataFrame(rows, columns=["gene", "grp", "col", "dauc", "pct", "type", "det"])
D = D.sort_values("dauc", ascending=True).reset_index(drop=True)

fig = plt.figure(figsize=(9.0, 8.6))
gs = GridSpec(2, 2, height_ratios=[1.0, 1.55], width_ratios=[1.15, 1.0],
              hspace=0.50, wspace=0.30,
              left=0.095, right=0.975, top=0.940, bottom=0.115)

# ---------- A ----------
ax = fig.add_subplot(gs[0, 0])
ax.hist(v, bins=90, color="#D3D1C7", edgecolor="white", linewidth=0.3)
ax.axvspan(P[95], v.max(), color="#E1F5EE", alpha=0.85, zorder=0)
for q, ls in ((50, ":"), (90, "--"), (95, "-")):
    ax.axvline(P[q], color="#5F5E5A", lw=0.9, ls=ls, zorder=3)
ax.set_xlim(-0.02, v.max() + 0.03)
ax.set_xlabel("任一基因的最大 ΔAUC（跨 13 类神经元）", fontsize=12, labelpad=15)
ax.set_ylabel("基因数", fontsize=12)
ax.tick_params(labelsize=10.5)
ax.set_xticks([0, 0.2, 0.4])
ax.set_xticks([P[50], P[90], P[95]], minor=True)
ax.set_xticklabels(["P50\n%.3f" % P[50], "P90\n%.3f" % P[90], "P95\n%.3f" % P[95]],
                   minor=True, fontsize=8.6, color="#5F5E5A")
for lb in ax.get_xticklabels(minor=True):
    lb.set_ha("center")
ax.tick_params(axis="x", which="minor", length=8, width=0.8, color="#5F5E5A")
ax.set_title("A  16 969 个基因的亚型选择性经验分布", fontsize=12.2, loc="left", pad=7)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# ---------- C (upper right): detection-rate matched ----------
axc = fig.add_subplot(gs[0, 1])
x = s.detect_pct.values
y = s.max_dAUC.values
axc.scatter(x, y, s=1.6, color="#B4B2A9", alpha=0.13, linewidths=0, rasterized=True)
edges = np.linspace(5, 100, 25)
ctr, med, p95 = [], [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (x >= a) & (x < b)
    if m.sum() >= 30:
        ctr.append((a + b) / 2); med.append(np.median(y[m])); p95.append(np.percentile(y[m], 95))
axc.plot(ctr, med, color="#5F5E5A", lw=1.2, ls="--", label="同检出率基因的中位数")
axc.plot(ctr, p95, color="#0F6E56", lw=1.2, ls="-", alpha=0.55, label="同检出率基因的 P95")
for _, r in D.iterrows():
    axc.scatter(r.det, r.dauc, s=52, color=r.col, zorder=4, linewidths=0.7,
                edgecolors="white")
for g in ("ADORA1", "OPRK1", "CNR1", "OPRM1", "SCN11A"):
    r = D[D.gene == g].iloc[0]
    dy = 0.030 if g in ("ADORA1", "CNR1") else -0.034
    axc.annotate(g, (r.det, r.dauc), xytext=(0, dy * 100), textcoords="offset points",
                 fontsize=10, color=r.col, ha="center",
                 va="bottom" if dy > 0 else "top")
axc.set_xlim(0, 105); axc.set_ylim(-0.02, 0.56)
axc.set_xlabel("该基因在 1 074 个神经元中的检出率 /%", fontsize=11.5)
axc.set_ylabel("最大 ΔAUC", fontsize=12)
axc.tick_params(labelsize=10.5)
axc.legend(fontsize=9.2, frameon=True, framealpha=0.88, edgecolor="none",
           facecolor="white", loc="lower right", handlelength=1.8,
           borderaxespad=0.35, labelspacing=0.25)
axc.set_title("C  检出率匹配后靶点仍居零分布之上", fontsize=12.2, loc="left", pad=7)
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)

# ---------- B ----------
ax2 = fig.add_subplot(gs[1, :])
for q, note in ((50, "P50"), (90, "P90"), (95, "P95")):
    ax2.axvline(P[q], color="#B4B2A9", lw=0.8, ls="--", zorder=0)
    ax2.text(P[q], len(D) + 0.12, note, fontsize=9.6, color="#5F5E5A",
             ha="center", va="bottom")
for i, r in D.iterrows():
    ax2.plot([0, r.dauc], [i, i], color=r.col, lw=1.5, alpha=0.42,
             solid_capstyle="round", zorder=2)
    ax2.scatter(r.dauc, i, s=56, color=r.col, zorder=3, linewidths=0)
    lab = "P>99.9" if r.pct >= 99.95 else "P%.1f" % r.pct
    ax2.text(r.dauc + 0.018, i, lab + " · " + r.type,
             fontsize=9.8, color=r.col, va="center", ha="left", zorder=4)
ax2.set_yticks(np.arange(len(D)))
ax2.set_yticklabels(D.gene, fontsize=11.6)
ax2.set_ylim(-0.75, len(D) + 0.95)
ax2.set_xlim(-0.02, 0.80)
ax2.set_xticks(np.arange(0, 0.71, 0.1))
ax2.set_xlabel("最大 ΔAUC（相对全基因零分布的富集量）", fontsize=11.5)
ax2.tick_params(axis="x", labelsize=10.5)
ax2.set_title("B  镇痛靶点在上述分布中的位置（右侧为全基因组百分位与最高富集细胞类型）",
              fontsize=12.2, loc="left", pad=7)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)

handles = [Line2D([], [], marker="o", ls="", color=c, label=g, markersize=6)
           for g, c in (("亚型标志（参照）", "#8C8C8C"), ("高选择性靶点", "#0F6E56"),
                        ("中等选择性", "#BA7517"), ("无亚型选择性", "#A32D2D"))]
fig.legend(handles=handles, fontsize=10.2, frameon=False, ncol=4,
           loc="lower center", bbox_to_anchor=(0.5, 0.018),
           handletextpad=0.3, columnspacing=1.8)

fig.savefig("fig2_genomewide_benchmark.png", dpi=400, facecolor="white")
print("saved")
print(D[["gene", "grp", "dauc", "pct", "type", "det"]].round(1).to_string(index=False))
