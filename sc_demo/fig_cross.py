# -*- coding: utf-8 -*-
"""跨物种图（重排版）：双热图 + 共享色条 + 同一把尺子的百分位散点。

面板 C 的口径
-------------
每个靶点取其在 5 个共有大类中的最高 ΔAUC，再分别用**本物种**全基因组 maxΔAUC
经验分布映射为百分位。横纵轴因此是同一把尺子上的两个读数，可直接对照；
不能再画 ΔAUC 对 ΔAUC——两条分布的尺度不同（人类中位 0.066、小鼠 0.065 虽近，
但尾部宽度不同），原始 ΔAUC 的散点会让「谁更靠右」取决于数据集的绝对幅度。
"""
import warnings
warnings.filterwarnings("ignore")
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
from matplotlib.patches import Ellipse

T = pd.read_csv("cross_species_final.csv")
GENES = ["SCN10A", "SCN11A", "P2RX3", "CALCA", "NMB", "MRGPRX1", "OPRM1", "OPRD1",
         "OPRK1", "CNR1", "ADORA1", "TRPM8", "ADORA2A", "CNR2", "GPR55"]
HC = ["NP", "PEP", "C-LTMR", "A-LTMR", "Aδ", "Propr"]
MC = ["NP", "PEP", "C-LTMR", "TH-C", "Aδ", "Propr"]
H = T.set_index("gene").reindex(GENES)[["H_" + c for c in HC]].rename(columns=lambda x: x[2:])
M = T.set_index("gene").reindex(GENES)[["M_" + c for c in MC]].rename(columns=lambda x: x[2:])

fig = plt.figure(figsize=(7.8, 7.15))
gs = fig.add_gridspec(3, 2, height_ratios=[1.62, 0.055, 1.32], hspace=0.44, wspace=0.34,
                      left=0.09, right=0.975, top=0.935, bottom=0.075)

vmax = 0.5
axes = []
for k, (D, ttl) in enumerate([(H, "A  人类  1074 个神经元  Smart-seq2（Yu 2024）"),
                              (M, "B  小鼠  433 个神经元  Smart-seq（GSE59739）")]):
    ax = fig.add_subplot(gs[0, k]); axes.append(ax)
    im = ax.imshow(D.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(D.shape[1]))
    ax.set_xticklabels(D.columns, rotation=40, ha="right", fontsize=10)
    ax.set_yticks(range(D.shape[0])); ax.set_yticklabels(D.index, fontsize=10)
    ax.set_title(ttl, fontsize=10.5, loc="left", pad=10)
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            v = D.values[i, j]
            if np.isnan(v):
                ax.text(j, i, "—", ha="center", va="center", fontsize=9, color="#aaaaaa")
            elif abs(v) >= 0.30:
                ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=7.6,
                        color="white" if abs(v) > 0.40 else "#222222")
    ax.set_xticks(np.arange(-.5, D.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, D.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

cax = fig.add_subplot(gs[1, :])
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("ΔAUC：相对该细胞类型全基因零分布的富集量（>0 表示该型富集）", fontsize=9)
cb.ax.tick_params(labelsize=8.5)

# ---------------- 面板 C：百分位对百分位 ----------------
SHARED = ["NP", "PEP", "C-LTMR", "Aδ", "Propr"]
HG = np.sort(pd.read_csv("human_collapsed_genomewide.csv", index_col=0)["max_dAUC"].values)
MG = np.sort(pd.read_csv("mouse_ref_genomewide_res1.2.csv", index_col=0)["max_dAUC"].values)
DD = pd.read_csv("cross_species_delta_auc.csv")


def pct(v, ref):
    return float((ref < v).mean() * 100)


rows = []
KEEP = ["SCN10A", "SCN11A", "P2RX3", "CALCA", "OPRM1", "TRPV1",
        "ADORA1", "CNR1", "MRGPRX1"]
for g in KEEP:
    s = DD[(DD.species == "human") & (DD.gene == g) & (DD.cls.isin(SHARED))]
    m = DD[(DD.species == "mouse") & (DD.gene == g) & (DD.cls.isin(SHARED))]
    if s.empty or m.empty:
        continue
    hv, mv = s.dAUC.max(), m.dAUC.max()
    rows.append(dict(gene=g, h=pct(hv, HG), m=pct(mv, MG)))
P = pd.DataFrame(rows)
P.to_csv("fig4c_percentiles.csv", index=False, encoding="utf-8-sig")

ax = fig.add_subplot(gs[2, :])
ax.set_xlim(-3, 104); ax.set_ylim(-4, 104)
ax.add_patch(plt.Rectangle((90, 90), 13, 14, fc="#F3F7FB", ec="#B9CBDD", lw=0.8, ls=(0, (3, 2)), zorder=0))
ax.plot([0, 100], [0, 100], ls="--", lw=1.1, color="#888888", zorder=2)
ax.axhline(90, color="#D5E0EA", lw=0.9, zorder=1)
ax.axvline(90, color="#D5E0EA", lw=0.9, zorder=1)

CAT = {"SCN10A": "#C0392B", "SCN11A": "#C0392B", "P2RX3": "#C0392B", "CALCA": "#2E7D32",
       "TRPV1": "#00838F", "MRGPRX1": "#8E24AA", "OPRM1": "#1565C0", "CNR1": "#E67E22",
       "ADORA1": "#6A1B9A"}
for _, r in P.iterrows():
    if r.gene in CAT:
        ax.scatter(r.h, r.m, s=52, c=CAT[r.gene], alpha=0.92, edgecolors="white",
                   linewidths=0.8, zorder=4)

ANN = {"ADORA1": (10, 4, "left"), "CNR1": (10, 4, "left"), "MRGPRX1": (-8, -21, "right"),
       "OPRM1": (9, 4, "left"), "TRPV1": (-9, -11, "right")}
for g, (dx, dy, ha) in ANN.items():
    s = P[P.gene == g]
    if s.empty:
        continue
    ax.annotate(g, (s.h.iloc[0], s.m.iloc[0]), textcoords="offset points", xytext=(dx, dy),
                fontsize=9.4, color=CAT[g], fontweight="bold", ha=ha, zorder=6)
ax.add_patch(Ellipse((99.6, 99.7), 7.6, 5.2, angle=0, fill=False, ec="#C0392B",
                     lw=1.0, ls=(0, (2, 2)), zorder=3))

ax.set_xlabel("人类：该靶点最大 ΔAUC 在本物种全基因组分布中的百分位", fontsize=10)
ax.set_ylabel("小鼠：同一量的百分位", fontsize=10)
ax.set_xticks([0, 25, 50, 75, 90, 100]); ax.set_yticks([0, 25, 50, 75, 90, 100])
ax.tick_params(labelsize=9)
ax.set_title("C  同一把尺子：两侧百分位对照（每点＝一个靶点，取其最高类型；n＝9）",
             fontsize=10.5, loc="left", pad=8)
ax.grid(alpha=0.18, lw=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.text(101, 33,
        "虚线＝两侧完全一致；浅蓝区＝两侧均居前十分之一。\n"
        "红虚线椭圆内为 Nav1.8、Nav1.9、P2X3 与 CGRP 四点，两侧均在 P99 以上。",
        fontsize=8.6, va="top", ha="right", color="#333333", linespacing=1.7,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#dddddd", alpha=0.94))
plt.savefig("fig4_cross_species.png", dpi=400, bbox_inches="tight", facecolor="white")
print("saved")
print(P.sort_values("h").to_string(index=False))
