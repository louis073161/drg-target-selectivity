#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""步骤 4：出图（图1 图谱与靶点分布；图2 ΔAUC 特异性热图；图3 共表达与特异性谱）"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, matplotlib
from matplotlib import font_manager
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.stats import rankdata

for f in ["C:/Windows/Fonts/msyh.ttc", "C:/Windows/Fonts/simhei.ttf", "C:/Windows/Fonts/simsun.ttc"]:
    try:
        font_manager.fontManager.addfont(f); break
    except Exception:
        pass
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["pdf.fonttype"] = 42

a = sc.read_h5ad("atlas_s3_final.h5ad")
per = pd.read_csv("atlas_target_stats_by_type.csv")
ann = pd.read_csv("atlas_cluster_annotation.csv")
co = pd.read_csv("atlas_coexpression_top20.csv", index_col=0)
dau = pd.read_csv("atlas_delta_auc_by_type.csv", index_col=0)
mean = pd.read_csv("atlas_mean_by_type.csv", index_col=0)

order = ann.sort_values("n_cells", ascending=False).cell_type.tolist()
dau = dau[order]; mean = mean[order]
CATS = a.obs.cell_type.cat.categories.tolist()
short = {"hPEP.TAC1+": "hPEP.TAC1+", "hPEP.NMB+": "hPEP.NMB+", "hPEP.KIT": "hPEP.KIT",
         "hPEP.CHRNA7": "hPEP.CHRNA7", "hSST": "hSST", "hNP1": "hNP1", "hC.LTMR": "hC.LTMR",
         "hAβ.LTMR": "hAβ.LTMR", "hAβ.LTMR.SPP1+": "hAβ.LTMR\n.SPP1+", "hAδ.LTMR": "hAδ.LTMR",
         "hPropr": "hPropr", "hTRPM8": "hTRPM8", "hATF3": "hATF3", "非神经元·巨噬细胞": "非神经元\n巨噬细胞"}
lab = [short[c] for c in order]
COL = dict(zip(CATS, plt.cm.tab20(np.linspace(0, 1, len(CATS)))))

# ============ 图 1 ============
fig = plt.figure(figsize=(7.4, 6.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1, 1], wspace=0.30, hspace=0.45)
ax = fig.add_subplot(gs[:, 0])
um = a.obsm["X_umap"]
for i, c in enumerate(order):
    m = (a.obs.cell_type == c).values
    ax.scatter(um[m, 0], um[m, 1], s=4, color=COL[c], linewidths=0, rasterized=True)

# 标签避让：先取质心，再在归一化坐标里做斥力松弛
cent = np.array([[np.median(um[(a.obs.cell_type == c).values, 0]),
                  np.median(um[(a.obs.cell_type == c).values, 1])] for c in order])
x0, x1 = um[:, 0].min(), um[:, 0].max(); y0, y1 = um[:, 1].min(), um[:, 1].max()
pad_x, pad_y = (x1 - x0) * 0.42, (y1 - y0) * 0.10
nx0, nx1 = x0 - pad_x, x1 + pad_x
ny0, ny1 = y0 - pad_y, y1 + pad_y
# 标签：按质心分左右两列，沿纵坐标贪心排布（最小间距），连线短且不交叉
mid = (x0 + x1) / 2
gap = (ny1 - ny0) * 0.058
UMAP_LAB = {"非神经元·巨噬细胞": "巨噬细胞"}


def pack(idxs):
    idxs = sorted(idxs, key=lambda i: -cent[i, 1])
    ys = [cent[i, 1] for i in idxs]
    for k in range(1, len(ys)):
        ys[k] = min(ys[k], ys[k - 1] - gap)
    for k in range(len(ys) - 2, -1, -1):          # 下边界回推
        ys[k] = min(ys[k], ny1 - (ny1 - ny0) * 0.02 - (len(ys) - 1 - k) * gap / 6)
    return idxs, ys


for side, idxs in (("L", [i for i in range(len(order)) if cent[i, 0] < mid]),
                   ("R", [i for i in range(len(order)) if cent[i, 0] >= mid])):
    idxs, ys = pack(idxs)
    tx = nx0 + (nx1 - nx0) * 0.01 if side == "L" else nx1 - (nx1 - nx0) * 0.01
    for k, i in enumerate(idxs):
        c = order[i]
        ax.annotate(UMAP_LAB.get(c, short[c].replace("\n", "")),
                    xy=(cent[i, 0], cent[i, 1]), xytext=(tx, ys[k]),
                    fontsize=5.8, ha="left" if side == "L" else "right", va="center",
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.9),
                    arrowprops=dict(arrowstyle="-", lw=0.35, color="#999999",
                                    shrinkA=1, shrinkB=2))
ax.set_xlim(nx0, nx1); ax.set_ylim(ny0, ny1)
ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values():
    s.set_visible(False)
ax.set_title("A  人类 DRG 1129 个细胞的 14 个簇", fontsize=8.5, loc="left")
ax.set_xlabel("UMAP-1", fontsize=7); ax.set_ylabel("UMAP-2", fontsize=7)
ax.tick_params(labelsize=6)

# 点图：颜色 = 组均值 z 分数；大小 = 该类型中表达量处于全局前 20% 的细胞比例
genes = ["OPRM1", "OPRD1", "OPRK1", "ADORA1", "CNR1", "CNR2", "SCN10A", "SCN11A", "TRPV1", "TRPM8"]
X = a.X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X, dtype=np.float32)
neu = (a.obs.cell_class != "非神经元").values
typ = a.obs.cell_type.values
frac = {}
for g in genes:
    j = a.var_names.get_loc(g)
    v = X[:, j]
    thr = np.quantile(v[neu], 0.8)
    hi = (v >= thr) & neu if thr > 1e-8 else np.zeros_like(v, dtype=bool)
    frac[g] = np.array([hi[typ == c].mean() * 100 if (typ == c).sum() else 0 for c in order])
frac = pd.DataFrame(frac, index=order).T

ax2 = fig.add_subplot(gs[0, 1])
gmean = {g: float(X[neu, a.var_names.get_loc(g)].mean()) for g in genes}
for j, g in enumerate(genes):
    row = mean.loc[g, order].astype(float)
    z = (row - row.mean()) / (row.std() + 1e-9)
    if gmean[g] < 0.03:                      # 全库近不检出：灰点表示
        ax2.scatter([j] * len(order), np.arange(len(order)), s=5, c="#b0b0b0", edgecolors="none")
        continue
    ax2.scatter([j] * len(order), np.arange(len(order)), s=np.clip(frac.loc[g].values * 1.5, 3, 70),
                c=z.values, cmap="RdBu_r", vmin=-1.5, vmax=1.5, edgecolors="none")
ax2.set_xticks(range(len(genes))); ax2.set_xticklabels(genes, rotation=90, fontsize=6)
ax2.set_yticks(range(len(order))); ax2.set_yticklabels(lab, fontsize=6)
ax2.set_xlim(-0.6, len(genes) - 0.4); ax2.set_ylim(-0.6, len(order) - 0.4)
ax2.invert_yaxis()
ax2.set_title("B  颜色=平均表达（z），大小=前20%细胞比例", fontsize=7.6, loc="left")
for s in ax2.spines.values():
    s.set_linewidth(0.5)

# 柱状：各类型细胞数
ax3 = fig.add_subplot(gs[1, 1])
n = [int((a.obs.cell_type == c).sum()) for c in order]
ax3.barh(range(len(order)), n, color=[COL[c] for c in order], height=0.68)
ax3.set_yticks(range(len(order))); ax3.set_yticklabels(lab, fontsize=6)
ax3.invert_yaxis(); ax3.set_xlabel("细胞数", fontsize=7)
ax3.tick_params(labelsize=6)
ax3.set_title("C  各类型细胞数", fontsize=8.5, loc="left")
for s in ["top", "right"]:
    ax3.spines[s].set_visible(False)
plt.savefig("fig1_atlas_targets.png", dpi=400, bbox_inches="tight")
plt.close()

# ============ 图 2：ΔAUC 热图 ============
ordg = ["OPRM1", "OPRD1", "OPRK1", "OPRL1", "ADORA1", "ADORA2A", "ADORA2B", "ADORA3",
        "CNR1", "CNR2", "GPR55", "SCN10A", "SCN9A", "SCN11A", "TRPV1", "TRPA1", "TRPM8", "P2RX3", "NTRK1"]
d = dau.loc[ordg]
fig, ax = plt.subplots(figsize=(7.4, 5.0))
norm = TwoSlopeNorm(vmin=-0.45, vcenter=0, vmax=0.55)
im = ax.imshow(d.values, cmap="RdBu_r", norm=norm, aspect="auto")
ax.set_xticks(range(len(order))); ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=6.2)
ax.set_yticks(range(len(ordg))); ax.set_yticklabels(ordg, fontsize=7)
for i in range(len(ordg)):
    for j in range(len(order)):
        v = d.values[i, j]
        if abs(v) >= 0.20:
            ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=5.2,
                    color="white" if abs(v) > 0.38 else "black")
ax.set_xticks(np.arange(-.5, len(order), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(ordg), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=0.6)
ax.tick_params(which="minor", length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.015)
cb.set_label("ΔAUC（相对经验零分布）", fontsize=7); cb.ax.tick_params(labelsize=6)
ax.set_title("镇痛靶点在人类 DRG 感觉神经元亚型中的富集程度", fontsize=9, loc="left")
plt.tight_layout()
plt.savefig("fig2_delta_auc_heatmap.png", dpi=400, bbox_inches="tight")
plt.close()

# ============ 图 3：共表达 + 特异性谱 ============
fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.5), gridspec_kw={"width_ratios": [1.15, 1]})
axc = axes[0]
cg = ["OPRM1", "OPRD1", "ADORA1", "CNR1", "SCN10A", "TRPV1", "CALCA", "MRGPRX1"]
cc = co.loc[cg, cg].astype(float)
im = axc.imshow(cc.values, cmap="YlOrRd", vmin=0, vmax=65)
axc.set_xticks(range(len(cg))); axc.set_xticklabels(cg, rotation=45, ha="right", fontsize=6.5)
axc.set_yticks(range(len(cg))); axc.set_yticklabels(cg, fontsize=6.5)
for i in range(len(cg)):
    for j in range(len(cg)):
        axc.text(j, i, "%.0f" % cc.values[i, j], ha="center", va="center", fontsize=5.6,
                 color="white" if cc.values[i, j] > 45 else "black")
axc.set_title("A  前 20% 高表达细胞的共表达率（%）\n独立基线 = 20%", fontsize=8, loc="left")
cb = fig.colorbar(im, ax=axc, fraction=0.03, pad=0.02); cb.ax.tick_params(labelsize=6)

axs = axes[1]
mx = dau.loc[ordg].max(1).sort_values()
axs.barh(range(len(mx)), mx.values, color=["#c0392b" if v >= 0.25 else "#7f8c8d" for v in mx.values], height=0.7)
axs.set_yticks(range(len(mx))); axs.set_yticklabels(mx.index, fontsize=6.5)
axs.axvline(0.10, ls="--", lw=0.8, color="k")
axs.set_xlabel("最大 ΔAUC（亚型选择性强度）", fontsize=7)
axs.tick_params(labelsize=6)
axs.set_title("B  靶点亚型选择性谱\n虚线=0.10（有意义的富集下限）", fontsize=8, loc="left")
for s in ["top", "right"]:
    axs.spines[s].set_visible(False)
plt.tight_layout()
plt.savefig("fig3_coexpression_specificity.png", dpi=400, bbox_inches="tight")
plt.close()
print("图已生成：fig1_atlas_targets.png / fig2_delta_auc_heatmap.png / fig3_coexpression_specificity.png")
