#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""图 7：镇痛通路的三层架构——受体、配体与效应器的错配。

核心命题：一个靶点能否做到亚型选择性，不只取决于受体表达在哪，
还取决于它的配体由谁提供、效应器是否同处一细胞。
通道类靶点"自给自足"（无需配体），GPCR 类则被配体分布绑死。

Panel A  三条通路的四格齐备度（配体合成酶 / 受体 / 效应器 / 三元同阳），
         末列为配体/受体比值——小于 1 说明本地配体不足以喂饱受体，只能靠外源或
         上位神经元释放（这正是阿片在 DRG 的处境）。
Panel B  可药性下限：即使药物 100% 只结合该靶点，在疼痛相关神经元中的覆盖率
         与在本体感觉神经元中的误伤率之比（覆盖比）。
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

CLASSMAP = {"肽能伤害感受器": "PEP", "非肽能伤害感受器": "NP", "Aδ-LTMR": "Aδ",
            "Aβ-LTMR": "Aβ", "本体感觉": "Propr", "C-LTMR": "C-LTMR",
            "SST/NP3": "SST", "TRPM8 冷感受": "TRPM8", "损伤应激": "ATF3"}
PAIN = ["PEP", "NP", "Aδ", "C-LTMR", "SST", "TRPM8", "ATF3"]

PATHWAYS = [
    ("腺苷—A1R", ["NT5E", "ENTPD1"], ["ADORA1"], ["KCNJ3", "KCNJ6"]),
    ("内源大麻素—CB1", ["DAGLB", "NAPEPLD"], ["CNR1"], ["KCNJ3", "KCNJ6"]),
    ("阿片肽—μ/δ/κ", ["PENK", "POMC", "PDYN"], ["OPRM1", "OPRD1", "OPRK1"], ["KCNJ3", "KCNJ6"]),
]
TARG = ["SCN9A", "SCN10A", "SCN11A", "P2RX3", "TRPV1", "OPRM1", "OPRD1", "OPRK1", "ADORA1", "CNR1"]
LABEL = {"SCN9A": "Nav1.7", "SCN10A": "Nav1.8", "SCN11A": "Nav1.9", "P2RX3": "P2X3",
         "TRPV1": "TRPV1", "OPRM1": "μ 阿片", "OPRD1": "δ 阿片", "OPRK1": "κ 阿片",
         "ADORA1": "A1R", "CNR1": "CB1"}


def load():
    ad = sc.read_h5ad(ROOT / "sc_demo" / "atlas_s3_final.h5ad")
    neu = ad[ad.obs.cell_class != "非神经元"].copy()
    cls = np.array([CLASSMAP.get(c, c) for c in neu.obs.cell_class])
    return neu, cls


def panel_a(neu, cls):
    genes = sorted({g for _, l, r, e in PATHWAYS for g in l + r + e})
    X = neu[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    B = X > 0
    idx = {g: k for k, g in enumerate(genes)}

    mat, cvs = [], []
    for name, lig, rec, eff in PATHWAYS:
        hasL = B[:, [idx[g] for g in lig]].any(1)
        hasR = B[:, [idx[g] for g in rec]].any(1)
        hasE = B[:, [idx[g] for g in eff]].any(1)
        mat.append([hasL.mean(), hasR.mean(), hasE.mean(), (hasL & hasR & hasE).mean()])
        cvs.append(hasL.mean() / max(hasR.mean(), 1e-9))
    return np.array(mat) * 100, np.array(cvs)


def panel_b(neu, cls):
    X = neu[:, TARG].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    B = X > 0
    pain = np.isin(cls, PAIN)
    propr = cls == "Propr"
    rows = []
    for k, g in enumerate(TARG):
        cp = B[pain, k].mean() * 100
        cr = B[propr, k].mean() * 100
        rows.append({"g": g, "pain": cp, "propr": cr, "win": cp / max(cr, 0.5)})
    return pd.DataFrame(rows)


def main():
    neu, cls = load()
    mat, cvs = panel_a(neu, cls)
    cov = panel_b(neu, cls)

    fig = plt.figure(figsize=(7.8, 4.35))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.18], wspace=0.42)

    # ---------------- A 齐备度 ----------------
    axA = fig.add_subplot(gs[0, 0])
    cols = ["配体合成酶", "受体", "效应器", "三元同阳"]
    cmap = LinearSegmentedColormap.from_list("b", ["#ffffff", "#cfe0f3", "#4a7fb5", "#1a3d66"])
    im = axA.imshow(mat, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    axA.set_xticks(range(4))
    axA.set_xticklabels(cols, fontsize=8.4)
    axA.set_yticks(range(len(PATHWAYS)))
    axA.set_yticklabels([p[0] for p in PATHWAYS], fontsize=8.4)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            axA.text(j, i, "%.0f" % v, ha="center", va="center", fontsize=8.2,
                     color="white" if v > 55 else "#2d3748")
    for i, cv in enumerate(cvs):
        col = "#c53030" if cv < 1 else "#2f855a"
        axA.text(4.05, i, "%.2f" % cv, ha="left", va="center", fontsize=8.4,
                 color=col, fontweight="bold")
    axA.text(4.05, -0.80, "配体/受体", ha="left", va="center", fontsize=7.6, color="#4a5568")
    axA.set_xlim(-0.5, 4.5)
    axA.set_ylim(len(PATHWAYS) - 0.5, -1.0)
    axA.set_title("A　通路的四格齐备度（神经元 %）", fontsize=10, loc="left")
    for s in axA.spines.values():
        s.set_visible(False)
    axA.tick_params(length=0)

    # ---------------- B 覆盖比 ----------------
    axB = fig.add_subplot(gs[0, 1])
    cov = cov.sort_values("win")
    y = np.arange(len(cov))
    axB.barh(y, cov["pain"], height=0.62, color="#4a7fb5", label="疼痛相关神经元")
    axB.barh(y, cov["propr"], height=0.28, color="#e0a26a", label="本体感觉神经元")
    for k, r in enumerate(cov.itertuples()):
        good = r.win >= 1.5
        bad = r.win < 1.0
        col = "#2f855a" if good else ("#c53030" if bad else "#4a5568")
        txt = ("%.2f" % r.win) if r.win < 1 else ("%.1f" % r.win)
        axB.text(max(r.pain, r.propr) + 1.6, k, txt, va="center",
                 fontsize=8.2, color=col, fontweight="bold" if (good or bad) else "normal")
    axB.set_yticks(y)
    axB.set_yticklabels([LABEL[g] for g in cov["g"]], fontsize=8.4)
    axB.set_xlim(0, 112)
    axB.set_xlabel("检出率（%）", fontsize=9.2)
    axB.set_title("B　覆盖与误伤：覆盖比＝疼痛／本体感觉", fontsize=10, loc="left")
    axB.axvline(50, ls=":", lw=0.8, c="#cbd5e0", zorder=0)
    axB.text(101, len(cov) - 0.35, "覆盖比", fontsize=7.6, color="#4a5568", ha="center")
    for s in ["top", "right"]:
        axB.spines[s].set_visible(False)
    axB.tick_params(labelsize=8.2)
    axB.legend(fontsize=7.4, frameon=False, loc="lower right", bbox_to_anchor=(1.0, -0.02))

    fig.savefig(ROOT / "sc_demo" / "fig7_pathway_architecture.png", dpi=400,
                facecolor="white", bbox_inches="tight")
    print("saved fig7")

    # 文本表
    t = pd.DataFrame(mat.round(1), index=[p[0] for p in PATHWAYS],
                     columns=["配体+", "受体+", "效应器+", "三元同阳"])
    t["配体/受体"] = cvs.round(2)
    t.to_csv(ROOT / "sc_demo" / "pathway_architecture.csv", encoding="utf-8-sig")
    print(t.to_string())
    print()
    print(cov.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
