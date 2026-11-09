#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
人类 DRG 图谱侧独立分析 —— 步骤 2：注释 + 靶点分布统计
输出：cluster 注释表、靶点×细胞类型指标表、统计检验表、类群汇总表
"""
import warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from sklearn.metrics import roc_auc_score

t0 = time.time()
np.random.seed(0)
a = sc.read_h5ad("atlas_s1_clustered.h5ad")
a.obs["cl"] = a.obs["leiden"].astype(str)

# ---------- 注释（依据 diag_panel_by_cluster.csv 的诊断基因证据） ----------
ANN = {
    "0":  ("hPEP.TAC1+",        "CALCA/CALCB/TAC1/NTRK1/RET/TRPV1", "肽能 C 类伤害感受器"),
    "2":  ("hPEP.NMB+",         "NMB/SERTM1/FABP3/NTRK1/SCN10A",     "肽能 C 类伤害感受器"),
    "7":  ("hPEP.KIT",          "KIT/CALCA/NMB/NTRK1/PTGIR",        "肽能 C 类（KIT+）"),
    "9":  ("hPEP.CHRNA7",       "CHRNA7/CD9/AKR1C2/SCN10A",         "肽能 A 类伤害感受器"),
    "12": ("hSST",              "SST/GFRA3/IL31RA/HTR3A/CHRNA3",    "SST+ 类 NP3"),
    "10": ("hNP1",              "MRGPRX1/OSMR/GFRA2/RET/SCN11A",    "非肽能 C 类伤害感受器"),
    "11": ("hC.LTMR",           "CASQ2/GFRA2/RET/POU4F2/TRPV1",     "C 类低阈值机械感受"),
    "5":  ("hAδ.LTMR",          "NTRK3/TPD52L1/CALCB",              "Aδ 低阈值机械感受"),
    "3":  ("hAβ.LTMR.SPP1+",    "SPP1/SFRP1/PDGFD/NEFH/NTRK2",      "有髓 A 类 LTMR（亚群）"),
    "8":  ("hAβ.LTMR.PIEZO2hi", "PIEZO2/QRFPR/CBLN2/PCDH7/KCNV1",   "有髓 A 类 LTMR（亚群）"),
    "1":  ("hPropr",            "PVALB/NEFH/NEFL/MCTP2/WLS",        "本体感觉"),
    "13": ("hTRPM8",            "TRPM8/STUM/ENTPD3",                "冷感受"),
    "6":  ("hATF3",             "ATF3/UCN/GADD45A/ADCYAP1",         "损伤/应激神经元"),
    "4":  ("非神经元·巨噬细胞",  "CD74/C1QC/C3/APOE/LAPTM5",          "免疫细胞（含环境 RNA）"),
}
CLASS = {
    "0": "肽能伤害感受器", "2": "肽能伤害感受器", "7": "肽能伤害感受器", "9": "肽能伤害感受器",
    "12": "SST/NP3", "10": "非肽能伤害感受器", "11": "C-LTMR",
    "5": "Aδ-LTMR", "3": "Aβ-LTMR", "8": "Aβ-LTMR", "1": "本体感觉",
    "13": "TRPM8 冷感受", "6": "损伤应激", "4": "非神经元",
}
a.obs["cell_type"] = a.obs.cl.map(lambda c: ANN[c][0]).astype("category")
a.obs["cell_class"] = a.obs.cl.map(lambda c: CLASS[c]).astype("category")

ann = pd.DataFrame([{"cluster": c, "cell_type": ANN[c][0], "cell_class": CLASS[c],
                     "marker_evidence": ANN[c][1], "note": ANN[c][2],
                     "n_cells": int((a.obs.cl == c).sum())} for c in ANN])
ann = ann.sort_values("n_cells", ascending=False)
ann.to_csv("atlas_cluster_annotation.csv", index=False, encoding="utf-8-sig")
print("[注释] 14 簇 -> %d 个类型 / %d 个类群" % (a.obs.cell_type.nunique(), a.obs.cell_class.nunique()))

# ---------- 靶点面板 ----------
PANEL = [
    ("OPRM1", "μ 阿片受体"), ("OPRD1", "δ 阿片受体"), ("OPRK1", "κ 阿片受体"), ("OPRL1", "痛敏肽受体"),
    ("ADORA1", "腺苷 A1"), ("ADORA2A", "腺苷 A2A"), ("ADORA2B", "腺苷 A2B"), ("ADORA3", "腺苷 A3"),
    ("CNR1", "大麻素 CB1"), ("CNR2", "大麻素 CB2"), ("GPR55", "GPR55"),
    ("SCN10A", "Nav1.8"), ("SCN9A", "Nav1.7"), ("SCN11A", "Nav1.9"),
    ("TRPV1", "TRPV1"), ("TRPA1", "TRPA1"), ("TRPM8", "TRPM8"),
    ("P2RX3", "P2X3"), ("NTRK1", "TrkA"),
]
PANEL = [(g, n) for g, n in PANEL if g in a.var_names]
genes = [g for g, _ in PANEL]
print("[靶点] 可用基因 %d/%d" % (len(genes), 19))

X = a[:, genes].X
X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
obs = a.obs[["cl", "cell_type", "cell_class", "donor"]].copy()

# ---------- 各细胞类型指标 ----------
rows = []
for gname in genes:
    v = X[:, genes.index(gname)].astype(float)
    glob = v.mean()
    for ct in a.obs.cell_type.cat.categories:
        m = (obs.cell_type == ct).values
        if m.sum() < 5:
            continue
        vv = v[m]
        rest = v[~m]
        # 单侧 Wilcoxon（该类型是否高于其余）
        try:
            p = stats.mannwhitneyu(vv, rest, alternative="greater").pvalue
        except Exception:
            p = np.nan
        # AUROC：区分该类型与其余细胞的能力
        y = m.astype(int)
        try:
            auc = roc_auc_score(y, v)
        except Exception:
            auc = np.nan
        # 供体一致性：该类型均值是否高于该供体全体均值
        cons, cons_n = 0, 0
        for d in obs.donor.unique():
            dm = (obs.donor == d).values
            if (dm & m).sum() >= 10:
                cons_n += 1
                if v[dm & m].mean() > v[dm].mean():
                    cons += 1
        rows.append({"gene": gname, "cell_type": ct, "n_cells": int(m.sum()),
                     "mean": vv.mean(), "pct_expr": (vv > 0).mean() * 100,
                     "enrich_ratio": vv.mean() / glob if glob > 0 else np.nan,
                     "auroc": auc, "p": p, "donor_n": cons_n, "donor_cons": cons})
per = pd.DataFrame(rows)

# BH 校正（每个基因内跨类型）
per["q"] = np.nan
for gname, idx in per.groupby("gene").groups.items():
    p = per.loc[idx, "p"].values
    order = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    per.loc[idx, "q"] = np.clip(q, 0, 1)

# 每个基因的整体差异检验（Kruskal-Wallis，仅神经元）
neu = (a.obs.cell_class != "非神经元").values
kw = {}
for gname in genes:
    v = X[neu, genes.index(gname)].astype(float)
    grp = [v[(a.obs.cell_class.values[neu] == c)] for c in a.obs.cell_class.cat.categories
           if c != "非神经元" and (a.obs.cell_class.values[neu] == c).sum() >= 5]
    kw[gname] = stats.kruskal(*grp).pvalue if len(grp) > 1 else np.nan
per["kw_p"] = per.gene.map(kw)
per.round(5).to_csv("atlas_target_stats_by_type.csv", index=False, encoding="utf-8-sig")

# ---------- 类群汇总 ----------
cls_rows = []
for gname in genes:
    v = X[:, genes.index(gname)].astype(float)
    for cl in a.obs.cell_class.cat.categories:
        m = (obs.cell_class == cl).values
        cls_rows.append({"gene": gname, "cell_class": cl, "n_cells": int(m.sum()),
                         "mean": v[m].mean(), "pct_expr": (v[m] > 0).mean() * 100})
byclass = pd.DataFrame(cls_rows).round(4)
byclass.to_csv("atlas_target_by_class.csv", index=False, encoding="utf-8-sig")

# ---------- 主表打印 ----------
piv_mean = per.pivot(index="gene", columns="cell_type", values="mean").reindex(genes)
piv_auc = per.pivot(index="gene", columns="cell_type", values="auroc").reindex(genes)
piv_q = per.pivot(index="gene", columns="cell_type", values="q").reindex(genes)
order = [ANN[c][0] for c in sorted(ANN, key=lambda x: -ann.set_index("cluster").loc[x, "n_cells"])]
order = [c for c in order if c in piv_mean.columns]
piv_mean = piv_mean[order]; piv_auc = piv_auc[order]; piv_q = piv_q[order]

print("\n[A] 平均表达 (log1p CP10K)，按细胞数排序")
print(piv_mean.round(2).to_string())
print("\n[B] AUROC（>0.7 视为该类型相对富集）")
print(piv_auc.round(2).to_string())

print("\n[C] 各基因最强富集类型")
for g in genes:
    r = per[per.gene == g].sort_values("auroc", ascending=False).iloc[0]
    print("  %-8s %-18s AUC=%.2f  mean=%.2f  富集%.1f×  q=%.2g  KW_p=%.2g  供体一致%d/%d"
          % (g, r.cell_type, r.auroc, r["mean"], r.enrich_ratio, r.q, r.kw_p, r.donor_cons, r.donor_n))

piv_mean.round(3).to_csv("atlas_mean_by_type.csv", encoding="utf-8-sig")
piv_auc.round(3).to_csv("atlas_auroc_by_type.csv", encoding="utf-8-sig")
a.write_h5ad("atlas_s2_annotated.h5ad")
print("\n[完成] %.1fs" % (time.time() - t0))
