#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
人类 DRG 图谱侧独立分析 —— 步骤 3（修正版）
1) 按原文判据纠正 A 类纤维命名
2) 经验零分布校正的富集统计（ΔAUC），消除"某类细胞整体表达偏高"的偏倚
3) 秩次共表达分析（20% 分位，独立基线 20%）
4) 基因丰度百分位（同一数据集内可比）
"""
import warnings, time
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from scipy.stats import rankdata

t0 = time.time()
np.random.seed(0)
a = sc.read_h5ad("atlas_s1_clustered.h5ad")
a.obs["cl"] = a.obs["leiden"].astype(str)

# ---------- 注释 v2：按原文（Yu et al., Nat Neurosci 2024）判据 ----------
ANN = {
    "0":  ("hPEP.TAC1+",       "CALCA/CALCB/TAC1/NTRK1/RET/TRPV1",  "肽能 C 类伤害感受器"),
    "2":  ("hPEP.NMB+",        "NMB/SERTM1/FABP3/NTRK1/SCN10A",      "肽能 C 类伤害感受器"),
    "7":  ("hPEP.KIT",         "KIT/CALCA/NMB/NTRK1/PTGIR",         "肽能 C 类（KIT+，偏内脏）"),
    "9":  ("hPEP.CHRNA7",      "CHRNA7/CD9/AKR1C2/SCN10A",          "肽能 A 类伤害感受器"),
    "12": ("hSST",             "SST/GFRA3/IL31RA/HTR3A/CHRNA3",     "SST+（对应啮齿类 NP3）"),
    "10": ("hNP1",             "MRGPRX1/OSMR/GFRA2/RET/SCN11A",     "非肽能 C 类伤害感受器"),
    "11": ("hC.LTMR",          "CASQ2/GFRA2/RET/POU4F2/TRPV1",      "C 类低阈值机械感受"),
    "5":  ("hAβ.LTMR",         "NTRK3 高 / NTRK2 低",                "Aβ 低阈值机械感受"),
    "8":  ("hAδ.LTMR",         "NTRK2 高 / PIEZO2 高 / KCNV1/低 NTRK3", "Aδ 低阈值机械感受"),
    "3":  ("hAβ.LTMR.SPP1+",   "SPP1/SFRP1/PDGFD/NEFH/NTRK3",       "Aβ-LTMR 亚群（分辨率细分）"),
    "1":  ("hPropr",           "PVALB/REEP5/NEFH/MCTP2/WLS",         "本体感觉"),
    "13": ("hTRPM8",           "TRPM8/STUM/ENTPD3",                  "冷感受"),
    "6":  ("hATF3",            "ATF3/UCN/GADD45A/ADCYAP1",           "损伤/应激（大直径）"),
    "4":  ("非神经元·巨噬细胞", "CD74/C1QC/C3/APOE/LAPTM5",           "免疫细胞（含环境 RNA）"),
}
CLASS = {
    "0": "肽能伤害感受器", "2": "肽能伤害感受器", "7": "肽能伤害感受器", "9": "肽能伤害感受器",
    "12": "SST/NP3", "10": "非肽能伤害感受器", "11": "C-LTMR",
    "5": "Aβ-LTMR", "3": "Aβ-LTMR", "8": "Aδ-LTMR", "1": "本体感觉",
    "13": "TRPM8 冷感受", "6": "损伤应激", "4": "非神经元",
}
a.obs["cell_type"] = a.obs.cl.map(lambda c: ANN[c][0]).astype("category")
a.obs["cell_class"] = a.obs.cl.map(lambda c: CLASS[c]).astype("category")

ann = pd.DataFrame([{"cluster": c, "cell_type": ANN[c][0], "cell_class": CLASS[c],
                     "marker_evidence": ANN[c][1], "note": ANN[c][2],
                     "n_cells": int((a.obs.cl == c).sum())} for c in ANN])
ann = ann.sort_values("n_cells", ascending=False)
ann.to_csv("atlas_cluster_annotation.csv", index=False, encoding="utf-8-sig")

PANEL = [("OPRM1", "μ 阿片受体"), ("OPRD1", "δ 阿片受体"), ("OPRK1", "κ 阿片受体"), ("OPRL1", "痛敏肽受体"),
         ("ADORA1", "腺苷 A1"), ("ADORA2A", "腺苷 A2A"), ("ADORA2B", "腺苷 A2B"), ("ADORA3", "腺苷 A3"),
         ("CNR1", "大麻素 CB1"), ("CNR2", "大麻素 CB2"), ("GPR55", "GPR55"),
         ("SCN10A", "Nav1.8"), ("SCN9A", "Nav1.7"), ("SCN11A", "Nav1.9"),
         ("TRPV1", "TRPV1"), ("TRPA1", "TRPA1"), ("TRPM8", "TRPM8"),
         ("P2RX3", "P2X3"), ("NTRK1", "TrkA")]
PANEL = [(g, n) for g, n in PANEL if g in a.var_names]
genes = [g for g, _ in PANEL]

# ---------- 全转录组秩变换（用于 ΔAUC 与共表达） ----------
Xall = a.X
Xall = np.asarray(Xall.todense()) if hasattr(Xall, "todense") else np.asarray(Xall, dtype=np.float32)
det = (Xall > 0).mean(0)
keep_idx = np.where(det >= 0.05)[0]           # 至少 5% 细胞检出
# 靶点基因无论检出率高低都保留（低表达靶点本身是结论的一部分）
pan_idx = [a.var_names.get_loc(g) for g, _ in PANEL]
keep_idx = np.union1d(keep_idx, np.array(pan_idx))
E = Xall[:, keep_idx].astype(np.float32)
E = np.log1p(E / (E.sum(1, keepdims=True) + 1e-9) * 1e4) if False else E  # X 已是 log1p CP10K
EG = a.var_names[keep_idx].tolist()
print("[秩变换] 检出率≥5%% 的基因 %d / %d" % (len(EG), a.n_vars))
R = rankdata(E, axis=0)                        # 每个基因在 1129 个细胞中的秩

ct = a.obs.cell_type.values
cats = [c for c in a.obs.cell_type.cat.categories]
n_all = len(ct)

def auc_matrix(Rm, mask):
    n1 = mask.sum(); n0 = n_all - n1
    mean_rank = Rm[mask].mean(0)
    return (mean_rank - (n1 + 1) / 2) / n0

# 经验零分布：所有基因的 AUC 中位数
base = {}
auc_all = {}
for c in cats:
    m = ct == c
    auc_all[c] = auc_matrix(R, m)
    base[c] = np.median(auc_all[c])
pd.Series(base).round(3).to_csv("atlas_auc_baseline.csv", encoding="utf-8-sig")
print("[零分布] 各类型 AUC 中位数:", {k: round(v, 3) for k, v in base.items()})

gidx = {g: i for i, g in enumerate(EG)}
rows = []
for g in genes:
    j = gidx[g]
    v = E[:, j]
    glob = v.mean()
    # 丰度百分位
    pct_rank = (E.mean(0) < glob).mean() * 100
    for c in cats:
        m = ct == c
        vv, rest = v[m], v[~m]
        auc = auc_all[c][j]
        p = stats.mannwhitneyu(vv, rest, alternative="greater").pvalue
        cons, cons_n = 0, 0
        for d in a.obs.donor.cat.categories:
            dm = (a.obs.donor == d).values
            if (dm & m).sum() >= 10:
                cons_n += 1
                cons += int(v[dm & m].mean() > v[dm].mean())
        rows.append(dict(gene=g, cell_type=c, n_cells=int(m.sum()), mean=vv.mean(),
                         pct_expr=float((vv > 0).mean() * 100), global_mean=glob,
                         enrich_ratio=vv.mean() / glob if glob > 0 else np.nan,
                         auroc=auc, delta_auc=auc - base[c], p=p,
                         donor_n=cons_n, donor_cons=cons))
per = pd.DataFrame(rows)

per["q"] = np.nan
for g, idx in per.groupby("gene").groups.items():
    p = per.loc[idx, "p"].values; n = len(p); order = np.argsort(p)
    q = np.empty(n); q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    per.loc[idx, "q"] = np.clip(q, 0, 1)

neu = (a.obs.cell_class != "非神经元").values
kw = {}
for g in genes:
    v = E[neu, gidx[g]]
    grp = [v[(ct[neu] == c)] for c in cats if (ct[neu] == c).sum() >= 5]
    kw[g] = stats.kruskal(*grp).pvalue if len(grp) > 1 else np.nan
per["kw_p"] = per.gene.map(kw)
per.round(5).to_csv("atlas_target_stats_by_type.csv", index=False, encoding="utf-8-sig")

# ---------- 类群汇总 ----------
cls_rows = []
for g in genes:
    v = E[:, gidx[g]]
    for cl in a.obs.cell_class.cat.categories:
        m = (a.obs.cell_class == cl).values
        cls_rows.append({"gene": g, "cell_class": cl, "n_cells": int(m.sum()),
                         "mean": v[m].mean(), "pct_expr": float((v[m] > 0).mean() * 100)})
byclass = pd.DataFrame(cls_rows).round(4)
byclass.to_csv("atlas_target_by_class.csv", index=False, encoding="utf-8-sig")

# ---------- 秩次共表达（前 20%） ----------
CO = genes + ["CALCA", "MRGPRX1", "SLC17A8", "TH", "NEFL"]
CO = [g for g in CO if g in gidx]
neu_mask = (a.obs.cell_class != "非神经元").values
top = {}
for g in CO:
    v = E[:, gidx[g]]
    thr = np.quantile(v[neu_mask], 0.8)
    top[g] = (v >= thr) & neu_mask
co = pd.DataFrame(index=CO, columns=CO, dtype=float)
for g1 in CO:
    for g2 in CO:
        co.loc[g1, g2] = (top[g2][top[g1]]).mean() * 100 if top[g1].sum() else np.nan
co.round(1).to_csv("atlas_coexpression_top20.csv", encoding="utf-8-sig")

# ---------- 打印 ----------
order = ann.sort_values("n_cells", ascending=False).cell_type.tolist()
piv = per.pivot(index="gene", columns="cell_type", values="mean").reindex(genes)[order]
dau = per.pivot(index="gene", columns="cell_type", values="delta_auc").reindex(genes)[order]
print("\n[A] 平均表达 (log1p CP10K)")
print(piv.round(2).to_string())
print("\n[B] ΔAUC（相对经验零分布；>0.10 视为有意义富集）")
print(dau.round(2).to_string())
print("\n[C] 最强富集类型（按 ΔAUC）")
abund = (E.mean(0))
for g in genes:
    r = per[per.gene == g].sort_values("delta_auc", ascending=False).iloc[0]
    tt = per[per.gene == g].nlargest(2, "delta_auc").cell_type.tolist()
    pra = (abund < E[:, gidx[g]].mean()).mean() * 100
    print("  %-8s %-16s ΔAUC=%.2f AUC=%.2f 均值%.2f 丰度P%.0f q=%.1e KW=%.1e 前二: %s"
          % (g, r.cell_type, r.delta_auc, r.auroc, r["mean"], pra, r.q, r.kw_p, "/".join(tt)))
print("\n[D] 共表达率（行→列，前20%重叠%，独立基线=20%）")
print(co.round(1).loc[["OPRM1", "OPRD1", "ADORA1", "CNR1", "SCN10A", "TRPV1"],
                      ["CALCA", "MRGPRX1", "TRPV1", "SCN10A", "OPRM1", "ADORA1", "CNR1"]].to_string())

piv.round(3).to_csv("atlas_mean_by_type.csv", encoding="utf-8-sig")
dau.round(3).to_csv("atlas_delta_auc_by_type.csv", encoding="utf-8-sig")
a.write_h5ad("atlas_s3_final.h5ad")
print("\n[完成] %.1fs" % (time.time() - t0))
