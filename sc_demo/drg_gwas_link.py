#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 DRG 亚型选择性（maxΔAUC）接到人类疼痛遗传学上。

问题：在 DRG 神经元中越是亚型选择性的基因，是否越富集人类慢性疼痛的遗传信号？

指标与模型
----------
基因层面指标来自 gene_level.py（基因体 ±10 kb，GRCh38）。
主统计量：maxΔAUC 与 -log10(P_min) 的 Spearman ρ。
零分布：分层置换——在「SNP 数十分位 × 表达量十分位」层内打乱 maxΔAUC，
        保留基因长度（SNP 密度）与表达水平两个主要混杂；10 000 次。
协变量模型：mlogp ~ maxΔAUC + log10(n_snp) + mean_expr，报告每 +0.1 ΔAUC 的效应。
阴性对照：抑郁（神经精神但非疼痛）、肥胖（代谢），同一流程比较。

输出：drg_gwas_rho.csv、drg_gwas_deciles.csv、drg_gwas_targets.csv
"""
import glob
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent.parent
GW = ROOT / "ext" / "gwas"
SEL = ROOT / "sc_demo" / "gene_selectivity_ranked.csv"
RNG_SEED = 20261109
NPERM = 10000

PAIN = {"G6_MIGRAINE": "偏头痛", "M13_FIBROMYALGIA": "纤维肌痛",
        "M13_ALGONEURO": "神经病理性疼痛", "M13_LOWBACKPAIN": "下背痛"}
CTRL = {"F5_DEPRESSIO": "抑郁", "E4_OBESITY": "肥胖"}
TARG = {"SCN11A": "Nav1.9", "SCN10A": "Nav1.8", "SCN9A": "Nav1.7", "P2RX3": "P2X3",
        "TRPV1": "TRPV1", "OPRM1": "μ 阿片", "OPRD1": "δ 阿片", "OPRK1": "κ 阿片",
        "ADORA1": "A1R", "CNR1": "CB1"}


def strata(df):
    a = pd.qcut(df.n_snp, 10, labels=False, duplicates="drop").astype(str)
    b = pd.qcut(df.mean_expr, 10, labels=False, duplicates="drop").astype(str)
    return (a + "_" + b).values


def perm_rho(x, y, strat, nperm=NPERM, seed=0):
    rng = np.random.default_rng(seed)
    obs = stats.spearmanr(x, y).statistic
    groups = [np.where(strat == s)[0] for s in np.unique(strat)]
    null = np.empty(nperm)
    xw = x.copy()
    for k in range(nperm):
        for ii in groups:
            xw[ii] = rng.permutation(x[ii])
        null[k] = stats.spearmanr(xw, y).statistic
    p = (np.sum(np.abs(null) >= abs(obs)) + 1) / (nperm + 1)
    return obs, p, null


def main():
    sel = pd.read_csv(SEL, index_col=0)[["max_dAUC", "detect_pct", "mean_expr"]]
    rows, dec_rows, targ_rows = [], [], []
    for f in sorted(glob.glob(str(GW / "genelvl_*.csv"))):
        tag = Path(f).stem.replace("genelvl_", "")
        if tag not in PAIN and tag not in CTRL:
            continue
        label = PAIN.get(tag) or CTRL.get(tag)
        g = pd.read_csv(f).set_index("gene")
        j = sel.join(g[["n_snp", "p_min", "mlogp_max", "n_gws", "z_max", "chi2_norm"]], how="inner")
        j = j[(j.n_snp >= 50) & j.p_min.notna() & (j.detect_pct >= 5)]
        if len(j) < 500:
            print("%-22s 可比基因不足" % tag, flush=True)
            continue

        st = strata(j)
        rho, p, _ = perm_rho(j.max_dAUC.values, j.mlogp_max.values, st,
                             seed=RNG_SEED + zlib.crc32(tag.encode()) % 997)
        # 协变量线性模型
        X = np.column_stack([np.ones(len(j)), j.max_dAUC.values,
                             np.log10(j.n_snp.values), j.mean_expr.values])
        beta, *_ = np.linalg.lstsq(X, j.mlogp_max.values, rcond=None)
        resid = j.mlogp_max.values - X @ beta
        dof = len(j) - X.shape[1]
        s2 = (resid ** 2).sum() / dof
        cov = s2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
        tval = beta[1] / se[1]
        p_lm = 2 * stats.t.sf(abs(tval), dof)

        print("%-22s %-8s n=%5d  ρ=%+.4f  置换P=%.4f  每+0.1ΔAUC=%+.3f log10P (P=%.3g)"
              % (tag, label, len(j), rho, p, beta[1] * 0.1, p_lm), flush=True)
        rows.append(dict(表型=label, tag=tag, n=len(j), rho=round(rho, 4),
                         p_perm=p, 每0p1ΔAUC=round(beta[1] * 0.1, 3), p_lm=p_lm,
                         显著基因=int((j.n_gws > 0).sum()),
                         最小P=float(j.p_min.min())))

        d = j.copy()
        d["dec"] = pd.qcut(d.max_dAUC, 10, labels=False) + 1
        agg = d.groupby("dec").agg(n=("p_min", "size"), 中位mlogp=("mlogp_max", "median"),
                                   显著比例=("n_gws", lambda s: (s > 0).mean() * 100))
        agg["表型"] = label
        dec_rows.append(agg.reset_index())

        for gname, short in TARG.items():
            if gname in j.index:
                r = j.loc[gname]
                targ_rows.append(dict(表型=label, 靶点=short, gene=gname,
                                      百分位=round(float(j.loc[gname, "max_dAUC"] /
                                                       j.max_dAUC.max() * 100), 1),
                                      maxΔAUC=round(float(r.max_dAUC), 3),
                                      n_snp=int(r.n_snp), p_min=float(r.p_min),
                                      mlogp=round(float(r.mlogp_max), 2)))

    pd.DataFrame(rows).to_csv(GW / "drg_gwas_rho.csv", index=False, encoding="utf-8-sig")
    pd.concat(dec_rows).to_csv(GW / "drg_gwas_deciles.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(targ_rows).to_csv(GW / "drg_gwas_targets.csv", index=False, encoding="utf-8-sig")
    print("\n-> drg_gwas_rho.csv / drg_gwas_deciles.csv / drg_gwas_targets.csv")


if __name__ == "__main__":
    main()
