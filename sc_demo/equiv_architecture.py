#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""3.10 节的统计底座：等效性检验、贝叶斯因子、供体重测信度。

口径与主分析完全一致
--------------------
ΔAUC(g,t) ＝ AUC(g,t) − median_{g'} AUC(g',t)        （中位数在同一簇内、全部
基因上取，不是分块取）
maxΔAUC(g) ＝ max_t ΔAUC(g,t)
全基因组经验分布由 16 969 个可评价基因（1 074 个神经元中检出率≥5%）构成，
与 selectivity_family.py 的 maxΔAUC 分布逐基因一致（P50＝0.112）。

等效性（TOST）
--------------
问的不是"有没有差异"，而是"能否证明它没有值得一提的选择性"。
等价界 δ 由数据自身标定：δ ＝ P75 − P50，即从全基因组中位数走到前四分之一
所需的距离。靶点相对中位基因的偏移落在 δ 内，即判为"与中位基因不可区分"。

贝叶斯因子
----------
Savage–Dickey 密度比：先验 Cauchy(0, δ)，后验取自细胞 bootstrap（固定各簇中位数）。
BF01 ＞ 3 记为"支持无选择性"的实质证据。
"""
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RNG = np.random.default_rng(20261108)

TARG = ["SCN11A", "SCN10A", "SCN9A", "P2RX3", "TRPV1", "OPRM1", "OPRD1", "OPRK1", "ADORA1", "CNR1"]


def auc_matrix(R, cl, labels):
    """R: 基因×细胞 的秩矩阵；返回 簇×基因 的 AUC。"""
    n = R.shape[1]
    out = {}
    for t in labels:
        sub = cl == t
        n1 = int(sub.sum())
        n2 = n - n1
        rs = R[:, sub].sum(1)
        out[t] = (rs - n1 * (n1 + 1) / 2) / (n1 * n2)
    return out


def main():
    ad = sc.read_h5ad(ROOT / "sc_demo" / "atlas_s3_final.h5ad")
    ad = ad[ad.obs["cell_class"] != "非神经元"].copy()   # 与主分析同：只留 1 074 个神经元
    cl = ad.obs["cl"].astype(str).values
    labels = np.unique(cl)

    # ---------- 1. 全基因组基准 ----------
    print("[1/4] 全基因组基准 ...")
    det = np.asarray((ad.X > 0).mean(0)).ravel()
    keep = det >= 0.05                                # 检出率≥5%，不强制纳入任何基因
    G = ad.var_names[keep]
    X = ad[:, G].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    R = stats.rankdata(X.T, method="average", axis=1)   # 基因×细胞
    aucs = auc_matrix(R, cl, labels)
    n_gene = len(G)

    dauc = np.empty((len(labels), n_gene))
    for i, t in enumerate(labels):
        med = np.median(aucs[t])
        dauc[i] = aucs[t] - med
    max_dauc = dauc.max(0)
    P50 = float(np.percentile(max_dauc, 50))
    P75 = float(np.percentile(max_dauc, 75))
    P95 = float(np.percentile(max_dauc, 95))
    delta_bound = P75 - P50
    print("    基因 %d；P50＝%.3f  P75＝%.3f  P95＝%.3f；δ＝%.3f" % (n_gene, P50, P75, P95, delta_bound))

    # 靶点实测
    obs = {}
    for g in TARG:
        k = list(G).index(g) if g in G else None
        if k is None:
            continue
        vals = dauc[:, k]
        t_best = labels[int(np.argmax(vals))]
        obs[g] = dict(best_t=t_best, maxv=float(vals.max()),
                      pct=float((max_dauc < vals.max()).mean() * 100))
    print("    靶点 maxΔAUC 与百分位：%s" % {g: (round(v["maxv"], 3), round(v["pct"], 1))
                                     for g, v in obs.items()})

    # ---------- 2. bootstrap ----------
    print("[2/4] Bootstrap（300 次） ...")
    Xt = ad[:, TARG].X
    Xt = Xt.toarray() if hasattr(Xt, "toarray") else np.asarray(Xt)
    n = ad.n_obs
    B = 300
    boot = {t: np.zeros((B, len(TARG))) for t in labels}
    for b in range(B):
        idx = RNG.integers(0, n, n)
        clb = cl[idx]
        Rt = stats.rankdata(Xt[idx].T, method="average", axis=1)
        for t in labels:
            sub = clb == t
            n1 = int(sub.sum()); n2 = n - n1
            if n1 < 5:
                boot[t][b] = np.nan; continue
            rs = Rt[:, sub].sum(1)
            boot[t][b] = (rs - n1 * (n1 + 1) / 2) / (n1 * n2) - np.median(aucs[t])
        if (b + 1) % 100 == 0:
            print("    %d/%d" % (b + 1, B))

    # ---------- 3. 等效性 ----------
    print("[3/4] 等效性检验 ...")
    rows = []
    for j, g in enumerate(TARG):
        if g not in obs:
            continue
        bt, mv = obs[g]["best_t"], obs[g]["maxv"]
        # 零假设是"与全基因组中位基因不可区分"，故中心化到 P50
        centre = mv - P50
        dd = boot[bt][:, j]
        dd = dd[~np.isnan(dd)] - P50
        se = dd.std(ddof=1)
        t_lo = (centre + delta_bound) / se
        t_hi = (delta_bound - centre) / se
        p_tost = max(1 - stats.norm.cdf(t_lo), 1 - stats.norm.cdf(t_hi))
        prior_at0 = 1.0 / (np.pi * delta_bound)
        post_at0 = float(stats.gaussian_kde(dd)(0.0)[0])
        bf01 = post_at0 / prior_at0
        lo, hi = np.percentile(dd, [2.5, 97.5])
        rows.append(dict(基因=g, 最富集簇=bt, maxΔAUC=round(mv, 3), 偏移=round(centre, 3),
                         百分位=round(obs[g]["pct"], 1),
                         CI低=round(lo, 3), CI高=round(hi, 3), SE=round(se, 4),
                         P_TOST=p_tost, BF01=bf01))
    eq = pd.DataFrame(rows)
    eq["判定"] = np.where(eq.P_TOST < 0.05, "等效", "非等效")
    print(eq.round(3).to_string(index=False))
    eq.to_csv(ROOT / "sc_demo" / "equivalence_tests.csv", index=False, encoding="utf-8-sig")

    # ---------- 4. ICC ----------
    print("[4/4] 供体重测信度 ...")
    donors = sorted(ad.obs["donor"].unique())
    icc = {}
    for g in TARG:
        k = list(G).index(g) if g in G else None
        if k is None:
            continue
        M = np.zeros((len(labels), len(donors)))
        for a, d in enumerate(donors):
            m = (ad.obs["donor"] == d).values
            Xd = Xt[m]
            Rd = stats.rankdata(Xd.T, method="average", axis=1)
            cd = cl[m]
            for i, t in enumerate(labels):
                sub = cd == t
                M[i, a] = dauc[i, k] if sub.sum() < 5 else (
                    (Rd[:, sub].sum(1)[list(TARG).index(g)] - sub.sum() * (sub.sum() + 1) / 2)
                    / (sub.sum() * (len(cd) - sub.sum())) - np.median(aucs[t]))
        gm = M.mean(); rm = M.mean(1, keepdims=True); cm = M.mean(0, keepdims=True)
        msr = len(donors) * ((rm - gm) ** 2).sum() / (len(labels) - 1)
        msc = len(labels) * ((cm - gm) ** 2).sum() / (len(donors) - 1)
        mse = ((M - rm - cm + gm) ** 2).sum() / ((len(labels) - 1) * (len(donors) - 1))
        icc[g] = (msr - mse) / (msr + (len(donors) - 1) * mse + len(donors) * (msc - mse) / len(labels))
    print({k: round(v, 3) for k, v in icc.items()})
    pd.DataFrame({"基因": list(icc), "ICC21": [round(v, 3) for v in icc.values()]}).to_csv(
        ROOT / "sc_demo" / "icc_by_target.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
