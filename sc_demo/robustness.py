"""稳健性三件套：①指标对照 ②聚类分辨率敏感性 ③标签置换
注意：全程限定在检出率≥5% 的 16 969 个基因上，与主分析同口径。
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARG = ["SCN11A", "SCN10A", "SCN9A", "P2RX3", "TRPV1", "OPRM1", "OPRD1", "OPRK1", "ADORA1", "CNR1"]
NAME = {"SCN11A": "Nav1.9", "SCN10A": "Nav1.8", "SCN9A": "Nav1.7", "P2RX3": "P2X3", "TRPV1": "TRPV1",
        "OPRM1": "μ", "OPRD1": "δ", "OPRK1": "κ", "ADORA1": "A1R", "CNR1": "CB1"}


def dAUC(X, cl):
    """X: cells×genes；cl: 类型。返回 genes×types 的 ΔAUC（减去该类型全基因中位数）"""
    R = stats.rankdata(X.T, method="average", axis=1)
    types = np.unique(cl); n = len(cl)
    out = np.empty((X.shape[1], len(types)), dtype=np.float32)
    for j, t in enumerate(types):
        sub = cl == t; n1 = int(sub.sum()); n2 = n - n1
        auc = (R[:, sub].sum(1) - n1 * (n1 + 1) / 2) / (n1 * n2)
        out[:, j] = auc - np.median(auc)
    return out, types


def main():
    A = sc.read_h5ad(ROOT / "sc_demo" / "atlas_s3_final.h5ad")
    A = A[A.obs.cell_class != "非神经元"].copy()
    Xall = A.X.toarray() if hasattr(A.X, "toarray") else np.asarray(A.X, dtype=np.float32)
    keep = (Xall > 0).mean(0) >= 0.05
    X = Xall[:, keep]; genes = A.var_names.to_numpy()[keep]
    cl = A.obs["cl"].astype(str).values
    print("细胞 %d，类型 %d，基因 %d" % (len(cl), len(np.unique(cl)), X.shape[1]))

    dA, types = dAUC(X, cl)
    mx = dA.max(1); base = np.median(mx)
    pct = stats.rankdata(mx) / len(mx) * 100
    print("maxΔAUC P50=%.4f P75=%.4f P95=%.4f（应与 0.112/0.178/0.320 一致）" %
          (base, *np.percentile(mx, [75, 95])))

    # ---------- ① 指标对照 ----------
    M = np.vstack([X[cl == t].mean(0) for t in types])
    mxm = M.max(0)
    tau = np.where(mxm > 1e-9, (1 - M / np.maximum(mxm, 1e-9)).sum(0) / (len(types) - 1), 0)
    best = dA.argmax(1)
    d = np.empty(X.shape[1])
    for i in range(X.shape[1]):
        a = X[cl == types[best[i]], i]; b = X[cl != types[best[i]], i]
        sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / max(len(a)+len(b)-2, 1))
        d[i] = (a.mean()-b.mean())/sp if sp > 1e-9 else 0

    print("\n=== 全基因组秩相关（n=%d）===" % X.shape[1])
    for nm, v in [("Tau 特异性指数", tau), ("Cohen's d", d)]:
        r = stats.spearmanr(mx, v)
        print("  maxΔAUC vs %-14s ρ=%+.3f (P=%.1e)" % (nm, r.statistic, r.pvalue))

    rows = []
    for g in TARG:
        if g not in genes: continue
        i = int(np.where(genes == g)[0][0])
        rows.append({"基因": g, "靶点": NAME[g], "maxΔAUC": round(float(mx[i]), 3),
                     "百分位": round(float(pct[i]), 1),
                     "Tau": round(float(tau[i]), 3),
                     "Tau百分位": round(float(stats.rankdata(tau)[i]/len(tau)*100), 1),
                     "Cohen_d": round(float(d[i]), 2),
                     "d百分位": round(float(stats.rankdata(d)[i]/len(d)*100), 1)})
    t1 = pd.DataFrame(rows).sort_values("百分位", ascending=False)
    print("\n=== 三种指标下十个靶点的位置 ===")
    print(t1.to_string(index=False))
    R = np.vstack([stats.rankdata(t1["百分位"]), stats.rankdata(t1["Tau百分位"]),
                   stats.rankdata(t1["d百分位"])])
    m, n = R.shape
    S = ((R.sum(0) - R.sum()/n)**2).sum(); W = 12*S/(m**2*(n**3-n))
    print("\nKendall W=%.3f，χ²=%.2f，P=%.4f" % (W, m*(n-1)*W, 1-stats.chi2.cdf(m*(n-1)*W, n-1)))
    hi, lo = ["Nav1.9", "Nav1.8", "P2X3", "TRPV1", "μ", "δ"], ["A1R", "CB1", "κ"]
    for col in ["百分位", "Tau百分位", "d百分位"]:
        a = t1[t1.靶点.isin(hi)][col]; b = t1[t1.靶点.isin(lo)][col]
        print("  %-9s 高选择性组均值 %5.1f vs 无选择性组 %5.1f，差 %5.1f，P=%.3f" %
              (col, a.mean(), b.mean(), a.mean()-b.mean(),
               stats.mannwhitneyu(a, b, alternative="greater")[1]))
    t1.to_csv(ROOT/"sc_demo"/"metric_comparison.csv", index=False, encoding="utf-8-sig")

    # ---------- ② 分辨率敏感性 ----------
    print("\n=== 聚类分辨率敏感性 ===")
    rows = []
    for res in [0.6, 0.8, 1.0, 1.2, 1.5]:
        B = A.copy()
        sc.tl.leiden(B, resolution=res, key_added="cl_tmp", random_state=0,
                     flavor="igraph", n_iterations=2, directed=False)
        c2 = B.obs["cl_tmp"].astype(str).values
        dA2, _ = dAUC(X, c2); mx2 = dA2.max(1)
        p2 = stats.rankdata(mx2)/len(mx2)*100
        r = {"分辨率": res, "簇数": len(np.unique(c2)), "P50": round(float(np.median(mx2)), 4)}
        for g in TARG:
            if g in genes: r[NAME[g]] = round(float(p2[np.where(genes == g)[0][0]]), 1)
        r["高组均值"] = round(np.mean([r[k] for k in ["Nav1.9", "Nav1.8", "P2X3", "TRPV1", "μ", "δ"]]), 1)
        r["低组均值"] = round(np.mean([r[k] for k in ["A1R", "CB1", "κ"]]), 1)
        rows.append(r)
    t2 = pd.DataFrame(rows)
    print(t2.to_string(index=False))
    t2.to_csv(ROOT/"sc_demo"/"resolution_sensitivity.csv", index=False, encoding="utf-8-sig")

    # ---------- ③ 标签置换 ----------
    print("\n=== 细胞类型标签置换（500 次）===")
    rng = np.random.default_rng(0)
    tgt = {g: int(np.where(genes == g)[0][0]) for g in TARG if g in genes}
    obs = {g: float(mx[i]) for g, i in tgt.items()}
    null = {g: [] for g in tgt}
    for _ in range(500):
        csh = rng.permutation(cl)
        parts = [np.where(csh == t)[0] for t in np.unique(csh)]
        for g, i in tgt.items():
            r = stats.rankdata(X[:, i], method="average")
            b = max((r[s].sum() - len(s)*(len(s)+1)/2)/(len(s)*(len(X)-len(s))) for s in parts)
            null[g].append(b - 0.5)
    out = []
    for g, i in tgt.items():
        nl = np.array(null[g]); p = (nl >= obs[g]).mean()
        out.append({"基因": g, "靶点": NAME[g], "观测": round(obs[g], 3),
                    "置换中位": round(float(np.median(nl)), 3),
                    "置换P": "<0.002" if p == 0 else "%.3f" % p,
                    "倍数": round(obs[g]/max(np.median(nl), 1e-6), 1)})
    t3 = pd.DataFrame(out).sort_values("观测", ascending=False)
    print(t3.to_string(index=False))
    t3.to_csv(ROOT/"sc_demo"/"label_permutation.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()
