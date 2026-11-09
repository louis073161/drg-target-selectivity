# -*- coding: utf-8 -*-
"""功效模拟：若 A1R（或任一低表达 GPCR）在某一亚型中确有 1.5～3 倍富集，
本流程能否检出？做法是保留 A1R 的检出率与阳性细胞表达量分布，
只改变阳性细胞在各类之间的归属（注入富集），按同一函数重算 maxΔAUC 与百分位。"""
import numpy as np, pandas as pd
import anndata as ad
from scipy.stats import rankdata

rng = np.random.default_rng(0)

a = ad.read_h5ad("atlas_s3_final.h5ad")
ct = a.obs["cell_type"].astype(str).values
keep = np.array(["巨噬" not in str(t) for t in ct])
a = a[keep].copy()
ct = ct[keep]
types = [t for t in pd.unique(ct)]
N = len(ct)

base = pd.read_csv("atlas_auc_baseline.csv", index_col=0).iloc[:, 0].to_dict()
gw = pd.read_csv("genomewide_delta_auc.csv", index_col=0).max(axis=1).values

gene = "ADORA1"
v = np.asarray(a[:, gene].X.todense()).ravel()
print("细胞数 %d，类型 %d，%s 阳性 %d（%.1f%%）"
      % (N, len(types), gene, (v > 0).sum(), (v > 0).mean() * 100))


def max_dauc(vec):
    rk = rankdata(vec)  # 平均秩处理并列
    best = -9
    for t in types:
        m = ct == t
        n = m.sum()
        R = rk[m].sum()
        auc = (R - n * (n + 1) / 2) / (n * (N - n))
        best = max(best, auc - base.get(t, 0.5))
    return best


obs = max_dauc(v)
print("观测 maxΔAUC = %.3f，百分位 = %.1f" % (obs, (gw < obs).mean() * 100))

pos_rate = {t: (v[ct == t] > 0).mean() for t in types}
print("\n各类型 A1R 阳性率：")
for t in types:
    print("  %-16s n=%3d  %.1f%%" % (t, (ct == t).sum(), pos_rate[t] * 100))

vals = v[v > 0]
K = len(vals)
TARGET = max(types, key=lambda t: (ct == t).sum())  # 最大的类型，作为"若富集"的落点
TARGET = "hNP1"
print("\n注入富集的落点：%s（n=%d，观测阳性率 %.1f%%）"
      % (TARGET, (ct == TARGET).sum(), pos_rate[TARGET] * 100))

rows = []
for fold in (1.5, 2.0, 3.0):
    res, pcts = [], []
    for _ in range(200):
        k_t = min(int(round(min(fold * pos_rate[TARGET], 1.0) * (ct == TARGET).sum())),
                  (ct == TARGET).sum())
        obs_t = int((v[ct == TARGET] > 0).sum())
        need_other = K - k_t
        obs_other = K - obs_t
        scale = need_other / obs_other if obs_other > 0 else 0
        chosen = []
        for t in types:
            n_t = (ct == t).sum()
            idx = np.where(ct == t)[0]
            if t == TARGET:
                k = k_t
            else:
                k = int(round((v[idx] > 0).sum() * scale))
            k = min(k, n_t)
            if k > 0:
                chosen.append(rng.choice(idx, size=k, replace=False))
        idx_new = np.concatenate(chosen)
        # 阳性细胞数不足 K 时补齐（随机从其余细胞取）
        if len(idx_new) < K:
            rest = np.setdiff1d(np.arange(N), idx_new)
            idx_new = np.concatenate([idx_new, rng.choice(rest, size=K - len(idx_new),
                                                          replace=False)])
        idx_new = idx_new[:K]
        sim = np.zeros(N)
        sim[idx_new] = rng.permutation(vals)
        d = max_dauc(sim)
        res.append(d)
        pcts.append((gw < d).mean() * 100)
    res = np.array(res); pcts = np.array(pcts)
    rows.append([fold, np.median(res), np.percentile(res, 5), np.percentile(res, 95),
                 np.median(pcts), (pcts >= 95).mean() * 100])

out = pd.DataFrame(rows, columns=["富集倍数", "maxΔAUC中位", "5%分位", "95%分位",
                                  "百分位中位", "判为选择性/%"])
out.to_csv("power_simulation.csv", index=False, encoding="utf-8-sig")
print("\n" + out.round(3).to_string(index=False))
