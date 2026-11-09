#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""第三套数据：GSE201586（Jung 等 Nat Commun 2023）18 个样本人 DRG snRNA-seq。

与前两套用完全相同的 ΔAUC 定义，但实现改为**按基因分块**：
18 个样本合并后细胞数上万，原始实现把 细胞×全基因 一次性 densify 会吃到数 GB，
分块后峰值内存降到 数百 MB，而结果在数学上完全等价（每个基因的秩只依赖该列）。

变换等价性说明：
  ΔAUC(g,c) = AUC(g,c) − median_{g'} AUC(g',c)
  每个基因的秩是在「全部细胞」上算的，与其它基因无关 → 可分块；
  簇内零分布取该簇在所有基因上的 AUC 中位数 → 先把 基因×簇 的 AUC 矩阵整体存下
  （约 2 万 × 20 个浮点，可忽略），循环结束后再按列取中位数，与一次算完完全一致。

用法：
  python run_201586.py                 # 全流程
  python run_201586.py --verify ext/GSE310000.h5ad   # 顺带核对分块实现与原始实现一致
"""
import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parent
EXT = ROOT.parent
sys.path.insert(0, str(EXT.parent / "sc_demo"))
from replicate_external import TARGETS, DISCOVERY, DISCOVERY_P50, delta_auc  # noqa: E402

H5 = ROOT / "h5v2" if (ROOT / "h5v2").exists() else ROOT / "h5"
NEURON = ["NEFL", "NEFH", "PRPH", "STMN2", "ELAVL2", "POU4F1", "SCN9A", "SCN10A", "TUBB3"]
GLIA = ["KCNJ10", "SLC1A3", "GFAP", "MPZ", "SOX10", "CD74", "COL1A1", "DCN"]


# ---------------------------------------------------------------- ΔAUC（分块）
def delta_auc_blocks(ad, targets, min_detect=0.05, block=1200):
    """与 replicate_external.delta_auc 数学等价，但按基因分块以控内存。"""
    X = ad.X
    genes = ad.var_names.to_numpy()
    n_cells, n_genes = X.shape
    cl = ad.obs["cl"].to_numpy()
    cats = sorted(pd.unique(cl), key=lambda c: -int((cl == c).sum()))
    masks = {c: (cl == c) for c in cats}
    n1 = {c: int(m.sum()) for c, m in masks.items()}

    det = np.asarray((X > 0).mean(0)).ravel() if hasattr(X, "tocsr") else (X > 0).mean(0)
    det = np.asarray(det).ravel()

    pan = [ad.var_names.get_loc(g) for g, _ in targets if g in ad.var_names]
    keep = np.where(det >= min_detect)[0]
    if pan:
        keep = np.union1d(keep, np.array(pan, dtype=int))

    auc = np.zeros((len(keep), len(cats)), dtype=np.float64)
    for b0 in range(0, len(keep), block):
        idx = keep[b0:b0 + block]
        E = X[:, idx]
        E = E.toarray() if hasattr(E, "toarray") else np.asarray(E)
        R = rankdata(E.astype(np.float32), axis=0)      # 每个基因在全部细胞上的秩
        for j, c in enumerate(cats):
            m = masks[c]
            auc[b0:b0 + len(idx), j] = (R[m].mean(0) - (n1[c] + 1) / 2) / (n_cells - n1[c])

    base = np.median(auc, axis=0)                        # 每个簇的全基因组 AUC 中位数
    d = auc - base[None, :]
    jmax = d.argmax(1)
    dmax = d[np.arange(len(keep)), jmax]
    df = pd.DataFrame({"gene": genes[keep], "max_dAUC": dmax,
                       "best_cluster": [cats[j] for j in jmax],
                       "detect_pct": det[keep] * 100})
    df["percentile"] = df.max_dAUC.rank(pct=True) * 100
    return df, {c: int(m.sum()) for c, m in masks.items()}, float(np.median(dmax))


# ---------------------------------------------------------------- 载入
def donor_map():
    """从 sample_meta.txt（GEO 样本元数据）取 GSM → 供体号。"""
    f = ROOT / "sample_meta.txt"
    m = {}
    if f.exists():
        for ln in f.read_text(encoding="utf-8").splitlines():
            if not ln.startswith("GSM"):
                continue
            gsm = ln.split()[0]
            for tok in ln.split("|"):
                if tok.strip().startswith("donor #:"):
                    m[gsm] = tok.split(":", 1)[1].strip()
    return m


def load_all():
    files = sorted(H5.glob("*_filtered_feature_bc_matrix.h5"))
    if not files:
        sys.exit("没有找到 h5 文件，请先下载：python rep/seg_batch.py --list ext/gse201586/manifest.tsv")
    dmap = donor_map()
    ads = []
    for f in files:
        a = sc.read_10x_h5(str(f))
        a.var_names_make_unique()
        gsm = f.name.split("_")[0]
        a.obs["sample"] = gsm
        a.obs["donor"] = dmap.get(gsm, gsm)
        a.obs["file"] = f.name
        ads.append(a)
        print("  读入 %-52s %6d 核 × %5d 基因  供体 %s"
              % (f.name, a.n_obs, a.n_vars, a.obs["donor"].iloc[0]))
    ad = sc.concat(ads, join="outer", index_unique="-", fill_value=0)
    ad.var_names_make_unique()
    print("[合并] %d 个细胞核 × %d 个基因，来自 %d 个样本"
          % (ad.n_obs, ad.n_vars, ad.obs["sample"].nunique()))

    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=10)
    print("[质控] %d 个细胞核 × %d 个基因" % (ad.shape))

    X = ad.X
    Xd = X[:50].toarray() if hasattr(X, "toarray") else np.asarray(X[:50])
    if np.allclose(Xd, np.round(Xd)) and Xd.max() > 30:
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        print("[预处理] 原始计数 → CP10K + log1p")
    return ad


def score(ad, genes):
    g = ad[:, [x for x in genes if x in ad.var_names]].X
    g = g.toarray() if hasattr(g, "toarray") else np.asarray(g)
    z = (g - g.mean(0)) / (g.std(0) + 1e-9)
    return z.mean(1)


def cluster(ad, res=1.0):
    sub = ad.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
    sub = sub[:, sub.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(sub, resolution=res, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False)
    ad.obs["cl"] = sub.obs["cl"].astype(str).values
    print("[聚类] leiden(res=%.2f) 得 %d 个簇" % (res, ad.obs.cl.nunique()))
    return ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=float, default=1.0)
    ap.add_argument("--verify", type=Path, default=None,
                    help="给一份 .h5ad，核对分块实现与原始 delta_auc 完全一致")
    ap.add_argument("--min-neuron", type=int, default=300, help="样本内神经元核达到此数才单独复算")
    a = ap.parse_args()

    if a.verify:
        v = sc.read_h5ad(a.verify)
        tg = [(g, n) for g, n in TARGETS if g in v.var_names]
        o1, _, p1 = delta_auc(v, tg)
        o2, _, p2 = delta_auc_blocks(v, tg)
        m = o1.merge(o2, on="gene", suffixes=("_ref", "_blk"))
        dd = (m.max_dAUC_ref - m.max_dAUC_blk).abs()
        print("[核对] %d 个基因，maxΔAUC 最大差 %.2e，P50 %.4f vs %.4f"
              % (len(m), dd.max(), p1, p2))
        print("[核对] 排名完全一致：%s" % bool((m.percentile_ref.round(6)
                                                == m.percentile_blk.round(6)).all()))

    ad = load_all()
    ad = cluster(ad, a.res)
    ad.obs["neu_score"] = score(ad, NEURON)
    ad.obs["glia_score"] = score(ad, GLIA)

    # 细胞层面的神经元判定：直接数标志基因的检出个数，比 z 分数均值稳。
    # 该数据深度浅（每核数百至三千 UMI），若只在簇层面比 neu/glia 均值，
    # 会把 S100B/MPZ 高表达的胶质簇误判为神经元——第一版就踩了这个坑。
    NPAN = ["NEFL", "NEFH", "PRPH", "STMN2", "ELAVL2", "POU4F1", "SCN9A",
            "SCN10A", "SCN11A", "TUBB3"]
    GPAN = ["S100B", "MPZ", "KCNJ10", "SOX10", "CD74", "COL1A1", "DCN"]
    def hits(genes):
        g = [x for x in genes if x in ad.var_names]
        M = ad[:, g].X
        M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
        return (M > 0).sum(1)
    nh, gh = hits(NPAN), hits(GPAN)
    ad.obs["neu_hits"] = nh
    ad.obs["glia_hits"] = gh
    mask = (nh >= 2) & (nh > gh)
    ad.obs["is_neu"] = mask
    print("[神经元筛] 标志规则：神经元标志检出≥2 个且多于胶质标志；保留 %d/%d 个细胞核"
          % (mask.sum(), ad.n_obs))

    prof = ad.obs.groupby("cl", observed=True).agg(
        n=("cl", "size"), n_neu=("is_neu", "sum"),
        neu=("neu_score", "mean"), glia=("glia_score", "mean"),
        n_sample=("sample", "nunique"))
    prof["神经元比例"] = prof.n_neu / prof.n
    prof.round(3).to_csv(ROOT / "gse201586_clusters.csv", encoding="utf-8-sig")
    print("\n[簇概况]\n%s" % prof.round(2).to_string())

    tg = [(g, n) for g, n in TARGETS if g in ad.var_names]
    print("\n[靶点] 命中 %d/%d：%s" % (len(tg), len(TARGETS), ", ".join(g for g, _ in tg)))

    out_all, sizes_all, p50_all = delta_auc_blocks(ad, tg)
    out_all.round(4).to_csv(ROOT / "gse201586_all_clusters.csv", index=False, encoding="utf-8-sig")
    print("[A] 全部 %d 簇 · 细胞核 %d · 全基因组 maxΔAUC 中位数 %.3f"
          % (len(sizes_all), ad.n_obs, p50_all))

    adn = ad[mask].copy()
    # 关键一步：在神经元核**内部重新分簇**。若沿用全局聚类的簇标签，
    # 落在胶质主导簇里的神经元核会被按胶质簇的表达谱去比，亚型含义就没了。
    sub = adn.copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat")
    sub = sub[:, sub.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=30, svd_solver="arpack")
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(sub, resolution=a.res, key_added="cl", flavor="igraph",
                 n_iterations=2, directed=False)
    adn.obs["cl"] = sub.obs["cl"].astype(str).values
    print("[神经元内分簇] %d 个簇" % adn.obs.cl.nunique())
    # 少于 200 个核的簇不参与：小簇的全基因组零分布不稳，会虚抬 ΔAUC
    keep_cl = [c for c, n in adn.obs["cl"].value_counts().items() if n >= 200]
    # 第二道：簇内泛神经元标志 Nav1.7（SCN9A）检出率须≥30%，否则说明仍有胶质残留
    s9 = adn[:, ["SCN9A"]].X
    s9 = s9.toarray().ravel() if hasattr(s9, "toarray") else np.asarray(s9).ravel()
    tmp = pd.Series(s9 > 0, index=adn.obs_names).groupby(adn.obs["cl"].values).mean()
    weak = [c for c in keep_cl if float(tmp.get(c, 0)) < 0.30]
    keep_cl = [c for c in keep_cl if c not in weak]
    print("[簇质控] 剔除神经元信号不足的簇：%s" % (", ".join(weak) if weak else "无"))
    adn = adn[adn.obs["cl"].isin(keep_cl)].copy()
    out_neu, sizes_neu, p50_neu = delta_auc_blocks(adn, tg)
    out_neu.round(4).to_csv(ROOT / "gse201586_neuron_clusters.csv", index=False, encoding="utf-8-sig")
    print("[B] 神经元簇 %d 个（≥200 核）· 细胞核 %d · 全基因组 maxΔAUC 中位数 %.3f"
          % (len(keep_cl), adn.n_obs, p50_neu))
    print("[B] 神经元簇大小：%s" % sizes_neu)

    # 与发现队列、与第二套外部数据三方对照
    ext2 = {}
    p2 = EXT / "ext_compare_neuron.csv"
    if p2.exists():
        d2 = pd.read_csv(p2)
        ext2 = dict(zip(d2["基因"], d2["外部百分位"]))

    rows = []
    for g, name in tg:
        r = out_neu[out_neu.gene == g].iloc[0]
        dd, dp = DISCOVERY.get(g, (np.nan, np.nan))
        rows.append({"基因": g, "靶点": name,
                     "外部maxΔAUC": r.max_dAUC, "外部百分位": r.percentile,
                     "外部最高簇": r.best_cluster, "外部检出%": r.detect_pct,
                     "发现百分位": dp, "GSE310000百分位": ext2.get(g, np.nan)})
    df = pd.DataFrame(rows).sort_values("发现百分位", ascending=False)
    df.round(3).to_csv(ROOT / "gse201586_compare.csv", index=False, encoding="utf-8-sig")

    ok = df.dropna(subset=["发现百分位"])
    rho, p = spearmanr(ok["外部百分位"], ok["发现百分位"])
    print("\n===== 与发现队列对照（全部 %d 个靶点）ρ=%.3f, P=%.4f =====" % (len(ok), rho, p))
    print(df.round(2).to_string(index=False))

    hi = ok[ok["发现百分位"] >= 95]
    lo = ok[ok["发现百分位"] <= 60]
    print("\n发现判为高选择性 %d 项 → 外部均值 %.1f（%s）"
          % (len(hi), hi["外部百分位"].mean(),
             "～".join("%.1f" % x for x in sorted(hi["外部百分位"]))))
    print("发现判为无选择性 %d 项 → 外部均值 %.1f（%s）"
          % (len(lo), lo["外部百分位"].mean(),
             "～".join("%.1f" % x for x in sorted(lo["外部百分位"]))))

    # 供体层面复算：只有这一层才是"独立个体"意义上的重复
    print("\n[供体复算] 每名供体内独立重算百分位（仅神经元核，≥%d 个者）" % a.min_neuron)
    rows = []
    for s, sub in adn.obs.groupby("donor", observed=True):
        if len(sub) < a.min_neuron:
            continue
        ads = adn[adn.obs["donor"] == s].copy()
        o, _, _ = delta_auc_blocks(ads, tg)
        rec = {"donor": s, "n": len(sub), "libraries": sub["sample"].nunique()}
        for g, _ in tg:
            r = o[o.gene == g]
            rec[g] = float(r.percentile.iloc[0]) if len(r) else np.nan
        rows.append(rec)
    dn = pd.DataFrame(rows)
    cols = [g for g in ["SCN11A", "SCN10A", "TRPV1", "P2RX3", "OPRM1", "OPRD1",
                        "OPRK1", "OPRL1", "CNR1", "ADORA1"] if g in dn]
    if not dn.empty:
        dn.round(1).to_csv(ROOT / "gse201586_per_donor.csv", index=False, encoding="utf-8-sig")
        print(dn.set_index("donor")[["n", "libraries"] + cols].round(1).to_string())
        rng = (dn[cols].max() - dn[cols].min()).round(1)
        sd = dn[cols].std().round(1)
        print("\n供体间极差 / 标准差：")
        print(pd.DataFrame({"极差": rng, "标准差": sd,
                            "最低": dn[cols].min().round(1),
                            "最高": dn[cols].max().round(1)}).to_string())

    # 文库层面复算：同一供体的两个文库可提示技术重复的一致性
    print("\n[文库复算] 每个文库内独立重算（仅神经元核，≥%d 个者）" % a.min_neuron)
    rows = []
    for s, sub in adn.obs.groupby("sample", observed=True):
        if len(sub) < a.min_neuron:
            continue
        ads = adn[adn.obs["sample"] == s].copy()
        o, _, _ = delta_auc_blocks(ads, tg)
        rec = {"sample": s, "donor": sub["donor"].iloc[0], "n": len(sub)}
        for g, _ in tg:
            r = o[o.gene == g]
            rec[g] = float(r.percentile.iloc[0]) if len(r) else np.nan
        rows.append(rec)
    dl = pd.DataFrame(rows)
    if not dl.empty:
        dl.round(1).to_csv(ROOT / "gse201586_per_library.csv", index=False, encoding="utf-8-sig")
        print(dl.set_index("sample")[["donor", "n"] + cols].round(1).to_string())
        # 同供体两文库之间的差（技术重复的一致性）
        gaps = []
        for dnr, g2 in dl.groupby("donor"):
            if len(g2) >= 2:
                gaps.append((g2[cols].max() - g2[cols].min()).rename(dnr))
        if gaps:
            gp = pd.DataFrame(gaps).round(1)
            print("\n同供体内两文库之间的百分位差（技术重复一致性）：")
            print(gp.to_string())

    # 与 GSE310000 的一致性
    if ext2:
        two = df.dropna(subset=["GSE310000百分位"])
        r2, pv2 = spearmanr(two["外部百分位"], two["GSE310000百分位"])
        print("\n[两套外部数据之间] 百分位 ρ=%.3f, P=%.4f（n=%d）" % (r2, pv2, len(two)))

    # ---- 供体配对统计：这才是本队列最有信息量的部分 ----
    if not dn.empty:
        from scipy.stats import wilcoxon
        rows = []
        contrasts = [("SCN11A", "ADORA1"), ("SCN10A", "ADORA1"), ("OPRM1", "ADORA1"),
                     ("TRPV1", "ADORA1"), ("SCN11A", "CNR1"), ("SCN11A", "OPRK1"),
                     ("SCN11A", "TRPV1")]
        for a_, b_ in contrasts:
            if a_ not in dn or b_ not in dn:
                continue
            diff = (dn[a_] - dn[b_]).dropna()
            rows.append({"对比": "%s − %s" % (a_, b_), "为正供体": "%d/%d" % ((diff > 0).sum(), len(diff)),
                         "中位差": round(diff.median(), 1),
                         "最小": round(diff.min(), 1), "最大": round(diff.max(), 1),
                         "Wilcoxon P": round(float(wilcoxon(diff).pvalue), 4)})
        hi = [g for g in ["SCN11A", "SCN10A", "OPRM1", "TRPV1"] if g in dn]
        lo = [g for g in ["ADORA1", "CNR1", "OPRK1"] if g in dn]
        if hi and lo:
            dm = (dn[hi].mean(axis=1) - dn[lo].mean(axis=1))
            rows.append({"对比": "高组均值 − 低组均值", "为正供体": "%d/%d" % ((dm > 0).sum(), len(dm)),
                         "中位差": round(dm.median(), 1), "最小": round(dm.min(), 1),
                         "最大": round(dm.max(), 1),
                         "Wilcoxon P": round(float(wilcoxon(dm).pvalue), 4)})
        pt = pd.DataFrame(rows)
        pt.to_csv(ROOT / "gse201586_paired_tests.csv", index=False, encoding="utf-8-sig")
        print("\n[供体配对] 正差 = 前者更高（n=%d 名供体）" % len(dn))
        print(pt.to_string(index=False))

    # ---- 簇标志核对：给表 7 的"最高簇"一个可读的身份 ----
    panel = ["SCN11A", "SCN10A", "SCN9A", "TRPV1", "TRPM8", "CALCA", "CALCB",
             "MRGPRX1", "MRGPRD", "P2RX3", "NTRK2", "NTRK3", "IL31RA", "SST",
             "KCNIP4", "S100B", "MPZ", "KCNJ10", "CD74", "COL1A1"]
    panel = [g for g in panel if g in adn.var_names]
    X = adn[:, panel].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    mk = pd.DataFrame(X > 0, columns=panel)
    mk["cl"] = adn.obs["cl"].to_numpy()
    tab = mk.groupby("cl", observed=True).mean() * 100
    tab["n"] = adn.obs["cl"].value_counts().reindex(tab.index).values
    tab.round(1).to_csv(ROOT / "gse201586_neuron_markers.csv", encoding="utf-8-sig")
    print("\n[神经元簇标志检出率 %%]\n%s" % tab.round(0).astype(int).to_string())

    ad.write(ROOT / "gse201586_annotated.h5ad")
    print("\n[saved] gse201586_annotated.h5ad")


if __name__ == "__main__":
    main()
