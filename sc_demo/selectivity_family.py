# -*- coding: utf-8 -*-
"""
全基因组 ΔAUC 扫描 + 靶点家族选择性对比
目的：检验"GPCR 类镇痛靶点的亚型选择性系统性低于离子通道类"这一可证伪假说
"""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
from scipy.stats import rankdata, mannwhitneyu

A = sc.read_h5ad("atlas_s3_final.h5ad")
A = A[A.obs.cell_class != "非神经元"].copy()          # 只留神经元
print("神经元细胞:", A.n_obs, "| 类型数:", A.obs.cell_type.nunique())

X = A.X.toarray() if sparse.issparse(A.X) else np.asarray(A.X, dtype=np.float32)
X = X.astype(np.float64)
genes = np.array(A.var_names)
det = (X > 0).mean(0)
keep = np.where(det >= 0.05)[0]                        # 至少 5% 细胞检出
X = X[:, keep]; genes = genes[keep]
print("基因数（检出>=5%）:", len(genes), "| 矩阵:", X.shape)

ct = A.obs.cell_type.astype(str).values
types = sorted(pd.unique(ct))
print("细胞类型:", types)

# ---- AUC（含并列校正）：AUC = (R1 - n1(n1+1)/2) / (n1*n2)
X = X.T                                                # 转成 基因 × 细胞
R = rankdata(X, axis=1)                                # 每个基因内跨细胞的中位秩
n_all = X.shape[1]
AUC = np.zeros((len(genes), len(types)), dtype=np.float64)
for j, t in enumerate(types):
    m1 = (ct == t); n1 = m1.sum(); n2 = n_all - n1
    R1 = R[:, m1].sum(1)
    AUC[:, j] = (R1 - n1 * (n1 + 1) / 2) / (n1 * n2)

# ---- ΔAUC：减去该类型"全基因 AUC 中位数"（经验零分布校正）
DAUC = AUC - np.median(AUC, axis=0, keepdims=True)
pd.DataFrame(DAUC, index=genes, columns=types).to_csv("genomewide_delta_auc.csv",
                                                      encoding="utf-8-sig")
print("ΔAUC 矩阵已存: genomewide_delta_auc.csv")

# ---- 每个基因的最大富集（跨类型）与丰度
sel = pd.DataFrame({"max_dAUC": DAUC.max(1), "best_type": np.array(types)[DAUC.argmax(1)],
                    "detect_pct": det[keep] * 100, "mean_expr": X.mean(1)}, index=genes)
sel["abundance_rank"] = sel.mean_expr.rank(pct=True) * 100

# ---- 家族分类（按基因名前缀 + IUPHAR 常用 GPCR 清单）
CHAN_PRE = ("SCN", "KCN", "CACNA", "CACNB", "CACNG", "TRPV", "TRPM", "TRPA", "TRPC",
            "P2RX", "ASIC", "HCN", "CLCN", "ANO", "PIEZO", "CATSPER", "TMEM63",
            "GRIN", "GRIA", "GABRA", "GABRB", "GABRG", "CHRNA", "CHRNB", "CHRND",
            "CHRNE", "CHRNG")
GPCR = """OPRM1 OPRD1 OPRK1 OPRL1 ADORA1 ADORA2A ADORA2B ADORA3 CNR1 CNR2 GPR55
MRGPRX1 MRGPRX2 MRGPRX3 MRGPRX4 MRGPRD MRGPRE MRGPRF MRGPRG HTR1B HTR1D HTR1F
HTR2A HTR2B HTR3A HTR4 HTR7 PTGER1 PTGER2 PTGER3 PTGER4 PTGIR PTGFR BDKRB1 BDKRB2
ADRB1 ADRB2 ADRB3 ADRA1A ADRA1B ADRA2A ADRA2B ADRA2C CHRM1 CHRM2 CHRM3 CHRM4 CHRM5
GABBR1 GABBR2 GRM1 GRM2 GRM3 GRM4 GRM5 GRM7 GRM8 NPY1R NPY2R NPY4R NPY5R SSTR1
SSTR2 SSTR3 SSTR4 SSTR5 GALR1 GALR2 GALR3 NTSR1 NTSR2 CCR2 CCR4 CCR5 CXCR4 CXCR5
CX3CR1 LPAR1 LPAR2 LPAR3 LPAR6 S1PR1 S1PR2 S1PR3 S1PR5 EDNRA EDNRB VIPR1 VIPR2
CRHR1 CRHR2 MC1R MC4R HCRTR1 HCRTR2 OXTR AVPR1A AVPR1B CALCRL CALCR RAMP1 ACKR1
ACKR3 C3AR1 C5AR1 F2R F2RL1 F2RL2 P2RY1 P2RY2 P2RY12 P2RY13 P2RY14 ADCYAP1R1
GIPR GLP1R SCTR NPSR1 PROKR2 QRFPR NPFFR1 NPFFR2 NPBWR1 GPR149 GPR83 GPR151
GPR174 GPR183 HCAR2 HTR5A HTR6 ADGRB1 ADGRG6 ADGRE5""".split()
GPCR = set(GPCR)

is_ch = np.array([g.startswith(CHAN_PRE) for g in genes])
is_gp = np.array([g in GPCR for g in genes])
print("\n家族基因数: 离子通道样 %d | GPCR(清单命中) %d" % (is_ch.sum(), is_gp.sum()))

sub = sel.copy(); sub["family"] = np.where(is_ch, "channel", np.where(is_gp, "GPCR", "other"))
for fam in ("channel", "GPCR"):
    s = sub[sub.family == fam]
    print("\n== %s (n=%d)" % (fam, len(s)))
    print("   max_dAUC: 中位 %.3f | 四分位 %.3f-%.3f | 最大 %.3f"
          % (s.max_dAUC.median(), s.max_dAUC.quantile(.25), s.max_dAUC.quantile(.75),
             s.max_dAUC.max()))
    print("   丰度中位百分位: %.0f" % s.abundance_rank.median())

a = sub[sub.family == "channel"].max_dAUC
b = sub[sub.family == "GPCR"].max_dAUC
u = mannwhitneyu(a, b, alternative="greater")
# 秩双列相关作为效应量
allv = np.r_[a, b]; rk = rankdata(allv)
rr = (rk[:len(a)].mean() - rk[len(a):].mean()) / len(allv)
print("\n>>> 假说检验 channel > GPCR : U=%.0f  p=%.3e  秩双列 r=%.3f" % (u.statistic, u.pvalue, rr))

# 丰度匹配后比较（排除"低表达导致 AUC 压缩"的混杂）
lo, hi = sub.abundance_rank.quantile([.25, .75])
m = sub[(sub.abundance_rank >= lo) & (sub.abundance_rank <= hi)]
a2 = m[m.family == "channel"].max_dAUC; b2 = m[m.family == "GPCR"].max_dAUC
if len(a2) > 2 and len(b2) > 2:
    u2 = mannwhitneyu(a2, b2, alternative="greater")
    print(">>> 丰度匹配(25-75百分位)后: channel n=%d 中位%.3f | GPCR n=%d 中位%.3f | p=%.3e"
          % (len(a2), a2.median(), len(b2), b2.median(), u2.pvalue))

print("\n--- 已知靶点核对（应与前稿一致）---")
chk = ["OPRM1", "OPRD1", "OPRK1", "ADORA1", "CNR1", "SCN10A", "SCN11A", "TRPV1", "P2RX3", "CALCA"]
print(sub.reindex([g for g in chk if g in sub.index])[["max_dAUC", "best_type", "detect_pct"]]
      .round(3).to_string())

sub.sort_values("max_dAUC", ascending=False).to_csv("gene_selectivity_ranked.csv",
                                                    encoding="utf-8-sig")
print("\n完成: gene_selectivity_ranked.csv")
