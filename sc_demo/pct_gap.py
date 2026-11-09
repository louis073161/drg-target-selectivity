# -*- coding: utf-8 -*-
"""表 4 增列「本侧百分位差」+ 图 4C 改百分位对百分位。

定义
----
本侧百分位差 = P(该侧首位类 ΔAUC) − P(该侧次位类 ΔAUC)，
两侧的 ΔAUC 各自通过本物种全基因组 maxΔAUC 经验分布映射为百分位（同一个单调标尺）。
它衡量「首位类型」赢得有多干脆：差值小说明首位与次位几乎并列，
此时「最高类型」这个判读本身不牢靠，跨物种「类型相同/不同」也随之失去意义。
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd

SHARED = ["NP", "PEP", "C-LTMR", "Aδ", "Propr"]
TARG = ["SCN10A", "SCN11A", "P2RX3", "CALCA", "OPRM1", "TRPV1",
        "ADORA1", "CNR1", "MRGPRX1", "OPRD1", "OPRK1", "TRPM8"]

D = pd.read_csv("cross_species_delta_auc.csv")
HG = pd.read_csv("human_collapsed_genomewide.csv", index_col=0)["max_dAUC"].values
MG = pd.read_csv("mouse_ref_genomewide_res1.2.csv", index_col=0)["max_dAUC"].values


def pct(v, ref):
    """v 在 ref 经验分布中的百分位（严格小于的比例）。"""
    return float((ref < v).mean() * 100)


rows = []
for sp, ref in (("human", HG), ("mouse", MG)):
    for g in TARG:
        s = D[(D.species == sp) & (D.gene == g) & (D.cls.isin(SHARED))]
        if s.empty:
            continue
        s = s.sort_values("dAUC", ascending=False).reset_index(drop=True)
        r1, r2 = s.iloc[0], s.iloc[1]
        rows.append(dict(species=sp, gene=g,
                         top_cls=r1.cls, top_dAUC=round(r1.dAUC, 3), top_pct=round(pct(r1.dAUC, ref), 1),
                         snd_cls=r2.cls, snd_dAUC=round(r2.dAUC, 3), snd_pct=round(pct(r2.dAUC, ref), 1),
                         gap=round(pct(r1.dAUC, ref) - pct(r2.dAUC, ref), 1)))
G = pd.DataFrame(rows)
G.to_csv("table4_pct_gap.csv", index=False, encoding="utf-8-sig")
print(G.pivot_table(index="gene", columns="species",
                    values=["top_cls", "top_pct", "snd_cls", "snd_pct", "gap"], aggfunc="first")
      .reindex(TARG).to_string())
print("\n===== 表 4 用的配对格式 =====")
for g in TARG:
    h = G[(G.species == "human") & (G.gene == g)]
    m = G[(G.species == "mouse") & (G.gene == g)]
    if h.empty or m.empty:
        continue
    print("%-9s 人类 %-7s P%-5.1f 次 %-7s P%-5.1f 差 %5.1f || 小鼠 %-7s P%-5.1f 次 %-7s P%-5.1f 差 %5.1f"
          % (g, h.top_cls.iloc[0], h.top_pct.iloc[0], h.snd_cls.iloc[0], h.snd_pct.iloc[0],
             h.gap.iloc[0], m.top_cls.iloc[0], m.top_pct.iloc[0], m.snd_cls.iloc[0],
             m.snd_pct.iloc[0], m.gap.iloc[0]))
