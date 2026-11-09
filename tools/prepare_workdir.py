# -*- coding: utf-8 -*-
"""把 data/ 下的派生表平铺到 sc_demo/，使分析脚本无需改动即可按原顺序运行。

脚本沿用写作时的工作目录约定：以 sc_demo/ 为当前目录，读取同目录下的
.csv 中间件，并向上取 ext/gwas、ext/prot。本仓库把中间件按主题收在
data/01…12 中，运行前用本脚本铺平即可。

用法（仓库根目录下）：
    python tools/prepare_workdir.py
"""
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DEST = ROOT / "sc_demo"

def main() -> int:
    if not DATA.is_dir():
        print("找不到 data/ 目录"); return 1

    copied, seen = 0, {}
    for f in sorted(DATA.rglob("*.csv")):
        tgt = DEST / f.name
        if f.name in seen:
            print(f"同名冲突，跳过：{f}（已由 {seen[f.name]} 提供）")
            continue
        seen[f.name] = f.relative_to(DATA)
        shutil.copy2(f, tgt)
        copied += 1

    print(f"已铺平 {copied} 个 CSV 到 sc_demo/")
    print("需要联网获取、未纳入仓库的输入：")
    print("  01/GSE249746_Expression_matrix_raw_counts.csv.gz  （GEO 下载，11 MB，"
          "atlas_s1_build.py 的起点）")
    print("  02/GSE59739_*                                    （GEO 下载，小鼠）")
    print("  08/genelvl_*.csv                                 （FinnGen R12，"
          "由 gene_level.py 从 sumstats 重算）")
    print("  11/ 大体积原始件                                  （见 data/00_说明.txt 的 URL）")
    return 0

if __name__ == "__main__":
    sys.exit(main())
