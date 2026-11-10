# -*- coding: utf-8 -*-
"""一键复核：重跑三个不依赖大文件的脚本，比对输出与稿件所用版本是否一致。

覆盖稿件中的表 10、3.13 节、3.11 节的匹配零分布与图 2 —— 即审稿人最可能
亲手复算的四处。

用法（仓库根目录下）：
    python tools/prepare_workdir.py
    python tools/quick_check.py
"""
from pathlib import Path
import hashlib
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "sc_demo"

# 与稿件 2026-11-09 稿所用输出逐字节一致的校验和
EXPECT = {
    "table10_proteome.csv":          "1639e0f770b63a6a3c55cd757b68f267",
    "adora1_by_class.csv":           "99173d5f4d50ce20031879701a227c77",
    "matched_null.csv":              "e8b523c10000a4699100faf0503eefb5",
    "fig2_genomewide_benchmark.png": "9c22a0d811b21841b60915a16f6daaae",
}
STEPS = [
    ("prot_check.py",   ["table10_proteome.csv", "adora1_by_class.csv"]),
    ("matched_null.py", ["matched_null.csv"]),
    ("fig_benchmark.py", ["fig2_genomewide_benchmark.png"]),
]


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> int:
    csv_only = "--csv-only" in sys.argv      # CI 用：PNG 字节不跨平台一致
    if not (WORK / "gene_selectivity_ranked.csv").exists():
        print("请先运行 python tools/prepare_workdir.py"); return 1

    bad = 0
    for script, outs in STEPS:
        if csv_only and script == "fig_benchmark.py":
            print("\n>>> fig_benchmark.py  跳过（--csv-only）")
            continue
        print(f"\n>>> {script}")
        r = subprocess.run([sys.executable, script], cwd=WORK,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("    失败：", (r.stderr or "").strip().splitlines()[-1:]); bad += 1; continue
        for o in outs:
            p = WORK / o
            if not p.exists():
                print(f"    缺少输出 {o}"); bad += 1; continue
            got = md5(p)
            ok = got == EXPECT[o]
            print(f"    {'一致' if ok else '不一致'}  {o}  {got[:8]}")
            bad += 0 if ok else 1

    print("\n" + ("全部一致，稿件数字可复现。" if bad == 0 else f"{bad} 项异常。"))
    print("power_sim.py（3.12 节注入模拟）需 atlas_s3_final.h5ad，"
          "该文件在 Zenodo 归档中，不在本仓库。")
    return 0 if bad == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
