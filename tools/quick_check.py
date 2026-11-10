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

# 与稿件 2026-11-09 稿所用输出一致的校验和（换行符已归一为 LF）
EXPECT = {
    "table10_proteome.csv":          "7b4dacb392d7549ea8d7925b169d4e96",
    "adora1_by_class.csv":           "7993deec62c1184ca9335de8ef47cf7a",
    "matched_null.csv":              "b5223738ec10524f429708cb7849f75b",
    "fig2_genomewide_benchmark.png": "9c22a0d811b21841b60915a16f6daaae",
}
STEPS = [
    ("prot_check.py",   ["table10_proteome.csv", "adora1_by_class.csv"]),
    ("matched_null.py", ["matched_null.csv"]),
    ("fig_benchmark.py", ["fig2_genomewide_benchmark.png"]),
]


def md5(p: Path) -> str:
    """按内容校验。CSV 先把 CRLF 归一成 LF，二进制文件按原样。

    同一份表格在 Windows 上由 pandas 写出为 \\r\\n、Linux 上为 \\n，字节不同但内容
    完全一致。归一化后校验和跨平台可比，CI 才有意义。PNG 不能这么处理——它的
    二进制里本来就可能含 \\r\\n，改写会破坏文件。
    """
    data = p.read_bytes()
    if p.suffix.lower() in (".csv", ".tsv", ".txt"):
        data = data.replace(b"\r\n", b"\n")
    return hashlib.md5(data).hexdigest()


def first_diff(a: Path, b: Path, n: int = 3) -> list[str]:
    """出错时给出前几行差异，省得下次还要翻日志猜。"""
    la = a.read_text(encoding="utf-8-sig", errors="replace").replace("\r\n", "\n").splitlines()
    lb = b.read_text(encoding="utf-8-sig", errors="replace").replace("\r\n", "\n").splitlines()
    out = []
    for i, (x, y) in enumerate(zip(la, lb), 1):
        if x != y:
            out.append(f"    第 {i} 行\n      期望 {y[:90]}\n      实得 {x[:90]}")
            if len(out) >= n:
                break
    if not out and len(la) != len(lb):
        out.append(f"    行数不同：{len(la)} vs {len(lb)}")
    return out


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
            if not ok:
                ref = next((q for q in (ROOT / "data").rglob(o)), None)
                if ref:
                    for line in first_diff(p, ref):
                        print(line)
            bad += 0 if ok else 1

    print("\n" + ("全部一致，稿件数字可复现。" if bad == 0 else f"{bad} 项异常。"))
    print("power_sim.py（3.12 节注入模拟）需 atlas_s3_final.h5ad，"
          "该文件在 Zenodo 归档中，不在本仓库。")
    return 0 if bad == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
