# Subtype selectivity of analgesic targets in human dorsal root ganglia

[![quick check](https://github.com/louis073161/drg-target-selectivity/actions/workflows/quick_check.yml/badge.svg)](https://github.com/louis073161/drg-target-selectivity/actions/workflows/quick_check.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24109653.svg)](https://doi.org/10.5281/zenodo.24109653)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Code, derived data and figure source files for the manuscript:

> **The adenosine A1 receptor lacks reproducible subtype selectivity in human dorsal root ganglia:
> a genome-wide benchmark across four cohorts**

The analysis asks a narrow, testable question. Targets repeatedly invoked to explain
acupuncture analgesia — the adenosine A1 receptor (ADORA1), cannabinoid CB1 (CNR1), the
opioid receptors, P2X3, TRPV1, and the voltage-gated sodium channels Nav1.7/1.8/1.9 — are
placed on one common yardstick: the distribution of ΔAUC over all expressed genes in the
same neurons. The question is not whether a target is "expressed in nociceptors" but whether
it is *more* restricted to one neuronal subtype than an arbitrary gene of similar abundance.

Answer, in one line: sodium channels, TRPV1, P2X3 and the μ/δ opioid receptors are subtype-restricted;
ADORA1, CNR1 and the κ opioid receptor sit at the genome-wide median. The split does not follow
molecular family and does not survive as a cross-family rule — see `docs/` for the caveats that
the manuscript states explicitly.

---

## What is in this repository

| Path | Contents |
|---|---|
| `sc_demo/` | 31 analysis scripts, in run order (discovery atlas → mouse → cross-species → two external cohorts → robustness/equivalence → pathway architecture → genetics → protein check → figures) |
| `data/01…12/` | Derived tables, grouped exactly as in the archived dataset: every number in every table and figure can be traced to one CSV here |
| `figures/` | The seven published figures, 400 dpi PNG |
| `ext/gwas/` | Gene-level GWAS scripts and GRCh38 gene coordinates |
| `ext/prot/` | Public human DRG proteome matrix (PRIDE PXD061129), used for Table 10 |
| `tools/` | `prepare_workdir.py` (flatten derived tables into the working directory), `quick_check.py` (re-run three scripts and verify checksums) |
| `docs/` | Data dictionary, run order, known limits, and the release guide (GitHub + Zenodo) |

Raw sequencing data are **not** here. They are public and are cited by accession:
GSE249746 (discovery, deep single-soma), GSE59739 (mouse), GSE310000 (2 donors),
GSE201586 (7 donors, 146 260 nuclei), FinnGen R12 (six endpoints), PRIDE PXD061129 (proteome).
The 3.1 GB archive of intermediates, including the `.h5ad` objects, is deposited on Zenodo
(DOI in `CITATION.cff`).

## Quick start

```bash
git clone https://github.com/louis073161/drg-target-selectivity.git
cd drg-target-selectivity
python -m pip install -r requirements.txt      # Python 3.13
python tools/prepare_workdir.py                # flatten data/**/*.csv into sc_demo/
python tools/quick_check.py                    # reproduce Table 10, 3.13, 3.11, Figure 2
```

`quick_check.py` re-runs `prot_check.py`, `matched_null.py` and `fig_benchmark.py` and compares
MD5 checksums against the versions used in the manuscript. It should print
`全部一致，稿件数字可复现。` — four outputs, byte-identical.

The three table-producing checks also run in continuous integration on Linux
(`--csv-only`, since PNG bytes are not platform-invariant); see the badge above and
`.github/workflows/quick_check.yml`. A green badge means the published numbers are
regenerable from the published tables by a third party, not merely by the authors.

Scripts that need the intermediate `.h5ad` objects (`atlas_s3_final.h5ad`, 216 MB;
`GSE310000.h5ad`; the mouse objects) require the Zenodo archive; unzip it into `archive/`
and re-run `prepare_workdir.py`. See `docs/run_order.md`.

## Three yardsticks, not one

The single most common way to misread this paper is to mix up its reference distributions.
`docs/data_dictionary.md` states which analysis uses which:

- **Main distribution** — 13 neuronal types, 16 969 genes, median 0.112. Used for Table 3,
  the two-tier structure in 3.4, and the equivalence margin δ = 0.178 − 0.112 = 0.066.
- **Human 7-class collapsed** — median 0.066. Used for Table 4 and Figure 4C.
- **Mouse 6-class** — 10 796 genes, median 0.065. Used for Table 4 and Figure 4C.

Two of these medians happen to equal 0.066. That is a coincidence, not a shared scale.

## Reproducibility notes

- Random seed fixed at 0 throughout; permutation seeds are `zlib.crc32` based, so they are
  reproducible across processes and Python builds (`abs(hash(...))` is not — it is salted per
  process, and an earlier version of the code was wrong on this point).
- `matched_null.py` matches each target to 800 nearest neighbours on detection rate and gene
  length; three matching recipes are reported and their rank correlation exceeds 0.97.
- `power_sim.py` reports what the negative result can and cannot exclude: enrichment of
  ≥3× is ruled out, ~2× is not.

## Citation

If you use this code or these tables, please cite the manuscript and the archived release:

> Lu Y. *Subtype selectivity of analgesic targets in human dorsal root ganglia* (v1.0.1).
> Zenodo. https://doi.org/10.5281/zenodo.24109653

`10.5281/zenodo.24109653` is the **concept DOI** — it always resolves to the newest archived
version. The DOI of the v1.0.0 snapshot is `10.5281/zenodo.24122624`. `CITATION.cff` carries
the metadata GitHub and Zenodo read automatically.

**Author.** Yi Lu, MD, PhD — The Affiliated Traditional Chinese Medicine Hospital,
Guangzhou Medical University. ORCID [0000-0002-0758-0648](https://orcid.org/0000-0002-0758-0648).

## License

Code: MIT (`LICENSE`). Derived tables and figures: CC BY 4.0.
Public data accessed through GEO, FinnGen and PRIDE remain under their own terms.

---

## 中文说明（给合作者）

本仓库是投稿用的代码与派生数据包，不含原始测序数据（原始数据按 GEO 编号公开，
3.1 GB 中间件另行存放 Zenodo，取得 DOI 后回填 `CITATION.cff` 与稿件"数据可用性"一节）。

三步复现：

```bash
python -m pip install -r requirements.txt
python tools/prepare_workdir.py
python tools/quick_check.py
```

`quick_check.py` 会把表 10、3.13 节、3.11 节匹配零分布与图 2 重算一遍，并与稿件所用
版本比对 MD5（CSV 先归一换行符，跨平台可比），全部一致才算复现成功。需要 `.h5ad` 中间件的脚本（含 3.12 节注入模拟）
要先解压 Zenodo 归档，见 `docs/run_order.md`。

发布步骤（GitHub 建库 → 打 tag → Zenodo 领 DOI → 回填）见 `docs/release_guide.md`。
