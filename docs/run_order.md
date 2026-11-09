# Run order and inputs

All scripts assume the working directory is `sc_demo/`, exactly as during the analysis.
Run `python tools/prepare_workdir.py` first — it copies every CSV from `data/**` into
`sc_demo/` so that the relative paths inside the scripts resolve.

## Tier 1 — runs immediately from this repository (no external download)

| Script | Produces | Verified |
|---|---|---|
| `prot_check.py` | `table10_proteome.csv`, `adora1_by_class.csv` (Table 10, section 3.13) | ✅ byte-identical |
| `matched_null.py` | `matched_null.csv` (section 3.11 matched null) | ✅ byte-identical |
| `fig_benchmark.py` | `fig2_genomewide_benchmark.png` (Figure 2) | ✅ byte-identical |

`python tools/quick_check.py` runs exactly these three and compares MD5s.

## Tier 2 — needs the Zenodo archive (`archive/`)

Unzip the Zenodo deposit into `archive/` at the repository root and re-run
`prepare_workdir.py`; the `.h5ad` objects are then linked into `sc_demo/`.

```
01_人类发现队列_GSE249746 → atlas_s1_build.py → atlas_s2_targets.py
                          → atlas_s3_stats.py  → atlas_s4_figs.py
```

| Script | Needs | Notes |
|---|---|---|
| `atlas_s1_build.py` | `GSE249746_Expression_matrix_raw_counts.csv.gz` | GEO download, 11 MB; the only true entry point |
| `atlas_s2_targets.py`, `atlas_s3_stats.py` | `atlas_s1_clustered.h5ad` | |
| `atlas_s4_figs.py`, `donor_replication.py`, `fig7_architecture.py`, `robustness.py`, `selectivity_family.py`, `power_sim.py`, `equiv_architecture.py` | `atlas_s3_final.h5ad` (216 MB) | |
| `mdrg_s1.py` → `mdrg_s2.py` → `mdrg_s3.py` | mouse FPKM matrix | mouse cohort |
| `mouse_ref_1_recluster.py` → `mouse_ref_2_dist.py` → `rebuild_cross.py` → `pct_gap.py` | `mdrg_reclustered_59739.h5ad` | cross-species; Table 4, Figure 4C |
| `prep_gse310000.py` → `replicate_external.py`, `ext_report.py` | `GSE310000.h5ad` | second human cohort |
| `run_201586.py` | third cohort h5ad | 7 donors, 146 260 nuclei |

## Tier 3 — needs public summary statistics

| Script | Input |
|---|---|
| `gene_level.py` | FinnGen R12 endpoints (`*.gz`), one run per endpoint |
| `drg_gwas_link.py` | `ext/gwas/gene_coords_grch38.tsv` + the six per-endpoint tables |
| `build_gene_table.py` | gene coordinates → `gene_len.pkl` |

The six endpoints used are `M13_ALGONEURO`, `M13_FIBROMYALGIA`, `M13_LOWBACKPAIN`,
`G6_MIGRAINE`, `F5_DEPRESSIO` (depression) and `E4_OBESITY` (the last two as positive
controls for the method, not as pain phenotypes).

## Manuscript build

`build_article.py` regenerates the .docx; it needs `figures/*.png` and a Word installation
for the PDF conversion. It is included for transparency, not required for data review.

## Determinism

Seed 0 throughout. Permutation seeds use `zlib.crc32(tag)`; an earlier version used
`abs(hash(tag))`, which is salted per process and therefore not reproducible. If you
find any remaining `hash(` in a seed path, that is a bug — please open an issue.
