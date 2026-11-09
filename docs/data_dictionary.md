# Data dictionary

Every number in the manuscript resolves to one file under `data/`. Folder numbering
follows the archived dataset (01–12). Folder 10 (scripts) is `sc_demo/` here, 11 (large raw
downloads) is not redistributed, and 09 (figure files) is `figures/` at the repository root.

## 1. Table and figure → file

| Item in manuscript | File |
|---|---|
| Table 1 — 14 clusters and marker genes | `data/01_…/cluster_annotation_final.csv`, `cluster_marker_profile.csv` |
| Table 3 — 13 targets in the main distribution | `data/01_…/atlas_target_by_class.csv` (numerator), `genomewide_delta_auc.csv` (denominator, 16 969 genes, median 0.112) |
| Table 3, detection rate and per-class ΔAUC | `data/01_…/target_detection_by_celltype.csv`, `atlas_delta_auc_by_type.csv` |
| Table 3, matched percentile column | `data/06_…/matched_null.csv`, column `匹配百分位` (detection-rate matched) and `匹配百分位_含长度` (detection rate + gene length) |
| Table 4 — cross-species comparison | `data/03_…/human_collapsed_targets.csv`, `table4_pct_gap.csv` |
| Figure 4C — percentile vs percentile | `data/03_…/fig4c_percentiles.csv` (n = 9) |
| Table 5 — per-donor recomputation | `data/01_…/donor_replication.csv` |
| Table 6 — external cohort GSE310000 | `data/04_…/ext_donor_split.csv` and siblings |
| Table 7 — third cohort GSE201586 | `data/05_…/gse201586_per_donor.csv` and siblings |
| Table 8 — three-layer pathway completeness | `data/07_…/pathway_architecture.csv`, `pathway_detection_by_type.csv`, `target_coverage.csv` |
| Table 9 — metric comparison and permutation P | `data/06_…/metric_comparison.csv`, `label_permutation.csv` |
| Table 10 — proteome detection and CV | `data/12_…/ST11-ExpressionMatrix.csv` → `table10_proteome.csv` |
| 3.11 — resolution sensitivity, TOST, ICC | `data/06_…/resolution_sensitivity.csv`, `equivalence_tests.csv`, `icc_by_target.csv` |
| 3.11 — matched null distribution | `data/06_…/matched_null.csv` |
| 3.12 — injection simulation | `data/06_…/power_simulation.csv` |
| 3.13 — ADORA1 detection by class | `data/12_…/adora1_by_class.csv` |
| 3.10 end, 4.6 item 9 — genetic test | `data/08_…/drg_gwas_rho.csv` (ρ, permutation P, covariate P), `drg_gwas_deciles.csv`, `drg_gwas_targets.csv` |
| Figures 1–7 | `figures/*.png`; scripts in `sc_demo/` |

## 2. The three reference distributions

Do not interchange them. Each is internally consistent; mixing them produces numbers that
appear plausible and are wrong.

| # | Distribution | Scope | Median | Used by |
|---|---|---|---|---|
| ① | Main | 13 neuronal types, 16 969 genes | 0.112 (P75 0.178, P95 0.320) | Table 3; two-tier structure in 3.4; equivalence margin δ = 0.178 − 0.112 = 0.066 |
| ② | Human collapsed | 5 shared classes after collapsing 13 → 7, max ΔAUC over 16 969 genes | 0.066 (P90 0.198, P95 0.258) | Table 4; Figure 4C human percentiles |
| ③ | Mouse | 6 classes, 10 796 genes | 0.065 (P95 0.217) | Table 4; Figure 4C mouse percentiles |

②'s median and ①'s equivalence margin are both 0.066 by coincidence. This is stated in
section 3.6 of the manuscript because it is the kind of coincidence a reviewer will notice.

## 3. Definition of ΔAUC

For a gene with expression vector over *n* cells of the target type and *N − n* cells of all
other types:

```
ΔAUC = [ R − n(n+1)/2 ] / [ n(N − n) ] − 0.5        # R = Mann–Whitney U statistic
```

i.e. AUROC of that gene as a classifier of the target type, minus 0.5, so that 0 means
"indistinguishable from a coin flip". The genome-wide percentile is the position of that
value within the same distribution recomputed for every expressed gene.

This ΔAUC is unrelated to the "ΔAUC" of two correlated ROC curves in diagnostic-test
literature. The manuscript says so explicitly (section 2.4).

## 4. Notation that is easy to misread

- **P96.5** — genome-wide percentile, not a p-value.
- **"N4 不同／N2 不同"** (Table 5, note) — that donor's top-enriched type differs from the other
  two donors and from the pooled analysis; it does **not** mean all three donors disagree.
- **"未检出"** (Table 10) — below the detection limit of a whole-ganglion DIA-PASEF run, not
  evidence that the protein is absent. ADORA1 protein in particular is not measurable in that
  dataset; the manuscript does not claim otherwise.
