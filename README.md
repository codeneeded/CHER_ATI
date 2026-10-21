# CHER ATI — PBMC 5′ scRNA-seq + TCR

Single-cell analysis of PBMC from the CHER trial around analytic treatment
interruption (ATI). Three participants sampled before interruption (B4) and
during interruption (Interrupt), plus one healthy control run in technical
duplicate.

Structure and conventions follow `codeneeded/CD8_Longitudinal`; the QC and
integration scripts are adapted from that repo's `CP018 - FACTseq` pipeline,
which used the same chemistry.

## Design

Eight samples across two 10x on-chip-multiplexed (OCM) GEM wells. Each well
carries four samples separated by OCM barcode `OB1`–`OB4`; one 5′ gene
expression library and one matched TCR library per well.

| Participant | B4 | Interrupt | Well |
|---|---|---|---|
| 149496 | OB2 | OB1 | 34349 |
| 150815 | OB3 | OB4 | 34349 |
| 154026 | OB1 | OB2 | 33445 |
| BC control | A = OB3, B = OB4 | — | 33445 |

The design is paired. Every participant has both timepoints in the same
well, so the within-participant B4-vs-Interrupt contrast is not confounded with
well. The two healthy-control replicates are both in well 33445, so
control-vs-ATI comparisons are confounded with well and must not be made
naively — see `Analysis/CHER_ATI_FINDINGS.md`.

With n = 3 pairs, the supportable analyses are cell-type composition shifts
within participant, pseudobulk per cell type with participant as a blocking
factor, and signature scoring. A per-gene differential expression test across
three pairs is not supportable and should not be presented as one.

## Data locations (outside this repo)

| What | Where |
|---|---|
| FASTQ | `/mnt/data/CHER_ATI/00_fastq/` |
| Cell Ranger configs and logs | `/mnt/data/CHER_ATI/01_configs/`, `logs/` |
| Cell Ranger output | `~/Desktop/Projects/cellranger_out/CHER_{34349,33445}_multi_vdj/` |
| Reference | `~/Desktop/Projects/refdata/human/{gex,vdj}/` |
| Object checkpoints | `saved_R_data/` (in-repo, gitignored) |

Alignment used Cell Ranger 10.1.0 `multi` against
`refdata-gex-GRCh38-2024-A` and
`refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0`.

## Run order

Every script sources `Scripts/config.R` first, which holds all paths, the sample
sheet, QC thresholds and shared helpers. Scripts are single-purpose and numbered
in run order so a failure can be isolated to one step.

```
Scripts/1. QC/01_build_objects.R     # H5 -> merged Seurat, QC metrics      -> QS_RAW
Scripts/1. QC/02_qc_plots_pre.R      # pre-filter diagnostics
Scripts/1. QC/03_filter.R            # apply thresholds, retention table    -> QS_FILTERED
Scripts/1. QC/04_cellcycle.R         # normalise, PCA/UMAP, phase scoring   -> QS_CELLCYCLE
Scripts/1. QC/05_doublets.R          # scDblFinder per library, singlets    -> QS_SINGLETS
Scripts/1. QC/06_integrate.R         # CCA + FastMNN, clustree sweep        -> QS_INTEGRATED
```

Run them with the project R environment:

```bash
cd ~/Desktop/Projects/Repositories/CHER_ATI
Rscript "Scripts/1. QC/01_build_objects.R"
```

## Output conventions

- Analysis output goes under `Analysis/`; nothing is written to `Manuscript/`
  (enforced by `assert_not_manuscript()` in `config.R`).
- Every figure ships its numbers — a `*_source_data.csv` beside each panel.
- Object checkpoints are `.qs2` (`qs2::qs_save` / `qs_read`) in `saved_R_data/`,
  which is gitignored: they are regenerable and large.
