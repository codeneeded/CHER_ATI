# CHER ATI — QC and integration findings

Covers alignment through integration. Every number here traces to a CSV under
`Analysis/`; where a judgement was made rather than measured, it says so.

Pipeline: `Scripts/1. QC/01` → `07`. Cell Ranger 10.1.0 `multi`,
`refdata-gex-GRCh38-2024-A`, R 4.4.3 / Seurat 5.5.1.

## 1. Cell accounting

| Stage | Cells | Note |
|---|---|---|
| Cell Ranger called | 24,233 | across 8 samples, two OCM wells |
| Loaded into Seurat | 23,835 | after `min.cells=3`, `min.features=200` |
| Passing QC filters | 21,076 (88.4%) | `CHER_QC_filtering_summary.csv` |
| Singlets | 20,050 | `CHER_doublet_summary_source_data.csv` |

Per-sample retention ranged 79.2–98.0%. The lowest were `149496_Interrupt`
(79.2%), `BC_control_A` (81.0%) and `BC_control_B` (82.8%); the highest
`150815_B4` (98.0%).

Genes: 25,921 loaded → 22,376 retained (detected in ≥10 cells).

## 2. Two inherited QC thresholds were changed

Thresholds started from CP018 (`Scripts/CP018 - FACTseq/config.R` in
`CD8_Longitudinal`) so the projects stay comparable. CP018 sorted **CD8 T
cells**; CHER is **whole PBMC**, and two thresholds did not transfer.
Per-criterion cost per sample is in `CHER_threshold_cost_per_sample.csv`.

**Ribosomal floor (`percent_ribo > 5`) — disabled. This is a judgement call.**
It would have removed 3,647 cells (15.3%), but 64.6% of those already fail a
direct quality criterion, so most of its work is redundant. The 1,291 cells
unique to it (5.4% of the dataset) are borderline rather than clearly bad:
median 1,488 genes against 2,362 in retained cells, mito 7.65% against 2.88% —
but complexity 0.899, *above* the retained median of 0.870. A ribosomal minimum
is also not a standard whole-PBMC criterion. Decision: let nFeature, mito and
complexity do the filtering directly rather than have an indirect proxy account
for 5.4% of removals. Re-enabling it is one line in `config.R`.

**Platelet metric — redefined from `PECAM1|PF4` to `PF4|PPBP`.** PECAM1 (CD31)
is broadly expressed across leukocytes, so it does not measure platelet
contamination; PPBP is the dominant platelet transcript and was absent from the
CP018 pattern. Net effect is to flag *more* cells, not fewer (4.54% at ≥2%
versus 1.84%), and those cells are genuinely platelet-dominated — median 356
genes against 2,298 in retained cells.

The haemoglobin pattern was also anchored to the actual globin genes; the
commonly-copied `^HB[^(P)]` idiom additionally matches HBEGF and HBS1L.

## 3. The duplicate healthy control is the noise floor

`BC_control_A` and `BC_control_B` are the same donor through the same pipeline,
so the gap between them bounds what counts as a meaningless difference.

| Metric | A | B | Difference |
|---|---|---|---|
| Cells called (Cell Ranger) | 4,118 | 4,030 | 2.2% |
| Median UMI | 7,025 | 6,842 | 2.6% |
| Median genes | 2,170 | 2,154 | 0.7% |
| QC retention | 81.0% | 82.8% | 1.8 pp |
| Doublet rate | 5.19% | 5.00% | 0.19 pp |
| Cells with paired TCR | 502 | 467 | 7.2% |

Tight on gene expression, looser on the TCR side. **A B4-versus-Interrupt
difference smaller than the corresponding row here is not interpretable.**

## 4. Doublets

`scDblFinder` run **per library, never pooled** — it simulates doublets by
combining cells from the same capture, so pooling would build doublets from
cells that were never in the same GEM well. Rates were 4.37–5.57% with no
outlier (`Analysis/QC/Doublets/`).

Note on OCM: Cell Ranger already removed cross-sample doublets identifiable by
conflicting OCM barcode. What this step targets is *within*-sample doublets,
which carry no barcode conflict and are invisible to the demultiplexer.

## 5. Integration

Batch is the **GEM well**. The wells differ systematically in depth (median UMI
6,842–10,279 in 33445 against 5,176–6,393 in 34349), so this is not a nominal
batch variable.

Both CCA and FastMNN were run on the same 3,000 variable features
(`integrated.cca.rna`, `integrated.mnn.rna`). FastMNN's own lost-variance
diagnostic is small — 0.58% and 0.55% per batch
(`CHER_fastMNN_lost_variance.csv`) — so the correction was gentle.

At resolution 0.8 (an evaluation point, **not** a chosen resolution):

| | Clusters | Median well skew | Clusters >90% one well | Clusters with all 4 donors |
|---|---|---|---|---|
| FastMNN | 21 | 14.3 pp | 4 | 21 of 21 |
| CCA | 20 | 11.8 pp | 3 | 19 of 20 |

**The residual well skew is donor identity, not uncorrected batch.** Across all
21 FastMNN clusters, the correlation between "% from well 33445" and "% healthy
control" is **r = 0.872**. Of the four well-skewed clusters (3,467 cells, 17.3%
of the dataset), two are essentially pure healthy control — cluster 3 is 99.4%
well 33445 and 99.2% control; cluster 9 is 95.9% and 95.7%. Those are a
different donor forming its own populations, which integration should not and
cannot remove.

The one exception is **cluster 14** (334 cells, 98.2% well 34349, 1.5% control)
— a population restricted to participants 149496 and 150815. Whether that is
donor biology or residual batch cannot be settled from these data, because
those two participants are the only occupants of that well. Flagged rather than
resolved.

## 6. Design constraints that bound the analysis

**The paired contrast is protected.** Every participant has both timepoints in
the *same* well, so B4-versus-Interrupt within a participant carries no well
effect. This is the comparison the design supports.

**Control-versus-ATI is confounded three ways** — with well, with donor, and
with the fact that both control replicates are the same person. Section 5 shows
the confound is real and large. The healthy control is a QC reference and a
noise floor; it is not a comparison arm.

**n = 3 pairs.** Supportable: cell-type composition shifts within participant,
pseudobulk per cell type with participant as a blocking factor, signature
scoring. Not supportable: a per-gene differential expression test across three
pairs. A volcano plot here would imply power the design does not have.

**Clonotype diversity must be rarefied before comparison.** Cell Ranger's
per-sample paired-clonotype diversity ranges 69.7–999.8, but the paired-cell
counts behind those values range 445–1,301. Diversity is depth-dependent;
comparing the raw values measures library size.

## 7. Open items

1. **Resolution is deliberately unset.** `MNN_RES_FINAL` and `CCA_RES_FINAL`
   remain `NA` in `config.R`. Choose from `CHER_clustree_MNN.png` /
   `_CCA.png` at annotation, not from this document.
2. **CCA or FastMNN?** CCA mixes marginally better (11.8 vs 14.3 pp median
   skew); FastMNN keeps all four donors in every cluster. Not yet decided.
3. **Cell cycle phase needs a sanity check.** ATI samples score ~40% S phase
   against ~20% in the controls. `CellCycleScoring` assigns phase by *relative*
   score, so in a mostly-resting PBMC population this can over-call S. Check
   against MKI67/TOP2A before reading it as proliferation.
4. **TCR is aligned but not yet attached** to the Seurat object.
5. **Azimuth is not installed**, so there are no reference labels yet.
   Annotation will be marker-based regardless.
