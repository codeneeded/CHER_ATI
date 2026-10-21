# CHER ATI — QC and integration findings

Covers alignment through integration. Every number here traces to a CSV under
`Analysis/`; where a judgement was made rather than measured, it says so.

Pipeline: `Scripts/1. QC/01` → `07`. Cell ranger 10.1.0 `multi`,
`refdata-gex-GRCh38-2024-A`, R 4.4.3 / Seurat 5.5.1.

## 1. Cell accounting

| Stage | Cells | Note |
|---|---|---|
| Cell ranger called | 24,233 | across 8 samples, two OCM wells |
| Loaded into seurat | 23,835 | after `min.cells=3`, `min.features=200` |
| Passing QC filters | 21,076 (88.4%) | `CHER_QC_filtering_summary.csv` |
| Singlets | 20,050 | `CHER_doublet_summary_source_data.csv` |

Per-sample retention ranged 79.2–98.0%. The lowest were `149496_Interrupt`
(79.2%), `BC_control_A` (81.0%) and `BC_control_B` (82.8%); the highest
`150815_B4` (98.0%).

Genes: 25,921 loaded → 22,376 retained (detected in ≥10 cells).

## 2. Two inherited QC thresholds were changed

Thresholds started from CP018 (`Scripts/CP018 - FACTseq/config.R` in
`CD8_Longitudinal`) so the projects stay comparable. CP018 sorted CD8 T
cells; CHER is whole PBMC, and two thresholds did not transfer.
Per-criterion cost per sample is in `CHER_threshold_cost_per_sample.csv`.

Ribosomal floor (`percent_ribo > 5`) — disabled. This is a judgement call.
It would have removed 3,647 cells (15.3%), but 64.6% of those already fail a
direct quality criterion, so most of its work is redundant. The 1,291 cells
unique to it (5.4% of the dataset) are borderline rather than clearly bad:
median 1,488 genes against 2,362 in retained cells, mito 7.65% against 2.88% —
but complexity 0.899, above the retained median of 0.870. A ribosomal minimum
is also not a standard whole-PBMC criterion. Decision: let nFeature, mito and
complexity do the filtering directly rather than have an indirect proxy account
for 5.4% of removals. Re-enabling it is one line in `config.R`.

Platelet metric — redefined from `PECAM1|PF4` to `PF4|PPBP`. PECAM1 (CD31)
is broadly expressed across leukocytes, so it does not measure platelet
contamination; PPBP is the dominant platelet transcript and was absent from the
CP018 pattern. Net effect is to flag more cells, not fewer (4.54% at ≥2%
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

Tight on gene expression, looser on the TCR side. A B4-versus-Interrupt
difference smaller than the corresponding row here is not interpretable.

## 4. Doublets

`scDblFinder` run per library, never pooled — it simulates doublets by
combining cells from the same capture, so pooling would build doublets from
cells that were never in the same GEM well. Rates were 4.37–5.57% with no
outlier (`Analysis/QC/Doublets/`).

Note on OCM: Cell ranger already removed cross-sample doublets identifiable by
conflicting OCM barcode. What this step targets is within-sample doublets,
which carry no barcode conflict and are invisible to the demultiplexer.

## 5. Integration

Batch is the GEM well. The wells differ systematically in depth (median UMI
6,842–10,279 in 33445 against 5,176–6,393 in 34349), so this is not a nominal
batch variable.

Both CCA and FastMNN were run on the same 3,000 variable features
(`integrated.cca.rna`, `integrated.mnn.rna`). FastMNN's own lost-variance
diagnostic is small — 0.58% and 0.55% per batch
(`CHER_fastMNN_lost_variance.csv`) — so the correction was gentle.

At resolution 0.8 (an evaluation point, not a chosen resolution):

| | Clusters | Median well skew | Clusters >90% one well | Clusters with all 4 donors |
|---|---|---|---|---|
| FastMNN | 21 | 14.3 pp | 4 | 21 of 21 |
| CCA | 20 | 11.8 pp | 3 | 19 of 20 |

The residual well skew is donor identity, not uncorrected batch. Across all
21 FastMNN clusters, the correlation between "% from well 33445" and "% healthy
control" is r = 0.872. Of the four well-skewed clusters (3,467 cells, 17.3%
of the dataset), two are essentially pure healthy control — cluster 3 is 99.4%
well 33445 and 99.2% control; cluster 9 is 95.9% and 95.7%. Those are a
different donor forming its own populations, which integration should not and
cannot remove.

The one exception is cluster 14 (334 cells, 98.2% well 34349, 1.5% control)
— a population restricted to participants 149496 and 150815. Whether that is
donor biology or residual batch cannot be settled from these data, because
those two participants are the only occupants of that well. Flagged rather than
resolved.

## 6. Design constraints that bound the analysis

The paired contrast is protected. Every participant has both timepoints in
the same well, so B4-versus-Interrupt within a participant carries no well
effect. This is the comparison the design supports.

Control-versus-ATI is confounded three ways — with well, with donor, and
with the fact that both control replicates are the same person. Section 5 shows
the confound is real and large. The healthy control is a QC reference and a
noise floor; it is not a comparison arm.

n = 3 pairs. Supportable: cell-type composition shifts within participant,
pseudobulk per cell type with participant as a blocking factor, signature
scoring. Not supportable: a per-gene differential expression test across three
pairs. A volcano plot here would imply power the design does not have.

Clonotype diversity must be rarefied before comparison. Cell Ranger's
per-sample paired-clonotype diversity ranges 69.7–999.8, but the paired-cell
counts behind those values range 445–1,301. Diversity is depth-dependent;
comparing the raw values measures library size.

## 7. Annotation

FastMNN was taken forward over CCA: every one of its clusters contains all four
donors, which matters more here than the 2.5 pp mixing advantage CCA showed.

Resolution 0.6 (16 clusters) was chosen on nameability — a judgement, not a
measurement. The marker-redundancy criterion carried over from miR29a did not
discriminate: it found zero redundant marker pairs anywhere from 0.2 to 1.0 and
would therefore have selected 1.0. Disclosed rather than presented as derived.

14 labels survive (`ANNO_ORDER` in `config.R`). Cluster 13 was dropped as
debris — PTPRC detected in 43% of its cells against 97% elsewhere, median 808
genes against 2,322, ~48% of UMIs in histone genes, no lineage marker, and
present in all 8 samples. Final: 19,826 cells.

Cluster 9 was resolved by TCR carriage, which is independent of the GEX used
to cluster it. It sat between T and NK. Measured receptor carriage split it
cleanly: subclusters 9_0 (60.2%) and 9_5 (59.3%) carry TCRs at T-cell rates and
became CD8 effector/TEMRA; 9_1–9_4 (0–1.5%) match the NK rate of 0.6% and
became NK cell.

No treg cluster is recoverable. FOXP3 is detected in 227 of 20,050 cells
(1.13%); 75 co-express IL2RA and 58 of those are IL7R-negative. Subclustering
C7/C8 and a dedicated CD4 re-embedding both failed — the best CD4 subcluster
reaches FOXP3 at 12.8% with IL7R still at 84.8%. Carried as a cell-level flag
(`Treg_like`, 39 cells strict), never as a cluster.

Audit (`10_annotation_audit.R`) forced two renames. Part 1 falsification:
14/14 pass. Part 2 discrimination: 12/14 unique, the two non-unique being the
expected programme-qualified pairs. Part 3: only CD8 effector T (cycling)
exceeds 30% programme markers (36%, all cycling), which is the intended
meaning of that label.
- Memory B → B cell (unresolved): failed the IGHD ≤25% ceiling at 41.5%,
  and CD27 reached only 3.2% — below naive B's 5.6%. The data do not support
  a memory call.
- CD4 memory → CD4 effector memory T: CCR4 41.9%, KLRB1 57.5%, S100A4
  69.6%. Thresholds for the CD4 CM/EM split were tightened post-hoc after the
  Part 2 failure; this is disclosed in the script.

## 8. Dissociation stress separates the controls from the ATI samples

A stress / immediate-early programme (FOS, FOSB, JUN, JUNB, HSPA1A/B, DUSP1,
PER1, GADD45G) accounts for 1.51% and 1.45% of UMIs in BC_control_A/B
against 0.24–0.41% in the ATI samples — roughly 5×. The duplicate agrees
with itself, so this is a property of how those two libraries were handled, not
noise. Within the ATI samples it is nearly flat (B4 0.33%, interrupt 0.28%).

This is the third independent reason control-versus-ATI is not a comparison
arm, alongside well and donor.

## 9. Cell cycle: the S-phase signal is largely artefactual

The ~40% vs ~20% S-phase difference flagged in section 6 does not survive
checking against proliferation transcripts.

Median proliferation-gene count is 0 in every phase, including S and G2M.
Only 11.5% of S-labelled cells detect ≥2 proliferation transcripts, against
2.6% of G1 cells. `CellCycleScoring` assigns phase by relative score, so in a
mostly-resting PBMC population it over-calls S.

Measured on transcripts rather than phase labels, real proliferation is ATI
median 8.81% against control 2.70% — not 43% vs 21%. It is also not
consistent per participant: 149496 +13.9 pp, 150815 +2.9 pp, 154026 −3.1 pp.

Internal controls behave on the transcript measure: plasma cells are 43.4%
proliferating against 0.72–2.01% across the three monocyte subsets
(`CHER_proliferation_by_annotation.csv`).

The phase label does not behave as cleanly even in monocytes, which is the
point of this section: S-phase label runs 0.62% (inflammatory classical) to
3.29% (classical) but reaches 11.34% in non-classical monocytes — a
population with 1.64% proliferating cells by transcript count. Quote the
proliferation column, not `pct_S`.

## 10. TCR repertoire

7,862 cells (39.7%) carry a receptor across 6,706 clonotypes. Coverage
validates the annotation independently of the GEX it was built from: T
clusters 61.6–92.2%, NK 0.7%, B / monocyte / pDC 0–1.2%.

Clonal expansion rises B4 → Interrupt in all three participants
(`CHER_clonal_expansion_paired.csv`; clone size counted within sample):

| Participant | B4 | Interrupt | Change |
|---|---|---|---|
| 149496 | 20.3% | 41.0% | +20.7 pp |
| 154026 | 10.0% | 28.8% | +18.8 pp |
| 150815 | 2.9% | 6.9% | +4.0 pp |
| BC_A vs BC_B (floor) | 15.6% | 14.0% | 1.6 pp |

The expansion is confined to CD8: effector/TEMRA 55.7% → 62.4% and effector T
(cycling) 33.3% → 62.5%, while every CD4 subset stays at 0–1.8%.

Diversity falls in the same three participants, rarefied to 593 TCR+ cells
over 100 draws (raw diversity is depth-dependent — section 6):

| Participant | Δ richness | Δ Shannon | Δ inverse simpson |
|---|---|---|---|
| 149496 | −122.6 | −0.68 | −142.7 |
| 154026 | −66.9 | −0.23 | −167.3 |
| 150815 | −9.8 | −0.03 | −27.0 |
| floor (BC_A vs BC_B) | 6.1 | 0.02 | 22.2 |

149496 and 154026 clear the floor decisively on every metric. 150815 clears it
only marginally (Shannon 0.03 vs 0.02; inverse simpson 27.0 vs 22.2) — the
direction agrees but that participant carries little weight.

## 10a. Repertoire structure (`18_tcr_repertoire_depth.R`)

Metrics computed from the metadata with dplyr rather than scRepertoire's
plotting functions, so each number traces to a CSV. Clone size counted within
sample throughout.

Clonal occupancy — share of TCR+ cells sitting in a clone of each size:

| Sample | Single | Small (2–5) | Medium (6–20) | Large (>20) |
|---|---|---|---|---|
| 149496_B4 | 79.66 | 5.79 | 4.24 | 10.30 |
| 149496_Interrupt | 59.02 | 11.47 | 11.47 | 18.04 |
| 150815_B4 | 97.14 | 2.86 | 0.00 | 0.00 |
| 150815_Interrupt | 93.15 | 5.79 | 1.06 | 0.00 |
| 154026_B4 | 90.02 | 7.40 | 2.58 | 0.00 |
| 154026_Interrupt | 71.20 | 13.31 | 15.48 | 0.00 |
| BC_control_A | 84.39 | 12.37 | 3.24 | 0.00 |
| BC_control_B | 85.97 | 11.94 | 2.10 | 0.00 |

Clonal proportion — share of TCR+ cells held by the top-N clones. The
duplicate is the floor and it is remarkably tight (25.04% vs 25.00% at N=100):

| Sample | top1 | top10 | top100 |
|---|---|---|---|
| 149496_B4 | 3.84 | 14.55 | 24.18 |
| 149496_Interrupt | 7.08 | 30.35 | 52.28 |
| 154026_B4 | 0.97 | 5.47 | 17.49 |
| 154026_Interrupt | 1.45 | 9.62 | 29.67 |
| 150815_B4 | 0.53 | 2.96 | 12.49 |
| 150815_Interrupt | 0.57 | 3.83 | 12.81 |
| BC_control_A / _B | 1.18 / 1.13 | 7.07 / 6.61 | 25.04 / 25.00 |

### Clonotype overlap, and why it cannot be read naively

Two libraries from the same blood draw share only 8.14% of their clonotypes
(Jaccard 0.0408). At this sequencing depth each library samples a largely
different slice of the repertoire. Low B4→Interrupt overlap is therefore the
detection floor, not turnover — all three ATI pairs have jaccard below the
duplicate (0.0171, 0.0062, 0.0146) and that means nothing on its own.

The metric that does carry information is what fraction of the expanded
(≥2-cell) Interrupt clones were detectable at B4, read against the same
ceiling from the duplicate:

| Participant | Expanded clones at interrupt | Seen at B4 | vs duplicate ceiling 65.6% |
|---|---|---|---|
| 149496 | 33 | 60.6% | at ceiling |
| 150815 | 27 | 22.2% | below |
| 154026 | 88 | 13.6% | well below |

So 149496 expanded clones that were already there, at essentially the
technical ceiling, whereas 154026's expanded clones were mostly not
detectable before interruption — on 88 clones, the largest expanded set of
the three. Two different routes to the same composition shift. Note the
asymmetry in sampling: 149496_Interrupt yielded only 383 clones against 1,226
at B4, so its "seen at B4" rate is favoured by the deeper B4 library; 154026 is
the reverse (869 → 1,072).

Negative controls behave. CDR3 length is flat across every sample (median
α 13–14, β 14–15), and expanded clones spanning more than one cell type run
11–27% against a duplicate range of 15.6–19.4% — no shift beyond the floor
except 150815_Interrupt at 44.4%, which rests on 27 clones.

150815 is the non-responder on every repertoire measure — occupancy flat,
top-100 share 12.49 → 12.81 against a 0.04 floor, diversity at the floor. The
expansion result therefore holds in 2 of 3 participants.

## 11. CoNGA

Pipeline: `14_conga_export.R` → `14_conga_run.sh` → `14_conga_interpret.py`.
4,769 paired-αβ clones; 303 hit rows over 209 unique clones.

The run is pooled across all 8 libraries. A graph-vs-graph hit is therefore
a dataset-level statement — a clone whose transcriptional neighbours and whose
receptor neighbours overlap more than chance — and is not a per-timepoint or
per-participant call. Library composition is reported only as a descriptive
diagnostic and was deliberately not split by timepoint.

The logo panel (`CHER_graph_vs_graph_logos.png`) renders via cairosvg. CoNGA
probes for a converter command with `shutil.which()`, so invoking its
interpreter by absolute path leaves the env's `bin/` off PATH and the logo step
fails while every table is written correctly — `14_conga_run.sh` now prepends it.

Clone-level QC was set to this project's thresholds, not CoNGA's defaults,
which would have imposed a second undocumented 3,500-gene ceiling selecting
against exactly the activated CD8 cells this analysis concerns.

The dominant hit signal is MAIT, which is a positive control. TRAV1-2 is in
29.2% of hit clones against 3.4% of background (8.6×); the full MAIT receptor
configuration is 17.2% against 1.3% (13×); 31 of 209 clones score
`mait_fraction` > 0.5. An invariant receptor paired with a distinct
transcriptional programme is precisely what graph-vs-graph is built to detect.

Two caveats on the hit set. 51.4% of hit cells come from `149496_B4` alone
(enrichment 2.85× over its share of clone-assigned cells), so the hits are not
evenly sourced. And the hit set is dominated by germline-segment structure
rather than convergent CDR3s.

### 11a. CoNGA surfaced a MAIT population the annotation had absorbed

Followed up in `16_mait_check.R`, testing the receptor (VDJ library)
against the programme (GEX library) — independent evidence.

81 cells carry a MAIT-configuration receptor. Receptor and programme agree in
the CD4 memory compartments and disagree in the CD8 ones:

| Label carrying the receptor | n | SLC4A10 (receptor+ vs −) | Verdict |
|---|---|---|---|
| CD4 effector memory T | 39 | 53.8% vs 7.2% (×7.5) | MAIT |
| CD4 central memory T | 10 | 70.0% vs 3.3% (×21.2) | MAIT |
| CD8 effector/TEMRA | 22 | 0.0% vs 0.5% | TRAV1-2+, not MAIT |
| CD8 naive/CM T | 8 | 0.0% vs 0.1% | TRAV1-2+, not MAIT |
| CD4 naive T | 2 | — (below n=5) | not assessed |

So ~49 genuine MAIT cells sit inside the CD4 memory labels, carrying KLRB1
(97.4% / 80.0%), NCR3, ZBTB16, CXCR6 and DPP4 well above their label-mates. At
0.25% of the dataset they will not form a cluster at any resolution. The
remaining 30 TRAV1-2+ CD8 cells are conventional T cells using that V segment,
not MAIT — the α-chain alone does not make the call.

Per library they run 0.12–0.56% of cells with no consistent timepoint pattern
(149496 0.38→0.54, 150815 0.26→0.12, 154026 0.56→0.52). The duplicate agrees
with itself (0.45% vs 0.48%).

These cells are not represented in `ANNO_ORDER`; the receptor and programme
flags in `CHER_mait_check_programme_by_label.csv` are the only record of them.

## 11b. Known-receptor matching (TCR Explorer)

`20_tcr_explorer_export.R` then `21_tcr_explorer_query.py`, run locally against
the harmonized VDJdb, IEDB, McPAS and TCR3d index built by
`tcr-explorer-refresh` (486,542 records; 329,746 human beta chains carrying a
CDR3). The export groups by beta and alpha chain, so it is collapsed to the
beta receptor before matching: 6,618 beta-alpha rows reduce to 6,467 unique
beta receptors per sample, 262 of them from expanded clones. Matching on a
beta-level key without that collapse discards the second alpha pairing of a
receptor and can drop an expanded clone.

Two levels are reported because they carry different weight: CDR3-only (same
CDR3 beta) and CDR3 plus V (same CDR3 beta and same TRBV). Convergent CDR3s
arise on different V genes with different specificity, so CDR3 plus V is the
defensible level.

| Match level | Group | Receptors | Hits | Rate |
|---|---|---|---|---|
| CDR3 + V | expanded | 262 | 3 | 1.15% |
| CDR3 + V | not expanded | 6,205 | 83 | 1.34% |
| CDR3 only | expanded | 262 | 16 | 6.11% |
| CDR3 only | not expanded | 6,205 | 442 | 7.12% |

Fisher exact on CDR3 plus V gives an odds ratio of 0.85 at p = 1. Expanded
clones therefore match known receptors no more often than the unexpanded
background. Were the expansion a bystander reactivation of common chronic-virus
memory, the expanded set would be expected to match known CMV, EBV and
influenza receptors more often than background. With 3 hits in the expanded set
the comparison excludes a large effect and not a moderate one.

The hit set reflects database composition. SARS-CoV-2 is the largest category
at 280 receptors, consistent with the post-2020 saturation of these databases,
followed by CMV at 141, M. tuberculosis at 68, influenza A at 56 and EBV at 45.
A hit rate here is not an unbiased estimate of what the repertoire recognises.

### HIV-annotated matches

The index carries 4,325 HIV-annotated beta records. Thirteen of our receptors
match an HIV record at CDR3-only level, across 20 cells. Two are from expanded
clones, both in 154026_Interrupt and both labelled CD8 effector/TEMRA:
TRBV11-3 CASSLTGELFF (3 cells) and TRBV7-3 CASSRLAGGTDTQYF (4 cells). The
second also matches at CDR3 plus V level.

Three observations bound how much these can carry:

- Both expanded receptors are promiscuous. CASSLTGELFF matches records
  annotated to 4 distinct organisms and CASSRLAGGTDTQYF to 3, so neither
  identifies a single specificity. Across all HIV-matching receptors the median
  is 2 organisms, against 1 for matching receptors generally, so HIV-annotated
  matches are more promiscuous than matches in general.
- The CDR3s are short and germline-proximal, median 12 aa against 15 aa across
  all receptors. Sequences such as CASSLSYEQYF and CASSYGGTEAFF are public
  across donors.
- Two of the 15 matching rows are in BC_control, who is HIV-negative. A match
  in an HIV-negative library is a false positive by construction, which bounds
  the specificity any CDR3-only hit can support.

HLA restriction was not checked; the participants' genotypes are not available
here, and public HIV-specific receptors are concentrated in a small number of
restrictions. The MHC fields are retained in `CHER_tcr_explorer_hits.csv`.

## 12. Composition shift, B4 → Interrupt

`17_composition_shift.R`. No p-values, deliberately: three pairs cannot support
a per-cell-type test, and shares are compositionally dependent (one subset
expanding forces the others down). What is reported is direction in each
participant, size, and whether it clears the BC duplicate floor.

Two shifts are consistent in all three participants and clear the floor:

| Cell type | 149496 | 150815 | 154026 | floor | min cells moved |
|---|---|---|---|---|---|
| CD8 effector/TEMRA | +3.20 | +4.25 | +13.46 | 0.68 | 114 |
| CD4 central memory T | −1.52 | −5.66 | −2.48 | 1.13 | 51 |

A floor test alone also passed plasma cell and CD8 effector T (cycling), on
changes of 1 and 8 cells in the weakest participant. A minimum 10-cell guard
was added and both are now reported as floor-only rather than consistent.

The plasma cell row is one participant, not a trend. 149496 goes from 0 to
323 plasma cells (21.7% of that library) between B4 and interrupt, while 150815
moves by 3 cells and 154026 by 1. It is a single observation and is not
supported by the other two participants.

Nothing else is consistent: NK cell moves +2.76/+3.03/−11.17, CD4 naive T
−20.42/+0.50/+0.46 (the mirror of 149496's plasma-cell burst).

## 13. Differential expression is null

`15_dge_pseudobulk.R`. Counts summed per participant × timepoint × cell type;
edgeR QL F-test with `~ Participant + Timepoint`, so the between-participant
offset is removed and the timepoint effect is tested within pair. Pseudobulk
rather than per-cell because the unit of replication is the participant —
a per-cell test would return p-values driven by cell count.

4 significant gene × cell-type calls at FDR < 0.05 and |logFC| > 0.5 across
the 10 testable cell types:

| Cell type | Gene | logFC | FDR | Direction |
|---|---|---|---|---|
| CD4 central memory T | HSP90AA1 | −1.17 | 0.033 | up in B4 |
| CD8 effector/TEMRA | S100B | −2.40 | 0.049 | up in B4 |
| CD8 effector/TEMRA | FOS | −2.20 | 0.049 | up in B4 |
| CD8 effector/TEMRA | SMIM43 | +5.02 | 0.049 | up in interrupt |

Three of the four are the handling confound from section 8. HSP90AA1 and
FOS are core members of the stress / immediate-early programme, and S100B
travels with it; all three are up in B4, matching the small residual stress
difference within the ATI samples (0.33% vs 0.28%). SMIM43 sits at logCPM 2.8,
where a handful of counts moves a 5-log-fold call.

Four cell types were not testable at all — 17 of 81 pseudobulk libraries fall
below the 20-cell floor (Plasma cell and classical monocyte have essentially no
cells in most libraries).

This is a negative result at low power rather than evidence of no change.
Three pairs is near the floor for pseudobulk testing.

## 14. Summary of supported conclusions

1. CD8 effector/TEMRA expand during interruption, in share of the library
   (section 12) and in clonal structure (section 10), consistently in all three
   participants and above the technical floor. These are the same cells by two
   independent measures — GEX composition and VDJ clonality.
2. Repertoire diversity contracts in the same direction in all three,
   decisively in two.
3. CD4 central memory contracts in all three.
4. Transcriptional differences are not detectable at n = 3 beyond a stress
   signature attributable to handling.
5. Control-versus-ATI remains uninterpretable — well, donor and
   dissociation stress, all three confounded.
