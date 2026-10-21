############################################################
# CHER ATI — central configuration
#
# source() this at the top of every CHER script so paths, the sample sheet and
# the QC thresholds live in ONE place.
#
#   source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")
#
# Adapted from CD8_Longitudinal "Scripts/CP018 - FACTseq/config.R". Same OCM
# chemistry and the same Cell Ranger 10.1.0, so the structure carries over; the
# sample sheet, the two-well layout and the pairing are CHER-specific.
############################################################

# ── Root paths ────────────────────────────────────────────────────────────────
PROJ_ROOT <- path.expand("~/Desktop/Projects")
REPO_DIR  <- file.path(PROJ_ROOT, "Repositories", "CHER_ATI")
SAVED_DIR <- file.path(REPO_DIR, "saved_R_data")            # gitignored
CR_ROOT   <- file.path(PROJ_ROOT, "cellranger_out")
DATADRIVE <- "/mnt/data/CHER_ATI"                            # FASTQ + CR configs

# Two OCM GEM wells = two cellranger multi runs. Sample-level output sits at
#   <run>/outs/per_sample_outs/<sample_id>/
RUN_DIRS <- c("34349" = file.path(CR_ROOT, "CHER_34349_multi_vdj"),
              "33445" = file.path(CR_ROOT, "CHER_33445_multi_vdj"))

# ── Project-local R library ──────────────────────────────────────────────────
# The conda env (cd8-seurat2) ships scDblFinder, batchelor and
# SingleCellExperiment but NOT GenomeInfoDbData, which GenomeInfoDb requires at
# load time. Without it every Bioconductor-dependent package fails with
# "there is no package called 'GenomeInfoDbData'" — which looks like the package
# is missing rather than its data dependency. Installing
# bioconductor-genomeinfodbdata via conda does NOT fix it (it is a post-link
# script package and the script does not populate the env's R library), so the
# working fix is a project-local library placed FIRST on the path.
#
# Do NOT set R_LIBS_SITE to /usr/lib/R/site-library when running this env: the
# system `sp` build fails with "undefined symbol: NO_REFERENCES" and kills any
# Seurat script at load.
.CHER_RLIB <- file.path(REPO_DIR, ".Rlib")
if (dir.exists(.CHER_RLIB)) .libPaths(c(.CHER_RLIB, .libPaths()))

# ── Analysis output roots ────────────────────────────────────────────────────
ANALYSIS_ROOT <- file.path(REPO_DIR, "Analysis")
QC_ROOT       <- file.path(ANALYSIS_ROOT, "QC")
INT_ROOT      <- file.path(ANALYSIS_ROOT, "Integration")
VDJ_ROOT      <- file.path(ANALYSIS_ROOT, "VDJ")

# ── Object checkpoints (qs2) ─────────────────────────────────────────────────
# One per pipeline stage, so a failure resumes from the previous step rather
# than from the H5 files.
QS_RAW        <- file.path(SAVED_DIR, "CHER_01_merged_raw.qs2")
QS_FILTERED   <- file.path(SAVED_DIR, "CHER_03_postQC_filtered.qs2")
QS_CELLCYCLE  <- file.path(SAVED_DIR, "CHER_04_cellcycle.qs2")
QS_SINGLETS   <- file.path(SAVED_DIR, "CHER_05_postQC_scDblFinder_singlets.qs2")
QS_INTEGRATED <- file.path(SAVED_DIR, "CHER_06_RNA_integrated_CCA_MNN.qs2")
QS_ANNOTATED  <- file.path(SAVED_DIR, "CHER_09_annotated.qs2")
QS_WITH_TCR   <- file.path(SAVED_DIR, "CHER_10_annotated_withTCR.qs2")

# ── Sample sheet ─────────────────────────────────────────────────────────────
# From L. de Armas, confirmed against both BaseSpace share pages 2026-10-19.
# Sample_ID is the cellranger `sample_id` and therefore the per_sample_outs
# folder name. Spaces are illegal in a cellranger sample_id, so "BC control_A"
# became "BC_control_A".
#
# PAIRING IS THE POINT: each participant has BOTH timepoints in the SAME well,
# so the within-participant contrast carries no well effect. The two healthy
# controls are both in 33445 — a control-vs-ATI contrast IS confounded with well.
CHER_SAMPLES <- data.frame(
  Sample_ID   = c("149496_B4", "149496_Interrupt",
                  "150815_B4", "150815_Interrupt",
                  "154026_B4", "154026_Interrupt",
                  "BC_control_A", "BC_control_B"),
  Participant = c("149496", "149496", "150815", "150815",
                  "154026", "154026", "BC_control", "BC_control"),
  Timepoint   = c("B4", "Interrupt", "B4", "Interrupt",
                  "B4", "Interrupt", "Control", "Control"),
  # Replicate applies only to the healthy control, which was run in duplicate.
  Replicate   = c(NA, NA, NA, NA, NA, NA, "A", "B"),
  OCM_Barcode = c("OB2", "OB1", "OB3", "OB4", "OB1", "OB2", "OB3", "OB4"),
  Tube_Label  = c("2", "1", "3", "4", "1", "2", "3", "4"),
  Well        = c("34349", "34349", "34349", "34349",
                  "33445", "33445", "33445", "33445"),
  stringsAsFactors = FALSE
)
# Group separates the three ATI participants from the healthy control donor.
CHER_SAMPLES$Group <- ifelse(CHER_SAMPLES$Participant == "BC_control",
                             "Healthy_control", "ATI")
# Participant_Timepoint is the per-sample label used on most plots.
CHER_SAMPLES$Participant_Timepoint <- paste0(CHER_SAMPLES$Participant, "_",
                                             CHER_SAMPLES$Timepoint)

# Resolve the filtered-matrix H5 for a sample, failing loudly rather than
# silently returning NA. Cell Ranger 10.x writes it directly in the sample dir;
# 9.x used a count/ subfolder — both are checked so this works if the project
# is ever re-run on an older aligner.
sample_h5 <- function(sample_id) {
  row <- CHER_SAMPLES[CHER_SAMPLES$Sample_ID == sample_id, ]
  if (nrow(row) != 1L) stop("Unknown Sample_ID: ", sample_id, call. = FALSE)
  base <- file.path(RUN_DIRS[[row$Well]], "outs", "per_sample_outs", sample_id)
  cand <- c(file.path(base, "sample_filtered_feature_bc_matrix.h5"),
            file.path(base, "count", "sample_filtered_feature_bc_matrix.h5"))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("No filtered matrix H5 for ", sample_id, ". Looked in:\n  ",
         paste(cand, collapse = "\n  "), call. = FALSE)
  }
  hit[1]
}

# Matched TCR contig annotations for a sample (used from script 07 onward).
sample_vdj <- function(sample_id) {
  row <- CHER_SAMPLES[CHER_SAMPLES$Sample_ID == sample_id, ]
  p <- file.path(RUN_DIRS[[row$Well]], "outs", "per_sample_outs", sample_id,
                 "vdj_t", "filtered_contig_annotations.csv")
  if (!file.exists(p)) stop("No VDJ contigs for ", sample_id, ":\n  ", p, call. = FALSE)
  p
}

# Cells carrying a productive TCR, as a named logical keyed by Seurat cell name.
# NOTE the `productive` field in Cell Ranger 10.1 contig CSVs is lowercase
# "true"; matching "True"/TRUE returns ZERO rows and reads as "no T cells
# anywhere" rather than as an error.
tcr_cells <- function(cell_names, paired_only = FALSE) {
  tc <- do.call(rbind, lapply(CHER_SAMPLES$Sample_ID, function(s) {
    d <- utils::read.csv(sample_vdj(s), stringsAsFactors = FALSE)
    d <- d[tolower(as.character(d$productive)) == "true" &
           tolower(as.character(d$high_confidence)) == "true" &
           d$chain %in% c("TRA", "TRB"), ]
    if (!nrow(d)) return(NULL)
    u   <- unique(d[, c("barcode", "chain")])
    agg <- stats::aggregate(chain ~ barcode, data = u,
                            FUN = function(x) paste(sort(unique(x)), collapse = "+"))
    data.frame(cell = paste0(s, "_", agg$barcode), chains = agg$chain,
               stringsAsFactors = FALSE)
  }))
  keep <- if (paired_only) tc$cell[tc$chains == "TRA+TRB"] else tc$cell
  stats::setNames(cell_names %in% keep, cell_names)
}

# ── Shared aesthetics ────────────────────────────────────────────────────────
# B4 / Interrupt are the contrast; the healthy control is deliberately neutral
# grey so it never reads as a third arm of the comparison.
TP_COLORS    <- c("B4" = "#2171B5", "Interrupt" = "#D95F02", "Control" = "#777777")
GROUP_COLORS <- c("ATI" = "#2171B5", "Healthy_control" = "#777777")
WELL_COLORS  <- c("34349" = "#52B788", "33445" = "#B07AA1")
PHASE_COLS   <- c("G1" = "#AEC6CF", "S" = "#FFD700", "G2M" = "#E63946")

# ── QC thresholds ────────────────────────────────────────────────────────────
# Starting values were CP018's. TWO WERE CHANGED for CHER, each because the
# measured cost showed the inherited value was wrong for whole PBMC rather than
# sorted CD8 T cells. Both measurements are in
# Analysis/QC/CHER_threshold_cost_per_sample.csv and the findings document.
QC_MIN_COUNT          <- 500
QC_MIN_FEATURE        <- 600
QC_MIN_COMPLEX        <- 0.80    # log10(genes)/log10(UMI)
QC_MAX_MITO           <- 15      # percent
QC_MAX_HB             <- 20      # percent
QC_MAX_PLAT           <- 2       # percent, now on PF4/PPBP only (see below)
QC_MIN_CELLS_PER_GENE <- 10

# CHANGE 1 — the ribosomal FLOOR is disabled (was: percent_ribo > 5).
# This is a JUDGEMENT CALL, not a defect in the CP018 value; the numbers are
# recorded so it can be overruled.
#   - The floor would remove 3,647 cells (15.3%), but 64.6% of those ALREADY
#     fail a direct quality criterion (counts/features/complexity/mito/hb/
#     platelet), so most of its work is redundant.
#   - The 1,291 cells unique to it (5.4% of the dataset) are BORDERLINE rather
#     than clearly bad: median 1,488 genes against 2,362 in retained cells, mito
#     7.65% against 2.88%, but complexity 0.899 — ABOVE the retained median of
#     0.870.
#   - A ribosomal minimum is also not a standard whole-PBMC criterion.
# Decision: let nFeature/mito/complexity carry the quality filtering directly
# rather than have an indirect proxy account for 5.4% of removals.
# NOTE an earlier version of this comment claimed the floor preferentially
# removed monocytes. That was measured on a partially-read counts matrix and is
# WRONG: monocyte-marker-high cells are 0.1% of removed and 0.1% of kept.
# percent_ribo is still COMPUTED and plotted; it is simply not a filter.
QC_MIN_RIBO           <- 0       # 0 = disabled

# CHANGE 2 — percent_plat is PF4/PPBP, not PECAM1|PF4.
# Two reasons. (a) PECAM1 (CD31) is broadly expressed across leukocytes, so it
# does not measure platelet contamination. (b) PPBP is the dominant platelet
# transcript and was absent from the CP018 pattern. Net effect is to flag MORE
# cells, not fewer: 4.54% at >=2% by PF4|PPBP against 1.84% by PECAM1|PF4. Those
# cells are genuinely platelet-dominated — median 356 genes against 2,298 in
# retained cells — so this is a real quality filter rather than a cell type.
PLATELET_PATTERN <- "^PF4$|^PPBP$"
RIBO_PATTERN     <- "^RP[SL]"
MITO_PATTERN     <- "^MT-"
# Haemoglobin: anchored to the actual globin genes. The commonly-copied
# "^HB[^(P)]" also matches HBEGF, HBS1L and similar non-globin genes.
HB_PATTERN       <- "^HB[ABDEGQZ][0-9]?$|^HBA[12]$"

# ── Clustering ───────────────────────────────────────────────────────────────
MIN_CLUSTER_SIZE <- 10
# FastMNN at 0.6, 16 clusters. NOTE the miR29a marker-redundancy criterion did
# NOT discriminate here: zero redundant pairs at every level 0.2-1.0, so it
# would have picked 1.0. 0.6 was chosen on nameability instead — every cluster
# gets a defensible lineage call, which higher levels start to lose. That is a
# judgement, recorded as such.
MNN_RES_FINAL    <- 0.6
CCA_RES_FINAL    <- NA_real_   # CCA embedding kept but not taken forward

# ── Annotation vocabulary ────────────────────────────────────────────────────
# Lineage identity first; a PROGRAMME (cycling, inflammatory) is only ever a
# bracketed qualifier, never the noun — naming a cluster after a programme hides
# what the cell is and makes composition analysis incoherent.
# Two labels were CORRECTED by the audit (script 10), not chosen up front:
#   "Memory B cell" -> "B cell (unresolved)". It failed its own falsification
#      criterion (IGHD <=25 required, observed 41.5) and CD27 — the memory B
#      marker — is detected in only 3.2% of its cells, LOWER than the naive B
#      cluster's 5.6%. It is a distinct population (1,351 DE genes vs naive B,
#      FCRL1 50.5% vs 23.9%, BANK1 81.2% vs 68.2%) but is neither memory nor
#      class-switched on the evidence, so it is named unresolved rather than
#      guessed.
#   "CD4 memory T" -> "CD4 effector memory T". Justified by CCR4 41.9%,
#      KLRB1 57.5%, S100A4 69.6%, AQP3 47.9% against 3.0 / 2.5 / 8.9 / 24.5 in
#      CD4 naive — an effector-memory phenotype, not generic memory.
ANNO_ORDER <- c(
  "CD4 naive T", "CD4 central memory T", "CD4 effector memory T",
  "CD8 naive/CM T", "CD8 effector/TEMRA", "CD8 effector T (cycling)",
  "NK cell",
  "Naive B cell", "B cell (unresolved)", "Plasma cell",
  "Classical monocyte", "Classical monocyte (inflammatory)",
  "Non-classical monocyte", "pDC")

ANNO_COLS <- c(
  "CD4 naive T"                      = "#2171B5",
  "CD4 central memory T"             = "#6BAED6",
  "CD4 effector memory T"            = "#9ECAE1",
  "CD8 naive/CM T"                   = "#1A9850",
  "CD8 effector/TEMRA"               = "#E63946",
  "CD8 effector T (cycling)"         = "#9D4EDD",
  "NK cell"                          = "#F4A261",
  "Naive B cell"                     = "#762A83",
  "B cell (unresolved)"              = "#C2A5CF",
  "Plasma cell"                      = "#40004B",
  "Classical monocyte"               = "#8C510A",
  "Classical monocyte (inflammatory)"= "#BF812D",
  "Non-classical monocyte"           = "#DFC27D",
  "pDC"                              = "#00838F")

# Clusters removed before annotation, with the reason. C13 is NOT a cell type:
# PTPRC (CD45) detected in only 43% of its cells against 97% elsewhere, median
# 808 genes against 2,322, replication-dependent histones at ~48%, no lineage
# marker at all, spread across all 8 samples. Debris the QC filters missed.
DROP_CLUSTERS <- c("13")

# ── Output location policy ───────────────────────────────────────────────────
# Everything written goes under Analysis/. Manuscript/ holds figures only once
# we decide to promote them; save_png()/save_csv() enforce this rather than
# relying on every script getting the path right.
MANUSCRIPT_DIR <- file.path(REPO_DIR, "Manuscript")

assert_not_manuscript <- function(path) {
  ap <- normalizePath(dirname(path), mustWork = FALSE)
  if (startsWith(ap, normalizePath(MANUSCRIPT_DIR, mustWork = FALSE))) {
    stop("REFUSING to write into Manuscript/:\n  ", path,
         "\nCHER outputs belong under Analysis/.", call. = FALSE)
  }
  invisible(TRUE)
}

require_file <- function(path, what = "input") {
  if (!file.exists(path)) {
    stop(what, " not found:\n  ", path,
         "\nCheck config.R paths, or that the upstream script has been run.",
         call. = FALSE)
  }
  invisible(path)
}

mkdirs <- function(...) invisible(lapply(c(...), dir.create,
                                        recursive = TRUE, showWarnings = FALSE))

# ── Namespace conflict guard ─────────────────────────────────────────────────
# Bioconductor packages (SingleCellExperiment/S4Vectors/IRanges, pulled in by
# scDblFinder and batchelor) mask several dplyr/tidyr verbs. Running the scripts
# sequentially in ONE session leaves those attachments on the search path, and
# dplyr::count then resolves to the wrong function. Binding the dplyr/tidyr
# versions into globalenv fixes it for every downstream script, because
# globalenv is searched before any attached package.
.prefer_verbs <- function(pkg, fns) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(invisible(NULL))
  ns <- asNamespace(pkg)
  for (f in fns) {
    if (exists(f, envir = ns, inherits = FALSE)) {
      assign(f, get(f, envir = ns), envir = globalenv())
    }
  }
}
.prefer_verbs("dplyr", c("count","filter","select","rename","slice","mutate",
                         "summarise","summarize","arrange","group_by","ungroup",
                         "first","last","desc","n","n_distinct","bind_rows",
                         "left_join","inner_join","full_join","transmute",
                         "recode","pull","distinct","across","all_of","case_when",
                         "if_else","lag","lead","between","setdiff","union",
                         "intersect","collapse"))
.prefer_verbs("tidyr", c("pivot_longer","pivot_wider","expand","separate","unite",
                         "replace_na","drop_na","fill","complete"))

# ── Optional presentation dependencies ───────────────────────────────────────
# SeuratExtend (DimPlot2/VlnPlot2) and Azimuth are presentation/annotation
# conveniences, not load-bearing for any result. Guard both so a missing
# package cannot fail the pipeline, and record which path was taken so a figure
# made with a fallback is not passed off as the real thing.
HAS_SEURATEXTEND <- requireNamespace("SeuratExtend", quietly = TRUE)
HAS_AZIMUTH      <- requireNamespace("Azimuth", quietly = TRUE)

if (HAS_SEURATEXTEND) {
  suppressPackageStartupMessages(library(SeuratExtend))
} else {
  message("SeuratExtend not installed - using base Seurat fallbacks for ",
          "DimPlot2 / VlnPlot2.")
  DimPlot2 <- function(object, features = NULL, group.by = NULL, reduction = NULL,
                       label = FALSE, cols = NULL, ncol = NULL, pt.size = NULL,
                       box = NULL, theme = NULL, ...) {
    if (!is.null(features)) {
      Seurat::FeaturePlot(object, features = features, reduction = reduction,
                          label = label, ncol = ncol, pt.size = pt.size, ...)
    } else {
      Seurat::DimPlot(object, group.by = group.by, reduction = reduction,
                      label = label, cols = cols, pt.size = pt.size, ...)
    }
  }
  VlnPlot2 <- function(object, features, group.by = NULL, show.mean = FALSE,
                       pt.size = 0, cols = NULL, ...) {
    p <- Seurat::VlnPlot(object, features = features, group.by = group.by,
                         pt.size = pt.size, cols = cols, ...)
    if (isTRUE(show.mean)) {
      p <- p + ggplot2::stat_summary(fun = mean, geom = "point",
                                     size = 2.4, colour = "black")
    }
    p + ggplot2::labs(x = NULL) + ggplot2::theme(legend.position = "none")
  }
}

# ── Plot saving ──────────────────────────────────────────────────────────────
FIG_DPI <- 350

theme_cher <- function(base_size = 16) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = base_size * 1.25, face = "bold",
                                               hjust = 0.5,
                                               margin = ggplot2::margin(b = 8)),
      plot.subtitle    = ggplot2::element_text(size = base_size * 0.95, hjust = 0.5,
                                               colour = "grey30"),
      axis.title       = ggplot2::element_text(size = base_size * 1.05, face = "bold"),
      axis.text        = ggplot2::element_text(size = base_size * 0.9, colour = "black"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5),
      strip.text       = ggplot2::element_text(size = base_size * 0.95, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.title     = ggplot2::element_text(size = base_size * 0.95, face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size * 0.9),
      plot.margin      = ggplot2::margin(12, 14, 10, 12)
    )
}

# width/height are INCHES.
save_png <- function(filename, plot_obj, width = 9, height = 6.5,
                     dpi = FIG_DPI, limitsize = FALSE) {
  assert_not_manuscript(filename)
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = filename, plot = plot_obj,
                  width = width, height = height, dpi = dpi,
                  bg = "white", limitsize = limitsize)
  invisible(filename)
}

# Size a multi-panel figure by its panel count so nothing gets crushed.
panel_size <- function(n_panels, ncol = 2, panel_w = 5.0, panel_h = 3.6) {
  nrow <- ceiling(n_panels / ncol)
  c(width = ncol * panel_w, height = nrow * panel_h)
}

# ── Filenames ────────────────────────────────────────────────────────────────
# Cell-type labels contain spaces, slashes and parentheses ("CD8 naive/CM T",
# "Classical monocyte (inflammatory)"). fs_safe turns one into a filename
# fragment that survives a shell and a Windows checkout.
fs_safe <- function(x) {
  x <- gsub("[^[:alnum:]]+", "_", as.character(x))
  gsub("(^_|_$)", "", x)
}

# ── CSV export ───────────────────────────────────────────────────────────────
# Every figure ships the numbers behind it (house rule). save_csv echoes the
# path so a skipped export is visible in the log.
save_csv <- function(x, path, row.names = FALSE) {
  assert_not_manuscript(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = row.names)
  message("  [csv] ", basename(path), "  (", nrow(x), " rows)")
  invisible(path)
}

# Cell-level metadata dump: the master table most downstream questions can be
# answered from without reloading the Seurat object.
export_metadata <- function(obj, path, extra_cols = NULL) {
  md <- obj@meta.data
  md$cell_barcode <- rownames(md)
  keep <- c("cell_barcode", intersect(
    c("Sample_ID","Participant","Timepoint","Group","Well","OCM_Barcode",
      "Replicate","Participant_Timepoint",
      "nCount_RNA","nFeature_RNA","percent_mito","percent_ribo","percent_hb",
      "percent_plat","log10GenesPerUMI",
      "S.Score","G2M.Score","Phase","scDblFinder.score","scDblFinder.class",
      "predicted.celltype.l1","predicted.celltype.l2",
      "cca_clusters_rna","mnn_clusters_rna","Annotation"),
    colnames(md)), extra_cols)
  save_csv(md[, unique(keep), drop = FALSE], path)
}

# Composition table: counts + within-group percentages for any two factors.
composition_table <- function(obj, group, split) {
  md <- obj@meta.data
  tb <- as.data.frame(table(group = md[[group]], split = md[[split]]))
  names(tb) <- c(group, split, "n_cells")
  tot <- tapply(tb$n_cells, tb[[split]], sum)
  tb$pct_of_sample <- round(100 * tb$n_cells / tot[as.character(tb[[split]])], 2)
  tb
}

rotate_x <- function(angle = 35, size = NULL) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(
    angle = angle, hjust = 1, vjust = 1, size = size))
}

# ── UMAP presentation ────────────────────────────────────────────────────────
# Seurat names UMAP axes after the reduction key ("umapmnnrna_1"). tidy_umap
# replaces that with readable labels; umap_labels adds haloed cluster labels
# that stay legible over dense points (DimPlot(label=TRUE) has no contrast
# control and is unreadable on dark clusters).
tidy_umap <- function(p, title = NULL, subtitle = NULL,
                      xlab = "UMAP 1", ylab = "UMAP 2", legend = TRUE) {
  p <- p + ggplot2::labs(x = xlab, y = ylab)
  if (!is.null(title))    p <- p + ggplot2::ggtitle(title)
  if (!is.null(subtitle)) p <- p + ggplot2::labs(subtitle = subtitle)
  p <- p + theme_cher() +
    ggplot2::theme(plot.title.position = "plot",
                   plot.title    = ggplot2::element_text(face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(hjust = 0),
                   axis.title    = ggplot2::element_text(face = "plain"),
                   legend.key    = ggplot2::element_blank())
  if (!legend) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# ── House UMAP idiom ─────────────────────────────────────────────────────────
# Ported from miR29a_PDAC_CD8/Scripts/config.R (umap_house + add_umap_arrows),
# which in turn takes its parameters verbatim from CD8_Longitudinal
# "Scripts/3. Manuscript/Figure_4.R". Boxed state names ON the embedding with
# leader lines, corner axis arrows instead of axes, square canvas.
#
# THE CANVAS MUST BE SQUARE. Each arrow arm is 0.12 of its OWN axis range, so
# the two arms only render equal on a square drawing area; on a wide panel the
# horizontal arm comes out nearly twice the vertical one.
add_umap_arrows <- function(p, df, x = "UMAP1", y = "UMAP2", size = 5) {
  xr <- diff(range(df[[x]], na.rm = TRUE)); yr <- diff(range(df[[y]], na.rm = TRUE))
  af <- 0.12
  xa <- min(df[[x]], na.rm = TRUE) - xr * 0.05
  ya <- min(df[[y]], na.rm = TRUE) - yr * 0.08
  p +
    ggplot2::annotate("segment", x = xa, xend = xa + xr * af, y = ya, yend = ya,
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
                      linewidth = 0.9, colour = "black") +
    ggplot2::annotate("text", x = xa + xr * af / 2, y = ya - yr * 0.04,
                      label = "UMAP1", size = size, fontface = "bold") +
    ggplot2::annotate("segment", x = xa, xend = xa, y = ya, yend = ya + yr * af,
                      arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm"), type = "closed"),
                      linewidth = 0.9, colour = "black") +
    ggplot2::annotate("text", x = xa - xr * 0.04, y = ya + yr * af / 2,
                      label = "UMAP2", size = size, fontface = "bold", angle = 90) +
    ggplot2::coord_cartesian(clip = "off")
}

umap_house <- function(obj, group = "Annotation", reduction = "umap.mnn",
                       cols = ANNO_COLS, title = NULL, pt.size = NULL,
                       label = TRUE, legend = FALSE, label.size = 5.2,
                       arrow.size = 5, split.by = NULL) {
  df <- as.data.frame(SeuratObject::Embeddings(obj, reduction))[, 1:2]
  colnames(df) <- c("UMAP1", "UMAP2")
  df$grp <- factor(obj@meta.data[[group]])
  if (!is.null(split.by)) df$split <- obj@meta.data[[split.by]]

  # Point size vs cell count: 1.2 at ~12k cells is the house value.
  if (is.null(pt.size)) pt.size <- max(0.8, min(1.6, 1.2 * sqrt(12000 / nrow(df))))

  p <- ggplot2::ggplot(df, ggplot2::aes(UMAP1, UMAP2, colour = grp)) +
    ggplot2::geom_point(size = pt.size, alpha = 0.7, stroke = 0) +
    ggplot2::expand_limits(
      x = range(df$UMAP1, na.rm = TRUE) + c(-1, 1) * 0.10 * diff(range(df$UMAP1, na.rm = TRUE)),
      y = range(df$UMAP2, na.rm = TRUE) + c(-1, 1) * 0.10 * diff(range(df$UMAP2, na.rm = TRUE))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      legend.text     = ggplot2::element_text(size = 11, face = "bold"),
      legend.title    = ggplot2::element_blank(),
      legend.position = if (legend) "bottom" else "none",
      strip.text      = ggplot2::element_text(size = 14, face = "bold"),
      plot.margin     = ggplot2::margin(10, 10, 20, 20))

  if (!is.null(cols)) {
    keep <- intersect(names(cols), levels(df$grp))
    if (length(keep) == nlevels(df$grp))
      p <- p + ggplot2::scale_colour_manual(values = cols, drop = FALSE)
  }
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  if (!is.null(split.by)) p <- p + ggplot2::facet_wrap(~ split)

  if (label && is.null(split.by)) {
    ctr <- stats::aggregate(cbind(UMAP1, UMAP2) ~ grp, df, stats::median)
    # House geom_label_repel parameters, verbatim.
    p <- p + ggrepel::geom_label_repel(
      data = ctr, ggplot2::aes(label = grp, colour = grp),
      fill = "white", size = label.size, fontface = "bold",
      label.size = 0.6, label.padding = ggplot2::unit(0.35, "lines"),
      box.padding = 1.5, point.padding = 0.5,
      segment.color = "grey40", segment.size = 0.6,
      min.segment.length = 0, max.overlaps = 30,
      force = 5, force_pull = 0.2, seed = 42, show.legend = FALSE)
  }
  add_umap_arrows(p, df, size = arrow.size)
}

umap_labels <- function(p, size = 4.2) {
  d   <- p$data
  xy  <- names(d)[1:2]
  grp <- if ("ident" %in% names(d)) "ident" else names(d)[3]
  cen <- stats::aggregate(d[, xy], by = list(lab = d[[grp]]), FUN = stats::median)
  p + ggrepel::geom_label_repel(
        data = cen, ggplot2::aes(x = .data[[xy[1]]], y = .data[[xy[2]]], label = lab),
        inherit.aes = FALSE, size = size, fontface = "bold",
        fill = grDevices::adjustcolor("white", alpha.f = 0.78),
        label.size = 0.15, box.padding = 0.35, seed = 1,
        min.segment.length = 0.2, segment.colour = "grey35")
}

# ── Shared QC plotting ───────────────────────────────────────────────────────
# Defined HERE rather than in the QC script so the pre-filter and post-filter
# figures are produced by identical code. If this lived in the script and were
# copied, the two panels could drift and the comparison would be meaningless.
QC_FEATS <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo",
              "percent_hb", "percent_plat", "log10GenesPerUMI")

# Thresholds keyed by metric, for the dashed reference line on each density.
QC_LINES <- c(nCount_RNA = QC_MIN_COUNT, nFeature_RNA = QC_MIN_FEATURE,
              log10GenesPerUMI = QC_MIN_COMPLEX, percent_mito = QC_MAX_MITO,
              percent_ribo = QC_MIN_RIBO, percent_hb = QC_MAX_HB,
              percent_plat = QC_MAX_PLAT)

run_qc_plots <- function(seu, group_col, out_dir) {
  stopifnot(group_col %in% colnames(seu@meta.data))
  mkdirs(out_dir)
  meta <- seu@meta.data

  # Cells per group, with the count printed so the figure ships its own number.
  p_cells <- ggplot2::ggplot(meta, ggplot2::aes(x = .data[[group_col]],
                                                fill = .data[[group_col]])) +
    ggplot2::geom_bar(colour = "black", linewidth = 0.3, width = 0.7) +
    ggplot2::geom_text(stat = "count",
                       ggplot2::aes(label = ggplot2::after_stat(count)),
                       vjust = -0.4, size = 4.5, fontface = "bold") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12)),
                                labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(x = NULL, y = "Cells", title = paste0("Cells by ", group_col)) +
    theme_cher() + rotate_x(30) +
    ggplot2::theme(legend.position = "none")
  save_png(file.path(out_dir, "Cells_per_group.png"), p_cells, width = 8, height = 5.5)
  save_csv(as.data.frame(table(meta[[group_col]]), responseName = "n_cells"),
           file.path(out_dir, "Cells_per_group_source_data.csv"))

  # pt.size = 0: with thousands of cells the jitter forms a solid mass that
  # hides the violin, and distribution shape is the point of this plot.
  p_vln <- Seurat::VlnPlot(seu, group.by = group_col, features = QC_FEATS,
                           pt.size = 0, ncol = 2) &
    theme_cher(base_size = 13) &
    ggplot2::labs(x = NULL) &
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  ps <- panel_size(length(QC_FEATS), ncol = 2, panel_w = 5.4, panel_h = 3.4)
  save_png(file.path(out_dir, "QC_features_grouped.png"), p_vln,
           width = ps[["width"]], height = ps[["height"]])

  dens <- function(var, vline, logx = TRUE) {
    p <- ggplot2::ggplot(meta, ggplot2::aes(colour = .data[[group_col]],
                                            x = .data[[var]],
                                            fill = .data[[group_col]])) +
      ggplot2::geom_density(alpha = 0.22, linewidth = 0.7) +
      ggplot2::geom_vline(xintercept = vline, linetype = "dashed",
                          colour = "grey25", linewidth = 0.7) +
      ggplot2::labs(x = var, y = "Cell density", colour = NULL, fill = NULL,
                    title = var,
                    subtitle = paste0("dashed line = QC threshold (", vline, ")")) +
      theme_cher(base_size = 14)
    if (logx) {
      p + ggplot2::scale_x_log10(labels = scales::label_number(big.mark = ",",
                                                              drop0trailing = TRUE))
    } else {
      p + ggplot2::scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE))
    }
  }

  for (v in names(QC_LINES)) {
    logx <- !(v %in% c("log10GenesPerUMI"))
    save_png(file.path(out_dir, paste0(v, ".png")),
             dens(v, QC_LINES[[v]], logx = logx), width = 8, height = 5.5)
  }

  # Per-group summary of every QC metric — the numbers behind the violins.
  summ <- do.call(rbind, lapply(split(meta, meta[[group_col]]), function(d) {
    data.frame(group = d[[group_col]][1], n_cells = nrow(d),
               t(vapply(QC_FEATS, function(f) stats::median(d[[f]], na.rm = TRUE),
                        numeric(1))), check.names = FALSE)
  }))
  names(summ)[1] <- group_col
  save_csv(summ, file.path(out_dir, "QC_metric_medians_source_data.csv"))

  invisible(out_dir)
}

# ── Canonical clonotype column ───────────────────────────────────────────────
# Ported from CD8_Longitudinal "Scripts/CP018 - FACTseq/config.R".
#
# Cell Ranger 10.1.0 writes an EMPTY D-gene field where 9.0.1 and older
# alignments write a literal "NA":
#   CR 9.0.1   TRBV6-1.NA.TRBJ2-5.TRBC2;TGTGCC...
#   CR 10.1.0  TRBV29-1..TRBJ1-2.TRBC1;TGCAGC...   <- ".."
# So a CHER clonotype string can NEVER match one from an older run on raw
# CTstrict, and the result is zero matches — indistinguishable from a real
# biological absence. normalise_ct() makes the conventions comparable and is
# idempotent, so re-stamping an object is harmless.
CT_CANON <- "CTstrict_norm"

normalise_ct <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA"] <- NA_character_
  for (i in 1:3) x <- gsub("\\.\\.", ".NA.", x, perl = TRUE)
  x
}

# Stamp CT_CANON ON the object, and record which convention it came from so
# provenance travels with the data rather than living in a script comment.
add_ct_canon <- function(obj, verbose = TRUE) {
  md <- obj@meta.data
  if (!"CTstrict" %in% colnames(md)) {
    warning("No CTstrict column - cannot add ", CT_CANON, "; TCR not attached?")
    return(obj)
  }
  raw <- as.character(md$CTstrict)
  nm  <- normalise_ct(raw)
  obj[[CT_CANON]] <- nm
  has <- !is.na(raw) & raw != "NA" & raw != ""
  n_changed <- sum(has & raw != nm, na.rm = TRUE)
  obj$ct_convention <- if (n_changed > 0) "CR10-style (empty D field, normalised)"
                       else                "CR9/legacy-style (.NA. D field)"
  if (verbose) cat(sprintf("  %s: %d clonotype strings | %d normalised (%.1f%%) | %s\n",
                           CT_CANON, sum(has), n_changed,
                           if (sum(has)) 100 * n_changed / sum(has) else 0,
                           obj$ct_convention[1]))
  obj
}

require_ct_canon <- function(obj, what = "object") {
  if (!CT_CANON %in% colnames(obj@meta.data)) {
    stop(what, " has no ", CT_CANON, " column. Run 12_tcr_attach.R first.\n",
         "Comparing raw CTstrict across Cell Ranger versions silently returns\n",
         "zero matches - see the normalise_ct() note in config.R.", call. = FALSE)
  }
  invisible(TRUE)
}

# ── Marker panel for whole PBMC ──────────────────────────────────────────────
# Lineage first, then state. Used by the resolution ladder, the dot plot and the
# annotation audit, so the panel is defined ONCE here rather than per script.
# Order matters: the dot plot reads as a lineage map top-to-bottom.
PBMC_PANEL <- list(
  `T cell`          = c("CD3D","CD3E","CD3G","TRAC"),
  `CD4`             = c("CD4","IL7R"),
  `CD8`             = c("CD8A","CD8B"),
  `Naive/CM`        = c("CCR7","SELL","LEF1","TCF7"),
  `Treg`            = c("FOXP3","IL2RA","CTLA4"),
  `MAIT`            = c("SLC4A10","KLRB1","TRAV1-2"),
  `gdT`             = c("TRDC","TRGC1"),
  `NK`              = c("GNLY","NKG7","KLRD1","KLRF1","NCAM1"),
  `Cytotoxic`       = c("GZMB","GZMK","PRF1","FGFBP2","CX3CR1"),
  `B`               = c("MS4A1","CD79A","TCL1A","IGHD"),
  `Plasma`          = c("MZB1","JCHAIN"),
  `Mono classical`  = c("CD14","LYZ","S100A8","S100A9","VCAN","FCN1"),
  `Mono non-class.` = c("FCGR3A","MS4A7","CDKN1C"),
  `cDC`             = c("FCER1A","CD1C","CLEC9A","CST3"),
  `pDC`             = c("LILRA4","CLEC4C","IRF7"),
  `Platelet`        = c("PPBP","PF4"),
  `HSPC`            = c("CD34"),
  `Exhaustion`      = c("PDCD1","TOX","LAG3","HAVCR2","TIGIT"),
  `Proliferation`   = c("MKI67","TOP2A","PCNA"),
  `Interferon`      = c("ISG15","IFI6","MX1","IFIT3")
)
PBMC_PANEL_FLAT <- unique(unlist(PBMC_PANEL, use.names = FALSE))

# Programme genes, for the identity-versus-programme audit. A cluster must never
# be NAMED after a programme (cycling / interferon / activation) — those are
# things a cell of any identity can do, and naming a cluster after one hides what
# the cell is and makes composition analysis incoherent. The audit rescores
# identity with these removed.
PROGRAMME_GENES <- list(
  cycling    = unique(c(Seurat::cc.genes.updated.2019$s.genes,
                        Seurat::cc.genes.updated.2019$g2m.genes)),
  interferon = c("ISG15","IFI6","MX1","MX2","IFIT1","IFIT2","IFIT3","OAS1","OAS2",
                 "OAS3","IFI44L","IFI44","RSAD2","STAT1","STAT2","IRF7","XAF1",
                 "EIF2AK2","HERC5","CMPK2","LY6E","SIGLEC1","EPSTI1","PARP9"),
  activation = c("CD69","FOS","FOSB","JUN","JUNB","EGR1","NR4A1","NR4A2","DUSP1",
                 "HSPA1A","HSPA1B","HSPB1","DNAJB1","ZFP36","IER2","KLF6")
)

# ── FastMNN integration without SeuratWrappers ───────────────────────────────
# CP018 used SeuratWrappers::FastMNNIntegration via IntegrateLayers.
# SeuratWrappers is GitHub-only and not installed here, but it is a thin wrapper
# around batchelor::fastMNN, which IS installed. Calling the engine directly
# avoids a source install and makes the actual computation explicit.
#
# Returns the object with a new DimReduc. `batch_col` is the technical grouping
# to absorb (library/well); biological grouping must NOT be passed here.
run_fastmnn <- function(seu, batch_col = "batch", features = NULL,
                        d = 50, assay = "RNA",
                        reduction.name = "integrated.mnn.rna",
                        reduction.key = "mnnRNA_") {
  stopifnot(batch_col %in% colnames(seu@meta.data))
  if (is.null(features)) features <- Seurat::VariableFeatures(seu)
  if (!length(features)) stop("No variable features; run FindVariableFeatures first.",
                              call. = FALSE)

  # fastMNN wants genes x cells log-normalised expression. Layers must be joined
  # so the matrix is a single block whose column order matches the batch vector.
  mat <- SeuratObject::LayerData(seu, assay = assay, layer = "data")
  mat <- mat[intersect(features, rownames(mat)), , drop = FALSE]

  batch <- as.character(seu@meta.data[[batch_col]])
  stopifnot(length(batch) == ncol(mat))

  mnn <- batchelor::fastMNN(mat, batch = batch, d = d)
  emb <- SingleCellExperiment::reducedDim(mnn, "corrected")
  rownames(emb) <- colnames(mat)
  colnames(emb) <- paste0(reduction.key, seq_len(ncol(emb)))

  seu[[reduction.name]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = reduction.key, assay = assay)

  # Per-gene "lost variance" per batch is fastMNN's own diagnostic for how much
  # signal the correction removed. Worth keeping rather than discarding.
  lv <- tryCatch(as.data.frame(S4Vectors::metadata(mnn)$merge.info$lost.var),
                 error = function(e) NULL)
  if (!is.null(lv)) attr(seu[[reduction.name]], "lost.var") <- lv
  seu
}

message("CHER config loaded | ", nrow(CHER_SAMPLES), " samples | ",
        length(RUN_DIRS), " OCM wells | SeuratExtend: ", HAS_SEURATEXTEND,
        " | Azimuth: ", HAS_AZIMUTH)
