############################################################
# CHER ATI — 12: attach TCR to the annotated object
#
# Single purpose: load the eight contig files, combine per sample, attach to the
# annotated object, and report coverage. Repertoire metrics are script 13.
#
# Two conventions, both from CD8_Longitudinal:
#
#   Clone size is counted within sample. proportion = FALSE and cloneSize in
#   absolute counts, so "expanded" means expanded in that library. Counting
#   within cluster x sample instead conflates expansion with cluster confinement
#   and can differ by tens of percentage points.
#
#   The canonical clonotype column is stamped once, on the object. Cell Ranger
#   10.1 writes an empty D-gene field where 9.0.1 wrote ".NA.", so raw CTstrict
#   strings never match across aligner versions and a cross-experiment
#   comparison silently returns zero. normalise_ct() is idempotent; downstream
#   scripts must read CT_CANON and never recompute it.
#
# Input : QS_ANNOTATED + per-sample filtered_contig_annotations.csv
# Output: QS_WITH_TCR + Analysis/VDJ/ coverage tables
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(scRepertoire); library(dplyr); library(ggplot2)
  library(qs2)
})

mkdirs(VDJ_ROOT)

require_file(QS_ANNOTATED, "annotated object (run 09_annotate.R)")
seu <- qs_read(QS_ANNOTATED)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)

# -----------------------------
# 1) Load and combine contigs
# -----------------------------
# loadContigs() is required — passing raw read.csv frames to combineTCR returns
# zero rows without an error.
contigs <- lapply(CHER_SAMPLES$Sample_ID, function(s)
  read.csv(sample_vdj(s), stringsAsFactors = FALSE))
names(contigs) <- CHER_SAMPLES$Sample_ID
contigs <- loadContigs(input = contigs, format = "10X")

# combineTCR is called without `samples`. Passing `samples` while leaving `ID`
# NULL makes it return zero rows silently (scRepertoire 2.2.0) — the same
# sample parsed correctly returns 1,606. Barcodes are prefixed here instead,
# which also guarantees the exact "<Sample_ID>_<barcode>" form that script 01
# gave the seurat cells.
combined <- combineTCR(contigs, removeNA = FALSE, removeMulti = FALSE,
                       filterMulti = TRUE)
names(combined) <- CHER_SAMPLES$Sample_ID
combined <- lapply(names(combined), function(s) {
  d <- combined[[s]]
  d$barcode <- paste0(s, "_", d$barcode)
  d$sample  <- s
  d
})
names(combined) <- CHER_SAMPLES$Sample_ID

# A key mismatch here attaches nothing and reads as absent TCR, so assert.
ov <- length(intersect(unlist(lapply(combined, function(d) d$barcode)), colnames(seu)))
cat("Barcode overlap with Seurat object: ", ov, "\n", sep = "")
stopifnot(ov > 1000)

# -----------------------------
# 2) Attach
# -----------------------------
seu <- combineExpression(combined, seu, cloneCall = "strict", chain = "both",
                         proportion = FALSE,
                         cloneSize = c(Single = 1, Small = 5, Medium = 20,
                                       Large = 100, Hyperexpanded = 500))

# -----------------------------
# 3) Canonical clonotype column
# -----------------------------
seu <- add_ct_canon(seu)

# -----------------------------
# 4) Coverage
# -----------------------------
seu$has_TCR <- !is.na(seu$CTstrict)

cov_sample <- seu@meta.data %>% group_by(Sample_ID) %>%
  summarise(n_cells = n(), n_TCR = sum(has_TCR),
            pct_TCR = round(100 * mean(has_TCR), 1),
            n_clonotypes = n_distinct(.data[[CT_CANON]][has_TCR]),
            .groups = "drop") %>% as.data.frame()
cov_sample <- merge(CHER_SAMPLES[, c("Sample_ID","Participant","Timepoint","Group","Well")],
                    cov_sample, by = "Sample_ID")
save_csv(cov_sample, file.path(VDJ_ROOT, "CHER_TCR_coverage_by_sample.csv"))
cat("\n=== TCR coverage by sample ===\n"); print(cov_sample, row.names = FALSE)

cov_ann <- seu@meta.data %>% group_by(Annotation) %>%
  summarise(n_cells = n(), n_TCR = sum(has_TCR),
            pct_TCR = round(100 * mean(has_TCR), 1), .groups = "drop") %>% as.data.frame()
save_csv(cov_ann, file.path(VDJ_ROOT, "CHER_TCR_coverage_by_annotation.csv"))
cat("\n=== TCR coverage by cell type ===\n"); print(cov_ann, row.names = FALSE)

cat("\nTotal cells with a receptor: ", sum(seu$has_TCR), " of ", ncol(seu),
    sprintf(" (%.1f%%)\n", 100 * mean(seu$has_TCR)), sep = "")
cat("Distinct clonotypes (", CT_CANON, "): ",
    dplyr::n_distinct(seu@meta.data[[CT_CANON]][seu$has_TCR]), "\n", sep = "")

qs_save(seu, file = QS_WITH_TCR)
message("Saved: ", QS_WITH_TCR)
