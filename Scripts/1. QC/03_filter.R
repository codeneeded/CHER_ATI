############################################################
# CHER ATI — 03: apply QC thresholds
#
# Single purpose: filter cells and genes, report what each sample lost, and
# redraw the diagnostics post-filter using the SAME code as 02 (run_qc_plots
# lives in config.R precisely so the two cannot drift).
#
# Thresholds are in config.R. Two differ from CP018 and the reason is recorded
# there with the measurement behind it: the ribosomal floor is disabled, and
# percent_plat is PF4/PPBP rather than PECAM1|PF4.
#
# Input : QS_RAW
# Output: QS_FILTERED + Analysis/QC/Post-QC/ + retention tables
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(ggplot2)
  library(patchwork); library(scales); library(qs2); library(scCustomize)
})

require_file(QS_RAW, "merged raw object (run 01_build_objects.R)")
seu <- qs_read(QS_RAW)

post_root <- file.path(QC_ROOT, "Post-QC")
mkdirs(post_root)

raw_counts <- read.csv(file.path(QC_ROOT, "CHER_cells_per_sample_raw.csv"),
                       stringsAsFactors = FALSE)

# -----------------------------
# 1) Cell filtering
# -----------------------------
# QC_MIN_RIBO is 0 (disabled), so the ribo term is a no-op; it is written out
# explicitly rather than dropped so that re-enabling the floor is a one-line
# change in config.R and the applied rule stays visible here.
keep <- (seu$nCount_RNA       >= QC_MIN_COUNT)   &
        (seu$nFeature_RNA     >= QC_MIN_FEATURE) &
        (seu$log10GenesPerUMI >  QC_MIN_COMPLEX) &
        (seu$percent_mito     <  QC_MAX_MITO)    &
        (seu$percent_ribo     >  QC_MIN_RIBO)    &
        (seu$percent_hb       <  QC_MAX_HB)      &
        (seu$percent_plat     <  QC_MAX_PLAT)

message("Cells passing: ", sum(keep), " / ", length(keep),
        " (", round(100 * mean(keep), 1), "%)")
filtered <- subset(seu, cells = colnames(seu)[keep])

# -----------------------------
# 2) Gene filtering
# -----------------------------
# Genes detected in fewer than QC_MIN_CELLS_PER_GENE cells cannot support a
# marker test and inflate the multiple-testing burden.
counts_mat <- LayerData(filtered, assay = "RNA", layer = "counts")
stopifnot(ncol(counts_mat) == ncol(filtered))        # guards layer misalignment
genes_keep <- Matrix::rowSums(counts_mat > 0) >= QC_MIN_CELLS_PER_GENE
message("Genes passing: ", sum(genes_keep), " / ", length(genes_keep))
filtered <- subset(filtered, features = names(genes_keep)[genes_keep])

# -----------------------------
# 3) Retention summary
# -----------------------------
post <- as.data.frame(table(filtered$Sample_ID), responseName = "n_cells_postQC")
names(post)[1] <- "Sample_ID"
filt_summary <- merge(raw_counts, post, by = "Sample_ID", all = TRUE)
filt_summary$n_cells_postQC[is.na(filt_summary$n_cells_postQC)] <- 0L
filt_summary$pct_retained <- round(100 * filt_summary$n_cells_postQC /
                                     filt_summary$n_cells_raw, 1)
save_csv(filt_summary, file.path(QC_ROOT, "CHER_QC_filtering_summary.csv"))
print(filt_summary[, c("Sample_ID","Group","Well","n_cells_raw",
                       "n_cells_postQC","pct_retained")], row.names = FALSE)

# The thresholds actually applied, written beside the results so a figure can
# never be read against the wrong cutoffs.
thresholds <- data.frame(
  parameter = c("QC_MIN_COUNT","QC_MIN_FEATURE","QC_MIN_COMPLEX","QC_MAX_MITO",
                "QC_MIN_RIBO","QC_MAX_HB","QC_MAX_PLAT","QC_MIN_CELLS_PER_GENE"),
  value = c(QC_MIN_COUNT, QC_MIN_FEATURE, QC_MIN_COMPLEX, QC_MAX_MITO,
            QC_MIN_RIBO, QC_MAX_HB, QC_MAX_PLAT, QC_MIN_CELLS_PER_GENE),
  note = c("", "", "", "",
           "0 = ribosomal floor DISABLED (see config.R)",
           "pattern anchored to globin genes only",
           "PF4/PPBP, not PECAM1|PF4", "")
)
save_csv(thresholds, file.path(QC_ROOT, "CHER_QC_thresholds_applied.csv"))

# -----------------------------
# 4) Post-filter diagnostics
# -----------------------------
for (g in c("Sample_ID", "Participant", "Timepoint", "Group", "Well")) {
  message("\n== post-QC plots by ", g, " ==")
  run_qc_plots(filtered, g, file.path(post_root, g))
}

# -----------------------------
# 5) Checkpoint
# -----------------------------
cat("\nPost-QC total:", ncol(filtered), "cells |", nrow(filtered), "genes\n")
qs_save(filtered, file = QS_FILTERED)
message("Saved: ", QS_FILTERED)
