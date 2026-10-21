############################################################
# CHER ATI — 02: pre-filter QC diagnostics
#
# Single purpose: draw the distributions that justify the thresholds, BEFORE
# any cell is removed. The plotting code lives in config.R (run_qc_plots) so
# that 03_filter.R renders the post-filter panels with identical code — a copy
# in each script could drift and make the comparison meaningless.
#
# Input : QS_RAW
# Output: Analysis/QC/Pre-QC/<grouping>/ figures + source CSVs
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(scales); library(qs2); library(scCustomize)
})

require_file(QS_RAW, "merged raw object (run 01_build_objects.R)")
seu <- qs_read(QS_RAW)

pre_root <- file.path(QC_ROOT, "Pre-QC")
mkdirs(pre_root)

# Sample_ID is the unit of QC. Participant/Timepoint/Group/Well are the
# groupings a threshold decision has to be fair across — in particular Well,
# because the two OCM wells differ systematically in sequencing depth and a
# single global threshold could remove cells from one well preferentially.
group_vars <- c("Sample_ID", "Participant", "Timepoint", "Group", "Well")

for (g in group_vars) {
  message("\n== pre-QC plots by ", g, " ==")
  run_qc_plots(seu, g, file.path(pre_root, g))
}

# -----------------------------
# Threshold cost, per sample
# -----------------------------
# How many cells EACH criterion would remove on its own. This is the table the
# threshold decision should be made from: a cutoff that is reasonable overall
# can still gut one sample, and that is invisible in a pooled density plot.
md <- seu@meta.data
crit <- list(
  min_count   = md$nCount_RNA       <  QC_MIN_COUNT,
  min_feature = md$nFeature_RNA     <  QC_MIN_FEATURE,
  complexity  = md$log10GenesPerUMI <= QC_MIN_COMPLEX,
  max_mito    = md$percent_mito     >= QC_MAX_MITO,
  min_ribo    = md$percent_ribo     <= QC_MIN_RIBO,
  max_hb      = md$percent_hb       >= QC_MAX_HB,
  max_plat    = md$percent_plat     >= QC_MAX_PLAT
)

cost <- do.call(rbind, lapply(names(crit), function(k) {
  tb <- tapply(crit[[k]], md$Sample_ID, sum)
  n  <- tapply(crit[[k]], md$Sample_ID, length)
  data.frame(criterion = k, Sample_ID = names(tb),
             n_failing = as.integer(tb),
             pct_failing = round(100 * as.integer(tb) / as.integer(n), 2),
             row.names = NULL)
}))
cost_wide <- reshape(cost[, c("criterion","Sample_ID","pct_failing")],
                     idvar = "Sample_ID", timevar = "criterion",
                     direction = "wide")
names(cost_wide) <- sub("^pct_failing\\.", "pct_", names(cost_wide))

# Cells failing ANY criterion — the actual retention if all are applied.
md$fails_any <- Reduce(`|`, crit)
anyfail <- as.data.frame(table(md$Sample_ID, md$fails_any))
names(anyfail) <- c("Sample_ID", "fails_any", "n_cells")

save_csv(cost,      file.path(QC_ROOT, "CHER_threshold_cost_long.csv"))
save_csv(cost_wide, file.path(QC_ROOT, "CHER_threshold_cost_per_sample.csv"))
save_csv(anyfail,   file.path(QC_ROOT, "CHER_threshold_fails_any.csv"))

cat("\n--- % of cells failing each criterion, per sample ---\n")
print(cost_wide, row.names = FALSE)
cat("\n--- cells failing at least one criterion ---\n")
print(round(100 * tapply(md$fails_any, md$Sample_ID, mean), 2))
