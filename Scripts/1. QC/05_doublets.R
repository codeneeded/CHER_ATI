############################################################
# CHER ATI — 05: doublet detection and removal
#
# Single purpose: call doublets with scDblFinder and keep singlets.
#
# PER LIBRARY, NEVER POOLED. scDblFinder simulates doublets by combining cells
# from the same capture; pooling the eight samples would let it build doublets
# from cells that were never in the same GEM well, which is not a doublet that
# can exist. The two OCM wells make this sharper still — a "doublet" simulated
# across wells is meaningless.
#
# NOTE on OCM: on-chip multiplexing means cellranger has already removed
# cross-sample doublets it could identify by conflicting OCM barcode. What
# remains, and what this script targets, is WITHIN-sample doublets, which carry
# no barcode conflict and are invisible to the demultiplexer.
#
# Input : QS_CELLCYCLE
# Output: QS_SINGLETS + Analysis/QC/Doublets/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(qs2)
  library(SingleCellExperiment); library(scDblFinder)
})

dbl_dir <- file.path(QC_ROOT, "Doublets")
mkdirs(dbl_dir)

require_file(QS_CELLCYCLE, "cell-cycle object (run 04_cellcycle.R)")
seu <- qs_read(QS_CELLCYCLE)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

set.seed(1)

# -----------------------------
# 1) Per-sample doublet calls
# -----------------------------
counts <- LayerData(seu, assay = "RNA", layer = "counts")
stopifnot(ncol(counts) == ncol(seu))

score <- rep(NA_real_, ncol(seu))
class <- rep(NA_character_, ncol(seu))
names(score) <- names(class) <- colnames(seu)

for (s in CHER_SAMPLES$Sample_ID) {
  idx <- which(seu$Sample_ID == s)
  sce <- scDblFinder(counts[, idx, drop = FALSE])
  score[idx] <- sce$scDblFinder.score
  class[idx] <- as.character(sce$scDblFinder.class)
  message("  ", s, ": ", sum(sce$scDblFinder.class == "doublet"), " / ",
          length(idx), " doublets")
}
stopifnot(!anyNA(class))

seu$scDblFinder.score <- score
seu$scDblFinder.class <- class

# -----------------------------
# 2) Summary + diagnostics
# -----------------------------
dbl_summary <- do.call(rbind, lapply(CHER_SAMPLES$Sample_ID, function(s) {
  i <- seu$Sample_ID == s
  data.frame(Sample_ID = s,
             Group     = seu$Group[i][1],
             Well      = seu$Well[i][1],
             n_cells   = sum(i),
             n_doublet = sum(seu$scDblFinder.class[i] == "doublet"),
             pct_doublet = round(100 * mean(seu$scDblFinder.class[i] == "doublet"), 2))
}))
save_csv(dbl_summary, file.path(dbl_dir, "CHER_doublet_summary_source_data.csv"))
print(dbl_summary, row.names = FALSE)

# The duplicate healthy control is the technical reference: BC_A and BC_B are
# the same donor through the same pipeline, so the gap between their doublet
# rates is the noise floor for any rate difference reported between samples.
bc <- dbl_summary[dbl_summary$Group == "Healthy_control", ]
if (nrow(bc) == 2L) {
  cat(sprintf("\nBC control duplicate doublet rate: %.2f%% vs %.2f%% (gap %.2f pp)\n",
              bc$pct_doublet[1], bc$pct_doublet[2],
              abs(diff(bc$pct_doublet))))
}

save_png(file.path(dbl_dir, "UMAP_pre_by_doublet_class.png"),
         tidy_umap(DimPlot(seu, reduction = "umap.pre",
                           group.by = "scDblFinder.class",
                           cols = c(singlet = "#BBBBBB", doublet = "#E63946"),
                           pt.size = 0.3),
                   title = "scDblFinder calls, pre-integration embedding"),
         width = 8.5, height = 6.5)

save_png(file.path(dbl_dir, "UMAP_pre_doublet_by_sample.png"),
         tidy_umap(DimPlot(seu, reduction = "umap.pre",
                           group.by = "scDblFinder.class", split.by = "Sample_ID",
                           ncol = 4,
                           cols = c(singlet = "#BBBBBB", doublet = "#E63946"),
                           pt.size = 0.3)),
         width = 18, height = 9)

p_rate <- ggplot(dbl_summary, aes(x = Sample_ID, y = pct_doublet, fill = Group)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", pct_doublet)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% doublets", title = "scDblFinder doublet rate per library") +
  theme_cher(base_size = 14) + rotate_x(30)
save_png(file.path(dbl_dir, "Doublet_rate_per_library.png"), p_rate,
         width = 9, height = 5.5)

# Doublets should be enriched for high UMI/gene content; if they are not, the
# calls are suspect rather than informative.
save_csv(
  data.frame(
    metric = c("median nCount_RNA", "median nFeature_RNA"),
    doublet = c(median(seu$nCount_RNA[seu$scDblFinder.class == "doublet"]),
                median(seu$nFeature_RNA[seu$scDblFinder.class == "doublet"])),
    singlet = c(median(seu$nCount_RNA[seu$scDblFinder.class == "singlet"]),
                median(seu$nFeature_RNA[seu$scDblFinder.class == "singlet"]))),
  file.path(dbl_dir, "CHER_doublet_vs_singlet_depth.csv"))

# -----------------------------
# 3) Keep singlets
# -----------------------------
singlets <- subset(seu, cells = colnames(seu)[seu$scDblFinder.class == "singlet"])
cat("\nSinglets per library:\n"); print(table(singlets$Sample_ID))
cat("Total singlets:", ncol(singlets), "of", ncol(seu), "\n")

export_metadata(singlets, file.path(QC_ROOT, "CHER_postQC_cell_metadata.csv"))

qs_save(singlets, file = QS_SINGLETS)
message("Saved: ", QS_SINGLETS)
