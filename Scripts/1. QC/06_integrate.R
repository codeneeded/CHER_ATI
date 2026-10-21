############################################################
# CHER ATI — 06: integration (CCA + FastMNN)
#
# Single purpose: produce the integrated embeddings. Clustering resolution and
# integration-quality assessment are script 07, so this one can be re-run
# without redoing the sweep and vice versa.
#
# What counts as batch here. The technical unit is the GEM Well, and it matters more
# than usual: the two wells differ systematically in depth (median UMI
# 6,842-10,279 in 33445 against 5,176-6,393 in 34349). Library and well are the
# same partition for GEX, so batch = Well.
#
# What must survive. Participant and timepoint are biological. Each
# participant's B4 and interrupt samples sit in the same well, so the paired
# contrast is not a batch contrast and integration cannot manufacture it.
# Healthy control vs ATI is confounded with well — integration reduces but does
# not remove that, and no amount of correction makes it a clean comparison.
#
# FastMNN is run through run_fastmnn() in config.R, which calls
# batchelor::fastMNN directly rather than SeuratWrappers::FastMNNIntegration
# (SeuratWrappers is GitHub-only and not installed; the engine is the same).
#
# Input : QS_SINGLETS
# Output: QS_INTEGRATED + Analysis/Integration/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(qs2); library(batchelor); library(SingleCellExperiment)
})

mkdirs(INT_ROOT)
set.seed(1)

require_file(QS_SINGLETS, "singlet object (run 05_doublets.R)")
seu <- qs_read(QS_SINGLETS)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

seu$batch <- seu$Well

# -----------------------------
# 1) Unintegrated baseline
# -----------------------------
# Split layers so IntegrateLayers has per-batch layers to work from. This is
# also what makes the unintegrated UMAP an honest baseline: same features, same
# dims, only the correction differs.
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$batch)

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000,
                            verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap.unintegrated", verbose = FALSE)

# shuffle = TRUE is required on any panel coloured by an overlapping grouping.
# DimPlot draws groups in factor order, so without it the last level is painted
# over the first: at a 45/55 well split the minority well reads as a scattered
# handful of cells and the figure implies integration failed when it did not.
p_un <- wrap_plots(
  tidy_umap(DimPlot(seu, reduction = "umap.unintegrated", group.by = "Well",
                    cols = WELL_COLORS, pt.size = 0.3,
                    shuffle = TRUE, seed = 1), title = "Well (batch)"),
  tidy_umap(DimPlot(seu, reduction = "umap.unintegrated", group.by = "Participant",
                    pt.size = 0.3, shuffle = TRUE, seed = 1), title = "Participant"),
  tidy_umap(DimPlot(seu, reduction = "umap.unintegrated", group.by = "Timepoint",
                    cols = TP_COLORS, pt.size = 0.3,
                    shuffle = TRUE, seed = 1), title = "Timepoint"),
  tidy_umap(DimPlot(seu, reduction = "umap.unintegrated", group.by = "Group",
                    cols = GROUP_COLORS, pt.size = 0.3,
                    shuffle = TRUE, seed = 1), title = "Group"),
  ncol = 2)
save_png(file.path(INT_ROOT, "CHER_RNA_Unintegrated_UMAPs.png"), p_un,
         width = 15, height = 11)

# -----------------------------
# 2) CCA integration
# -----------------------------
seu <- IntegrateLayers(object = seu, method = CCAIntegration,
                       orig.reduction = "pca", assay = "RNA",
                       new.reduction = "integrated.cca.rna", verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "integrated.cca.rna",
               reduction.name = "umap.cca", verbose = FALSE)

# -----------------------------
# 3) FastMNN integration
# -----------------------------
# fastMNN needs one joined matrix whose column order matches the batch vector,
# so layers are rejoined first. Variable features are preserved across the
# rejoin, so both integrations use the same 3,000 features.
hvg <- VariableFeatures(seu)
seu <- JoinLayers(seu)
VariableFeatures(seu) <- hvg

seu <- run_fastmnn(seu, batch_col = "batch", features = hvg, d = 50,
                   reduction.name = "integrated.mnn.rna",
                   reduction.key = "mnnRNA_")
seu <- RunUMAP(seu, dims = 1:30, reduction = "integrated.mnn.rna",
               reduction.name = "umap.mnn", verbose = FALSE)

lv <- attr(seu[["integrated.mnn.rna"]], "lost.var")
if (!is.null(lv)) {
  save_csv(lv, file.path(INT_ROOT, "CHER_fastMNN_lost_variance.csv"),
           row.names = TRUE)
  cat("\nfastMNN lost variance per merge step:\n"); print(lv)
}

# -----------------------------
# 4) Side-by-side comparison
# -----------------------------
for (grp in c("Well", "Participant", "Timepoint", "Group")) {
  cols <- switch(grp, Well = WELL_COLORS, Timepoint = TP_COLORS,
                 Group = GROUP_COLORS, NULL)
  p <- wrap_plots(
    tidy_umap(DimPlot(seu, reduction = "umap.unintegrated", group.by = grp,
                      cols = cols, pt.size = 0.3, shuffle = TRUE, seed = 1),
              title = "Unintegrated"),
    tidy_umap(DimPlot(seu, reduction = "umap.cca", group.by = grp,
                      cols = cols, pt.size = 0.3, shuffle = TRUE, seed = 1),
              title = "CCA"),
    tidy_umap(DimPlot(seu, reduction = "umap.mnn", group.by = grp,
                      cols = cols, pt.size = 0.3, shuffle = TRUE, seed = 1),
              title = "FastMNN"),
    ncol = 3) +
    patchwork::plot_annotation(
      title = paste0("Integration comparison, coloured by ", grp))
  save_png(file.path(INT_ROOT, paste0("CHER_Integration_comparison_", grp, ".png")),
           p, width = 19, height = 6.5)
}

# -----------------------------
# 5) Checkpoint
# -----------------------------
cat("\nReductions on object:", paste(Reductions(seu), collapse = ", "), "\n")
cat("Cells:", ncol(seu), "| genes:", nrow(seu), "\n")
qs_save(seu, file = QS_INTEGRATED)
message("Saved: ", QS_INTEGRATED)
