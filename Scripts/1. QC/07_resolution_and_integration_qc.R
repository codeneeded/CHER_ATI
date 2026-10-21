############################################################
# CHER ATI — 07: resolution sweep and integration quality
#
# Single purpose: cluster across a resolution grid on both integrated graphs,
# render clustree, and measure whether integration absorbed the well effect
# Without erasing biology. No resolution is chosen here — MNN_RES_FINAL and
# CCA_RES_FINAL stay NA in config.R until these figures have been reviewed.
# Picking a resolution is an annotation decision, not a QC one.
#
# Input : QS_INTEGRATED
# Output: Analysis/Integration/ clustree + mixing tables (object re-saved with
#         the clustering columns attached)
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(qs2); library(clustree)
})

mkdirs(INT_ROOT)
set.seed(1)

require_file(QS_INTEGRATED, "integrated object (run 06_integrate.R)")
seu <- qs_read(QS_INTEGRATED)
DefaultAssay(seu) <- "RNA"

RES_GRID <- seq(0.2, 1.6, by = 0.2)

# -----------------------------
# 1) Clustering grids + clustree
# -----------------------------
for (spec in list(c(red = "integrated.mnn.rna", graph = "mnn_snn", tag = "MNN"),
                  c(red = "integrated.cca.rna", graph = "cca_snn", tag = "CCA"))) {
  seu <- FindNeighbors(seu, reduction = spec[["red"]], dims = 1:30,
                       graph.name = spec[["graph"]], verbose = FALSE)
  for (r in RES_GRID) {
    seu <- FindClusters(seu, graph.name = spec[["graph"]], resolution = r,
                        algorithm = 1, verbose = FALSE)
  }
  # The clustree legend stack runs tall; the explicit bottom margin gives it
  # room rather than growing the whole canvas.
  save_png(file.path(INT_ROOT, paste0("CHER_clustree_", spec[["tag"]], ".png")),
           clustree(seu, prefix = paste0(spec[["graph"]], "_res.")) +
             ggplot2::theme(legend.box.margin = ggplot2::margin(b = 30)),
           width = 11, height = 13)

  n_tab <- data.frame(
    resolution = RES_GRID,
    n_clusters = vapply(RES_GRID, function(r)
      nlevels(factor(seu@meta.data[[paste0(spec[["graph"]], "_res.", r)]])),
      integer(1)))
  save_csv(n_tab, file.path(INT_ROOT,
           paste0("CHER_resolution_to_cluster_count_", spec[["tag"]], ".csv")))
  cat("\n--", spec[["tag"]], "resolution -> clusters --\n"); print(n_tab, row.names = FALSE)
}

# -----------------------------
# 2) Integration quality: did the well effect get absorbed?
# -----------------------------
# Two things have to be true at once. (a) Mixing: clusters should not be
# dominated by one well beyond what cell numbers imply. (b) Preservation:
# participant and timepoint structure must still exist. A correction that
# achieves (a) by destroying (b) is worse than no correction.
#
# Measured at a mid-grid resolution, reported for both embeddings so the choice
# between them is made on evidence.
RES_EVAL <- 0.8

well_frac <- prop.table(table(seu$Well))

mixing <- do.call(rbind, lapply(c("mnn_snn", "cca_snn"), function(g) {
  cl <- seu@meta.data[[paste0(g, "_res.", RES_EVAL)]]
  tb <- table(cl, seu$Well)
  pr <- prop.table(tb, 1)
  data.frame(
    embedding   = ifelse(g == "mnn_snn", "FastMNN", "CCA"),
    cluster     = rownames(tb),
    n_cells     = as.integer(rowSums(tb)),
    pct_34349   = round(100 * pr[, "34349"], 1),
    pct_33445   = round(100 * pr[, "33445"], 1),
    # Deviation from the dataset-wide well split: 0 = perfectly mixed.
    well_skew   = round(100 * abs(pr[, "34349"] - well_frac[["34349"]]), 1),
    n_participants = as.integer(apply(table(cl, seu$Participant) > 0, 1, sum))[
                       match(rownames(tb), rownames(table(cl, seu$Participant)))],
    row.names = NULL)
}))
save_csv(mixing, file.path(INT_ROOT,
         paste0("CHER_cluster_well_mixing_res", RES_EVAL, "_source_data.csv")))

cat("\n=== well mixing at resolution ", RES_EVAL, " ===\n", sep = "")
cat(sprintf("dataset-wide well split: 34349 %.1f%% / 33445 %.1f%%\n",
            100 * well_frac[["34349"]], 100 * well_frac[["33445"]]))
for (e in unique(mixing$embedding)) {
  m <- mixing[mixing$embedding == e, ]
  cat(sprintf("\n%-8s clusters=%2d | median well skew %.1f pp | max %.1f pp | clusters >90%% one well: %d\n",
              e, nrow(m), median(m$well_skew), max(m$well_skew),
              sum(m$pct_34349 > 90 | m$pct_33445 > 90)))
  cat(sprintf("         clusters containing all 4 donors: %d of %d\n",
              sum(m$n_participants == 4), nrow(m)))
}

# Composition tables that downstream composition analysis will start from.
for (v in c("Participant", "Timepoint", "Group", "Well", "Phase")) {
  save_csv(composition_table(seu, paste0("mnn_snn_res.", RES_EVAL), v),
           file.path(INT_ROOT, paste0("CHER_MNN_res", RES_EVAL,
                                      "_composition_by_", v, ".csv")))
}

# -----------------------------
# 3) Cluster UMAPs at the evaluation resolution
# -----------------------------
for (spec in list(c(red = "umap.mnn", g = "mnn_snn", tag = "MNN"),
                  c(red = "umap.cca", g = "cca_snn", tag = "CCA"))) {
  save_png(file.path(INT_ROOT, paste0("CHER_", spec[["tag"]], "_clusters_res",
                                      RES_EVAL, ".png")),
           umap_labels(tidy_umap(
             DimPlot(seu, reduction = spec[["red"]],
                     group.by = paste0(spec[["g"]], "_res.", RES_EVAL),
                     pt.size = 0.3),
             title = paste0(spec[["tag"]], " clusters, resolution ", RES_EVAL),
             legend = FALSE)),
           width = 9, height = 7.5)
}

# -----------------------------
# 4) Re-save with clustering attached
# -----------------------------
export_metadata(seu, file.path(INT_ROOT, "CHER_integrated_cell_metadata.csv"),
                extra_cols = grep("_res\\.", colnames(seu@meta.data), value = TRUE))
qs_save(seu, file = QS_INTEGRATED)
message("Re-saved with clustering: ", QS_INTEGRATED)
