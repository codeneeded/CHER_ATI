############################################################
# CHER ATI — 09: annotation
#
# Single purpose: turn FastMNN res-0.6 cluster numbers into names, split the
# one mixed cluster, drop the one non-cell cluster, and flag Treg-like cells.
#
# Three decisions, each with the measurement behind it:
#
# 1. C9 is mixed, and is split. The matched TCR settles it: NK cells have no
#    receptor. C9 carries a productive TCR in 22.8% of cells, against >79% in
#    every true T cluster and 0.6% in C5 (NK). Subclustering separates it
#    cleanly — subclusters 9_0 and 9_5 are 59-60% TCR+ (CD8A 67-69%, CD8B
#    41-52%) and go to CD8 effector; 9_1/9_2/9_3/9_4 are 0-1.5% TCR+ with
#    KLRF1 up to 97% and go to NK.
#
# 2. C13 is dropped, not named (see DROP_CLUSTERS in config.R).
#
# 3. There is no Treg cluster, and this is a data limit rather than a
#    clustering failure. FOXP3 is detected in 227 of 20,050 cells (1.13%); only
#    75 co-express IL2RA and 58 are IL7R-negative. Subclustering C7/C8 at 0.4
#    and a dedicated CD4 re-embedding at 0.8 both failed to isolate a
#    FOXP3-high / IL7R-low population — the best CD4 subcluster reaches FOXP3
#    12.8% while keeping IL7R at 84.8%, which is not a treg. Tregs are therefore
#    carried as a Cell-Level flag, never as a cluster, and any treg claim must
#    be made on that flag with its count stated.
#
# Input : QS_INTEGRATED
# Output: QS_ANNOTATED + Analysis/Annotation/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(qs2); library(Matrix)
})

ANNO_ROOT <- file.path(ANALYSIS_ROOT, "Annotation")
mkdirs(ANNO_ROOT)
set.seed(1)

require_file(QS_INTEGRATED, "integrated object (run 06/07)")
seu <- qs_read(QS_INTEGRATED)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

res_col <- paste0("mnn_snn_res.", MNN_RES_FINAL)
stopifnot(res_col %in% colnames(seu@meta.data))
seu$cluster_base <- droplevels(factor(seu@meta.data[[res_col]]))

# -----------------------------
# 1) Split the mixed NK/T cluster
# -----------------------------
# FindSubCluster operates on the active idents, not on the column you name.
# Without this the saved idents (the res-1.6 grid from script 07) are used and
# "cluster 9" is a different 897 cells than res-0.6 cluster 9.
Idents(seu) <- "cluster_base"
seu <- FindSubCluster(seu, cluster = "9", graph.name = "mnn_snn",
                      subcluster.name = "sub_9", resolution = 0.4, algorithm = 1)

# Assign each subcluster by measured TCR carriage rather than by index — the
# subcluster numbering is not stable across runs. The two populations are far
# apart (true T clusters >79% TCR+, C5 NK 0.6%), so any threshold between them
# gives the same split; 40% is used.
TCR_SPLIT_THRESHOLD <- 40
has_tcr <- tcr_cells(colnames(seu))
sub9    <- seu$sub_9[seu$cluster_base == "9"]
frac9   <- tapply(has_tcr[seu$cluster_base == "9"], droplevels(factor(sub9)),
                  function(x) 100 * mean(x))
cat("\nC9 subcluster TCR carriage (%):\n"); print(round(frac9, 1))

C9_TO_T  <- names(frac9)[frac9 >= TCR_SPLIT_THRESHOLD]
C9_TO_NK <- names(frac9)[frac9 <  TCR_SPLIT_THRESHOLD]
cat("  -> CD8 effector: ", paste(C9_TO_T, collapse = ", "),
    "\n  -> NK: ",         paste(C9_TO_NK, collapse = ", "), "\n", sep = "")
stopifnot(length(C9_TO_T) > 0, length(C9_TO_NK) > 0)

# -----------------------------
# 2) Apply labels
# -----------------------------
MAP <- c(
  "0"  = "Naive B cell",
  "1"  = "CD4 naive T",
  "2"  = "CD8 effector/TEMRA",
  "3"  = "Classical monocyte",
  "4"  = "CD8 naive/CM T",
  "5"  = "NK cell",
  "6"  = "Classical monocyte (inflammatory)",
  "7"  = "CD4 effector memory T",
  "8"  = "CD4 central memory T",
  "10" = "Non-classical monocyte",
  "11" = "B cell (unresolved)",
  "12" = "Plasma cell",
  "14" = "CD8 effector T (cycling)",
  "15" = "pDC")

ann <- rep(NA_character_, ncol(seu))
base <- as.character(seu$cluster_base)
for (k in names(MAP)) ann[base == k] <- MAP[[k]]
ann[seu$sub_9 %in% C9_TO_T]  <- "CD8 effector/TEMRA"
ann[seu$sub_9 %in% C9_TO_NK] <- "NK cell"
ann[base %in% DROP_CLUSTERS] <- "DROP"

stopifnot(!anyNA(ann))
seu$Annotation <- ann

# -----------------------------
# 3) Treg-like cell flag (Not a cluster)
# -----------------------------
cts <- LayerData(seu, assay = "RNA", layer = "counts")
seu$FOXP3_pos  <- cts["FOXP3", ] > 0
seu$Treg_like  <- seu$FOXP3_pos &
                  (cts["IL2RA", ] > 0 | cts["CTLA4", ] > 0) &
                  cts["IL7R", ] == 0
cat("FOXP3+ cells: ", sum(seu$FOXP3_pos),
    " | Treg-like (FOXP3+ & (IL2RA+|CTLA4+) & IL7R-): ", sum(seu$Treg_like), "\n", sep = "")

# -----------------------------
# 4) Drop the non-cell cluster
# -----------------------------
n_before <- ncol(seu)
seu <- subset(seu, cells = colnames(seu)[seu$Annotation != "DROP"])
seu$Annotation <- factor(seu$Annotation, levels = ANNO_ORDER)
stopifnot(!anyNA(seu$Annotation))
cat("Dropped ", n_before - ncol(seu), " cells (cluster ",
    paste(DROP_CLUSTERS, collapse = ","), "); retained ", ncol(seu), "\n", sep = "")

# -----------------------------
# 5) Tables + figures
# -----------------------------
comp <- as.data.frame(table(Annotation = seu$Annotation))
comp$pct <- round(100 * comp$Freq / sum(comp$Freq), 2)
names(comp)[2] <- "n_cells"
save_csv(comp, file.path(ANNO_ROOT, "CHER_annotation_cell_counts.csv"))
print(comp, row.names = FALSE)

for (v in c("Sample_ID", "Participant", "Timepoint", "Group", "Well", "Phase")) {
  save_csv(composition_table(seu, "Annotation", v),
           file.path(ANNO_ROOT, paste0("CHER_annotation_composition_by_", v, ".csv")))
}

# House idiom: boxed labels on the embedding, corner arrows, square canvas (the
# arrow arms only render equal on a square drawing area).
save_png(file.path(ANNO_ROOT, "CHER_UMAP_annotated.png"),
         umap_house(seu, group = "Annotation", reduction = "umap.mnn",
                    cols = ANNO_COLS, title = "CHER ATI PBMC"),
         width = 11, height = 11)

# Split panels carry a legend instead of per-panel boxes: repelling 14 boxed
# labels inside each of three narrow facets collides no matter the force.
save_png(file.path(ANNO_ROOT, "CHER_UMAP_annotated_by_Timepoint.png"),
         umap_house(seu, group = "Annotation", reduction = "umap.mnn",
                    cols = ANNO_COLS, split.by = "Timepoint",
                    label = FALSE, legend = TRUE,
                    title = "Annotation by timepoint"),
         width = 16, height = 7)

panel <- intersect(PBMC_PANEL_FLAT, rownames(seu))
save_png(file.path(ANNO_ROOT, "CHER_DotPlot_annotated.png"),
         DotPlot(seu, features = panel, group.by = "Annotation", dot.scale = 5) +
           scale_colour_gradient2(low = "#2166AC", mid = "grey94", high = "#B2182B",
                                  midpoint = 0, name = "scaled\nexpression") +
           labs(x = NULL, y = NULL) + theme_cher(base_size = 11) +
           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                            size = 8, face = "italic")),
         width = 17, height = 7)

# Per-gene feature plots and grouped violins for the lineage-Defining genes —
# the panels a reviewer uses to check a call cell-by-cell rather than trusting a
# cluster mean. Source CSV per figure, per house convention.
LINEAGE_GENES <- list(
  `T cell`   = c("CD3E","CD4","CD8A","CD8B"),
  `T state`  = c("CCR7","SELL","IL7R","CCR4"),
  `NK`       = c("KLRD1","KLRF1","GNLY","NCAM1"),
  `B`        = c("MS4A1","CD79A","IGHD","FCRL1"),
  `Plasma`   = c("MZB1","JCHAIN","CD38","TNFRSF13B"),
  `Myeloid`  = c("CD14","LYZ","FCGR3A","CDKN1C"),
  `DC`       = c("LILRA4","CLEC4C","IRF7","FCER1A"),
  `Programme`= c("MKI67","TOP2A","CCL2","IDO1"))

fp_dir <- file.path(ANNO_ROOT, "feature_plots")
vl_dir <- file.path(ANNO_ROOT, "violins")
mkdirs(fp_dir, vl_dir)

for (grp in names(LINEAGE_GENES)) {
  gs <- intersect(LINEAGE_GENES[[grp]], rownames(seu))
  if (!length(gs)) next
  tag <- gsub("[^A-Za-z0-9]+", "_", grp)

  save_png(file.path(fp_dir, paste0("FeaturePlot_", tag, ".png")),
           FeaturePlot(seu, features = gs, reduction = "umap.mnn", ncol = 2,
                       pt.size = 0.25, order = TRUE) &
             theme_cher(base_size = 11) &
             ggplot2::labs(x = "UMAP 1", y = "UMAP 2"),
           width = 10, height = 4.4 * ceiling(length(gs) / 2))

  save_png(file.path(vl_dir, paste0("Vln_", tag, ".png")),
           VlnPlot(seu, features = gs, group.by = "Annotation", pt.size = 0,
                   cols = ANNO_COLS, ncol = 1) &
             theme_cher(base_size = 11) & rotate_x(35) &
             ggplot2::labs(x = NULL) &
             ggplot2::theme(legend.position = "none"),
           width = 10, height = 3.0 * length(gs))

  # Source data: percent detected and mean normalised expression per label.
  dat <- LayerData(seu, assay = "RNA", layer = "data")[gs, , drop = FALSE]
  cnt <- LayerData(seu, assay = "RNA", layer = "counts")[gs, , drop = FALSE]
  src <- do.call(rbind, lapply(levels(seu$Annotation), function(k) {
    i <- which(seu$Annotation == k)
    data.frame(Annotation = k, gene = gs,
               pct_detected = round(100 * Matrix::rowMeans(cnt[, i, drop = FALSE] > 0), 2),
               mean_expression = round(Matrix::rowMeans(dat[, i, drop = FALSE]), 4),
               row.names = NULL)
  }))
  save_csv(src, file.path(fp_dir, paste0("FeaturePlot_", tag, "_source_data.csv")))
}

export_metadata(seu, file.path(ANNO_ROOT, "CHER_annotated_cell_metadata.csv"),
                extra_cols = c("Annotation", "cluster_base", "Treg_like", "FOXP3_pos"))

qs_save(seu, file = QS_ANNOTATED)
message("Saved: ", QS_ANNOTATED)
