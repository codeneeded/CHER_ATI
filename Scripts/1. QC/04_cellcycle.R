############################################################
# CHER ATI — 04: normalisation, pre-integration embedding, cell cycle
#
# Single purpose: get a first embedding and score cell cycle phase. The UMAP
# produced here (umap.pre) is pre-Integration and exists so that doublet calls
# in script 05 can be shown in a 2-D space, and so the batch structure is
# visible before correction. It is not the embedding any result should use.
#
# Azimuth reference annotation runs only if the package is installed; it is a
# convenience for orientation, not a result, so a missing package must not fail
# the pipeline (see HAS_AZIMUTH in config.R).
#
# Input : QS_FILTERED
# Output: QS_CELLCYCLE + Analysis/QC/Cell_Cycle/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(qs2)
})

cc_dir <- file.path(QC_ROOT, "Cell_Cycle")
mkdirs(cc_dir)

require_file(QS_FILTERED, "post-QC object (run 03_filter.R)")
seu <- qs_read(QS_FILTERED)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

set.seed(1)

# -----------------------------
# 1) Standard pre-integration processing
# -----------------------------
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000,
                            verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, reduction = "pca",
               reduction.name = "umap.pre", verbose = FALSE)

# -----------------------------
# 2) Cell cycle scoring
# -----------------------------
s.use   <- intersect(Seurat::cc.genes.updated.2019$s.genes,   rownames(seu))
g2m.use <- intersect(Seurat::cc.genes.updated.2019$g2m.genes, rownames(seu))
message("Cell-cycle genes found: S ", length(s.use), " / G2M ", length(g2m.use))

seu <- CellCycleScoring(seu, s.features = s.use, g2m.features = g2m.use,
                        set.ident = FALSE)
seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))

# -----------------------------
# 3) Diagnostics
# -----------------------------
save_png(file.path(cc_dir, "UMAP_pre_by_Phase.png"),
         tidy_umap(DimPlot(seu, reduction = "umap.pre", group.by = "Phase",
                           cols = PHASE_COLS, pt.size = 0.3),
                   title = "Cell cycle phase, pre-integration embedding"),
         width = 8, height = 6.5)

for (g in c("Sample_ID", "Timepoint", "Group", "Well")) {
  save_png(file.path(cc_dir, paste0("UMAP_pre_by_", g, ".png")),
           tidy_umap(DimPlot(seu, reduction = "umap.pre", group.by = g,
                             pt.size = 0.3),
                     title = paste0("Pre-integration embedding by ", g)),
           width = 8.5, height = 6.5)
  save_png(file.path(cc_dir, paste0("CellCycle_Score_Vln_by_", g, ".png")),
           VlnPlot(seu, features = c("S.Score", "G2M.Score"), group.by = g,
                   pt.size = 0) &
             theme_cher(base_size = 13) & rotate_x(30) &
             ggplot2::theme(legend.position = "none"),
           width = 11, height = 5)
}

# Phase composition — a real readout, not only a QC metric: proliferating
# fraction can genuinely differ between B4 and interrupt.
phase_tab <- composition_table(seu, "Phase", "Sample_ID")
save_csv(phase_tab, file.path(cc_dir, "CHER_phase_by_sample_source_data.csv"))
print(as.data.frame.matrix(table(seu$Sample_ID, seu$Phase)))

p_phase <- ggplot(phase_tab, aes(x = Sample_ID, y = pct_of_sample, fill = Phase)) +
  geom_col(colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = PHASE_COLS) +
  labs(x = NULL, y = "% of sample", fill = "Phase",
       title = "Cell cycle phase composition") +
  theme_cher(base_size = 14) + rotate_x(30)
save_png(file.path(cc_dir, "Phase_composition.png"), p_phase, width = 9, height = 5.5)

ridge_feats <- intersect(c("PCNA", "TOP2A", "MCM6", "MKI67"), rownames(seu))
if (length(ridge_feats)) {
  save_png(file.path(cc_dir, "Ridge_CellCycle_Markers.png"),
           RidgePlot(seu, features = ridge_feats, ncol = 2) &
             theme_cher(base_size = 12),
           width = 10, height = 7)
}

# -----------------------------
# 4) Optional: Azimuth reference annotation
# -----------------------------
if (HAS_AZIMUTH) {
  message("Running Azimuth (pbmcref) for orientation labels")
  seu <- Azimuth::RunAzimuth(seu, reference = "pbmcref")
  save_png(file.path(cc_dir, "UMAP_pre_Azimuth_l2.png"),
           umap_labels(tidy_umap(DimPlot(seu, reduction = "umap.pre",
                                         group.by = "predicted.celltype.l2"),
                                 title = "Azimuth pbmcref, predicted.celltype.l2",
                                 legend = FALSE)),
           width = 9, height = 7.5)
  save_csv(composition_table(seu, "predicted.celltype.l2", "Sample_ID"),
           file.path(cc_dir, "CHER_azimuth_l2_by_sample_source_data.csv"))
} else {
  message("Azimuth not installed - skipping reference annotation. ",
          "Cluster-based annotation happens after integration regardless; ",
          "this step only provides early orientation labels.")
}

# -----------------------------
# 5) Checkpoint
# -----------------------------
qs_save(seu, file = QS_CELLCYCLE)
message("Saved: ", QS_CELLCYCLE)
