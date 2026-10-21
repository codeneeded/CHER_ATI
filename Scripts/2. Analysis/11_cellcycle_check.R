############################################################
# CHER ATI — 11: is the S-phase difference real proliferation?
#
# Single purpose: settle one question. CellCycleScoring put ~40% of ATI cells in
# S phase against ~20% of control cells, which would be a striking result if
# true. But CellCycleScoring assigns phase by the relative score of the S and
# G2M modules within the dataset — every cell gets a phase, including a resting
# one, so in a mostly-quiescent population the S label can be assigned to cells
# that are not cycling at all.
#
# The test: a genuinely proliferating cell transcribes MKI67, TOP2A, pcna and
# the mcm helicases. If S-labelled cells do not carry those above G1-labelled
# cells, the phase call is a scoring artefact and must not be reported as
# proliferation.
#
# Input : QS_ANNOTATED
# Output: Analysis/QC/Cell_Cycle/ verdict tables + figure
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(qs2); library(Matrix)
})

cc_dir <- file.path(QC_ROOT, "Cell_Cycle")
mkdirs(cc_dir)

require_file(QS_ANNOTATED, "annotated object (run 09_annotate.R)")
seu <- qs_read(QS_ANNOTATED)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)
cts <- LayerData(seu, assay = "RNA", layer = "counts")

PROLIF <- intersect(c("MKI67","TOP2A","PCNA","MCM2","MCM6","TYMS","CCNA2","BIRC5","UBE2C"),
                    rownames(cts))
# A cell is proliferating if it detects at least two independent proliferation
# transcripts — one can be ambient, two together is a cycling cell.
n_prolif <- Matrix::colSums(cts[PROLIF, , drop = FALSE] > 0)
seu$prolif_pos <- n_prolif >= 2

# -----------------------------
# 1) Does the phase label track proliferation at all?
# -----------------------------
byphase <- data.frame(
  Phase = levels(seu$Phase),
  n_cells = as.integer(table(seu$Phase)),
  pct_prolif_pos = round(100 * tapply(seu$prolif_pos, seu$Phase, mean), 2),
  median_n_prolif_genes = tapply(n_prolif, seu$Phase, median),
  MKI67_pct = round(100 * tapply(cts["MKI67", ] > 0, seu$Phase, mean), 2),
  TOP2A_pct = round(100 * tapply(cts["TOP2A", ] > 0, seu$Phase, mean), 2))
save_csv(byphase, file.path(cc_dir, "CHER_phase_vs_proliferation_markers.csv"))
cat("\n=== does the phase label track proliferation? ===\n")
print(byphase, row.names = FALSE)

# -----------------------------
# 2) The actual comparison: ATI vs control
# -----------------------------
bysample <- data.frame(
  Sample_ID = CHER_SAMPLES$Sample_ID,
  Group     = CHER_SAMPLES$Group,
  Timepoint = CHER_SAMPLES$Timepoint)
bysample$n_cells     <- as.integer(table(seu$Sample_ID)[bysample$Sample_ID])
bysample$pct_S       <- round(100 * tapply(seu$Phase == "S",   seu$Sample_ID, mean)[bysample$Sample_ID], 2)
bysample$pct_G2M     <- round(100 * tapply(seu$Phase == "G2M", seu$Sample_ID, mean)[bysample$Sample_ID], 2)
bysample$pct_prolif  <- round(100 * tapply(seu$prolif_pos,     seu$Sample_ID, mean)[bysample$Sample_ID], 2)
bysample$MKI67_pct   <- round(100 * tapply(cts["MKI67", ] > 0, seu$Sample_ID, mean)[bysample$Sample_ID], 2)
save_csv(bysample, file.path(cc_dir, "CHER_proliferation_by_sample_source_data.csv"))
cat("\n=== per sample: phase label vs measured proliferation ===\n")
print(bysample, row.names = FALSE)

cat("\nS-phase LABEL      ATI median ", round(median(bysample$pct_S[bysample$Group == "ATI"]), 2),
    " vs control ", round(median(bysample$pct_S[bysample$Group != "ATI"]), 2), "\n", sep = "")
cat("PROLIFERATING cells ATI median ", round(median(bysample$pct_prolif[bysample$Group == "ATI"]), 2),
    " vs control ", round(median(bysample$pct_prolif[bysample$Group != "ATI"]), 2), "\n", sep = "")

# -----------------------------
# 3) Where the cycling cells actually are
# -----------------------------
byann <- as.data.frame(
  seu@meta.data %>% group_by(Annotation) %>%
    summarise(n_cells = n(),
              pct_S = round(100 * mean(Phase == "S"), 2),
              pct_prolif = round(100 * mean(prolif_pos), 2), .groups = "drop"))
save_csv(byann, file.path(cc_dir, "CHER_proliferation_by_annotation.csv"))
cat("\n=== by cell type ===\n"); print(byann, row.names = FALSE)

cat("\nTotal truly proliferating cells: ", sum(seu$prolif_pos), " of ", ncol(seu),
    sprintf(" (%.2f%%)\n", 100 * mean(seu$prolif_pos)), sep = "")
cat("Of those, in 'CD8 effector T (cycling)': ",
    sum(seu$prolif_pos & seu$Annotation == "CD8 effector T (cycling)"),
    sprintf(" (%.1f%%)\n", 100 * mean(seu$Annotation[seu$prolif_pos] ==
                                      "CD8 effector T (cycling)")), sep = "")

p <- ggplot(bysample, aes(x = pct_S, y = pct_prolif, colour = Group)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(aes(label = Sample_ID), size = 3.4, show.legend = FALSE) +
  scale_colour_manual(values = GROUP_COLORS) +
  labs(x = "% of cells LABELLED S phase", y = "% of cells detecting >=2 proliferation genes",
       title = "Phase label against measured proliferation") +
  theme_cher(base_size = 13)
save_png(file.path(cc_dir, "CHER_phase_label_vs_proliferation.png"), p, width = 8.5, height = 6)
