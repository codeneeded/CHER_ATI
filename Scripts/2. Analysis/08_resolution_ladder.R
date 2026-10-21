############################################################
# CHER ATI — 08: resolution ladder
#
# Single purpose: choose the clustering resolution FROM THE DATA rather than
# inheriting one. Criterion follows miR29a_PDAC_CD8 Scripts/02_analysis/
# 23_resolution_ladder.R: at each level, two clusters are MARKER-REDUNDANT if
# they share >=5 of their top 10 positive markers. A level with redundant pairs
# is over-split for this panel; the chosen level is the last one with zero.
#
# The dot-plot ladder is the other half of the decision — a level where a known
# population is buried inside another is under-split, and that is visible on the
# panel but not in the redundancy count.
#
# Input : QS_INTEGRATED (already carries the clustering grid from script 07)
# Output: Analysis/Integration/resolution_ladder/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
  library(qs2)
})

lad_dir <- file.path(INT_ROOT, "resolution_ladder")
mkdirs(lad_dir)
set.seed(1)

require_file(QS_INTEGRATED, "integrated object (run 06/07)")
seu <- qs_read(QS_INTEGRATED)
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# FastMNN is the embedding taken forward: at res 0.8 all 21 of its clusters
# contained all four donors against 19 of 20 for CCA, and its residual well
# skew is donor identity rather than batch (r = 0.872 with % healthy control).
GRAPH <- "mnn_snn"
RES_LADDER <- c(0.2, 0.4, 0.6, 0.8, 1.0)

panel_present <- intersect(PBMC_PANEL_FLAT, rownames(seu))
message("Panel genes present: ", length(panel_present), " / ", length(PBMC_PANEL_FLAT))
missing <- setdiff(PBMC_PANEL_FLAT, rownames(seu))
if (length(missing)) message("  absent: ", paste(missing, collapse = ", "))

# -----------------------------
# Marker redundancy per level
# -----------------------------
redundant <- function(col) {
  Idents(seu) <- droplevels(factor(seu@meta.data[[col]]))
  mk <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.5,
                       min.pct = 0.3, verbose = FALSE)
  top <- mk %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>%
    slice_max(avg_log2FC, n = 10) %>% summarise(g = list(gene), .groups = "drop")
  n <- nrow(top); pairs <- 0; which_pairs <- character(0)
  if (n > 1) for (i in 1:(n - 1)) for (j in (i + 1):n) {
    shared <- intersect(top$g[[i]], top$g[[j]])
    if (length(shared) >= 5) {
      pairs <- pairs + 1
      which_pairs <- c(which_pairs,
                       sprintf("%s~%s(%d)", top$cluster[i], top$cluster[j],
                               length(shared)))
    }
  }
  list(pairs = pairs, detail = which_pairs, markers = mk)
}

rows <- list()
for (r in RES_LADDER) {
  col <- paste0(GRAPH, "_res.", r)
  stopifnot(col %in% colnames(seu@meta.data))
  cl <- droplevels(factor(seu@meta.data[[col]]))
  res <- redundant(col)

  rows[[as.character(r)]] <- data.frame(
    resolution      = r,
    n_clusters      = nlevels(cl),
    redundant_pairs = res$pairs,
    smallest        = min(table(cl)),
    n_below_min     = sum(table(cl) < MIN_CLUSTER_SIZE),
    detail          = paste(res$detail, collapse = "; "))
  save_csv(res$markers, file.path(lad_dir, paste0("markers_res", r, ".csv")))

  seu$.lad <- cl
  p <- DotPlot(seu, features = panel_present, group.by = ".lad", dot.scale = 5) +
    scale_colour_gradient2(low = "#2166AC", mid = "grey94", high = "#B2182B",
                           midpoint = 0, name = "scaled\nexpression") +
    labs(x = NULL, y = NULL,
         title = sprintf("%d clusters (resolution %s) | redundant pairs: %d",
                         nlevels(cl), r, res$pairs)) +
    theme_cher(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 8, face = "italic"),
          plot.title = element_text(face = "bold", size = 12, hjust = 0))
  save_png(file.path(lad_dir, paste0("DotPlot_res", r, ".png")), p,
           width = 17, height = 1.6 + 0.30 * nlevels(cl))

  message("res ", r, ": ", nlevels(cl), " clusters, ", res$pairs,
          " redundant pairs, smallest ", min(table(cl)))
}

ladder <- do.call(rbind, rows)
save_csv(ladder, file.path(lad_dir, "CHER_resolution_ladder_summary.csv"))
cat("\n=== resolution ladder ===\n")
print(ladder[, c("resolution","n_clusters","redundant_pairs","smallest","n_below_min")],
      row.names = FALSE)

zero <- ladder$resolution[ladder$redundant_pairs == 0]
cat("\nResolutions with ZERO marker-redundant pairs: ",
    if (length(zero)) paste(zero, collapse = ", ") else "(none)", "\n")
if (length(zero)) cat("Last such resolution (the criterion's pick): ", max(zero), "\n")
cat("\nRedundant pair detail per level:\n")
for (i in seq_len(nrow(ladder))) {
  if (nzchar(ladder$detail[i])) cat("  res ", ladder$resolution[i], ": ",
                                    ladder$detail[i], "\n", sep = "")
}
