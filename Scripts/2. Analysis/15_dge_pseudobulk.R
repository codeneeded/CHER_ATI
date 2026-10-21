############################################################
# CHER ATI — 15: pseudobulk differential expression, B4 vs Interrupt
#
# Single purpose: per cell type, sum counts to one library per
# participant x timepoint and test the timepoint effect with participant as a
# blocking factor.
#
# WHY PSEUDOBULK AND NOT PER-CELL. There are three biological replicates, not
# 19,826. A per-cell Wilcoxon treats every cell as independent and returns
# p-values driven by cell count, so it would call hundreds of genes in the
# large clusters and nothing in the small ones regardless of the effect. The
# unit of replication here is the participant.
#
# WHY edgeR WITH A BLOCKING FACTOR. The design is PAIRED — both timepoints come
# from the same participant and, for 149496/150815, the same GEM well. Fitting
# ~ Participant + Timepoint removes the between-participant offset, which is
# much larger than the timepoint effect. The QL F-test is used rather than the
# LRT because with 3 pairs the dispersion is poorly estimated and the QL
# approach is the more conservative of the two.
#
# CONTROLS ARE NOT AN ARM. BC_control_A/B are excluded from testing: they
# differ from the ATI samples by donor, by GEM well and by dissociation stress
# (Analysis/CHER_ATI_FINDINGS.md). They are carried only as a technical floor
# elsewhere.
#
# Input : QS_WITH_TCR
# Output: Analysis/DGE/CHER_dge_<celltype>.csv, summary table + volcano
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(edgeR); library(ggplot2); library(ggrepel)
})

DGE_DIR <- file.path(ANALYSIS_ROOT, "DGE"); mkdirs(DGE_DIR)

# Minimum cells contributing to one pseudobulk library. Below this the summed
# profile is dominated by a handful of cells and the library is dropped rather
# than down-weighted.
MIN_CELLS_PER_PSEUDOBULK <- 20
FDR_CUT <- 0.05
LFC_CUT <- 0.5

require_file(QS_WITH_TCR, "object with TCR (run 12_tcr_attach.R)")
seu <- qs_read(QS_WITH_TCR)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)
require_ct_canon(seu)

# ATI participants only. Subsetting by an explicit cell vector rather than
# subset(subset = Group == "ATI"): Seurat's non-standard evaluation looks the
# column up in the calling frame, which fails inside Rscript.
seu <- subset(seu, cells = colnames(seu)[seu$Group == "ATI"])
cat("ATI cells: ", ncol(seu), " in ", dplyr::n_distinct(seu$Sample_ID),
    " libraries\n", sep = "")

# -----------------------------
# 1) Aggregate to pseudobulk
# -----------------------------
seu$pb_group <- paste(seu$Annotation, seu$Sample_ID, sep = "||")
agg <- AggregateExpression(seu, assays = "RNA", group.by = "pb_group",
                           slot = "counts", verbose = FALSE)$RNA
# AggregateExpression sanitises group names; rebuild the key from the object.
key <- seu@meta.data |>
  dplyr::count(Annotation, Sample_ID, Participant, Timepoint, name = "n_cells") |>
  dplyr::mutate(pb_group = paste(Annotation, Sample_ID, sep = "||"),
                col = gsub("_", "-", pb_group))
stopifnot(all(key$col %in% colnames(agg)))
agg <- agg[, key$col, drop = FALSE]
colnames(agg) <- key$pb_group

save_csv(key |> dplyr::select(Annotation, Sample_ID, Participant, Timepoint, n_cells),
         file.path(DGE_DIR, "CHER_dge_pseudobulk_cell_counts.csv"))

# -----------------------------
# 2) Test each cell type
# -----------------------------
results <- list(); summ <- list()

for (ct in ANNO_ORDER) {
  k <- key[key$Annotation == ct, ]
  k_ok <- k[k$n_cells >= MIN_CELLS_PER_PSEUDOBULK, ]
  # Require a complete pairing: all three participants at both timepoints.
  complete <- k_ok |> dplyr::count(Participant) |> dplyr::filter(n == 2)
  k_ok <- k_ok[k_ok$Participant %in% complete$Participant, ]
  n_pairs <- dplyr::n_distinct(k_ok$Participant)

  if (n_pairs < 2) {
    summ[[ct]] <- data.frame(Annotation = ct, n_pairs = n_pairs,
                             n_libraries = nrow(k_ok), n_genes_tested = NA_integer_,
                             n_up = NA_integer_, n_down = NA_integer_,
                             status = sprintf("not tested (%d complete pairs at >=%d cells)",
                                              n_pairs, MIN_CELLS_PER_PSEUDOBULK))
    next
  }

  m <- agg[, k_ok$pb_group, drop = FALSE]
  grp <- k_ok

  y <- DGEList(counts = m,
               samples = data.frame(Participant = factor(grp$Participant),
                                    Timepoint   = factor(grp$Timepoint,
                                                         levels = c("B4", "Interrupt"))))
  design <- model.matrix(~ Participant + Timepoint, data = y$samples)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "TimepointInterrupt")

  tt <- topTags(qlf, n = Inf)$table |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(Annotation = ct,
                  direction = dplyr::case_when(
                    FDR < FDR_CUT & logFC >=  LFC_CUT ~ "up in Interrupt",
                    FDR < FDR_CUT & logFC <= -LFC_CUT ~ "up in B4",
                    TRUE ~ "ns")) |>
    dplyr::arrange(FDR, dplyr::desc(abs(logFC)))

  results[[ct]] <- tt
  save_csv(tt, file.path(DGE_DIR, sprintf("CHER_dge_%s.csv", fs_safe(ct))))

  summ[[ct]] <- data.frame(
    Annotation = ct, n_pairs = n_pairs, n_libraries = nrow(grp),
    n_genes_tested = nrow(tt),
    n_up   = sum(tt$direction == "up in Interrupt"),
    n_down = sum(tt$direction == "up in B4"),
    status = "tested")
  cat(sprintf("  %-32s pairs=%d libs=%2d genes=%5d  up=%4d down=%4d\n",
              ct, n_pairs, nrow(grp), nrow(tt),
              sum(tt$direction == "up in Interrupt"),
              sum(tt$direction == "up in B4")))
}

summary_tbl <- dplyr::bind_rows(summ) |>
  dplyr::mutate(Annotation = factor(Annotation, levels = ANNO_ORDER)) |>
  dplyr::arrange(Annotation)
save_csv(summary_tbl, file.path(DGE_DIR, "CHER_dge_summary.csv"))

all_res <- dplyr::bind_rows(results)
save_csv(all_res |> dplyr::filter(direction != "ns"),
         file.path(DGE_DIR, "CHER_dge_significant_all_celltypes.csv"))

cat("\nTested ", sum(summary_tbl$status == "tested"), " of ", length(ANNO_ORDER),
    " cell types; ", nrow(dplyr::filter(all_res, direction != "ns")),
    " significant gene x cell-type calls at FDR<", FDR_CUT,
    " and |logFC|>", LFC_CUT, "\n", sep = "")
cat("15 OK\n")
