############################################################
# CHER ATI — 16: MAIT check, triggered by the CoNGA result
#
# Single purpose: CoNGA's strongest graph-vs-graph hits are canonical MAIT
# receptors (TRAV1-2 with TRAJ33/12/20, mait_fraction 1.0), and those clones
# carry CD4 naive T / CD4 memory / CD8 naive-CM labels from the res 0.6
# annotation. This script asks whether those cells also carry the MAIT
# TRANSCRIPTIONAL programme — i.e. whether the annotation missed a MAIT
# population, or whether they are TRAV1-2+ cells without the programme.
#
# The receptor and the programme are independent evidence here: the receptor
# comes from the VDJ library, the programme from the GEX library, so agreement
# between them is not circular.
#
# Input : QS_WITH_TCR + Analysis/CoNGA/CHER_clones.tsv (+ barcode mapping)
# Output: Analysis/CoNGA/CHER_mait_check_*.csv
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
})

CONGA_ROOT <- file.path(ANALYSIS_ROOT, "CoNGA")

# MAIT programme. SLC4A10 is the most specific; KLRB1/IL18R1/CCR6/ZBTB16/NCR3
# are the supporting set; DPP4 (CD26) is high on MAIT.
MAIT_GENES <- c("SLC4A10", "KLRB1", "IL18R1", "CCR6", "ZBTB16", "NCR3",
                "DPP4", "RORC", "CXCR6")

require_file(QS_WITH_TCR, "object with TCR")
seu <- qs_read(QS_WITH_TCR)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)

# -----------------------------
# 1) Flag MAIT-configuration receptors from the clones file
# -----------------------------
clones <- read.delim(file.path(CONGA_ROOT, "CHER_clones.tsv"), stringsAsFactors = FALSE)
bcmap  <- read.delim(file.path(CONGA_ROOT, "CHER_clones.tsv.barcode_mapping.tsv"),
                     stringsAsFactors = FALSE)
map <- bcmap |>
  tidyr::separate_rows(barcodes, sep = ",") |>
  dplyr::mutate(barcode = trimws(barcodes)) |>
  dplyr::select(clone_id, barcode) |>
  dplyr::left_join(clones |> dplyr::select(clone_id, va_gene, ja_gene, cdr3a),
                   by = "clone_id")

map$mait_receptor <- grepl("TRAV1-2", map$va_gene) &
                     grepl("TRAJ33|TRAJ12|TRAJ20", map$ja_gene)

seu$mait_receptor <- FALSE
seu$mait_receptor[colnames(seu) %in% map$barcode[map$mait_receptor]] <- TRUE
seu$trav1_2 <- FALSE
seu$trav1_2[colnames(seu) %in% map$barcode[grepl("TRAV1-2", map$va_gene)]] <- TRUE

cat("Cells with MAIT-configuration receptor: ", sum(seu$mait_receptor), "\n",
    "Cells with TRAV1-2 (any J):            ", sum(seu$trav1_2), "\n", sep = "")

# -----------------------------
# 2) Programme detection: MAIT-receptor cells vs their own label
# -----------------------------
genes <- intersect(MAIT_GENES, rownames(seu))
expr  <- LayerData(seu, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
det   <- as.matrix(expr > 0)

md <- data.frame(barcode = colnames(seu),
                 Annotation = as.character(seu$Annotation),
                 mait_receptor = seu$mait_receptor,
                 Sample_ID = seu$Sample_ID,
                 stringsAsFactors = FALSE)

# Compare receptor-positive cells against receptor-negative cells carrying the
# SAME label, so the contrast is not confounded by cell type.
labs <- md |> dplyr::filter(mait_receptor) |> dplyr::count(Annotation) |>
  dplyr::filter(n >= 5) |> dplyr::pull(Annotation)

rows <- list()
for (lb in labs) {
  i_pos <- which(md$Annotation == lb &  md$mait_receptor)
  i_neg <- which(md$Annotation == lb & !md$mait_receptor)
  for (g in genes) {
    rows[[length(rows) + 1]] <- data.frame(
      Annotation = lb, gene = g,
      n_receptor_pos = length(i_pos), n_receptor_neg = length(i_neg),
      pct_pos = round(100 * mean(det[g, i_pos]), 1),
      pct_neg = round(100 * mean(det[g, i_neg]), 1))
  }
}
tab <- dplyr::bind_rows(rows) |> dplyr::mutate(fold = round(pct_pos / pmax(pct_neg, 0.1), 1))
save_csv(tab, file.path(CONGA_ROOT, "CHER_mait_check_programme_by_label.csv"))

cat("\n-- MAIT programme detection, receptor+ vs receptor- within the same label --\n")
for (lb in labs) {
  s <- tab[tab$Annotation == lb, ]
  cat(sprintf("\n%s  (receptor+ n=%d, receptor- n=%d)\n", lb,
              s$n_receptor_pos[1], s$n_receptor_neg[1]))
  for (k in seq_len(nrow(s)))
    cat(sprintf("    %-8s  %5.1f%%  vs %5.1f%%   x%.1f\n",
                s$gene[k], s$pct_pos[k], s$pct_neg[k], s$fold[k]))
}

# -----------------------------
# 3) Where would a MAIT label sit?
# -----------------------------
loc <- md |> dplyr::filter(mait_receptor) |> dplyr::count(Annotation, name = "cells") |>
  dplyr::arrange(dplyr::desc(cells))
save_csv(loc, file.path(CONGA_ROOT, "CHER_mait_check_location.csv"))

# NOTE: compute pct from the ORIGINAL logical column before shadowing its name.
# summarise() evaluates sequentially, so naming the sum `mait_receptor` and
# then calling mean(mait_receptor) averages the already-summed scalar.
bysample <- md |> dplyr::group_by(Sample_ID) |>
  dplyr::summarise(cells = dplyr::n(),
                   n_mait_receptor = sum(mait_receptor),
                   pct_of_cells = round(100 * mean(mait_receptor), 3),
                   .groups = "drop")
save_csv(bysample, file.path(CONGA_ROOT, "CHER_mait_check_by_sample.csv"))
cat("\n-- MAIT-receptor cells per library --\n")
print(as.data.frame(bysample), row.names = FALSE)

cat("\n16 OK\n")
