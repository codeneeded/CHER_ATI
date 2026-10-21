############################################################
# CHER ATI — 18: repertoire structure in depth
#
# 13_tcr_repertoire.R established THAT clonal expansion rises B4 -> Interrupt.
# This script asks what the expanded repertoire is made of, following the
# miR29a 26_tcr_repertoire_seven.R pattern: metrics computed from the metadata
# with dplyr rather than from scRepertoire's plotting functions, so every
# number lands in a CSV that can be traced.
#
# Clone size is counted WITHIN SAMPLE throughout (as in miR29a line 197) —
# pooling B4 and Interrupt would let a clone's own expansion inflate its size
# at the earlier timepoint.
#
#   A) Clonal overlap B4 <-> Interrupt within participant  <- the ATI question:
#      are the clones that expand during interruption already there before it?
#   B) Clonal occupancy by cell type and timepoint (size bins)
#   C) Clonal homeostasis — share of the repertoire held by each size bin
#   D) Clonal proportion — share held by the top-N ranked clones
#   E) Clones spanning more than one cell type
#   F) V-gene usage shift
#   G) CDR3 length
#
# Input : QS_WITH_TCR
# Output: Analysis/VDJ/CHER_rep_*.csv + figures
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr); library(ggplot2)
})

VDJ <- file.path(ANALYSIS_ROOT, "VDJ"); mkdirs(VDJ)

require_file(QS_WITH_TCR, "object with TCR")
seu <- qs_read(QS_WITH_TCR)
require_ct_canon(seu)

md <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::select(barcode, Sample_ID, Participant, Timepoint, Group, Annotation,
                clone = dplyr::all_of(CT_CANON), CTaa, CTgene) |>
  dplyr::filter(!is.na(clone), clone != "", clone != "NA")

cat("Cells with a clonotype: ", nrow(md), "\n", sep = "")

# Clone size within sample, and the size bins used throughout.
md <- md |>
  dplyr::group_by(Sample_ID, clone) |>
  dplyr::mutate(clone_n = dplyr::n()) |>
  dplyr::ungroup() |>
  dplyr::mutate(size_bin = cut(clone_n, breaks = c(0, 1, 5, 20, Inf),
                               labels = c("Single (1)", "Small (2-5)",
                                          "Medium (6-20)", "Large (>20)")))

# ─────────────────────────────────────────────────────────────────────────────
# A) Overlap between timepoints within participant
# ─────────────────────────────────────────────────────────────────────────────
# The question the design was built for. If the Interrupt repertoire is
# dominated by clones ALREADY detectable at B4, interruption expanded a
# pre-existing population; if its large clones are mostly absent at B4, it
# recruited new ones. Detection is depth-limited, so "absent at B4" means
# "not sampled at B4" — stated, not glossed.
ov <- list()
for (p in unique(md$Participant[md$Group == "ATI"])) {
  b4  <- md |> dplyr::filter(Participant == p, Timepoint == "B4")
  itr <- md |> dplyr::filter(Participant == p, Timepoint == "Interrupt")
  s_b4 <- unique(b4$clone); s_it <- unique(itr$clone)
  shared <- intersect(s_b4, s_it)

  # Of the Interrupt CELLS, how many sit in a clone that was seen at B4?
  itr_cells_shared <- sum(itr$clone %in% s_b4)
  # Restrict to the clones driving the expansion signal (>=2 cells).
  exp_cl <- unique(itr$clone[itr$clone_n >= 2])
  ov[[p]] <- data.frame(
    Participant = p,
    clones_B4 = length(s_b4), clones_Interrupt = length(s_it),
    clones_shared = length(shared),
    jaccard = round(length(shared) / length(union(s_b4, s_it)), 4),
    pct_Interrupt_clones_seen_at_B4 = round(100 * length(shared) / length(s_it), 2),
    pct_Interrupt_cells_in_shared_clone = round(100 * itr_cells_shared / nrow(itr), 2),
    n_expanded_clones_Interrupt = length(exp_cl),
    pct_expanded_clones_seen_at_B4 =
      round(100 * sum(exp_cl %in% s_b4) / max(length(exp_cl), 1), 2))
}
# Same computation across the control duplicate: two libraries from ONE blood
# draw, so its overlap is the ceiling this metric can reach when the repertoire
# has genuinely not changed. Any B4/Interrupt pair must be read against it.
bc_a <- md |> dplyr::filter(Sample_ID == "BC_control_A")
bc_b <- md |> dplyr::filter(Sample_ID == "BC_control_B")
sa <- unique(bc_a$clone); sb <- unique(bc_b$clone)
exp_b <- unique(bc_b$clone[bc_b$clone_n >= 2])
ov[["BC_control"]] <- data.frame(
  Participant = "BC_control (duplicate ceiling)",
  clones_B4 = length(sa), clones_Interrupt = length(sb),
  clones_shared = length(intersect(sa, sb)),
  jaccard = round(length(intersect(sa, sb)) / length(union(sa, sb)), 4),
  pct_Interrupt_clones_seen_at_B4 = round(100 * length(intersect(sa, sb)) / length(sb), 2),
  pct_Interrupt_cells_in_shared_clone = round(100 * sum(bc_b$clone %in% sa) / nrow(bc_b), 2),
  n_expanded_clones_Interrupt = length(exp_b),
  pct_expanded_clones_seen_at_B4 =
    round(100 * sum(exp_b %in% sa) / max(length(exp_b), 1), 2))

overlap_tbl <- dplyr::bind_rows(ov)
save_csv(overlap_tbl, file.path(VDJ, "CHER_rep_overlap_timepoint.csv"))
cat("\n-- A) overlap B4 <-> Interrupt --\n")
print(as.data.frame(overlap_tbl), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# B) Clonal occupancy: share of CELLS in each size bin
# ─────────────────────────────────────────────────────────────────────────────
occ_ct <- md |>
  dplyr::count(Annotation, Timepoint, size_bin, name = "cells") |>
  dplyr::group_by(Annotation, Timepoint) |>
  dplyr::mutate(pct = round(100 * cells / sum(cells), 2),
                n_cells_total = sum(cells)) |>
  dplyr::ungroup()
save_csv(occ_ct, file.path(VDJ, "CHER_rep_occupancy_by_celltype.csv"))

occ_s <- md |>
  dplyr::count(Sample_ID, Participant, Timepoint, size_bin, name = "cells") |>
  dplyr::group_by(Sample_ID) |>
  dplyr::mutate(pct = round(100 * cells / sum(cells), 2)) |>
  dplyr::ungroup()
save_csv(occ_s, file.path(VDJ, "CHER_rep_occupancy_by_sample.csv"))

cat("\n-- B) occupancy by sample (% of TCR+ cells) --\n")
print(as.data.frame(occ_s |> dplyr::select(Sample_ID, size_bin, pct) |>
        tidyr::pivot_wider(names_from = size_bin, values_from = pct,
                           values_fill = 0)), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# C) Clonal homeostasis: share of the REPERTOIRE (clones, not cells) per bin
# ─────────────────────────────────────────────────────────────────────────────
homeo <- md |>
  dplyr::distinct(Sample_ID, Participant, Timepoint, clone, clone_n, size_bin) |>
  dplyr::count(Sample_ID, Participant, Timepoint, size_bin, name = "clones") |>
  dplyr::group_by(Sample_ID) |>
  dplyr::mutate(pct_of_clones = round(100 * clones / sum(clones), 2)) |>
  dplyr::ungroup()
save_csv(homeo, file.path(VDJ, "CHER_rep_homeostasis.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# D) Clonal proportion: share of cells held by the top-N ranked clones
# ─────────────────────────────────────────────────────────────────────────────
prop_rows <- list()
for (s in unique(md$Sample_ID)) {
  d <- md |> dplyr::filter(Sample_ID == s) |>
    dplyr::distinct(clone, clone_n) |> dplyr::arrange(dplyr::desc(clone_n))
  tot <- sum(d$clone_n)
  for (k in c(1, 5, 10, 20, 50, 100)) {
    prop_rows[[length(prop_rows) + 1]] <- data.frame(
      Sample_ID = s, top_n = k,
      pct_of_cells = round(100 * sum(head(d$clone_n, k)) / tot, 2))
  }
}
prop <- dplyr::bind_rows(prop_rows) |>
  dplyr::left_join(dplyr::distinct(md, Sample_ID, Participant, Timepoint),
                   by = "Sample_ID")
save_csv(prop, file.path(VDJ, "CHER_rep_proportion_topN.csv"))
cat("\n-- D) share of TCR+ cells held by the top N clones --\n")
print(as.data.frame(prop |> dplyr::select(Sample_ID, top_n, pct_of_cells) |>
        tidyr::pivot_wider(names_from = top_n, values_from = pct_of_cells,
                           names_prefix = "top")), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# E) Clones spanning more than one cell type
# ─────────────────────────────────────────────────────────────────────────────
span <- md |>
  dplyr::group_by(Sample_ID, Participant, Timepoint, clone) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   n_types = dplyr::n_distinct(Annotation),
                   types = paste(sort(unique(Annotation)), collapse = ";"),
                   .groups = "drop")
span_summary <- span |>
  dplyr::filter(n_cells >= 2) |>
  dplyr::group_by(Sample_ID, Participant, Timepoint) |>
  dplyr::summarise(expanded_clones = dplyr::n(),
                   pct_spanning_2plus = round(100 * mean(n_types >= 2), 2),
                   .groups = "drop")
save_csv(span, file.path(VDJ, "CHER_rep_clone_celltype_span.csv"))
save_csv(span_summary, file.path(VDJ, "CHER_rep_clone_span_summary.csv"))
cat("\n-- E) expanded clones spanning >1 cell type --\n")
print(as.data.frame(span_summary), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# F) V-gene usage
# ─────────────────────────────────────────────────────────────────────────────
# CTgene is "TRAV.TRAJ.TRAC_TRBV.TRBD.TRBJ.TRBC"; take the first field of each
# chain. Counted once per CLONE, not per cell, so one large clone does not
# dominate the usage profile.
vgene <- md |>
  dplyr::distinct(Sample_ID, Participant, Timepoint, clone, CTgene) |>
  dplyr::mutate(
    chains = strsplit(as.character(CTgene), "_", fixed = TRUE),
    TRAV = vapply(chains, function(z) if (length(z) >= 1) strsplit(z[1], "\\.")[[1]][1] else NA_character_, character(1)),
    TRBV = vapply(chains, function(z) if (length(z) >= 2) strsplit(z[2], "\\.")[[1]][1] else NA_character_, character(1))) |>
  dplyr::select(-chains)

trbv <- vgene |> dplyr::filter(!is.na(TRBV), TRBV != "NA") |>
  dplyr::count(Sample_ID, Participant, Timepoint, TRBV, name = "clones") |>
  dplyr::group_by(Sample_ID) |>
  dplyr::mutate(pct = round(100 * clones / sum(clones), 2)) |>
  dplyr::ungroup()
save_csv(trbv, file.path(VDJ, "CHER_rep_TRBV_usage.csv"))

# Paired shift per gene, against the duplicate ceiling.
trbv_shift <- trbv |>
  dplyr::filter(Timepoint %in% c("B4", "Interrupt")) |>
  dplyr::select(Participant, Timepoint, TRBV, pct) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = pct, values_fill = 0) |>
  dplyr::mutate(delta_pp = Interrupt - B4)
bc_floor_v <- trbv |>
  dplyr::filter(Sample_ID %in% c("BC_control_A", "BC_control_B")) |>
  dplyr::select(Sample_ID, TRBV, pct) |>
  tidyr::pivot_wider(names_from = Sample_ID, values_from = pct, values_fill = 0)
if (all(c("BC_control_A", "BC_control_B") %in% colnames(bc_floor_v))) {
  bc_floor_v$floor_pp <- abs(bc_floor_v$BC_control_A - bc_floor_v$BC_control_B)
  trbv_shift <- trbv_shift |>
    dplyr::left_join(bc_floor_v |> dplyr::select(TRBV, floor_pp), by = "TRBV")
}
save_csv(trbv_shift, file.path(VDJ, "CHER_rep_TRBV_shift.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# G) CDR3 length
# ─────────────────────────────────────────────────────────────────────────────
cdr3 <- md |>
  dplyr::distinct(Sample_ID, Participant, Timepoint, clone, CTaa) |>
  dplyr::mutate(
    parts = strsplit(as.character(CTaa), "_", fixed = TRUE),
    cdr3a = vapply(parts, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1)),
    cdr3b = vapply(parts, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1)),
    len_a = nchar(gsub("NA", "", cdr3a)), len_b = nchar(gsub("NA", "", cdr3b))) |>
  dplyr::select(-parts)
cdr3_summary <- cdr3 |>
  dplyr::group_by(Sample_ID, Participant, Timepoint) |>
  dplyr::summarise(n_clones = dplyr::n(),
                   median_cdr3a = median(len_a[len_a > 0], na.rm = TRUE),
                   median_cdr3b = median(len_b[len_b > 0], na.rm = TRUE),
                   .groups = "drop")
save_csv(cdr3_summary, file.path(VDJ, "CHER_rep_cdr3_length.csv"))
cat("\n-- G) median CDR3 length per sample --\n")
print(as.data.frame(cdr3_summary), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────
occ_plot <- occ_s |>
  dplyr::mutate(Sample_ID = factor(Sample_ID, levels = CHER_SAMPLES$Sample_ID))
p1 <- ggplot(occ_plot, aes(x = Sample_ID, y = pct, fill = size_bin)) +
  geom_col(width = 0.75) +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  labs(x = NULL, y = "% of TCR+ cells", fill = "Clone size",
       title = "Clonal occupancy",
       subtitle = "Clone size counted within sample; BC_control_A/B are one blood draw split in two") +
  theme_cher() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_png(file.path(VDJ, "CHER_rep_occupancy.png"), p1, width = 9.5, height = 5.4)

prop_plot <- prop |>
  dplyr::mutate(Sample_ID = factor(Sample_ID, levels = CHER_SAMPLES$Sample_ID))
p2 <- ggplot(prop_plot, aes(x = top_n, y = pct_of_cells,
                            colour = Participant, group = Sample_ID,
                            linetype = Timepoint)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_x_log10(breaks = c(1, 5, 10, 20, 50, 100)) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Top N clones (log scale)", y = "% of TCR+ cells held",
       title = "Clonal proportion",
       subtitle = "Higher curve = repertoire concentrated in fewer clones") +
  theme_cher()
save_png(file.path(VDJ, "CHER_rep_proportion.png"), p2, width = 8.5, height = 5.4)

cat("\n18 OK\n")
