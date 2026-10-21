############################################################
# CHER ATI — 17: composition shift within participant, B4 vs interrupt
#
# Single purpose: per cell type, the paired change in its share of the library,
# judged against the BC_control_A/B technical floor.
#
# Why the floor matters. BC_control_A and _B are the same blood split across two
# GEM wells, so |a - B| for a cell type is the change that this protocol
# produces when nothing biological has happened. A paired B4->Interrupt change
# smaller than that is not interpretable, regardless of its direction.
#
# No p-values are computed here, deliberately. Three pairs cannot support a composition
# test per cell type, and compositional data are not independent across types
# (shares sum to 1, so one subset expanding forces the others down). What is
# reported is the direction in each participant, the size of the change, and
# whether it clears the technical floor.
#
# Input : Analysis/Annotation/CHER_annotated_cell_metadata.csv
# Output: Analysis/Annotation/CHER_composition_shift_*.csv + figure
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

ANN <- file.path(ANALYSIS_ROOT, "Annotation")
md <- read.csv(file.path(ANN, "CHER_annotated_cell_metadata.csv"), stringsAsFactors = FALSE)
stopifnot(all(c("Sample_ID", "Annotation", "Participant", "Timepoint", "Group") %in% colnames(md)))

# Share of each library held by each cell type.
comp <- md |>
  dplyr::count(Sample_ID, Participant, Timepoint, Group, Annotation, name = "cells") |>
  dplyr::group_by(Sample_ID) |>
  dplyr::mutate(pct = 100 * cells / sum(cells)) |>
  dplyr::ungroup() |>
  tidyr::complete(tidyr::nesting(Sample_ID, Participant, Timepoint, Group),
                  Annotation = ANNO_ORDER,
                  fill = list(cells = 0, pct = 0))
save_csv(comp, file.path(ANN, "CHER_composition_shift_long.csv"))

# Technical floor from the duplicate.
floor_tbl <- comp |>
  dplyr::filter(Group != "ATI") |>
  dplyr::select(Sample_ID, Annotation, pct) |>
  tidyr::pivot_wider(names_from = Sample_ID, values_from = pct) |>
  dplyr::mutate(floor_pp = abs(BC_control_A - BC_control_B))

# Paired change per participant.
paired <- comp |>
  dplyr::filter(Group == "ATI") |>
  dplyr::select(Participant, Timepoint, Annotation, pct, cells) |>
  tidyr::pivot_wider(names_from = Timepoint, values_from = c(pct, cells)) |>
  dplyr::mutate(delta_pp    = pct_Interrupt - pct_B4,
                delta_cells = cells_Interrupt - cells_B4) |>
  dplyr::rename(B4 = pct_B4, Interrupt = pct_Interrupt)

# Smallest absolute Cell-Count change across the three participants. A
# percentage-point change can clear the duplicate floor while resting on one or
# two cells in two of the three participants, which is a single-participant
# event dressed as a consistent one. This column makes that visible.
cells_guard <- paired |>
  dplyr::group_by(Annotation) |>
  dplyr::summarise(min_abs_cells = min(abs(delta_cells)), .groups = "drop")

wide <- paired |>
  dplyr::select(Participant, Annotation, delta_pp) |>
  tidyr::pivot_wider(names_from = Participant, values_from = delta_pp,
                     names_prefix = "d_") |>
  dplyr::left_join(floor_tbl |> dplyr::select(Annotation, floor_pp), by = "Annotation") |>
  dplyr::left_join(cells_guard, by = "Annotation")

dcols <- grep("^d_", colnames(wide), value = TRUE)
wide <- wide |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_up        = sum(dplyr::c_across(dplyr::all_of(dcols)) > 0),
    n_down      = sum(dplyr::c_across(dplyr::all_of(dcols)) < 0),
    consistent  = (n_up == length(dcols)) | (n_down == length(dcols)),
    min_abs_pp  = min(abs(dplyr::c_across(dplyr::all_of(dcols)))),
    median_pp   = median(dplyr::c_across(dplyr::all_of(dcols)))) |>
  dplyr::ungroup() |>
  # "clears_floor" requires that every participant moved the same way and the
  # smallest of the three changes still exceeds what the duplicate produced.
  dplyr::mutate(
    clears_floor = consistent & (min_abs_pp > floor_pp),
    # The reportable call. clears_floor alone can pass on a 1-cell change in two
    # of three participants (Plasma cell does exactly that), so a shift is only
    # called consistent when the weakest participant also moved at least 10
    # cells — below that the "agreement" carries no information.
    consistent_call = clears_floor & (min_abs_cells >= 10)) |>
  dplyr::mutate(Annotation = factor(Annotation, levels = ANNO_ORDER)) |>
  dplyr::arrange(dplyr::desc(consistent_call), dplyr::desc(clears_floor),
                 dplyr::desc(abs(median_pp)))

save_csv(wide, file.path(ANN, "CHER_composition_shift_summary.csv"))

cat("\n-- paired composition change, B4 -> Interrupt (percentage points) --\n")
cat(sprintf("%-34s %8s %8s %8s %8s %8s %7s %s\n", "Cell type", "149496", "150815",
            "154026", "median", "floor", "minCell", "call"))
for (i in seq_len(nrow(wide))) {
  r <- wide[i, ]
  cat(sprintf("%-34s %+8.2f %+8.2f %+8.2f %+8.2f %8.2f %7d %s\n",
              as.character(r$Annotation), r$d_149496, r$d_150815, r$d_154026,
              r$median_pp, r$floor_pp, r$min_abs_cells,
              if (r$consistent_call) "CONSISTENT"
              else if (r$clears_floor) "floor only (too few cells)" else ""))
}

# ── Figure ───────────────────────────────────────────────────────────────────
plot_df <- paired |>
  dplyr::left_join(floor_tbl |> dplyr::select(Annotation, floor_pp), by = "Annotation") |>
  dplyr::mutate(Annotation = factor(Annotation, levels = rev(ANNO_ORDER)))

p <- ggplot(plot_df, aes(x = delta_pp, y = Annotation)) +
  geom_rect(data = dplyr::distinct(plot_df, Annotation, floor_pp),
            aes(xmin = -floor_pp, xmax = floor_pp,
                ymin = as.numeric(Annotation) - 0.45,
                ymax = as.numeric(Annotation) + 0.45),
            inherit.aes = FALSE, fill = "grey85", alpha = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_point(aes(colour = Participant), size = 2.4, alpha = 0.9) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Change in share of library, B4 \u2192 Interrupt (percentage points)",
       y = NULL, colour = "Participant",
       title = "Composition change during treatment interruption",
       subtitle = "Grey band = technical floor from the BC control duplicate") +
  theme_cher()

save_png(file.path(ANN, "CHER_composition_shift.png"), p, width = 10.5, height = 5.6)
cat("\n17 OK\n")
