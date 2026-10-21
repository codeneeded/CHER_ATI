############################################################
# CHER ATI — 19: clone tracking and clonal structure on the UMAP
#
# The three analyses from CD8_Longitudinal Scripts/2. Analysis/1. V(D)JscRepitoire.R
# that transfer to this design:
#
#   A) clonalCompare alluvial  — CD8 tracked the top clones across
#      CP003_entry -> V12 -> V24. Here the ordered series is B4 -> Interrupt
#      within each participant. This is the picture behind the overlap table in
#      18_tcr_repertoire_depth.R: which ribbons are already present at B4 and
#      which appear only during interruption.
#   B) clonalOverlay            — expanded clones as a density layer on the
#      integrated UMAP, faceted by timepoint. Shows where the expansion sits
#      relative to the annotation, independently of the composition table.
#   C) clonalNetwork            — clonal sharing between cell types.
#
# The BC Duplicate gets its own alluvial panel. BC_control_A -> BC_control_B is
# one blood draw split across two GEM wells, so that panel is what "nothing
# happened" looks like in this assay. Read the participant panels against it.
#
# API Note: the CD8 script was written against scRepertoire 1.x and passes
# `order.by` and `proportion`. Neither exists in 2.2.0 — ordering is taken from
# the order of `samples`, and the plot is always proportional. That default is
# the right one here anyway: the paired libraries differ markedly in recovered
# TCR+ cells (149496: 1,039 at B4 vs 403 at interrupt), so raw counts would
# show a "contraction" that is only library depth.
#
# Input : QS_WITH_TCR + VDJ contigs
# Output: Analysis/VDJ/tracking/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(scRepertoire)
  library(ggplot2); library(ggalluvial)
})

TRK <- file.path(ANALYSIS_ROOT, "VDJ", "tracking"); mkdirs(TRK)

require_file(QS_WITH_TCR, "object with TCR")
seu <- qs_read(QS_WITH_TCR)

# Rebuild the combined list exactly as 12_tcr_attach.R does (scRepertoire 2.2.0
# needs loadContigs first, and passing `samples` without `ID` silently returns
# zero rows, so barcodes are prefixed by hand).
contigs <- loadContigs(lapply(CHER_SAMPLES$Sample_ID, function(s)
                         read.csv(sample_vdj(s), stringsAsFactors = FALSE)),
                       format = "10X")
names(contigs) <- CHER_SAMPLES$Sample_ID
combined <- combineTCR(contigs, removeNA = FALSE, removeMulti = FALSE,
                       filterMulti = TRUE)
names(combined) <- CHER_SAMPLES$Sample_ID
combined <- lapply(names(combined), function(s) {
  d <- combined[[s]]; d$barcode <- paste0(s, "_", sub("^.*_", "", d$barcode)); d
})
names(combined) <- CHER_SAMPLES$Sample_ID
cat("combined rebuilt: ", sum(vapply(combined, nrow, integer(1))), " rows\n", sep = "")

# ─────────────────────────────────────────────────────────────────────────────
# A) Alluvial clone tracking
# ─────────────────────────────────────────────────────────────────────────────
pairs <- list(
  "149496"    = c("149496_B4",    "149496_Interrupt"),
  "150815"    = c("150815_B4",    "150815_Interrupt"),
  "154026"    = c("154026_B4",    "154026_Interrupt"),
  "BC_control" = c("BC_control_A", "BC_control_B"))

for (nm in names(pairs)) {
  sm <- pairs[[nm]]
  ttl <- if (nm == "BC_control")
    "BC control duplicate - one blood draw, two GEM wells (technical floor)"
  else paste0("Participant ", nm, ": B4 \u2192 Interrupt")

  p <- clonalCompare(combined, top.clones = 20, samples = sm,
                     cloneCall = "strict", relabel.clones = TRUE,
                     graph = "alluvial") +
    labs(title = ttl,
         subtitle = "Top 20 clones per sample, as a proportion of that sample's receptors")
  save_png(file.path(TRK, sprintf("CHER_alluvial_%s.png", nm)), p,
           width = 9, height = 7)

  # The numbers behind the ribbons.
  dat <- clonalCompare(combined, top.clones = 20, samples = sm,
                       cloneCall = "strict", relabel.clones = TRUE,
                       exportTable = TRUE)
  save_csv(dat, file.path(TRK, sprintf("CHER_alluvial_%s_source_data.csv", nm)))
}

# ─────────────────────────────────────────────────────────────────────────────
# B) Clonal overlay on the integrated UMAP
# ─────────────────────────────────────────────────────────────────────────────
red <- grep("umap", Reductions(seu), value = TRUE)
red <- if (any(grepl("mnn", red))) grep("mnn", red, value = TRUE)[1] else red[1]
cat("Using reduction: ", red, "\n", sep = "")

ov <- clonalOverlay(seu, reduction = red, freq.cutpoint = 5, bins = 25,
                    facet.by = "Timepoint") +
  guides(colour = "none") +
  labs(title = "Expanded clones on the integrated UMAP",
       subtitle = "Contours = density of cells in clones larger than 5; freq.cutpoint 5 rather than CD8's 10 because these libraries are shallower")
save_png(file.path(TRK, "CHER_clonal_overlay_by_Timepoint.png"), ov,
         width = 15, height = 5.6)

ov2 <- clonalOverlay(seu, reduction = red, freq.cutpoint = 5, bins = 25,
                     facet.by = "Sample_ID") +
  guides(colour = "none") +
  labs(title = "Expanded clones by library",
       subtitle = "BC_control_A vs _B shows the technical reproducibility of this view")
save_png(file.path(TRK, "CHER_clonal_overlay_by_Sample.png"), ov2,
         width = 16, height = 8)

# ─────────────────────────────────────────────────────────────────────────────
# C) Clonal sharing between cell types
# ─────────────────────────────────────────────────────────────────────────────
Idents(seu) <- "Annotation"

# The table is the result; the ggraph rendering is presentation. ggraph 2.2.2
# asks ggplot2 for an `edge_colourbar` guide that ggplot2 4.x no longer
# registers, so the network plot can fail on a version pairing that has nothing
# to do with the data. Export the table first, then attempt the plot.
net_tbl <- clonalNetwork(seu, reduction = red, group.by = "Annotation",
                         cloneCall = "strict", exportTable = TRUE)

# The weight is a ratio with the receiving cluster's receptor count in the
# Denominator, so a cluster with almost no TCR+ cells produces a huge weight off
# a handful of cells — and clonalNetwork draws that as a thick arrow. In this
# dataset NK (13 TCR+ cells of 1,909) and inflammatory classical monocytes (8 of
# 1,294) are exactly that case: their inbound edges are doublet / ambient
# artefacts, not clonal sharing. Attach the counts so the table carries its own
# caveat, and flag edges whose receiving cluster is below a usable floor.
tcr_n <- seu@meta.data |>
  dplyr::mutate(has_tcr = !is.na(.data[[CT_CANON]]) & .data[[CT_CANON]] != "" &
                          .data[[CT_CANON]] != "NA") |>
  dplyr::group_by(Annotation) |>
  dplyr::summarise(n_cells = dplyr::n(), n_TCR = sum(has_tcr), .groups = "drop") |>
  dplyr::mutate(Annotation = as.character(Annotation))

net_out <- as.data.frame(net_tbl) |>
  dplyr::mutate(to = as.character(to), from = as.character(from)) |>
  dplyr::left_join(tcr_n |> dplyr::rename(to = Annotation,
                                          to_n_TCR = n_TCR, to_n_cells = n_cells),
                   by = "to") |>
  dplyr::left_join(tcr_n |> dplyr::rename(from = Annotation,
                                          from_n_TCR = n_TCR, from_n_cells = n_cells),
                   by = "from") |>
  dplyr::mutate(interpretable = to_n_TCR >= 50 & from_n_TCR >= 50) |>
  dplyr::arrange(dplyr::desc(weight))
save_csv(net_out, file.path(TRK, "CHER_clonal_network_source_data.csv"))

cat("\n-- clonal network edges, strongest first --\n")
cat(sprintf("%-34s %-34s %8s %8s %s\n", "to", "from", "weight", "to_TCR", "usable"))
for (i in seq_len(min(nrow(net_out), 8))) {
  r <- net_out[i, ]
  cat(sprintf("%-34s %-34s %8.3f %8d %s\n", r$to, r$from, r$weight, r$to_n_TCR,
              ifelse(r$interpretable, "yes", "NO - low receptor count")))
}

ok <- tryCatch({
  net <- clonalNetwork(seu, reduction = red, group.by = "Annotation",
                       filter.clones = NULL, filter.identity = NULL,
                       cloneCall = "strict") +
    ggraph::scale_edge_colour_continuous(guide = "none")
  save_png(file.path(TRK, "CHER_clonal_network.png"), net, width = 10, height = 8)
  TRUE
}, error = function(e) { message("  clonalNetwork plot failed: ", conditionMessage(e)); FALSE })

if (!ok) {
  # Same information as a heatmap, in plain ggplot2: what share of each cell
  # type's clones is also seen in each other cell type.
  unlink(file.path(TRK, "CHER_clonal_network.png"))
  nt <- as.data.frame(net_tbl)
  stopifnot(all(c("to", "from", "weight") %in% colnames(nt)))
  hm <- ggplot(nt, aes(x = factor(from, levels = ANNO_ORDER),
                       y = factor(to,   levels = rev(ANNO_ORDER)),
                       fill = weight)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "grey95") +
    labs(x = "Clone also found in", y = "Cell type", fill = "Shared\nclone weight",
         title = "Clonal sharing between cell types",
         subtitle = "Rendered as a heatmap: ggraph 2.2.2 and ggplot2 4.x disagree on the edge guide") +
    theme_cher() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  save_png(file.path(TRK, "CHER_clonal_sharing_heatmap.png"), hm,
           width = 9.5, height = 7)
}

cat("\n19 OK\n")
