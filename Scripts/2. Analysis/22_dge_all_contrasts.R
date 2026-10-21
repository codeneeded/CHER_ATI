############################################################
# CHER ATI - 22: pseudobulk differential expression, all three contrasts
#
# Extends 15_dge_pseudobulk.R, which tested B4 against Interrupt only, to the
# two control contrasts, and writes one folder per contrast with one file per
# cell type:
#
#   Analysis/DGE/B4_vs_Interrupt/      paired, participant as blocking factor
#   Analysis/DGE/B4_vs_Control/        unpaired
#   Analysis/DGE/Interrupt_vs_Control/ unpaired
#
# The two control contrasts carry limitations the paired contrast does not, and
# the output columns are built so those limitations travel with the gene lists
# rather than sitting only in this header.
#
# 1. The control arm is one donor. BC_control_A and BC_control_B are two GEM
#    wells loaded from a single blood draw, so they are technical replicates,
#    not biological ones. The contrast is therefore 3 donors against 1, and
#    edgeR estimates dispersion as though there were 2 control replicates, which
#    understates biological variance and makes the p-values anticonservative.
#    They are reported, not withheld, but they are not a 3-against-2 comparison.
#
# 2. Donor, GEM well and dissociation stress are all confounded with group.
#    154026 and both controls sit in well 33445; the other two participants sit
#    in well 34349. The stress and immediate-early programme is roughly five
#    times higher in the controls (1.51% and 1.45% of UMI) than in the ATI
#    samples (0.24-0.41%), and the duplicate agrees with itself, so that
#    difference is a reproducible property of how those libraries were handled.
#
# Two columns address this per gene:
#   tech_floor_logFC  the absolute log fold change between BC_control_A and
#                     BC_control_B for that gene in that cell type. The two
#                     libraries are the same blood, so this is the change the
#                     protocol produces when nothing differs. A control-contrast
#                     logFC below it is not interpretable.
#   stress_programme  TRUE if the gene is in the stress and immediate-early set
#                     in config.R.
#
# Input : QS_WITH_TCR
# Output: Analysis/DGE/<contrast>/CHER_dge_<celltype>.csv + per-contrast summary
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
  library(edgeR)
})

DGE_DIR <- file.path(ANALYSIS_ROOT, "DGE")
# Donor_null is a calibration contrast, not a biological question. 149496 and
# 150815 are both ATI participants in the same GEM well, so a test between them
# differs by DONOR ONLY - same condition, same well, same handling. Each arm
# pools that participant's two timepoints, giving 2 libraries per arm, which
# mirrors the structure of the control contrasts. Its gene count is therefore
# the number of hits this design returns when nothing but donor identity
# differs, and the control contrasts should be read against it. The control
# duplicate cannot serve this purpose: BC_control_A and _B are the same blood,
# so their agreement measures technical reproducibility, not donor variation.
CONTRASTS <- c("B4_vs_Interrupt", "B4_vs_Control", "Interrupt_vs_Control",
               "Donor_null_149496_vs_150815")
for (cc in CONTRASTS) mkdirs(file.path(DGE_DIR, cc))

MIN_CELLS_PER_PSEUDOBULK <- 20
FDR_CUT <- 0.05
LFC_CUT <- 0.5

# Stress and immediate-early genes, for the flag column.
# PROGRAMME_GENES$activation holds the immediate-early and heat-shock genes.
# PER1 and GADD45G are added so the flag matches the set quantified in section 8
# of Analysis/CHER_ATI_FINDINGS.md.
STRESS <- unique(c(PROGRAMME_GENES$activation, "PER1", "GADD45G"))
stopifnot(all(c("FOS", "JUN", "HSPA1A", "DUSP1") %in% STRESS))

require_file(QS_WITH_TCR, "object with TCR (run 12_tcr_attach.R)")
seu <- qs_read(QS_WITH_TCR)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)

# -----------------------------
# 1) Pseudobulk, all 8 libraries
# -----------------------------
seu$pb_group <- paste(seu$Annotation, seu$Sample_ID, sep = "||")
agg <- AggregateExpression(seu, assays = "RNA", group.by = "pb_group",
                           slot = "counts", verbose = FALSE)$RNA
key <- seu@meta.data |>
  dplyr::count(Annotation, Sample_ID, Participant, Timepoint, Group, name = "n_cells") |>
  dplyr::mutate(pb_group = paste(Annotation, Sample_ID, sep = "||"),
                col = gsub("_", "-", pb_group))
stopifnot(all(key$col %in% colnames(agg)))
agg <- agg[, key$col, drop = FALSE]
colnames(agg) <- key$pb_group
cat("pseudobulk libraries: ", ncol(agg), " over ", nrow(agg), " genes\n", sep = "")

# -----------------------------
# 2) Per-gene technical floor from the control duplicate
# -----------------------------
# logCPM difference between the two control libraries of the same cell type.
tech_floor <- function(ct) {
  cols <- key$pb_group[key$Annotation == ct &
                       key$Sample_ID %in% c("BC_control_A", "BC_control_B")]
  if (length(cols) != 2) return(NULL)
  n <- key$n_cells[match(cols, key$pb_group)]
  if (any(n < MIN_CELLS_PER_PSEUDOBULK)) return(NULL)
  m <- agg[, cols, drop = FALSE]
  cpm_m <- edgeR::cpm(m, log = TRUE, prior.count = 1)
  data.frame(gene = rownames(cpm_m),
             tech_floor_logFC = abs(cpm_m[, 1] - cpm_m[, 2]),
             stringsAsFactors = FALSE)
}

# -----------------------------
# 3) Test
# -----------------------------
run_contrast <- function(ct, contrast) {
  k <- key[key$Annotation == ct & key$n_cells >= MIN_CELLS_PER_PSEUDOBULK, ]

  if (contrast == "B4_vs_Interrupt") {
    k <- k[k$Group == "ATI", ]
    complete <- k |> dplyr::count(Participant) |> dplyr::filter(n == 2)
    k <- k[k$Participant %in% complete$Participant, ]
    if (dplyr::n_distinct(k$Participant) < 2) return(NULL)
    lvl <- c("B4", "Interrupt"); fac <- factor(k$Timepoint, levels = lvl)
    design <- model.matrix(~ factor(k$Participant) + fac)
    coef_name <- paste0("fac", lvl[2])
    colnames(design) <- make.names(colnames(design))
    coef_name <- make.names(coef_name)
  } else if (contrast == "Donor_null_149496_vs_150815") {
    k <- k[k$Participant %in% c("149496", "150815"), ]
    if (sum(k$Participant == "149496") < 2 || sum(k$Participant == "150815") < 2)
      return(NULL)
    fac <- factor(k$Participant, levels = c("149496", "150815"))
    design <- model.matrix(~ fac)
    colnames(design) <- make.names(colnames(design))
    coef_name <- make.names("fac150815")
  } else {
    tp <- sub("_vs_Control", "", contrast)
    k <- k[(k$Group == "ATI" & k$Timepoint == tp) | k$Group != "ATI", ]
    # Need both arms present; controls have no biological replication.
    if (sum(k$Group == "ATI") < 2 || sum(k$Group != "ATI") < 2) return(NULL)
    fac <- factor(ifelse(k$Group == "ATI", tp, "Control"),
                  levels = c("Control", tp))
    design <- model.matrix(~ fac)
    colnames(design) <- make.names(colnames(design))
    coef_name <- make.names(paste0("fac", tp))
  }
  if (!coef_name %in% colnames(design)) return(NULL)

  m <- agg[, k$pb_group, drop = FALSE]
  y <- DGEList(counts = m)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = coef_name)

  up_lab <- switch(contrast,
    "B4_vs_Interrupt"             = "up in Interrupt",
    "Donor_null_149496_vs_150815" = "up in 150815",
    paste("up in", sub("_vs_Control", "", contrast)))
  dn_lab <- switch(contrast,
    "B4_vs_Interrupt"             = "up in B4",
    "Donor_null_149496_vs_150815" = "up in 149496",
    "up in Control")

  # Computed before the mutate: inside it, `contrast` would resolve to the
  # column being created rather than to this argument.
  n_ctrl_donors <- if (grepl("Control", contrast)) 1L else 0L

  tt <- topTags(qlf, n = Inf)$table |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(Annotation = ct, contrast = contrast,
                  n_libraries = nrow(k),
                  n_ATI_libraries = sum(k$Group == "ATI"),
                  n_control_donors = n_ctrl_donors,
                  stress_programme = gene %in% STRESS,
                  direction = dplyr::case_when(
                    FDR < FDR_CUT & logFC >=  LFC_CUT ~ up_lab,
                    FDR < FDR_CUT & logFC <= -LFC_CUT ~ dn_lab,
                    TRUE ~ "ns"))

  tf <- tech_floor(ct)
  tt <- if (is.null(tf)) dplyr::mutate(tt, tech_floor_logFC = NA_real_)
        else dplyr::left_join(tt, tf, by = "gene")
  tt |>
    dplyr::mutate(exceeds_tech_floor = !is.na(tech_floor_logFC) &
                                       abs(logFC) > tech_floor_logFC) |>
    dplyr::arrange(FDR, dplyr::desc(abs(logFC)))
}

summ <- list()
for (contrast in CONTRASTS) {
  cat("\n== ", contrast, " ==\n", sep = "")
  for (ct in ANNO_ORDER) {
    tt <- run_contrast(ct, contrast)
    if (is.null(tt)) {
      summ[[paste(contrast, ct)]] <- data.frame(
        contrast = contrast, Annotation = ct, n_genes_tested = NA_integer_,
        n_up = NA_integer_, n_down = NA_integer_, n_sig_over_floor = NA_integer_,
        n_sig_stress = NA_integer_,
        status = sprintf("not tested (<2 libraries per arm at >=%d cells)",
                         MIN_CELLS_PER_PSEUDOBULK))
      next
    }
    save_csv(tt, file.path(DGE_DIR, contrast,
                           sprintf("CHER_dge_%s.csv", fs_safe(ct))))
    sig <- tt[tt$direction != "ns", ]
    summ[[paste(contrast, ct)]] <- data.frame(
      contrast = contrast, Annotation = ct, n_genes_tested = nrow(tt),
      n_up   = sum(grepl("^up in", sig$direction) & sig$logFC > 0),
      n_down = sum(grepl("^up in", sig$direction) & sig$logFC < 0),
      n_sig_over_floor = sum(sig$exceeds_tech_floor, na.rm = TRUE),
      n_sig_stress = sum(sig$stress_programme),
      status = "tested")
    cat(sprintf("  %-34s genes=%5d sig=%4d over_floor=%4d stress=%3d\n",
                ct, nrow(tt), nrow(sig),
                sum(sig$exceeds_tech_floor, na.rm = TRUE),
                sum(sig$stress_programme)))
  }
}

summary_tbl <- dplyr::bind_rows(summ) |>
  dplyr::mutate(Annotation = factor(Annotation, levels = ANNO_ORDER),
                contrast = factor(contrast, levels = CONTRASTS)) |>
  dplyr::arrange(contrast, Annotation)
save_csv(summary_tbl, file.path(DGE_DIR, "CHER_dge_all_contrasts_summary.csv"))

readme <- c(
  "Pseudobulk differential expression, one folder per contrast.",
  "",
  "B4_vs_Interrupt/      paired within participant, participant as a blocking",
  "                      factor. This is the contrast the design supports.",
  "B4_vs_Control/        unpaired.",
  "Interrupt_vs_Control/ unpaired.",
  "",
  "The control arm is a single donor. BC_control_A and BC_control_B are two GEM",
  "wells from one blood draw, so they are technical replicates. Both control",
  "contrasts are 3 donors against 1, with dispersion estimated as if there were",
  "2 control replicates; their p-values are anticonservative.",
  "",
  "Donor, GEM well and dissociation stress are confounded with group. The",
  "stress and immediate-early programme is about five times higher in the",
  "control libraries than in the ATI libraries.",
  "",
  "Columns added for this reason:",
  "  tech_floor_logFC    absolute logFC between BC_control_A and BC_control_B",
  "                      for the same gene and cell type. Same blood, so this is",
  "                      the change the protocol produces with nothing differing.",
  "  exceeds_tech_floor  abs(logFC) greater than tech_floor_logFC.",
  "  stress_programme    gene is in the stress and immediate-early set.",
  "  n_control_donors    1 for the control contrasts, recorded per row.")
writeLines(readme, file.path(DGE_DIR, "README.txt"))

cat("\n22 OK\n")
