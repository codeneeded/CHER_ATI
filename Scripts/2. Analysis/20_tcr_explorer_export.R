############################################################
# CHER ATI — 20: export the beta-chain clone table for TCR Explorer
#
# TCR Explorer (Packages/tcr-explorer-main) runs LOCALLY against a harmonized
# index of VDJdb + IEDB + McPAS + TCR3d (486,542 records after
# `tcr-explorer-refresh`). This script writes the one table that
# 21_tcr_explorer_query.py queries: unique beta receptors with the metadata
# needed to join any hit back onto cell type, participant and timepoint.
#
# Beta chain is the queryable unit for the records index. Alpha is carried
# alongside so paired tcrdist can be attempted where both chains exist.
# Allele suffixes are stripped (TRBV20-1*01 -> TRBV20-1): the index and the
# tcrdist germline table are keyed on IMGT gene-level names.
#
# WHAT A HIT CAN AND CANNOT SAY. A match to a known CMV/EBV/influenza receptor
# argues the clone is bystander memory rather than interruption-driven. A
# NON-match is NOT evidence of HIV specificity — these databases are dominated
# by common viral and cancer epitopes, and this dataset was aligned to a stock
# GRCh38 reference with no HIV genome, so there is no viral readout to
# corroborate against. Read a negative as "not a known receptor", nothing more.
#
# Input : QS_WITH_TCR
# Output: Analysis/VDJ/tcr_explorer/CHER_clone_beta_table.csv
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(dplyr); library(tidyr)
})

TX <- file.path(ANALYSIS_ROOT, "VDJ", "tcr_explorer"); mkdirs(TX)

require_file(QS_WITH_TCR, "object with TCR")
seu <- qs_read(QS_WITH_TCR)
require_ct_canon(seu)

md <- seu@meta.data |>
  tibble::rownames_to_column("barcode") |>
  dplyr::select(barcode, Sample_ID, Participant, Timepoint, Group, Annotation,
                clone = dplyr::all_of(CT_CANON), CTaa, CTgene) |>
  dplyr::filter(!is.na(clone), clone != "", clone != "NA")

# CTgene is "TRAV.TRAJ.TRAC_TRBV.TRBD.TRBJ.TRBC"; CTaa is "CDR3a_CDR3b".
parse_beta <- function(ctgene, ctaa) {
  g <- strsplit(as.character(ctgene), "_", fixed = TRUE)
  a <- strsplit(as.character(ctaa),   "_", fixed = TRUE)
  beta_g <- vapply(g, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  beta_a <- vapply(a, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  bg <- strsplit(beta_g, ".", fixed = TRUE)
  data.frame(
    TRBV = sub("\\*.*$", "", vapply(bg, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))),
    TRBJ = sub("\\*.*$", "", vapply(bg, function(z) if (length(z) >= 3) z[3] else NA_character_, character(1))),
    CDR3_beta = beta_a, stringsAsFactors = FALSE)
}

md <- dplyr::bind_cols(md, parse_beta(md$CTgene, md$CTaa))

# Clone size within sample, as everywhere else in this project.
md <- md |>
  dplyr::group_by(Sample_ID, clone) |>
  dplyr::mutate(clone_n = dplyr::n()) |>
  dplyr::ungroup()

usable <- md |>
  dplyr::filter(!is.na(TRBV), !is.na(TRBJ), !is.na(CDR3_beta),
                TRBV != "NA", TRBJ != "NA", CDR3_beta != "NA", CDR3_beta != "",
                grepl("^TRBV", TRBV), grepl("^TRBJ", TRBJ),
                # TCRex filters non-canonical CDR3s; do it here so the upload
                # is not silently truncated at their end.
                grepl("^C[A-Z]+F$|^C[A-Z]+$", CDR3_beta), nchar(CDR3_beta) >= 8)

cat("Cells with a clonotype:        ", nrow(md), "\n",
    "Cells with a usable beta chain:", nrow(usable), " (",
    round(100 * nrow(usable) / nrow(md), 1), "%)\n", sep = "")

# -----------------------------
# Clone table: one row per (sample, beta receptor)
# -----------------------------
# Alpha chain kept for the paired tcrdist attempt. cell_types/clone_ids let a
# hit be attributed without reloading the Seurat object.
parse_alpha <- function(ctgene, ctaa) {
  g <- strsplit(as.character(ctgene), "_", fixed = TRUE)
  a <- strsplit(as.character(ctaa),   "_", fixed = TRUE)
  ag <- strsplit(vapply(g, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1)),
                 ".", fixed = TRUE)
  data.frame(
    TRAV = sub("\\*.*$", "", vapply(ag, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))),
    CDR3_alpha = vapply(a, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1)),
    stringsAsFactors = FALSE)
}
usable <- dplyr::bind_cols(usable, parse_alpha(usable$CTgene, usable$CTaa))

key <- usable |>
  dplyr::group_by(Sample_ID, Participant, Timepoint, Group,
                  TRBV, TRBJ, CDR3_beta, TRAV, CDR3_alpha) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   clone_size_in_sample = max(clone_n),
                   expanded = max(clone_n) >= 2,
                   cell_types = paste(sort(unique(Annotation)), collapse = ";"),
                   clone_ids = paste(sort(unique(clone)), collapse = ";"),
                   .groups = "drop")
save_csv(key, file.path(TX, "CHER_clone_beta_table.csv"))

cat("\nUnique beta receptors: ", nrow(key), " (",
    sum(key$expanded), " from expanded clones)\n", sep = "")
cat("With a paired alpha:   ",
    sum(!is.na(key$CDR3_alpha) & key$CDR3_alpha != "NA"), "\n", sep = "")

cat("\n20 OK\n")
