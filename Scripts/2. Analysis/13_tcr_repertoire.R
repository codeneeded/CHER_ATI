############################################################
# CHER ATI — 13: repertoire metrics
#
# Single purpose: clonal expansion and diversity, B4 vs interrupt, Within
# participant.
#
# Three conventions that change the answer, stated rather than assumed:
#
#   Clone size is counted within sample. "Expanded" means the clone is expanded
#   in that library. Counting within cell type x sample instead asks a different
#   question (is it expanded *here*) and the two differ by tens of percentage
#   points — a 50-cell clone spread over five cell types reads as unexpanded in
#   each.
#
#   Expansion is a property of the clonotype, not the cell. Clone size is
#   computed over all TCR+ cells of a sample and then joined back, so a lineage
#   is never truncated where its count in one cell type falls below threshold.
#
#   Diversity is rarefied. Raw diversity is depth-dependent and TCR+ cell counts
#   range 593-1,485 across these samples, so unrarefied values would largely
#   measure library size. Each sample is subsampled to the smallest TCR+ count,
#   100 times, and the mean reported.
#
# Input : QS_WITH_TCR
# Output: Analysis/VDJ/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(qs2); library(tidyr)
})

mkdirs(VDJ_ROOT)
set.seed(1)

require_file(QS_WITH_TCR, "object with TCR (run 12_tcr_attach.R)")
seu <- qs_read(QS_WITH_TCR)
require_ct_canon(seu, "TCR object")

md <- seu@meta.data
md$cell <- rownames(md)
tcr <- md[!is.na(md[[CT_CANON]]), ]
cat("TCR+ cells: ", nrow(tcr), " | clonotypes: ",
    dplyr::n_distinct(tcr[[CT_CANON]]), "\n", sep = "")

# -----------------------------
# 1) Clone size, within sample
# -----------------------------
tcr <- tcr %>% group_by(Sample_ID) %>%
  mutate(clone_n = ave(seq_along(.data[[CT_CANON]]),
                       .data[[CT_CANON]], FUN = length)) %>% ungroup()
tcr$expanded_2 <- tcr$clone_n >= 2
tcr$expanded_5 <- tcr$clone_n >= 5

exp_sample <- tcr %>% group_by(Sample_ID) %>%
  summarise(n_TCR = n(), n_clonotypes = n_distinct(.data[[CT_CANON]]),
            pct_in_clone_ge2 = round(100 * mean(expanded_2), 1),
            pct_in_clone_ge5 = round(100 * mean(expanded_5), 1),
            largest_clone = max(clone_n), .groups = "drop") %>% as.data.frame()
exp_sample <- merge(CHER_SAMPLES[, c("Sample_ID","Participant","Timepoint","Group")],
                    exp_sample, by = "Sample_ID")
save_csv(exp_sample, file.path(VDJ_ROOT, "CHER_clonal_expansion_by_sample_source_data.csv"))
cat("\n=== clonal expansion by sample (clone size within sample) ===\n")
print(exp_sample, row.names = FALSE)

# Paired contrast, the comparison the design supports.
paired <- exp_sample %>% filter(Group == "ATI") %>%
  select(Participant, Timepoint, pct_in_clone_ge2, pct_in_clone_ge5, n_clonotypes) %>%
  pivot_wider(names_from = Timepoint,
              values_from = c(pct_in_clone_ge2, pct_in_clone_ge5, n_clonotypes)) %>%
  as.data.frame()
paired$delta_ge2 <- round(paired$pct_in_clone_ge2_Interrupt - paired$pct_in_clone_ge2_B4, 1)
paired$delta_ge5 <- round(paired$pct_in_clone_ge5_Interrupt - paired$pct_in_clone_ge5_B4, 1)
save_csv(paired, file.path(VDJ_ROOT, "CHER_clonal_expansion_paired.csv"))
cat("\n=== paired B4 -> Interrupt, within participant ===\n")
print(paired[, c("Participant","pct_in_clone_ge2_B4","pct_in_clone_ge2_Interrupt",
                 "delta_ge2","delta_ge5")], row.names = FALSE)

bc <- exp_sample[exp_sample$Group != "ATI", ]
cat("\nBC control duplicate (technical floor): ", bc$pct_in_clone_ge2[1], "% vs ",
    bc$pct_in_clone_ge2[2], "% -> gap ",
    abs(diff(bc$pct_in_clone_ge2)), " pp\n", sep = "")

# -----------------------------
# 2) Expansion by cell type
# -----------------------------
exp_ann <- tcr %>% group_by(Annotation, Sample_ID) %>%
  summarise(n_TCR = n(), pct_in_clone_ge2 = round(100 * mean(expanded_2), 1),
            .groups = "drop") %>% as.data.frame()
save_csv(exp_ann, file.path(VDJ_ROOT, "CHER_clonal_expansion_by_celltype_source_data.csv"))

exp_ann_tp <- tcr %>% filter(Group == "ATI") %>%
  group_by(Annotation, Participant, Timepoint) %>%
  summarise(n_TCR = n(), pct_ge2 = round(100 * mean(expanded_2), 1), .groups = "drop") %>%
  as.data.frame()
save_csv(exp_ann_tp, file.path(VDJ_ROOT, "CHER_clonal_expansion_celltype_paired.csv"))

cat("\n=== expansion by cell type, ATI samples pooled by timepoint ===\n")
print(tcr %>% filter(Group == "ATI") %>% group_by(Annotation, Timepoint) %>%
      summarise(n_TCR = n(), pct_ge2 = round(100 * mean(expanded_2), 1), .groups = "drop") %>%
      pivot_wider(names_from = Timepoint, values_from = c(n_TCR, pct_ge2)) %>%
      as.data.frame(), row.names = FALSE)

# -----------------------------
# 3) Rarefied diversity
# -----------------------------
div_metrics <- function(v) {
  tb <- table(v); p <- tb / sum(tb)
  c(richness = length(tb),
    shannon  = -sum(p * log(p)),
    inv_simpson = 1 / sum(p^2),
    evenness = if (length(tb) > 1) -sum(p * log(p)) / log(length(tb)) else NA_real_)
}

n_min <- min(table(tcr$Sample_ID))
cat("\nRarefying every sample to ", n_min, " TCR+ cells, 100 draws\n", sep = "")

rare <- do.call(rbind, lapply(CHER_SAMPLES$Sample_ID, function(s) {
  v <- tcr[[CT_CANON]][tcr$Sample_ID == s]
  m <- replicate(100, div_metrics(sample(v, n_min, replace = FALSE)))
  data.frame(Sample_ID = s, metric = rownames(m),
             mean = round(rowMeans(m), 2),
             sd = round(apply(m, 1, stats::sd), 2), row.names = NULL)
}))
rare_w <- pivot_wider(rare[, c("Sample_ID","metric","mean")],
                      names_from = metric, values_from = mean) %>% as.data.frame()
rare_w <- merge(CHER_SAMPLES[, c("Sample_ID","Participant","Timepoint","Group")],
                rare_w, by = "Sample_ID")
save_csv(rare,   file.path(VDJ_ROOT, "CHER_rarefied_diversity_long.csv"))
save_csv(rare_w, file.path(VDJ_ROOT, "CHER_rarefied_diversity_source_data.csv"))
cat("\n=== rarefied diversity (mean of 100 draws at n=", n_min, ") ===\n", sep = "")
print(rare_w, row.names = FALSE)

# Unrarefied, for the contrast — this is what not to report.
unrare <- do.call(rbind, lapply(CHER_SAMPLES$Sample_ID, function(s) {
  v <- tcr[[CT_CANON]][tcr$Sample_ID == s]
  data.frame(Sample_ID = s, n_TCR = length(v), t(round(div_metrics(v), 2)))
}))
save_csv(unrare, file.path(VDJ_ROOT, "CHER_unrarefied_diversity_for_contrast.csv"))
cat("\n=== unrarefied, for contrast (NOT comparable across samples) ===\n")
print(unrare, row.names = FALSE)

# -----------------------------
# 4) Figures
# -----------------------------
pe <- ggplot(exp_sample, aes(x = Sample_ID, y = pct_in_clone_ge2, fill = Timepoint)) +
  geom_col(colour = "black", linewidth = 0.25) +
  geom_text(aes(label = pct_in_clone_ge2), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = TP_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of TCR+ cells in a clone of >=2",
       title = "Clonal expansion per library (clone size counted within sample)") +
  theme_cher(base_size = 13) + rotate_x(30)
save_png(file.path(VDJ_ROOT, "CHER_clonal_expansion_per_sample.png"), pe,
         width = 9.5, height = 5.5)

pd <- ggplot(rare_w, aes(x = Timepoint, y = inv_simpson, group = Participant,
                         colour = Participant)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3.5) +
  labs(x = NULL, y = "Inverse Simpson (rarefied)",
       title = paste0("Rarefied clonal diversity at n=", n_min, " TCR+ cells")) +
  theme_cher(base_size = 13)
save_png(file.path(VDJ_ROOT, "CHER_rarefied_diversity_paired.png"), pd,
         width = 8, height = 5.5)
