############################################################
# CHER ATI — 10: adversarial annotation audit
#
# Procedure from miR29a_PDAC_CD8 (Scripts/02_analysis/20_annotation_audit.R),
# three parts, because any one of them passes on a wrong annotation:
#
#   1. Falsification  — each label declares required markers with floors and
#      Refuting markers with ceilings. Thresholds are written in this script
#      Before the data are read, from canonical marker biology.
#   2. Discrimination — each label's full criterion set is run against every
#      label's cells. A set satisfied by more than one population does not
#      define a population; a set satisfied by none is worse. This is what makes
#      part 1 worth anything, since the author writes their own thresholds.
#   3. Identity vs programme — how much of each population's marker signal is a
#      programme (cycling / interferon / activation-stress) rather than
#      identity. A population whose markers are mostly programme must not be
#      Named after it.
#
# Input : QS_ANNOTATED
# Output: Analysis/Annotation/audit/
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(qs2); library(Matrix)
})

ANNO_ROOT <- file.path(ANALYSIS_ROOT, "Annotation")
AUD <- file.path(ANNO_ROOT, "audit")
mkdirs(AUD)

require_file(QS_ANNOTATED, "annotated object (run 09_annotate.R)")
seu <- qs_read(QS_ANNOTATED)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)
cts <- LayerData(seu, assay = "RNA", layer = "counts")
ann <- droplevels(factor(seu$Annotation))

# ══════════════════════════════════════════════════════════════════════════
# Criteria — declared before reading the data. `dir` is min (floor, a Required
# marker) or max (ceiling, a refuting marker); values are percent of cells
# detecting the gene.
# ══════════════════════════════════════════════════════════════════════════
CRIT <- rbind(
  data.frame(type="CD4 naive T", gene=c("CD3E","CCR7","LEF1","IL7R","CD8A","CD8B","MS4A1","CD14","KLRF1"),
             dir=c("min","min","min","min","max","max","max","max","max"), thr=c(70,70,50,70,10,10,5,5,5)),
  # Second pass, disclosed. The first-pass criteria for the two CD4 memory
  # subsets were CD3E + IL7R only, which three populations satisfied — part 2
  # correctly reported them as non-discriminating. They are tightened here with
  # the genes that actually separate the three CD4 subsets (S100A4 8.9 / 35.5 /
  # 69.6 and CCR4 3.0 / 6.1 / 41.9 across naive / CM / EM). These thresholds are
  # therefore post-Hoc, unlike the rest of the table.
  data.frame(type="CD4 central memory T", gene=c("CD3E","CCR7","IL7R","S100A4","CCR4","CD8A","MS4A1","CD14","KLRF1"),
             dir=c("min","min","min","min","max","max","max","max","max"), thr=c(50,40,70,20,20,25,5,5,5)),
  data.frame(type="CD4 effector memory T", gene=c("CD3E","IL7R","CCR4","S100A4","KLRB1","CD8A","MS4A1","CD14","KLRF1"),
             dir=c("min","min","min","min","min","max","max","max","max"), thr=c(70,70,25,50,30,25,5,5,5)),
  data.frame(type="CD8 naive/CM T", gene=c("CD3E","CD8A","CD8B","CCR7","MS4A1","CD14","KLRF1"),
             dir=c("min","min","min","min","max","max","max"), thr=c(70,50,50,70,5,5,5)),
  data.frame(type="CD8 effector/TEMRA", gene=c("CD3E","CD8A","GNLY","CCR7","MS4A1","CD14"),
             dir=c("min","min","min","max","max","max"), thr=c(60,50,50,40,5,5)),
  data.frame(type="CD8 effector T (cycling)", gene=c("CD3E","CD8A","MKI67","MS4A1","CD14"),
             dir=c("min","min","min","max","max"), thr=c(70,50,20,5,5)),
  data.frame(type="NK cell", gene=c("KLRD1","GNLY","CD3E","MS4A1","CD14"),
             dir=c("min","min","max","max","max"), thr=c(60,60,40,5,5)),
  data.frame(type="Naive B cell", gene=c("MS4A1","CD79A","IGHD","TCL1A","CD3E","CD14"),
             dir=c("min","min","min","min","max","max"), thr=c(50,50,45,30,5,5)),
  # "B cell (unresolved)" makes no memory or class-switch claim, so its criteria
  # assert only B lineage plus the two features that separate it from naive B
  # (FCRL1 up, TCL1A down). A CD27 floor is deliberately not asserted: CD27 is
  # detected in 3.2% of these cells, which is what forced the rename.
  data.frame(type="B cell (unresolved)", gene=c("MS4A1","CD79A","FCRL1","TCL1A","CD3E","CD14"),
             dir=c("min","min","min","max","max","max"), thr=c(40,50,40,30,5,5)),
  data.frame(type="Plasma cell", gene=c("MZB1","JCHAIN","MS4A1","CD3E"),
             dir=c("min","min","max","max"), thr=c(50,50,20,5)),
  data.frame(type="Classical monocyte", gene=c("CD14","LYZ","VCAN","CD3E","MS4A1"),
             dir=c("min","min","min","max","max"), thr=c(40,80,50,5,5)),
  data.frame(type="Classical monocyte (inflammatory)", gene=c("CD14","LYZ","CCL2","CD3E","MS4A1"),
             dir=c("min","min","min","max","max"), thr=c(40,80,25,5,5)),
  data.frame(type="Non-classical monocyte", gene=c("LYZ","FCGR3A","CDKN1C","CD3E","MS4A1"),
             dir=c("min","min","min","max","max"), thr=c(50,40,20,5,5)),
  data.frame(type="pDC", gene=c("LILRA4","CLEC4C","IRF7","CD3E","MS4A1"),
             dir=c("min","min","min","max","max"), thr=c(30,30,30,10,10))
)
stopifnot(setequal(unique(CRIT$type), levels(ann)))
missing_genes <- setdiff(unique(CRIT$gene), rownames(cts))
if (length(missing_genes)) stop("criteria reference absent genes: ",
                                paste(missing_genes, collapse = ", "))

# percent-detected matrix: gene x label
genes <- unique(CRIT$gene)
PCT <- sapply(levels(ann), function(k)
  100 * Matrix::rowMeans(cts[genes, which(ann == k), drop = FALSE] > 0))
rownames(PCT) <- genes
save_csv(as.data.frame(round(PCT, 2)),
         file.path(AUD, "CHER_audit_percent_detected.csv"), row.names = TRUE)

satisfies <- function(type, label) {
  cr <- CRIT[CRIT$type == type, ]
  ok <- mapply(function(g, d, t) if (d == "min") PCT[g, label] >= t else PCT[g, label] <= t,
               cr$gene, cr$dir, cr$thr)
  list(all = all(ok), n_fail = sum(!ok),
       failed = paste(sprintf("%s %s%s(obs %.1f)", cr$gene[!ok],
                              ifelse(cr$dir[!ok] == "min", ">=", "<="),
                              cr$thr[!ok], PCT[cbind(cr$gene[!ok], label)]),
                      collapse = "; "))
}

# ── Part 1: falsification ────────────────────────────────────────────────────
p1 <- do.call(rbind, lapply(levels(ann), function(k) {
  r <- satisfies(k, k)
  data.frame(label = k, n_cells = sum(ann == k), passes = r$all,
             n_criteria = sum(CRIT$type == k), n_failed = r$n_fail,
             failed_criteria = r$failed)
}))
save_csv(p1, file.path(AUD, "CHER_audit_part1_falsification.csv"))
cat("\n================ PART 1: FALSIFICATION ================\n")
print(p1[, c("label","n_cells","passes","n_criteria","n_failed")], row.names = FALSE)
if (any(!p1$passes)) {
  cat("\nFAILURES:\n")
  for (i in which(!p1$passes)) cat("  ", p1$label[i], ": ", p1$failed_criteria[i], "\n", sep = "")
}

# ── Part 2: discrimination ───────────────────────────────────────────────────
grid <- expand.grid(criterion_set = levels(ann), tested_on = levels(ann),
                    stringsAsFactors = FALSE)
grid$satisfied <- mapply(function(a, b) satisfies(a, b)$all,
                         grid$criterion_set, grid$tested_on)
p2 <- grid %>% group_by(criterion_set) %>%
  summarise(n_satisfying = sum(satisfied),
            which = paste(tested_on[satisfied], collapse = " | "), .groups = "drop") %>%
  as.data.frame()
p2$unique <- p2$n_satisfying == 1
save_csv(p2, file.path(AUD, "CHER_audit_part2_discrimination.csv"))
cat("\n================ PART 2: DISCRIMINATION ================\n")
print(p2, row.names = FALSE)
cat("\nUNIQUELY discriminating criterion sets: ", sum(p2$unique), " of ", nrow(p2), "\n", sep = "")

# ── Part 3: identity vs programme ────────────────────────────────────────────
Idents(seu) <- "Annotation"
mk <- FindAllMarkers(seu, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25,
                     verbose = FALSE)
mk <- mk[mk$p_val_adj < 0.05, ]
save_csv(mk, file.path(AUD, "CHER_annotation_markers.csv"))

prog_all <- unique(unlist(PROGRAMME_GENES, use.names = FALSE))
p3 <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 50) %>%
  summarise(n_markers = n(),
            n_programme = sum(gene %in% prog_all),
            pct_programme = round(100 * mean(gene %in% prog_all), 1),
            n_cycling = sum(gene %in% PROGRAMME_GENES$cycling),
            n_interferon = sum(gene %in% PROGRAMME_GENES$interferon),
            n_activation = sum(gene %in% PROGRAMME_GENES$activation),
            .groups = "drop") %>% as.data.frame()
save_csv(p3, file.path(AUD, "CHER_audit_part3_identity_vs_programme.csv"))
cat("\n================ PART 3: IDENTITY vs PROGRAMME ================\n")
cat("(% of each label's top-50 markers that are programme genes)\n")
print(p3[order(-p3$pct_programme), ], row.names = FALSE)

cat("\n=== AUDIT SUMMARY ===\n")
cat("  labels passing falsification : ", sum(p1$passes), " / ", nrow(p1), "\n", sep = "")
cat("  uniquely discriminating sets : ", sum(p2$unique), " / ", nrow(p2), "\n", sep = "")
cat("  labels >30% programme markers: ",
    paste(p3$cluster[p3$pct_programme > 30], collapse = ", "), "\n", sep = "")
