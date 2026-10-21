############################################################
# CHER ATI — 01: build the merged Seurat object
#
# Single purpose: read the 8 per-sample filtered matrices, stamp metadata,
# merge, and compute QC metrics. NO filtering happens here — thresholds are
# applied in 03_filter.R so that the pre-filter distributions can be inspected
# first, and so a threshold change does not require re-reading the H5 files.
#
# Input : <RUN>/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix.h5
# Output: QS_RAW  +  Analysis/QC/CHER_cells_per_sample_raw.csv
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(qs2)
})

mkdirs(QC_ROOT, SAVED_DIR)

# -----------------------------
# 1) Read each H5 -> Seurat
# -----------------------------
objs <- list()

for (i in seq_len(nrow(CHER_SAMPLES))) {
  samp    <- CHER_SAMPLES[i, ]
  h5_path <- sample_h5(samp$Sample_ID)        # stops loudly if absent

  x <- Read10X_h5(h5_path)

  # A multi run with a VDJ library returns a LIST from Read10X_h5. Taking x
  # directly would hand a list to CreateSeuratObject; take Gene Expression.
  if (is.list(x)) {
    nm <- intersect(c("Gene Expression", "RNA"), names(x))
    if (!length(nm)) {
      stop("H5 list has no 'Gene Expression'/'RNA': ", h5_path,
           "\nNames: ", paste(names(x), collapse = ", "), call. = FALSE)
    }
    x_rna <- x[[nm[1]]]
  } else {
    x_rna <- x
  }

  obj <- CreateSeuratObject(counts = x_rna, project = "CHER_ATI", assay = "RNA",
                            min.cells = 3, min.features = 200)

  obj$Sample_ID             <- samp$Sample_ID
  obj$Participant           <- samp$Participant
  obj$Timepoint             <- samp$Timepoint
  obj$Group                 <- samp$Group
  obj$Well                  <- samp$Well
  obj$OCM_Barcode           <- samp$OCM_Barcode
  obj$Tube_Label            <- samp$Tube_Label
  obj$Replicate             <- samp$Replicate
  obj$Participant_Timepoint <- samp$Participant_Timepoint

  # Cell barcodes are only unique WITHIN a GEM well, and the two wells share the
  # same 10x whitelist — without a per-sample prefix, merge() would silently
  # collide barcodes between wells.
  obj <- RenameCells(obj, add.cell.id = samp$Sample_ID)

  objs[[samp$Sample_ID]] <- obj
  message("Loaded: ", samp$Sample_ID,
          " | cells=", ncol(obj), " | features=", nrow(obj))
}

# -----------------------------
# 2) Merge
# -----------------------------
# Merge in ONE call, passing the rest of the list as y. Folding pairwise with
# Reduce(merge) instead produces eight separate counts layers with mangled
# names ("counts.CHER_ATI.SeuratProject.SeuratProject..."), and LayerData() then
# returns ONE of them rather than all cells — so every PercentageFeatureSet
# metric is computed on a partial matrix and silently misaligned.
seu <- merge(x = objs[[1]], y = objs[-1])
DefaultAssay(seu) <- "RNA"

# ncol() on a Seurat object returns a double, so vapply must expect numeric.
stopifnot(ncol(seu) == sum(vapply(objs, ncol, numeric(1))))
stopifnot(!anyDuplicated(colnames(seu)))

# Join the per-sample layers into one counts matrix BEFORE computing any
# metric. Everything downstream (QC fractions, filtering, normalisation) assumes
# a single block whose column order matches the metadata.
seu <- JoinLayers(seu)
stopifnot(identical(Layers(seu[["RNA"]]), "counts"))

# -----------------------------
# 3) QC metrics
# -----------------------------
# log10GenesPerUMI is library complexity: a low value means many reads over few
# genes, the signature of a dying or low-quality cell.
seu$log10GenesPerUMI <- log10(seu$nFeature_RNA) / log10(seu$nCount_RNA)

# Patterns come from config.R so the metric DEFINITION is versioned in one
# place — percent_plat in particular is PF4/PPBP here, not the PECAM1|PF4 used
# by CP018 (PECAM1 is pan-leukocyte; see the note in config.R).
seu <- PercentageFeatureSet(seu, pattern = MITO_PATTERN,     col.name = "percent_mito")
seu <- PercentageFeatureSet(seu, pattern = RIBO_PATTERN,     col.name = "percent_ribo")
seu <- PercentageFeatureSet(seu, pattern = HB_PATTERN,       col.name = "percent_hb")
seu <- PercentageFeatureSet(seu, pattern = PLATELET_PATTERN, col.name = "percent_plat")

# Verify the fractions against an independent recomputation. A multi-layer
# object silently yields wrong values here, so this assertion is the guard that
# the JoinLayers above actually took effect.
local({
  cts <- SeuratObject::LayerData(seu, assay = "RNA", layer = "counts")
  stopifnot(ncol(cts) == ncol(seu))
  nc <- Matrix::colSums(cts)
  stopifnot(max(abs(nc - seu$nCount_RNA)) < 1e-6)
  for (p in list(c("percent_mito", MITO_PATTERN), c("percent_ribo", RIBO_PATTERN),
                 c("percent_hb", HB_PATTERN),     c("percent_plat", PLATELET_PATTERN))) {
    g <- grep(p[2], rownames(cts), value = TRUE)
    man <- if (length(g)) Matrix::colSums(cts[g, , drop = FALSE]) / nc * 100 else rep(0, ncol(cts))
    if (max(abs(man - seu@meta.data[[p[1]]])) > 1e-6) {
      stop("QC metric ", p[1], " disagrees with recomputation — layer misalignment.",
           call. = FALSE)
    }
    message("  verified ", p[1], " (", length(g), " genes)")
  }
})

# -----------------------------
# 4) Raw per-sample counts
# -----------------------------
# Needed by 03_filter.R to report retention, and it is the number every later
# figure's denominator traces back to.
raw_counts <- as.data.frame(table(seu$Sample_ID), responseName = "n_cells_raw")
names(raw_counts)[1] <- "Sample_ID"
raw_counts <- merge(CHER_SAMPLES[, c("Sample_ID","Participant","Timepoint",
                                     "Group","Well","OCM_Barcode")],
                    raw_counts, by = "Sample_ID")
save_csv(raw_counts, file.path(QC_ROOT, "CHER_cells_per_sample_raw.csv"))
print(raw_counts)

cat("\nTotal cells loaded:", ncol(seu), "| genes:", nrow(seu), "\n")

# -----------------------------
# 5) Checkpoint
# -----------------------------
qs_save(seu, file = QS_RAW)
message("Saved: ", QS_RAW)
