############################################################
# CHER ATI — 14: export GEX + TCR for CoNGA
#
# Single purpose: write the 10x mtx triplet and the merged contig file that
# CoNGA's setup script expects, with barcodes in the object's naming
# ("<Sample_ID>_<barcode>") so the GEX and TCR sides join.
#
# Follows miR29a_PDAC_CD8 Scripts/02_analysis/26_conga_export.R and the stage
# layout in Analysis_Exploratory/CoNGA.
#
# Input : QS_WITH_TCR
# Output: Analysis/CoNGA/gex/{matrix.mtx,barcodes.tsv,features.tsv}
#         Analysis/CoNGA/CHER_filtered_contig_annotations.csv
#         Analysis/CoNGA/CHER_conga_cellmeta.csv
############################################################

source("~/Desktop/Projects/Repositories/CHER_ATI/Scripts/config.R")

suppressPackageStartupMessages({
  library(Seurat); library(qs2); library(Matrix); library(dplyr)
})

CONGA_ROOT <- file.path(ANALYSIS_ROOT, "CoNGA")
GEX <- file.path(CONGA_ROOT, "gex")
mkdirs(CONGA_ROOT, GEX)

require_file(QS_WITH_TCR, "object with TCR (run 12_tcr_attach.R)")
seu <- qs_read(QS_WITH_TCR)
DefaultAssay(seu) <- "RNA"; seu <- JoinLayers(seu)

cnt <- LayerData(seu, assay = "RNA", layer = "counts")
stopifnot(inherits(cnt, "sparseMatrix"), ncol(cnt) == ncol(seu))

# -----------------------------
# 1) GEX triplet
# -----------------------------
Matrix::writeMM(cnt, file.path(GEX, "matrix.mtx"))
write.table(colnames(cnt), file.path(GEX, "barcodes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# features.tsv needs three columns for 10x_mtx format: id, name, type.
write.table(data.frame(id = rownames(cnt), name = rownames(cnt),
                       type = "Gene Expression"),
            file.path(GEX, "features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# scanpy's read_10x_mtx looks for the gzipped triplet first and only falls back
# to the plain names on an older layout — with uncompressed files present it
# raises FileNotFoundError on matrix.mtx.gz rather than reading what is there.
for (f in c("matrix.mtx", "barcodes.tsv", "features.tsv")) {
  p <- file.path(GEX, f)
  if (file.exists(p)) {
    inc  <- file(p, "rb"); outc <- gzfile(paste0(p, ".gz"), "wb")
    repeat { b <- readBin(inc, "raw", n = 8e6); if (!length(b)) break; writeBin(b, outc) }
    close(inc); close(outc); unlink(p)
  }
}
cat("GEX written: ", nrow(cnt), " genes x ", ncol(cnt), " cells\n", sep = "")

# -----------------------------
# 2) Merged contig file
# -----------------------------
# Barcodes are prefixed to the object's naming, and the lowercase "true" that
# Cell ranger 10.1 writes is normalised to "True" — CoNGA's setup script, like
# scRepertoire, filters on the capitalised forms and silently keeps nothing
# otherwise.
contigs <- do.call(rbind, lapply(CHER_SAMPLES$Sample_ID, function(s) {
  d <- read.csv(sample_vdj(s), stringsAsFactors = FALSE)
  d$barcode <- paste0(s, "_", d$barcode)
  d$sample_id <- s
  d
}))
for (cc in intersect(c("is_cell","high_confidence","full_length","productive"),
                     colnames(contigs))) {
  contigs[[cc]] <- ifelse(tolower(as.character(contigs[[cc]])) == "true", "True", "False")
}
# Keep only contigs whose cell survived QC/annotation, so CoNGA's clone set and
# the GEX matrix describe the same cells.
contigs <- contigs[contigs$barcode %in% colnames(seu), ]
write.csv(contigs, file.path(CONGA_ROOT, "CHER_filtered_contig_annotations.csv"),
          row.names = FALSE, quote = TRUE)
cat("Contigs written: ", nrow(contigs), " rows, ",
    dplyr::n_distinct(contigs$barcode), " cells\n", sep = "")

# -----------------------------
# 3) Cell metadata for interpretation
# -----------------------------
md <- data.frame(
  barcode     = colnames(seu),
  Sample_ID   = seu$Sample_ID,
  Participant = seu$Participant,
  Timepoint   = seu$Timepoint,
  Group       = seu$Group,
  Well        = seu$Well,
  Annotation  = as.character(seu$Annotation),
  Phase       = as.character(seu$Phase),
  clone_id    = seu@meta.data[[CT_CANON]],
  stringsAsFactors = FALSE)
save_csv(md, file.path(CONGA_ROOT, "CHER_conga_cellmeta.csv"))

save_csv(as.data.frame(table(Sample = md$Sample_ID, Annotation = md$Annotation)),
         file.path(CONGA_ROOT, "CHER_conga_export_composition.csv"))
cat("Metadata written for ", nrow(md), " cells\n", sep = "")
