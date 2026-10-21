#!/usr/bin/env bash
# CHER ATI — CoNGA graph-vs-graph, run Post-Annotation.
#
# Prerequisite: Scripts/2. Analysis/14_conga_export.R has written the gex/
# triplet, the merged contig CSV and the cell metadata into Analysis/CoNGA/.
#
#   stage 1  setup_10x_for_conga.py   clones file + TCRdist kernel PCs
#   stage 2  run_conga.py             the graph-vs-graph test
#
# CoNGA is not in the conda env's site-packages (read-only); it lives under
# Packages/conga_build and is reached via pythonpath. Its tcrdist_cpp binaries
# are prebuilt there.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
STAGE="$REPO/Analysis/CoNGA"
CONGA_BUILD="/home/ammasakshay/Desktop/Projects/Packages/conga_build"
PY="$HOME/.claude-science/conda/envs/conga/bin/python"
export PYTHONPATH="$CONGA_BUILD/site-packages"
PREFIX="CHER"

# The logo figure needs A Converter ON *PATH*, not just AN Importable module.
# conga renders the graph-vs-graph logo panel to SVG and then shells out via
# conga.convert_svg_to_png, which probes for the `convert`, `inkscape`, `rsvg`,
# `cairosvg` and `magick` Commands with shutil.which(). The conga env ships a
# cairosvg console script, but calling $PY by absolute path never puts the env's
# bin/ on PATH, so which() found nothing and the logo step failed while every
# table was written correctly. Prepending the interpreter's own bin/ fixes it.
export PATH="$(dirname "$PY"):$PATH"
export CONGA_PNG_TO_SVG_UTILITY=cairosvg

cd "$STAGE" || exit 1

echo "== stage 1: clones file + TCRdist kernel PCs =="
"$PY" "$CONGA_BUILD/scripts/setup_10x_for_conga.py" \
  --filtered_contig_annotations_csvfile "${PREFIX}_filtered_contig_annotations.csv" \
  --output_clones_file "${PREFIX}_clones.tsv" \
  --organism human \
  > "${PREFIX}_setup_10x_log.txt" 2>&1
rc1=$?
[ -s "${PREFIX}_clones.tsv" ] || { echo "FAILED: no clones file (exit $rc1)";
                                   tail -25 "${PREFIX}_setup_10x_log.txt"; exit 1; }
echo "  clones: $(( $(wc -l < "${PREFIX}_clones.tsv") - 1 ))"

echo "== stage 2: graph vs graph =="
# Cell filters restate this Project'S QC, not CONGA'S Defaults.
# conga re-filters on read at max_genes_per_cell = 3500 and max_percent_mito =
# 0.1. The object handed to it has already been through Scripts/1. QC, whose
# thresholds were justified for this dataset (mito ceiling 15%, no upper gene
# bound — deliberately, because a gene-count ceiling removes the most
# transcriptionally active cells, which here are the expanded CD8 effectors the
# analysis is about). Leaving conga's defaults would impose a second,
# stricter, undocumented QC pass selecting against exactly those cells.
#
# Not fail-fast on exit status, deliberately: after writing the results tables
# run_conga tries to draw the TCR logo figure via convert_svg_to_png, which
# needs ImageMagick or inkscape. Neither is installed, so it exits non-zero
# After the analysis is complete. The file checks below are the real gate.
#
# nbr_fracs are left at the published defaults — thresholds are not tuned here.
"$PY" "$CONGA_BUILD/scripts/run_conga.py" \
  --graph_vs_graph --graph_vs_graph_stats \
  --gex_data gex --gex_data_type 10x_mtx \
  --clones_file "${PREFIX}_clones.tsv" \
  --organism human \
  --clustering_method louvain \
  --max_genes_per_cell 100000 \
  --max_percent_mito 0.15 \
  --checkpoint \
  --outfile_prefix "$PREFIX" \
  > "${PREFIX}_run_conga_log.txt" 2>&1
rc2=$?

if grep -q "convert_svg_to_png failed" "${PREFIX}_run_conga_log.txt" 2>/dev/null; then
  echo "  note: TCR logo PNGs not drawn (no ImageMagick/Inkscape); tables unaffected"
fi
for f in "${PREFIX}_graph_vs_graph.tsv" "${PREFIX}_graph_vs_graph_stats.tsv"; do
  [ -s "$f" ] || { echo "FAILED: run_conga did not write $f (exit $rc2)";
                   tail -30 "${PREFIX}_run_conga_log.txt"; exit 1; }
done
echo "  graph_vs_graph rows: $(( $(wc -l < "${PREFIX}_graph_vs_graph.tsv") - 1 ))"
echo "== done; outputs in $STAGE =="
