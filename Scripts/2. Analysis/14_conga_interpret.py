#!/usr/bin/env python
"""CHER ATI - CoNGA stage 4: accounting, confound diagnostics, annotation cross-tabs.

Reads the graph-vs-graph output and maps every hit clone back to the cells it
came from, so the hits can be read against the 14 annotation labels rather than
against CoNGA's own GEX cluster numbers.

The run is pooled. All 19,826 cells from all 8 libraries go into one
GEX graph and one TCR graph. A CoNGA hit is therefore a statement about the
dataset as a whole - a clone whose transcriptional neighbours and whose
receptor neighbours overlap more than chance - and is not a per-timepoint or
per-participant call. The sample composition of the hits is reported below as a
Descriptive diagnostic only: it says which libraries the hit cells happened to
come from, which is needed to spot a hit set driven by one library, but it is
not a test of timepoint and must not be presented as one. Splitting these hits
by timepoint and comparing would be a post-hoc partition of a pooled statistic.

Outputs (Analysis/CoNGA/):
  CHER_conga_hits_annotated.csv            one row per hit clone, with labels
  CHER_conga_hits_by_annotation.csv        hit cells vs TCR+ background
  CHER_conga_hits_by_library.csv           library composition diagnostic
  CHER_conga_hits_by_overlap_type.csv      overlap-type accounting
  CHER_conga_gex_cluster_vs_annotation.csv CoNGA GEX clusters vs our labels
"""
import os
import pandas as pd

STAGE = os.path.dirname(os.path.abspath(__file__))
STAGE = os.path.join(STAGE, "..", "..", "Analysis", "CoNGA")
STAGE = os.path.abspath(STAGE)
PREFIX = "CHER"


def p(name):
    return os.path.join(STAGE, name)


gvg = pd.read_csv(p(f"{PREFIX}_graph_vs_graph.tsv"), sep="\t")
clones = pd.read_csv(p(f"{PREFIX}_clones.tsv"), sep="\t")
bcmap = pd.read_csv(p(f"{PREFIX}_clones.tsv.barcode_mapping.tsv"), sep="\t")
meta = pd.read_csv(p(f"{PREFIX}_conga_cellmeta.csv"))
# The exported metadata carries scRepertoire's clonotype string under
# "clone_id" too. Rename it so it cannot collide with CoNGA's own clone_id,
# which is the key everything below joins on.
meta = meta.rename(columns={"clone_id": "screp_clonotype"})

# clone_index indexes into the clones file as read by conga.
clones = clones.reset_index(drop=True)
clones["clone_index"] = clones.index

# barcode_mapping gives clone_id -> comma-separated barcodes.
bc_col = [c for c in bcmap.columns if "barcode" in c.lower()][-1]
bcmap = bcmap[["clone_id", bc_col]].rename(columns={bc_col: "barcodes"})
long = (bcmap.assign(barcode=bcmap["barcodes"].str.split(","))
             .explode("barcode")
             .drop(columns="barcodes"))
long["barcode"] = long["barcode"].str.strip()
cell = long.merge(meta, on="barcode", how="inner")

n_clone_cells = len(cell)
print(f"clones: {len(clones)}  cells mapped to clones: {n_clone_cells}")

# -----------------------------
# 1) Hit accounting
# -----------------------------
hits = gvg.merge(clones[["clone_index", "clone_id"]], on="clone_index", how="left")
n_hit_rows = len(hits)
n_hit_clones = hits["clone_id"].nunique()
print(f"graph_vs_graph rows: {n_hit_rows}  unique hit clones: {n_hit_clones}")

ot = (hits.groupby(["graph_overlap_type", "nbr_frac"])
          .agg(rows=("clone_id", "size"), clones=("clone_id", "nunique"))
          .reset_index())
ot.to_csv(p(f"{PREFIX}_conga_hits_by_overlap_type.csv"), index=False)
print("\n-- overlap type x nbr_frac --")
print(ot.to_string(index=False))

# -----------------------------
# 2) Hits against the annotation
# -----------------------------
hit_clone_ids = set(hits["clone_id"].dropna())
cell["is_hit"] = cell["clone_id"].isin(hit_clone_ids)

bg = cell.groupby("Annotation").size().rename("clone_cells")
hc = cell[cell.is_hit].groupby("Annotation").size().rename("hit_cells")
by_anno = pd.concat([bg, hc], axis=1).fillna(0).astype(int)
by_anno["pct_of_hit_cells"] = (100 * by_anno.hit_cells / by_anno.hit_cells.sum()).round(2)
by_anno["pct_of_clone_cells"] = (100 * by_anno.clone_cells / by_anno.clone_cells.sum()).round(2)
by_anno["hit_rate_within_label"] = (100 * by_anno.hit_cells / by_anno.clone_cells).round(2)
by_anno = by_anno.sort_values("hit_cells", ascending=False).reset_index()
by_anno.to_csv(p(f"{PREFIX}_conga_hits_by_annotation.csv"), index=False)
print("\n-- hit cells by annotation (background = cells assigned to a clone) --")
print(by_anno.to_string(index=False))

# -----------------------------
# 3) Library composition - Descriptive confound diagnostic, not a Test
# -----------------------------
bgl = cell.groupby("Sample_ID").size().rename("clone_cells")
hcl = cell[cell.is_hit].groupby("Sample_ID").size().rename("hit_cells")
by_lib = pd.concat([bgl, hcl], axis=1).fillna(0).astype(int)
by_lib["pct_of_hit_cells"] = (100 * by_lib.hit_cells / by_lib.hit_cells.sum()).round(2)
by_lib["pct_of_clone_cells"] = (100 * by_lib.clone_cells / by_lib.clone_cells.sum()).round(2)
by_lib["enrichment"] = (by_lib.pct_of_hit_cells / by_lib.pct_of_clone_cells).round(2)
by_lib = by_lib.reset_index()
by_lib.to_csv(p(f"{PREFIX}_conga_hits_by_library.csv"), index=False)
print("\n-- library composition of hits (DESCRIPTIVE; pooled run, not a timepoint test) --")
print(by_lib.to_string(index=False))

# -----------------------------
# 4) CoNGA GEX clusters vs our labels
# -----------------------------
if "gex_cluster" in hits.columns:
    gx = (hits[["clone_id", "gex_cluster"]].drop_duplicates()
          .merge(cell[["clone_id", "Annotation"]].drop_duplicates(), on="clone_id"))
    ct = pd.crosstab(gx.gex_cluster, gx.Annotation)
    ct.to_csv(p(f"{PREFIX}_conga_gex_cluster_vs_annotation.csv"))
    print(f"\n-- CoNGA GEX cluster vs annotation: {ct.shape[0]} clusters x {ct.shape[1]} labels --")

# -----------------------------
# 5) Per-hit table with labels
# -----------------------------
lab = (cell.groupby("clone_id")
           .agg(n_cells=("barcode", "size"),
                annotations=("Annotation", lambda s: ";".join(sorted(set(s)))),
                libraries=("Sample_ID", lambda s: ";".join(sorted(set(s)))),
                participants=("Participant", lambda s: ";".join(sorted(set(s)))),
                timepoints=("Timepoint", lambda s: ";".join(sorted(set(s)))))
           .reset_index())
out = (hits.merge(lab, on="clone_id", how="left")
            .sort_values("conga_score")
            [["clone_id", "conga_score", "graph_overlap_type", "nbr_frac",
              "num_neighbors_gex", "num_neighbors_tcr", "overlap",
              "mait_fraction", "gex_cluster", "tcr_cluster",
              "va", "ja", "cdr3a", "vb", "jb", "cdr3b",
              "n_cells", "annotations", "libraries", "participants", "timepoints"]])
out.to_csv(p(f"{PREFIX}_conga_hits_annotated.csv"), index=False)

print("\n-- strongest 10 hits --")
print(out.head(10)[["conga_score", "graph_overlap_type", "cdr3a", "cdr3b",
                    "n_cells", "annotations", "participants"]].to_string(index=False))

# MAIT check: CoNGA scores MAIT-like receptors readily and they are an
# invariant-TCR artefact rather than an ATI finding, so the fraction is logged.
if "mait_fraction" in out.columns:
    mf = out.drop_duplicates("clone_id")["mait_fraction"]
    print(f"\nmait_fraction over unique hit clones: median {mf.median():.3f}, "
          f"clones with >0.5: {(mf > 0.5).sum()} of {len(mf)}")

print("\n14_conga_interpret OK")
