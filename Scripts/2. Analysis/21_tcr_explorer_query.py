#!/usr/bin/env python
"""CHER ATI - 21: match the repertoire against known TCR records (TCR Explorer).

Runs locally against the harmonized VDJdb + IEDB + McPAS + TCR3d index built by
`tcr-explorer-refresh` (486,542 records). No network call at query time.

The question. The composition and repertoire results say CD8 effector/TEMRA
clones expand during interruption. This asks whether those clones are known
receptors - i.e. whether the expansion looks like bystander memory against
common chronic viruses (CMV, EBV, influenza) rather than something specific to
the interruption.

Two levels of evidence, reported separately:
  CDR3-only match   - same CDR3 beta amino acid sequence. Weaker: convergent
                      CDR3s arise on different V genes with different
                      specificity.
  CDR3 + V match    - same CDR3 beta and same TRBV gene. The defensible call.

What this cannot do.
  - These databases are dominated by common viral and cancer epitopes, so a hit
    is biased toward exactly those specificities. The hit rate is not an
    unbiased estimate of what the repertoire recognises.
  - HLA restriction is not checked against the participants' genotypes (unknown
    here). A CDR3 match whose record carries an HLA the participant may not have
    is weaker still; the MHC field is carried through so this can be revisited.
  - a Non-match is not evidence of HIV specificity. It means "not a known
    receptor in these four databases".

Output: Analysis/VDJ/tcr_explorer/
"""
import os
import sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "Analysis", "VDJ", "tcr_explorer"))
PKG = "/home/ammasakshay/Desktop/Projects/Packages/tcr-explorer-main"
os.environ.setdefault("TCR_EXPLORER_DATA", os.path.join(PKG, "data_local"))

IDX = os.path.join(os.environ["TCR_EXPLORER_DATA"], "records_index.parquet")
if not os.path.isfile(IDX):
    sys.exit(f"records index missing at {IDX}\nRun: TCR_EXPLORER_DATA={os.environ['TCR_EXPLORER_DATA']} "
             f"{PKG}/.venv/bin/tcr-explorer-refresh")

clones_raw = pd.read_csv(os.path.join(OUT, "CHER_clone_beta_table.csv"))
rec = pd.read_parquet(IDX)

# Collapse to the beta receptor before matching.
# The export groups by beta and alpha, so one beta receptor paired with two
# different alpha chains occupies two rows. Matching is done on beta only, so
# those rows are the same query and must be merged first: otherwise a later
# de-duplication on a beta-level key keeps whichever row happens to come first
# and silently discards the other, which can drop an expanded clone and
# understate both the cell count and the expanded count.
clones = (clones_raw
          .groupby(["Sample_ID", "Participant", "Timepoint", "Group",
                    "TRBV", "TRBJ", "CDR3_beta"], as_index=False)
          .agg(n_cells=("n_cells", "sum"),
               clone_size_in_sample=("clone_size_in_sample", "max"),
               expanded=("expanded", "max"),
               n_alpha_pairings=("CDR3_alpha", "nunique"),
               cell_types=("cell_types", lambda s: ";".join(
                   sorted({t for v in s.dropna() for t in str(v).split(";") if t}))),
               clone_ids=("clone_ids", lambda s: ";".join(
                   sorted({t for v in s.dropna() for t in str(v).split(";") if t}))))
          )
clones["expanded"] = clones.expanded.astype(bool)
print(f"clone table: {len(clones_raw)} beta-alpha rows -> {len(clones)} unique beta receptors "
      f"per sample ({int(clones.expanded.sum())} expanded)")
assert clones.duplicated(["Sample_ID", "TRBV", "TRBJ", "CDR3_beta"]).sum() == 0

# Human beta records with a CDR3 and an organism annotation.
beta = rec[(rec.species.str.lower().isin(["human", "homosapiens", "homo sapiens"]))
           & (rec.chain.str.upper().str.contains("B", na=False))].copy()
beta["v_gene_clean"] = beta.v_gene.astype(str).str.replace(r"\*.*$", "", regex=True)
beta = beta[beta.cdr3_aa.notna() & (beta.cdr3_aa != "")]
print(f"index: {len(rec)} records -> {len(beta)} human beta with CDR3")

# -----------------------------
# 1) Exact matching
# -----------------------------
# Collapse the index to one row per (cdr3, v_gene) carrying the distinct
# organisms/epitopes, so a receptor seen in many records does not multiply rows.
agg = {"antigen_organism": lambda s: ";".join(sorted({x for x in s.dropna() if x})),
       "epitope_aa":       lambda s: ";".join(sorted({x for x in s.dropna() if x})[:5]),
       "antigen":          lambda s: ";".join(sorted({x for x in s.dropna() if x})[:5]),
       "mhc_a":            lambda s: ";".join(sorted({x for x in s.dropna() if x})[:5]),
       "source":           lambda s: ";".join(sorted(set(s.dropna()))),
       "source_record_id": "count"}

by_cdr3 = (beta.groupby("cdr3_aa").agg(agg)
           .rename(columns={"source_record_id": "n_records"}).reset_index())
by_cdr3_v = (beta.groupby(["cdr3_aa", "v_gene_clean"]).agg(agg)
             .rename(columns={"source_record_id": "n_records"}).reset_index())

m1 = clones.merge(by_cdr3, left_on="CDR3_beta", right_on="cdr3_aa", how="left")
m1["hit_cdr3"] = m1.cdr3_aa.notna()
m2 = clones.merge(by_cdr3_v, left_on=["CDR3_beta", "TRBV"],
                  right_on=["cdr3_aa", "v_gene_clean"], how="left")
m1["hit_cdr3_and_v"] = m2.cdr3_aa.notna().values
m1["matched_v_gene"] = m2.v_gene_clean.values

hits = m1[m1.hit_cdr3].copy()
hits.to_csv(os.path.join(OUT, "CHER_tcr_explorer_hits.csv"), index=False)

n_c = int(m1.hit_cdr3.sum()); n_cv = int(m1.hit_cdr3_and_v.sum())
print(f"\nCDR3-only matches: {n_c} of {len(m1)} ({100*n_c/len(m1):.2f}%)")
print(f"CDR3 + V matches : {n_cv} of {len(m1)} ({100*n_cv/len(m1):.2f}%)")

# -----------------------------
# 2) Expanded vs non-expanded
# -----------------------------
# The comparison that matters: if the expanding clones were bystander memory,
# they should match known records more often than the unexpanded background.
rows = []
for lvl, col in [("CDR3 only", "hit_cdr3"), ("CDR3 + V", "hit_cdr3_and_v")]:
    for grp, sub in [("expanded", m1[m1.expanded]), ("not expanded", m1[~m1.expanded])]:
        rows.append({"match_level": lvl, "group": grp, "n_receptors": len(sub),
                     "n_hits": int(sub[col].sum()),
                     "pct_hits": round(100 * sub[col].mean(), 2)})
enr = pd.DataFrame(rows)
enr.to_csv(os.path.join(OUT, "CHER_tcr_explorer_enrichment.csv"), index=False)
print("\n-- hit rate, expanded vs background --")
print(enr.to_string(index=False))

# Fisher exact on the defensible level, if scipy is available.
try:
    from scipy.stats import fisher_exact
    a = int(m1[m1.expanded].hit_cdr3_and_v.sum()); b = int((~m1[m1.expanded].hit_cdr3_and_v).sum())
    c = int(m1[~m1.expanded].hit_cdr3_and_v.sum()); d = int((~m1[~m1.expanded].hit_cdr3_and_v).sum())
    orr, pv = fisher_exact([[a, b], [c, d]])
    print(f"\nFisher exact (CDR3+V, expanded vs not): OR={orr:.2f}, p={pv:.3g}  [{a}/{a+b} vs {c}/{c+d}]")
except ImportError:
    print("\n(scipy absent - Fisher test skipped)")

# -----------------------------
# 3) What are the hits against?
# -----------------------------
# The four sources do not agree on organism naming - CMV appears as both "CMV"
# and "Human herpesvirus 5 (Human cytomegalovirus)", SARS-CoV-2 as "SARS-CoV2"
# and "SARS-CoV-2", influenza three different ways. Counting the raw strings
# splits one organism across several rows and understates every one of them.
ORG_MAP = [
    ("cytomegalo|herpesvirus 5|^cmv$",            "CMV"),
    ("epstein|herpesvirus 4|^ebv$",               "EBV"),
    ("influenza",                                 "Influenza A"),
    ("sars-?cov-?2|covid",                        "SARS-CoV-2"),
    ("tuberculosis|mycobacter",                   "M. tuberculosis"),
    ("hepatitis b",                               "Hepatitis B"),
    ("hepatitis c",                               "Hepatitis C"),
    ("immunodeficiency|^hiv",                     "HIV"),
    ("homo ?sapiens|^self$",                      "Human (self / tumour)"),
    ("yellow fever|^yfv$",                        "Yellow fever"),
    ("sars-?cov-?1|coronavirus tor2|^sars$",      "SARS-CoV-1"),
    ("herpes simplex|herpesvirus 1|herpesvirus 2", "HSV"),
]


def norm_org(x):
    s = str(x).strip().lower()
    for pat, lab in ORG_MAP:
        if pd.Series([s]).str.contains(pat, regex=True).iloc[0]:
            return lab
    return str(x).strip()


org = (hits.assign(organism=hits.antigen_organism.str.split(";"))
           .explode("organism"))
org = org[org.organism.notna() & (org.organism != "")]
org["organism"] = org.organism.map(norm_org)
# One receptor can carry several records for the same organism under different
# spellings; collapse so n_receptors counts receptors, not label variants.
org = org.drop_duplicates(subset=["Sample_ID", "TRBV", "TRBJ", "CDR3_beta", "organism"])
org_tab = (org.groupby("organism")
             .agg(n_receptors=("CDR3_beta", "nunique"),
                  n_cells=("n_cells", "sum"),
                  n_expanded=("expanded", "sum"))
             .sort_values("n_receptors", ascending=False).reset_index())
org_tab.to_csv(os.path.join(OUT, "CHER_tcr_explorer_organisms.csv"), index=False)
print("\n-- organisms behind the CDR3-only hits (top 12) --")
print(org_tab.head(12).to_string(index=False))

# -----------------------------
# 4) Attribution: cell type and timepoint
# -----------------------------
attr = (m1.assign(ct=m1.cell_types.str.split(";")).explode("ct")
          .groupby(["ct", "Timepoint"])
          .agg(n_receptors=("CDR3_beta", "nunique"),
               n_hits_cdr3_v=("hit_cdr3_and_v", "sum"))
          .reset_index())
attr["pct_hits"] = (100 * attr.n_hits_cdr3_v / attr.n_receptors).round(2)
attr.to_csv(os.path.join(OUT, "CHER_tcr_explorer_by_celltype.csv"), index=False)
print("\n-- CDR3+V hit rate by cell type and timepoint (receptors >= 20) --")
print(attr[attr.n_receptors >= 20].sort_values("pct_hits", ascending=False)
          .head(14).to_string(index=False))

# -----------------------------
# 5) HIV specifically
# -----------------------------
# The index does contain HIV-annotated receptors, so the absence of HIV matches
# among expanded clones is a measured negative rather than an untested one -
# though still bounded by HLA (unchecked) and by how sparsely HIV is covered.
n_hiv_idx = int(beta.antigen_organism.astype(str)
                .str.contains("immunodeficiency|^HIV", case=False, regex=True).sum())
hiv = org[org.organism == "HIV"].copy()

# How many distinct organisms does each matching receptor hit? a CDR3 matching
# records annotated to CMV, HIV, tuberculosis and SARS-CoV-2 at once is a
# public, promiscuous sequence rather than evidence for any one specificity.
n_org_per_rec = (org.groupby(["Sample_ID", "TRBV", "TRBJ", "CDR3_beta"])
                   .organism.nunique().rename("n_organisms").reset_index())
hiv = hiv.merge(n_org_per_rec, on=["Sample_ID", "TRBV", "TRBJ", "CDR3_beta"],
                how="left")

print(f"\n-- HIV --")
print(f"  HIV-annotated beta records in the index: {n_hiv_idx}")
print(f"  our receptors matching an HIV record (CDR3-only): {hiv.CDR3_beta.nunique()}"
      f" across {int(hiv.n_cells.sum())} cells")
print(f"  of those, from an expanded clone: {int(hiv.expanded.sum())}")
if len(hiv):
    hiv.to_csv(os.path.join(OUT, "CHER_tcr_explorer_hiv_hits.csv"), index=False)
    print(hiv[["Sample_ID", "Timepoint", "TRBV", "CDR3_beta", "n_cells",
               "expanded", "n_organisms", "cell_types"]].to_string(index=False))
    print(f"\n  median distinct organisms matched per HIV-matching receptor: "
          f"{hiv.n_organisms.median():.1f}")
    print(f"  of all matching receptors, median organisms: "
          f"{n_org_per_rec.n_organisms.median():.1f}")

# -----------------------------
# 6) The healthy control as a false-positive control
# -----------------------------
# BC_control is HIV-negative. Any "HIV" match in that library is by definition a
# false positive, which bounds how much weight a CDR3-only match can carry.
ctrl_hiv = hiv[hiv.Group != "ATI"] if len(hiv) else hiv
print("\n-- CDR3-only matching: false-positive bound from the HIV-negative control --")
if len(hiv):
    print(f"  HIV matches in the HIV-NEGATIVE control library: {len(ctrl_hiv)}")
    print(f"  median CDR3 length among HIV matches: {int(hiv.CDR3_beta.str.len().median())} aa"
          f"  (all receptors: {int(clones.CDR3_beta.str.len().median())} aa)")
    print("  Short, germline-proximal CDR3s of this kind are public across donors;"
          " a CDR3-only hit on them is not specificity evidence.")

print("\n21 OK")
