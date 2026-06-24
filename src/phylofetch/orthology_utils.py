"""
orthology_utils.py — orthology / paralog sanity check (D-036)
-------------------------------------------------------------
After extracting loci from the isolates and importing comparison tips (amplicon references),
align the two together **per locus** and surface divergence outliers — sequences that don't sit
with the rest. A paralog or artifact doesn't announce itself; it just lands on a long branch,
far from everything else. This is deliberately **source-blind**: it indicts an imported GenBank
amplicon reference exactly as readily as one of Exonerate's extractions, because the QC behind a
published barcode is usually opaque (terse methods, no voucher of orthology).

The science: build a pairwise p-distance matrix over the alignment, then flag a sequence when its
**median distance to all the others** is a robust (MAD-based) high-side outlier — this catches
both a lone paralog (no close relative) and a small paralog *cluster* (a couple of sequences that
agree with each other but are far from the majority). A neighbour-joining tree (Biopython, no new
dependency) gives the visual: outliers show as long branches. The automated flag is a guide; the
table + tree are there for the expert eye, which is better than any threshold.

Robustness (D-036 review): pairwise distances are computed only over comparable columns, NaN when
two sequences barely overlap — a partial amplicon tip that covers only a sub-region of a full-gene
extraction is **excluded** from the median (`MIN_OVERLAP`) rather than scored as maximally distant,
which would otherwise fabricate a paralog signal for a perfectly good short reference. The
non-comparable set is alphabet-aware (protein `N` is asparagine, not "any base"). Duplicate record
ids are de-duplicated (with a surfaced warning) so a double-imported accession can't crash the
Biopython matrix — it's exactly the import slip this QC exists to catch.

`analyze_alignment` is pure (operates on already-aligned records) and is the unit-tested core;
`orthology_check` runs MAFFT first and then calls it.
"""

from __future__ import annotations

import io
import math
import os
import re
import statistics
import tempfile

from Bio import Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from .alignment.mafft import run_mafft

# Below this many comparable columns a pairwise p-distance is too noisy to trust (common when an
# amplicon tip overlaps only a sub-region of a full-gene extraction). Such pairs are dropped from
# the median-distance flag metric rather than scored 1.0.
MIN_OVERLAP = 30
_GAPS = set("-.?~")
_ORG_RE = re.compile(r"\[organism=([^\]]+)\]")


def _noncomp_set(protein: bool) -> set:
    """Columns that can't be compared: gaps always; `X`/`x` (unknown residue/base) always; and for
    NUCLEOTIDE `N`/`n` (any base). For PROTEIN, `N` is asparagine — a real residue — so it is NOT
    excluded (D-036 review: excluding it silently dropped real signal on the protein substrate)."""
    s = set(_GAPS) | set("Xx")
    if not protein:
        s |= set("Nn")
    return s


def organism_from_desc(desc: str) -> str:
    """Pull `[organism=…]` from a rich FASTA header; "" when absent."""
    m = _ORG_RE.search(desc or "")
    return m.group(1).strip() if m else ""


def is_isolate(rec_id: str, isolate_ids) -> bool:
    """An extracted-isolate record is `{strain}_{locus}_{product}`; a tip is an accession. Source
    is decided by strain-id prefix membership, so the check stays source-blind but can *label* each
    row (the user must see whether an outlier is theirs or a reference)."""
    return any(rec_id == sid or rec_id.startswith(f"{sid}_") for sid in isolate_ids)


def _pdist_overlap(s1: str, s2: str, noncomp: set) -> tuple[float, int]:
    """(p-distance, comparable-column-count). NaN distance when nothing overlaps."""
    comp = diff = 0
    for a, b in zip(s1, s2):
        if a in noncomp or b in noncomp:
            continue
        comp += 1
        if a.upper() != b.upper():
            diff += 1
    return ((diff / comp) if comp else math.nan), comp


def pairwise_pdistance(s1: str, s2: str, *, protein: bool = False) -> float:
    """Uncorrected p-distance over comparable columns (case-insensitive); NaN when nothing overlaps.
    Pass `protein=True` on the protein substrate so `N` (asparagine) is treated as a residue."""
    return _pdist_overlap(s1, s2, _noncomp_set(protein))[0]


def _pairwise(seqs: list[str], noncomp: set) -> dict:
    pd = {}
    for i in range(len(seqs)):
        for j in range(i):
            pd[(i, j)] = _pdist_overlap(seqs[i], seqs[j], noncomp)
    return pd


def _pair(pd: dict, i: int, j: int) -> tuple[float, int]:
    return pd[(i, j)] if i > j else pd[(j, i)]


def distance_matrix(aln_records, *, protein: bool = False,
                    min_overlap: int = MIN_OVERLAP) -> DistanceMatrix:
    """Lower-triangular p-distance matrix for the NJ tree. Pairs that share fewer than `min_overlap`
    comparable columns (or none) are treated as missing and filled with the max *trustworthy*
    distance so a tree can still be built — the **same** overlap gate the flag metric uses, so the
    tree the UI tells the user to trust isn't distorted by a one-column noise distance (D-036
    review). Display aid only. Requires unique record ids."""
    names = [r.id for r in aln_records]
    seqs = [str(r.seq) for r in aln_records]
    pd = _pairwise(seqs, _noncomp_set(protein))
    trusted = [d for d, c in pd.values() if not math.isnan(d) and c >= min_overlap]
    fill = max(trusted) if trusted else 1.0
    matrix = []
    for i in range(len(seqs)):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                d, c = _pair(pd, i, j)
                row.append(d if (not math.isnan(d) and c >= min_overlap) else fill)
        matrix.append(row)
    return DistanceMatrix(names, matrix)


def robust_high_outliers(values: dict, k: float = 3.5) -> set:
    """High-side robust (median/MAD) outliers, conservatively. Needs ≥4 values to judge; on a tight
    distribution (MAD≈0) only clearly-separated values are flagged, so a homogeneous set flags
    nothing."""
    vals = list(values.values())
    if len(vals) < 4:
        return set()
    med = statistics.median(vals)
    mad = statistics.median([abs(v - med) for v in vals])
    if mad > 0:
        thr = med + k * 1.4826 * mad
    else:
        bigger = [v for v in vals if v > med]
        thr = (med + max(0.05, (max(vals) - med) * 0.5)) if bigger else float("inf")
    return {key for key, v in values.items() if v > thr and v > med * 1.3}


def analyze_alignment(aln_records, meta: dict, k: float = 3.5, *, protein: bool = False,
                      min_overlap: int = MIN_OVERLAP) -> dict:
    """
    Pure analysis of an already-aligned set (unit-tested core). `meta` maps record id →
    {"source": "isolate"|"tip", "organism": str}. Each row is flagged when its **median distance to
    all sufficiently-overlapping others** is a robust outlier; a sequence that overlaps no one by
    ≥`min_overlap` columns is reported as "insufficient overlap" and never flagged (can't assess).
    Returns {rows, ascii_tree, newick, n_seqs, n_flagged, n_assessed, median_nn, n_low_overlap}.
    """
    # Safety net: the Biopython DistanceMatrix rejects duplicate names with a ValueError, so the
    # public core must never be handed dup ids even if a caller skips `load_combined_records` (which
    # dedupes with user-facing warnings). Silent rename here is strictly better than a crash (D-036
    # review); the page path has already deduped, so this is a no-op there.
    aln_records, _ = _dedupe_ids(list(aln_records))
    names = [r.id for r in aln_records]
    seqs = [str(r.seq) for r in aln_records]
    pd = _pairwise(seqs, _noncomp_set(protein))
    n = len(seqs)

    nn, med, low_overlap = {}, {}, set()
    for i, name in enumerate(names):
        good = []
        for j in range(n):
            if j == i:
                continue
            d, c = _pair(pd, i, j)
            if math.isnan(d) or c < min_overlap:
                continue
            good.append((d, names[j]))
        if good:
            d_min, who = min(good, key=lambda x: x[0])
            nn[name] = (who, d_min)
            med[name] = statistics.median([d for d, _ in good])
        else:
            nn[name] = ("", math.nan)
            med[name] = math.nan
            low_overlap.add(name)

    branch, ascii_tree, newick = {}, "", ""
    if n >= 3:
        tree = DistanceTreeConstructor().nj(
            distance_matrix(aln_records, protein=protein, min_overlap=min_overlap))
        # NJ can emit tiny negative branch lengths on non-additive (real) distances. Clamp them on
        # the tree object itself so the persisted Newick / ascii — which the UI hands to a tree
        # viewer and tells the user to trust — never carries a negative branch (D-036 review).
        for clade in tree.find_clades():
            if clade.branch_length and clade.branch_length < 0:
                clade.branch_length = 0.0
        branch = {t.name: float(t.branch_length or 0.0) for t in tree.get_terminals()}
        buf = io.StringIO()
        Phylo.draw_ascii(tree, file=buf)
        ascii_tree = buf.getvalue()
        nbuf = io.StringIO()
        Phylo.write(tree, nbuf, "newick")
        newick = nbuf.getvalue().strip()

    assessable = {nm: v for nm, v in med.items() if not math.isnan(v)}
    flagged_ids = robust_high_outliers(assessable, k=k)
    nn_floor = statistics.median([d for (_w, d) in nn.values() if not math.isnan(d)]) \
        if assessable else 0.0

    rows = []
    for name in names:
        m = meta.get(name, {})
        who, nnd = nn[name]
        if name in low_overlap:
            flagged = False
            reason = f"insufficient overlap (<{min_overlap} shared cols) — not assessed"
        elif name in flagged_ids:
            flagged = True
            if nnd > max(0.05, nn_floor * 2):
                reason = "divergence outlier — no close relative (possible paralog/artifact)"
            else:
                reason = (f"divergence outlier — clusters with {who} "
                          f"(shared divergence; both suspect)")
        else:
            flagged, reason = False, ""
        rows.append({
            "id": name, "source": m.get("source", "?"), "organism": m.get("organism", ""),
            "median_dist": None if math.isnan(med[name]) else round(med[name], 4),
            "nn_id": who, "nn_dist": None if math.isnan(nnd) else round(nnd, 4),
            "branch_len": round(branch.get(name, 0.0), 4) if branch else None,
            "flagged": flagged, "reason": reason,
        })
    rows.sort(key=lambda x: (x["median_dist"] is not None, x["median_dist"] or 0.0), reverse=True)
    return {
        "rows": rows, "ascii_tree": ascii_tree, "newick": newick,
        "n_seqs": n, "n_flagged": len(flagged_ids), "n_assessed": len(assessable),
        "median_nn": round(nn_floor, 4), "n_low_overlap": len(low_overlap),
    }


def _dedupe_ids(records) -> tuple[list, list]:
    """Make record ids unique (Biopython's DistanceMatrix rejects duplicates), surfacing each
    rename as a warning — a duplicate id usually means a double-imported accession, which is itself
    worth flagging to the user."""
    seen, warnings = {}, []
    for r in records:
        if r.id in seen:
            seen[r.id] += 1
            new = f"{r.id}__dup{seen[r.id]}"
            warnings.append(f"duplicate sequence id '{r.id}' in matrix → renamed '{new}' "
                            "(likely a double import — check it)")
            r.id, r.name = new, new
        else:
            seen[r.id] = 1
    return records, warnings


def load_combined_records(matrix_fasta: str, isolate_ids):
    """Read a per-locus combined matrix → (records, meta, warnings) with source + organism, ids
    de-duplicated."""
    records, warnings = _dedupe_ids(list(SeqIO.parse(matrix_fasta, "fasta")))
    meta = {r.id: {"source": "isolate" if is_isolate(r.id, isolate_ids) else "tip",
                   "organism": organism_from_desc(r.description)} for r in records}
    return records, meta, warnings


def orthology_check(matrix_fasta: str, isolate_ids, *, mafft_bin: str = "mafft",
                    mode: str = "auto", threads: int = 4, manager=None, k: float = 3.5,
                    protein: bool = False, min_overlap: int = MIN_OVERLAP) -> dict:
    """
    Align a per-locus combined matrix (isolates + tips) with MAFFT, then analyze it. Returns the
    `analyze_alignment` report (with a `warnings` list), or `{"error": …, "warnings": […]}` if the
    matrix is too small or MAFFT fails.
    """
    records, meta, warnings = load_combined_records(matrix_fasta, isolate_ids)
    if len(records) < 3:
        return {"error": f"need ≥3 sequences to compare ({len(records)} in matrix)",
                "rows": [], "n_seqs": len(records), "n_flagged": 0, "warnings": warnings}
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "in.fasta")
        outp = os.path.join(d, "aln.fasta")
        SeqIO.write(records, inp, "fasta")
        rc, txt = run_mafft(inp, outp, mode=mode, threads=threads, mafft_bin=mafft_bin,
                            run_manager=manager)
        if rc != 0 or not os.path.exists(outp):
            return {"error": f"MAFFT failed: {str(txt)[:200]}", "rows": [],
                    "n_seqs": len(records), "n_flagged": 0, "warnings": warnings}
        aln = list(SeqIO.parse(outp, "fasta"))
    report = analyze_alignment(aln, meta, k=k, protein=protein, min_overlap=min_overlap)
    report["warnings"] = warnings
    return report
