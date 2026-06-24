"""
Tests for the orthology / paralog sanity check engine (D-036). Targets the pure
`analyze_alignment` core (no MAFFT) with hand-built alignments so paralog detection is
deterministic, plus the `orthology_check` MAFFT wrapper via a stub aligner.
"""

import math

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phylofetch import orthology_utils as ou
from phylofetch.orthology_utils import (
    analyze_alignment,
    distance_matrix,
    is_isolate,
    load_combined_records,
    organism_from_desc,
    pairwise_pdistance,
    robust_high_outliers,
)

# 40 bp so the default MIN_OVERLAP (30) is satisfied for full-length records.
BASE = "AAAAACCCCCGGGGGTTTTT" * 2
PARA = "TTTTTGGGGGCCCCCAAAAA" * 2          # every column differs from BASE → p-dist 1.0


def _rec(rid, seq, desc=""):
    return SeqRecord(Seq(seq), id=rid, description=desc)


def _ortho(rid, mut_index):
    s = list(BASE)
    s[mut_index] = "C" if s[mut_index] == "A" else "A"
    return _rec(rid, "".join(s))


def _near(rid, base_seq, n_subs):
    """A sequence close to base_seq: flip the first n_subs positions."""
    s = list(base_seq)
    for i in range(n_subs):
        s[i] = "A" if s[i] != "A" else "C"
    return _rec(rid, "".join(s))


class TestPrimitives:
    def test_pdistance_identical_and_disjoint(self):
        assert pairwise_pdistance("ACGT", "ACGT") == 0.0
        assert pairwise_pdistance("AAAA", "TTTT") == 1.0

    def test_pdistance_zero_overlap_is_nan(self):
        assert math.isnan(pairwise_pdistance("AAAA", "----"))

    def test_pdistance_ignores_gaps_and_is_case_insensitive(self):
        assert pairwise_pdistance("AC-T", "ac-a") == pytest.approx(1 / 3)

    def test_protein_N_is_a_residue_not_a_wildcard(self):
        # protein=True: N (Asn) is compared → N vs K is a real mismatch (1/3 over 3 residues).
        assert pairwise_pdistance("NND", "NKD", protein=True) == pytest.approx(1 / 3)
        # protein=False (nucleotide): N is "any base" → skipped, leaving only the matching D → 0.0.
        assert pairwise_pdistance("NND", "NKD", protein=False) == 0.0

    def test_protein_X_is_excluded(self):
        assert pairwise_pdistance("XD", "KD", protein=True) == 0.0   # X skipped, D matches

    def test_is_isolate_prefix_and_exact(self):
        ids = {"NS26-3-C2", "S9-1B-A2"}
        assert is_isolate("NS26-3-C2_RPB2_genomic", ids)      # prefix
        assert is_isolate("S9-1B-A2", ids)                    # exact
        assert not is_isolate("MH123456.1", ids)

    def test_organism_from_desc(self):
        assert organism_from_desc("x [organism=Alternaria alternata] [framed=no]") == \
            "Alternaria alternata"
        assert organism_from_desc("no brackets") == ""

    def test_robust_outliers_homogeneous_flags_nothing(self):
        assert robust_high_outliers({"a": 0.02, "b": 0.02, "c": 0.03, "d": 0.02}) == set()

    def test_robust_outliers_mad_positive_path(self):
        # Spread orthologs so MAD > 0 → exercises the `med + k·1.4826·mad` branch (NOT the MAD==0
        # fallback). median=0.30, MAD=0.10 (>0), thr≈0.82 → only the 1.5 outlier clears it.
        assert robust_high_outliers(
            {"a": 0.10, "b": 0.20, "c": 0.30, "d": 0.40, "p": 1.5}) == {"p"}

    def test_robust_outliers_mad_zero_path(self):
        # The realistic "N identical orthologs + 1 paralog" shape (MAD==0) — the fallback branch.
        assert robust_high_outliers({"a": 0.0, "b": 0.0, "c": 0.0, "d": 0.0, "p": 0.6}) == {"p"}


class TestAnalyzeAlignment:
    def _meta(self, ids_sources):
        return {i: {"source": s, "organism": ""} for i, s in ids_sources.items()}

    def test_lone_paralog_is_flagged_others_are_not(self):
        recs = [_rec("o1", BASE), _ortho("o2", 0), _ortho("o3", 1),
                _ortho("o4", 2), _rec("p", PARA)]
        meta = self._meta({"o1": "isolate", "o2": "isolate", "o3": "tip",
                           "o4": "tip", "p": "tip"})
        rep = analyze_alignment(recs, meta)
        flagged = {r["id"] for r in rep["rows"] if r["flagged"]}
        assert flagged == {"p"}
        assert rep["rows"][0]["id"] == "p"            # most-divergent first
        assert rep["rows"][0]["source"] == "tip"
        assert "paralog" in rep["rows"][0]["reason"]

    def test_flags_a_tip_paralog_not_just_isolates(self):
        recs = [_rec("iso1", BASE), _ortho("iso2", 0), _ortho("iso3", 1),
                _rec("REF_BAD", PARA)]
        meta = self._meta({"iso1": "isolate", "iso2": "isolate", "iso3": "isolate",
                           "REF_BAD": "tip"})
        rep = analyze_alignment(recs, meta)
        bad = next(r for r in rep["rows"] if r["id"] == "REF_BAD")
        assert bad["flagged"] and bad["source"] == "tip"

    def test_flag_follows_median_not_nearest_neighbour(self):
        # P and Q form a divergent PAIR: each other's nearest neighbour (SMALL nn_dist), but far
        # from the 3-ortholog majority (LARGE median). A nearest-neighbour rule would MISS them;
        # the median-distance metric must catch them. This pins the central scientific claim.
        Q = _near("Q", PARA, 2)                       # ~0.05 from P(PARA), ~1.0 from orthologs
        recs = [_rec("o1", BASE), _ortho("o2", 0), _ortho("o3", 1), _rec("P", PARA), Q]
        meta = self._meta({k: "tip" for k in ("o1", "o2", "o3", "P", "Q")})
        rep = analyze_alignment(recs, meta)
        flagged = {r["id"] for r in rep["rows"] if r["flagged"]}
        # Exact match: the divergent pair is flagged AND the 3-ortholog majority is NOT — pins
        # against the documented majority-inversion failure mode (over-flagging the true orthologs).
        assert flagged == {"P", "Q"}
        prow = next(r for r in rep["rows"] if r["id"] == "P")
        assert prow["nn_dist"] < 0.2                  # has a CLOSE neighbour (Q) …
        assert prow["median_dist"] > 0.4              # … yet flagged on high median
        assert "clusters with" in prow["reason"]

    def test_low_overlap_sequence_is_not_assessed_not_flagged(self):
        # A short tip overlapping the others by <MIN_OVERLAP columns must NOT be fabricated into a
        # paralog (D-036 review H2) — it's reported "insufficient overlap" and left unflagged.
        short = _rec("shorttip", "ACGTACGTAC" + "-" * 30)   # only 10 non-gap columns
        recs = [_rec("o1", BASE), _ortho("o2", 0), _ortho("o3", 1), short]
        meta = self._meta({k: "isolate" for k in ("o1", "o2", "o3")} | {"shorttip": "tip"})
        rep = analyze_alignment(recs, meta)
        srow = next(r for r in rep["rows"] if r["id"] == "shorttip")
        assert not srow["flagged"]
        assert srow["median_dist"] is None
        assert "insufficient overlap" in srow["reason"]
        assert rep["n_low_overlap"] == 1

    def test_homogeneous_set_flags_nothing_and_builds_tree(self):
        recs = [_rec("o1", BASE), _ortho("o2", 0), _ortho("o3", 1), _ortho("o4", 2)]
        meta = self._meta({k: "isolate" for k in ("o1", "o2", "o3", "o4")})
        rep = analyze_alignment(recs, meta)
        assert rep["n_flagged"] == 0
        assert rep["newick"]

    def test_tiny_set_does_not_crash(self):
        recs = [_rec("a", BASE), _rec("b", PARA)]
        rep = analyze_alignment(recs, self._meta({"a": "isolate", "b": "tip"}))
        assert rep["n_seqs"] == 2 and rep["n_flagged"] == 0

    def test_duplicate_ids_do_not_crash_the_core(self):
        # The Biopython DistanceMatrix raises ValueError on duplicate names; the public core must
        # defensively dedupe even when a caller skips load_combined_records (D-036 review). A double-
        # imported accession is exactly the slip this QC exists to catch — it must not be a crash.
        recs = [_rec("x", BASE), _rec("x", PARA), _ortho("y", 0), _ortho("z", 1)]
        rep = analyze_alignment(recs, self._meta({"x": "tip", "y": "isolate", "z": "isolate"}))
        ids = [r["id"] for r in rep["rows"]]
        assert len(ids) == len(set(ids)) == 4          # renamed, not crashed, none lost

    def test_fewer_than_four_assessable_never_flags(self):
        # 3 sequences (one wildly divergent) — the robust test needs ≥4 assessable values, so it
        # reports n_assessed and flags nothing rather than guessing from too little data (D-036
        # review: the page uses n_assessed to avoid a misleading "no outliers").
        recs = [_rec("o1", BASE), _ortho("o2", 0), _rec("p", PARA)]
        rep = analyze_alignment(recs, self._meta({"o1": "tip", "o2": "tip", "p": "tip"}))
        assert rep["n_flagged"] == 0
        assert rep["n_assessed"] == 3                  # all overlap, but < 4 → can't judge

    def test_persisted_newick_has_no_negative_branch_lengths(self):
        # NJ can emit tiny negative branches on real (non-additive) distances; the tree object is
        # clamped so the persisted Newick/ascii the UI hands a viewer never carries one (D-036).
        # This exact 6-record set provably yields a negative branch WITHOUT the clamp (verified),
        # so the assertion genuinely guards the clamp rather than passing vacuously.
        recs = [_rec("o1", BASE), _near("n0", BASE, 0), _near("n1", BASE, 4),
                _near("n2", BASE, 0), _near("n3", BASE, 8), _rec("p", PARA)]
        rep = analyze_alignment(recs, self._meta({k: "tip" for k in
                                                  ("o1", "n0", "n1", "n2", "n3", "p")}))
        assert ":-" not in rep["newick"]
        assert all(r["branch_len"] is None or r["branch_len"] >= 0 for r in rep["rows"])

    def test_distance_matrix_gates_low_overlap_pairs_like_the_flag(self):
        # The NJ tree (which the UI tells the user to trust) must use the SAME overlap gate as the
        # flag, so a one-column-noise distance can't pull a partial sequence misleadingly close
        # (D-036 review). A=BASE; B differs by 0.5 over the full 40 cols (trusted); C overlaps A in
        # only 5 identical cols (raw d=0 — pure noise). With the gate, C's cells become the max
        # trusted distance (0.5), not the seductive 0.0.
        A = _rec("A", BASE)
        B = _near("B", BASE, 20)                       # 20/40 differ → d(A,B)=0.5 over 40 cols
        C = _rec("C", "AAAAA" + "-" * 35)              # 5 cols overlap A, identical → raw d=0, c=5
        dm = distance_matrix([A, B, C])
        assert dm["A", "B"] == pytest.approx(0.5)      # trusted pair kept as-is
        assert dm["C", "A"] == pytest.approx(0.5)      # low-overlap → filled, NOT the raw 0.0
        assert dm["C", "B"] == pytest.approx(0.5)


class TestLoadAndCheck:
    def test_load_combined_records_dedupes_and_labels(self, tmp_path):
        fa = tmp_path / "RPB2_genomic_combined.fasta"
        fa.write_text(
            ">NS26_RPB2_genomic\nACGT\n"
            ">MH1.1 [organism=Alternaria alternata]\nACGT\n"
            ">MH1.1 [organism=Alternaria alternata]\nACGA\n"      # duplicate id
        )
        recs, meta, warnings = load_combined_records(str(fa), {"NS26"})
        ids = [r.id for r in recs]
        assert len(set(ids)) == 3                       # duplicate disambiguated
        assert any("duplicate" in w for w in warnings)
        assert meta["NS26_RPB2_genomic"]["source"] == "isolate"
        assert meta["MH1.1"]["source"] == "tip"
        assert meta["MH1.1"]["organism"] == "Alternaria alternata"
        # The RENAMED record must carry correctly-aligned metadata too — meta is keyed by the new
        # id, so a desync would drop its source/organism (D-036 review: half-covered before).
        assert "MH1.1__dup2" in ids
        assert meta["MH1.1__dup2"]["source"] == "tip"
        assert meta["MH1.1__dup2"]["organism"] == "Alternaria alternata"

    def test_orthology_check_preserves_ids_through_stub_mafft(self, tmp_path, monkeypatch):
        # Stub MAFFT: copy input → output unchanged (records are already equal length). This locks
        # the contract the page depends on: ids survive the align round-trip so source labels are
        # correct, and the report is well-formed.
        fa = tmp_path / "RPB2_genomic_combined.fasta"
        fa.write_text(f">NS26_RPB2_genomic\n{BASE}\n>NS27_RPB2_genomic\n{_ortho('x',0).seq}\n"
                      f">MH1.1 [organism=A. b]\n{PARA}\n")

        def _stub_mafft(inp, outp, **k):
            __import__("shutil").copyfile(inp, outp)
            return 0, "ok"
        monkeypatch.setattr(ou, "run_mafft", _stub_mafft)

        rep = ou.orthology_check(str(fa), {"NS26", "NS27"})
        assert "error" not in rep
        ids = {r["id"] for r in rep["rows"]}
        assert ids == {"NS26_RPB2_genomic", "NS27_RPB2_genomic", "MH1.1"}
        src = {r["id"]: r["source"] for r in rep["rows"]}
        assert src["NS26_RPB2_genomic"] == "isolate" and src["MH1.1"] == "tip"

    def test_orthology_check_surfaces_mafft_failure(self, tmp_path, monkeypatch):
        fa = tmp_path / "RPB2_genomic_combined.fasta"
        fa.write_text(f">a\n{BASE}\n>b\n{PARA}\n>c\n{_ortho('x',0).seq}\n")
        monkeypatch.setattr(ou, "run_mafft", lambda inp, outp, **k: (1, "boom"))
        rep = ou.orthology_check(str(fa), set())
        assert "error" in rep and "MAFFT" in rep["error"]

    def test_orthology_check_too_few_sequences(self, tmp_path):
        fa = tmp_path / "RPB2_genomic_combined.fasta"
        fa.write_text(f">a\n{BASE}\n>b\n{PARA}\n")
        rep = ou.orthology_check(str(fa), set())
        assert "error" in rep and "≥3" in rep["error"]
