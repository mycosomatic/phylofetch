"""
Tests for src/phylofetch/exonerate_utils.py

Covers the RM-002 / D-008 work: frame-safe CDS extraction via Exonerate spliced
alignment. Parser and result-building tests run **offline** against verbatim
Exonerate 2.4.0 output fixtures (tests/fixtures/exo_*.gff) generated from a
deterministic synthetic 2-exon gene, so CI without exonerate still exercises the
parsing/coordinate logic. The end-to-end class is guarded by shutil.which and
runs only where the exonerate binary is installed.
"""

import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from Bio.Seq import Seq

from phylofetch.exonerate_utils import (
    MODEL_FOR_QUERYTYPE,
    build_result_from_model,
    extract_locus_exonerate,
    parse_exonerate_gff,
    select_best_model,
    soft_mask_genomic,
    validate_cds,
)

FIX = Path(__file__).parent / "fixtures"
HAVE_EXONERATE = shutil.which("exonerate") is not None


def _fix(name: str) -> str:
    return (FIX / name).read_text()


@pytest.fixture(scope="module")
def expected_cds() -> str:
    return _fix("expected_cds.txt").strip()


@pytest.fixture(scope="module")
def contig_plus() -> str:
    rec = next(iter(_read_fasta(FIX / "contig_plus.fasta")))
    return rec


@pytest.fixture(scope="module")
def contig_minus() -> str:
    return next(iter(_read_fasta(FIX / "contig_minus.fasta")))


def _read_fasta(path: Path) -> list[str]:
    """Tiny FASTA reader → list of sequence strings (no Bio dependency churn)."""
    seqs, cur = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


# ── validate_cds (the RM-002 QC primitive) ─────────────────────────────────

class TestValidateCds:
    def test_clean_cds(self):
        v = validate_cds("ATGAAACTTGTTGCT")          # M K L V A
        assert v["clean"] is True
        assert v["len_mod3"] == 0
        assert v["n_internal_stops"] == 0
        assert v["protein"] == "MKLVA"
        assert v["starts_with_atg"] is True

    def test_internal_stop_detected(self):
        v = validate_cds("ATGTAAGGGCCC")             # M * G P
        assert v["n_internal_stops"] == 1
        assert v["clean"] is False

    def test_trailing_stop_is_not_internal(self):
        v = validate_cds("ATGAAATAA")                # M K *
        assert v["ends_with_stop"] is True
        assert v["n_internal_stops"] == 0

    def test_frameshift_length(self):
        v = validate_cds("ATGAAACTTGT")              # 11 bp
        assert v["len_mod3"] == 2
        assert v["clean"] is False

    def test_empty(self):
        v = validate_cds("")
        assert v["length"] == 0 and v["protein"] == ""


# ── Model selection from query alphabet ────────────────────────────────────

def test_model_for_querytype():
    assert MODEL_FOR_QUERYTYPE["protein"] == "protein2genome"
    assert MODEL_FOR_QUERYTYPE["nucleotide"] == "coding2genome"


# ── GFF + RYO parsing (offline, verbatim exonerate 2.4.0 fixtures) ─────────

class TestParseExonerateGff:
    def test_parse_plus(self):
        models = parse_exonerate_gff(_fix("exo_protein2genome_plus.gff"))
        assert len(models) == 1
        m = models[0]
        assert m["contig"] == "contig_plus"
        assert m["strand"] == "+"
        assert m["model"].startswith("exonerate:protein2genome")
        assert m["cds"] == [(46, 195), (254, 403)]
        assert m["introns"] == [(196, 253)]
        assert m["gene_start"] == 46 and m["gene_end"] == 403
        assert m["score"] == pytest.approx(507.0)
        assert m["pident"] == pytest.approx(100.0)
        assert m["tcs"].startswith("ATG") and len(m["tcs"]) == 300

    def test_parse_minus_strand(self):
        models = parse_exonerate_gff(_fix("exo_protein2genome_minus.gff"))
        assert len(models) == 1
        m = models[0]
        assert m["strand"] == "-"
        # both exons present regardless of emission order
        assert sorted(m["cds"]) == [(46, 195), (254, 403)]
        assert m["introns"] == [(196, 253)]

    def test_parse_coding2genome(self):
        models = parse_exonerate_gff(_fix("exo_coding2genome_plus.gff"))
        assert len(models) == 1
        assert models[0]["model"].startswith("exonerate:coding2genome")
        assert models[0]["cds"] == [(46, 195), (254, 403)]

    def test_parse_multi_results(self):
        models = parse_exonerate_gff(_fix("exo_protein2genome_multi.gff"))
        # tandem-duplicated gene → at least the two full-length copies
        assert len(models) >= 2
        top_scores = sorted((m["score"] for m in models), reverse=True)[:2]
        assert all(s == pytest.approx(507.0) for s in top_scores)

    def test_parse_empty_text(self):
        assert parse_exonerate_gff("") == []
        assert parse_exonerate_gff("no gff here\n-- completed exonerate analysis\n") == []


# ── select_best_model (paralog awareness) ──────────────────────────────────

class TestSelectBestModel:
    def test_picks_highest_score_and_returns_others(self):
        models = parse_exonerate_gff(_fix("exo_protein2genome_multi.gff"))
        best, others = select_best_model(models)
        assert best is not None
        assert best["score"] == max(m["score"] for m in models)
        assert len(others) == len(models) - 1

    def test_empty(self):
        assert select_best_model([]) == (None, [])

    def test_min_score_filters_all(self):
        models = parse_exonerate_gff(_fix("exo_protein2genome_plus.gff"))
        best, others = select_best_model(models, min_score=10_000)
        assert best is None and others == []


# ── build_result_from_model: coordinate / strand / QC correctness ──────────

class TestBuildResultFromModel:
    def test_plus_strand_rebuilds_expected_cds(self, contig_plus, expected_cds):
        models = parse_exonerate_gff(_fix("exo_protein2genome_plus.gff"))
        res = build_result_from_model(models[0], contig_plus, "ST1", "TEF1")
        assert res["cds_seq"] == expected_cds
        assert res["strand"] == "+"
        assert res["n_exons"] == 2
        assert res["n_introns"] == 1
        assert res["cds_length"] == 300
        assert res["len_mod3"] == 0
        assert res["n_internal_stops"] == 0
        assert res["tcs_matches"] is True
        # explicit intron splice sites are canonical GT-AG
        assert res["introns"][0]["splice_5"] == "GT"
        assert res["introns"][0]["splice_3"] == "AG"
        # genomic amplicon spans first..last exon
        assert res["amplicon_start"] == 46 and res["amplicon_end"] == 403

    def test_minus_strand_rebuilds_same_cds(self, contig_minus, expected_cds):
        models = parse_exonerate_gff(_fix("exo_protein2genome_minus.gff"))
        res = build_result_from_model(models[0], contig_minus, "ST1", "TEF1")
        assert res["strand"] == "-"
        # minus-strand exons reverse-complemented in transcript order → same CDS
        assert res["cds_seq"] == expected_cds
        assert res["tcs_matches"] is True
        assert res["introns"][0]["splice_5"] == "GT"
        assert res["introns"][0]["splice_3"] == "AG"

    def test_protein_translation_present(self, contig_plus):
        models = parse_exonerate_gff(_fix("exo_protein2genome_plus.gff"))
        res = build_result_from_model(models[0], contig_plus, "ST1", "TEF1")
        expected_prot = _fix("expected_protein.txt").strip()
        assert res["protein"] == expected_prot

    def test_query_coverage_computed(self, contig_plus):
        models = parse_exonerate_gff(_fix("exo_protein2genome_plus.gff"))
        res = build_result_from_model(models[0], contig_plus, "ST1", "TEF1")
        assert res["query_coverage"] == pytest.approx(100.0)

    @pytest.mark.parametrize("tag,contig_fix", [("plus", "contig_plus"),
                                                ("minus", "contig_minus")])
    def test_genomic_amplicon_is_soft_masked(self, tag, contig_fix, request):
        # The genomic product is exon-UPPER / intron-lower so isolates and tips carry the same
        # boundary annotation (D-022). Bases/length are unchanged; only case differs.
        contig = request.getfixturevalue(contig_fix)
        res = build_result_from_model(
            parse_exonerate_gff(_fix(f"exo_protein2genome_{tag}.gff"))[0], contig, "ST1", "TEF1")
        amp = res["amplicon_seq"]
        assert len(amp) == res["amplicon_length"]
        assert sum(c.islower() for c in amp) == res["introns"][0]["length"]   # introns masked
        assert sum(c.isupper() for c in amp) == res["cds_length"]             # exons unmasked


class TestSoftMaskGenomic:
    def test_plus_lowercases_intron_span(self):
        # intron at 1-based positions 4..6 of a 12-mer
        assert soft_mask_genomic("AAACCCGGGTTT", 1, 12, [(4, 6)], False) == "AAAcccGGGTTT"

    def test_minus_strand_preserves_case_through_revcomp(self):
        out = soft_mask_genomic("AAACCCGGGTTT", 1, 12, [(4, 6)], True)
        assert sum(c.islower() for c in out) == 3
        assert out.upper() == str(Seq("AAACCCGGGTTT").reverse_complement())

    def test_empty_span_returns_empty(self):
        assert soft_mask_genomic("ACGT", 0, 0, [], False) == ""


# ── End-to-end (requires the exonerate binary) ─────────────────────────────

@pytest.mark.skipif(not HAVE_EXONERATE, reason="exonerate not installed")
class TestExonerateIntegration:
    def _run(self, tmp_path, target, ref, **kw):
        return extract_locus_exonerate(
            assembly_fasta=str(FIX / target),
            reference_fasta=str(FIX / ref),
            output_dir=str(tmp_path), strain_id="ST1", locus_name="TEF1",
            narrow=False, **kw,
        )

    def test_protein2genome_plus_end_to_end(self, tmp_path, expected_cds):
        res, status = self._run(tmp_path, "contig_plus.fasta", "ref_protein.fasta")
        assert res is not None, status
        assert status.startswith("OK")
        assert res["cds_seq"] == expected_cds
        assert res["n_internal_stops"] == 0 and res["len_mod3"] == 0
        assert res["introns"][0]["splice_5"] == "GT"
        assert res["introns"][0]["splice_3"] == "AG"
        # output files written
        for fname in ["TEF1_CDS.fasta", "TEF1_protein.fasta", "TEF1_genomic.fasta",
                      "TEF1_introns.fasta", "TEF1.gff3", "TEF1_partition.nex",
                      "TEF1_extraction.log"]:
            assert (tmp_path / fname).exists(), f"missing {fname}"

    def test_protein2genome_minus_end_to_end(self, tmp_path, expected_cds):
        res, status = self._run(tmp_path, "contig_minus.fasta", "ref_protein.fasta")
        assert res is not None, status
        assert res["strand"] == "-"
        assert res["cds_seq"] == expected_cds

    def test_coding2genome_nucleotide_query(self, tmp_path, expected_cds):
        res, status = self._run(tmp_path, "contig_plus.fasta", "ref_cds.fasta")
        assert res is not None, status
        assert res["blast_type"].startswith("exonerate:coding2genome")
        assert res["cds_seq"] == expected_cds

    def test_hybrid_narrowing_selects_correct_contig(self, tmp_path, expected_cds):
        # narrow=True: tblastn must pick the gene contig out of a 2-contig assembly.
        # Build a 2-contig assembly on the fly.
        decoy = "ACGTTGCA" * 200
        gene = next(iter(_read_fasta(FIX / "contig_plus.fasta")))
        asm = tmp_path / "asm.fasta"
        asm.write_text(f">decoy\n{decoy}\n>gene_contig\n{gene}\n")
        res, status = extract_locus_exonerate(
            assembly_fasta=str(asm), reference_fasta=str(FIX / "ref_protein.fasta"),
            output_dir=str(tmp_path), strain_id="ST1", locus_name="TEF1", narrow=True,
        )
        assert res is not None, status
        assert res["contig_id"] == "gene_contig"
        assert res["cds_seq"] == expected_cds
        assert (tmp_path / "_exonerate_target.fasta").exists()

    def test_missing_binary_returns_error(self, tmp_path):
        res, status = extract_locus_exonerate(
            assembly_fasta=str(FIX / "contig_plus.fasta"),
            reference_fasta=str(FIX / "ref_protein.fasta"),
            output_dir=str(tmp_path), strain_id="ST1", locus_name="TEF1",
            exonerate_bin="exonerate_not_installed_xyz",
        )
        assert res is None
        assert "not found" in status.lower()
