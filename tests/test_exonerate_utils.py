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

from phylofetch.blast_loci_utils import write_region_gff3
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


class TestRegionGff3:
    """D-029: the region-relative GFF3 annotates the extracted LOCUS_genomic.fasta directly —
    its exon ranges must equal the UPPERCASE (exon) runs of the soft-masked genomic, on both
    strands, and the seqid must match the genomic FASTA record id."""

    def _exon_ranges(self, gff_text):
        out = []
        for line in gff_text.splitlines():
            f = line.split("\t")
            if len(f) > 4 and f[2] == "exon":
                out.append((int(f[3]), int(f[4])))
        return out

    @pytest.mark.parametrize("tag,contig_fix", [("plus", "contig_plus"),
                                                ("minus", "contig_minus")])
    def test_exon_ranges_match_uppercase_runs(self, tag, contig_fix, request, tmp_path):
        import re
        contig = request.getfixturevalue(contig_fix)
        res = build_result_from_model(
            parse_exonerate_gff(_fix(f"exo_protein2genome_{tag}.gff"))[0], contig, "ST1", "TEF1")
        path = write_region_gff3(res, str(tmp_path))
        text = Path(path).read_text()
        amp = res["amplicon_seq"]
        upper_runs = [(m.start() + 1, m.end()) for m in re.finditer(r"[ACGT]+", amp)]
        assert self._exon_ranges(text) == upper_runs            # exon track lines up with the FASTA
        assert f"##sequence-region ST1_TEF1_genomic 1 {len(amp)}" in text
        assert Path(path).name == "TEF1_genomic.gff3"
        # everything is on '+' in the region frame (FASTA is coding-oriented)
        feat_strands = {l.split("\t")[6] for l in text.splitlines() if len(l.split("\t")) > 6}
        assert feat_strands == {"+"}

    def test_introns_are_lowercase_in_fasta(self, contig_minus, tmp_path):
        import re
        res = build_result_from_model(
            parse_exonerate_gff(_fix("exo_protein2genome_minus.gff"))[0], contig_minus, "ST1", "TEF1")
        text = Path(write_region_gff3(res, str(tmp_path))).read_text()
        amp = res["amplicon_seq"]
        for line in text.splitlines():
            f = line.split("\t")
            if len(f) > 4 and f[2] == "intron":
                assert amp[int(f[3]) - 1:int(f[4])].islower()


class TestRefineEscalation:
    """D-030: escalate --refine when the CDS is frameshifted; keep the cleanest, never drop."""

    def _setup(self, tmp_path, monkeypatch, clean_at):
        """Stub the Exonerate run/build chain so refine level `clean_at` (and higher) yields a
        clean CDS and lower levels frameshift. Returns (calls, ref, asm)."""
        import phylofetch.exonerate_utils as eu
        ref = tmp_path / "ref.fasta"; ref.write_text(">r\nMKLVAAAA\n")     # protein → protein2genome
        asm = tmp_path / "asm.fasta"; asm.write_text(">c\nACGTACGT\n")
        calls = []
        monkeypatch.setattr(eu.shutil, "which", lambda b: "/usr/bin/" + b)
        monkeypatch.setattr(eu, "run_exonerate",
                            lambda *a, refine, **k: (calls.append(refine), (0, refine, "p"))[1])
        monkeypatch.setattr(eu, "parse_exonerate_gff", lambda raw: [{"raw": raw, "contig": "c"}])
        monkeypatch.setattr(eu, "select_best_model", lambda ms, **k: (ms[0], []))
        monkeypatch.setattr(eu, "_load_contig", lambda *a: "A" * 60)
        monkeypatch.setattr(eu, "_probe_exonerate_version", lambda b: "x")
        ORDER = ["none", "region", "full"]

        def fake_build(model, *a, **k):
            lvl = model["raw"]
            clean = clean_at in ORDER and ORDER.index(lvl) >= ORDER.index(clean_at)
            return {"strain_id": "S", "locus_name": "RPB2", "cds_length": 300,
                    "n_internal_stops": 0 if clean else 9, "len_mod3": 0 if clean else 1}
        monkeypatch.setattr(eu, "build_result_from_model", fake_build)
        for w in ("write_exonerate_fastas", "write_gff3", "write_region_gff3",
                  "write_codon_partition", "write_exonerate_log"):
            monkeypatch.setattr(eu, w, lambda *a, **k: "")
        return calls, str(ref), str(asm), eu

    def test_escalates_until_clean(self, tmp_path, monkeypatch):
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="region")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none")
        assert res["n_internal_stops"] == 0 and res["len_mod3"] == 0      # recovered clean CDS
        assert calls == ["none", "region"]                               # stopped at first clean
        assert "boundary-refined: refine=region" in status

    def test_clean_first_pass_no_escalation(self, tmp_path, monkeypatch):
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="none")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none")
        assert calls == ["none"] and "boundary-refined" not in status     # single fast pass

    def test_unfixable_keeps_best_and_flags_not_dropped(self, tmp_path, monkeypatch):
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="never")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none")
        assert res is not None                                            # written, not dropped
        assert res["n_internal_stops"] == 9 and "QC: review" in status
        assert calls == ["none", "region", "full"]                        # tried everything

    def test_disabled_does_not_escalate(self, tmp_path, monkeypatch):
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="region")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none", escalate_refine=False)
        assert calls == ["none"] and res["n_internal_stops"] == 9         # stayed frameshifted


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
