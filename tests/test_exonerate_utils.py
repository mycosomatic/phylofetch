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

    def test_prefers_authoritative_tcs_over_frameshifted_coord_rebuild(self):
        # D-039: protein2genome trims the split codon that straddles an intron, so its authoritative
        # %tcs is in-frame — but the naive coordinate rebuild keeps those few bp and comes out
        # frameshifted with bogus internal stops (observed on real RPB2/TUB2: a 98.9%-identity model
        # whose %tcs is clean, written out as a 40-stop CDS). We must emit %tcs, not the rebuild.
        clean = "ATGAAACCCGGGTTTAAACCCGGGTTTGGG"        # 30 bp, in frame, no stops
        contig = clean + "GG"                            # +2 bp the rebuild wrongly keeps -> 32 bp

        def _model(tcs):
            return {"strand": "+", "cds": [(1, 32)], "exons": [], "introns": [],
                    "gene_start": 1, "gene_end": 32, "query_len": 10, "q_aln_begin": 0,
                    "q_aln_end": 10, "tcs": tcs, "model": "exonerate:protein2genome:local",
                    "query_id": "REF1", "contig": "c1", "score": 100.0, "pident": 99.0}

        # %tcs present and disagreeing with the rebuild -> emit the clean %tcs, flag the mismatch.
        res = build_result_from_model(_model(clean), contig, "ST1", "RPB2")
        assert res["cds_seq"] == clean
        assert res["cds_source"] == "tcs"
        assert res["tcs_matches"] is False
        assert res["len_mod3"] == 0 and res["n_internal_stops"] == 0

        # No %tcs (defensive fallback) -> the coordinate rebuild, which here is the frameshifted
        # 32 bp -> proves the rebuild really was the broken sequence we were shipping before.
        res2 = build_result_from_model(_model(""), contig, "ST1", "RPB2")
        assert res2["cds_seq"] == contig.upper()
        assert res2["cds_source"] == "coords"
        assert res2["len_mod3"] == 2

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
        # D-035: default escalation now caps at `region` (full is opt-in), so an unfixable CDS
        # tries none→region only — and is still kept-and-flagged, never dropped.
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="never")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none")
        assert res is not None                                            # written, not dropped
        assert res["n_internal_stops"] == 9 and "QC: review" in status
        assert calls == ["none", "region"]                                # capped at region (D-035)

    def test_default_ceiling_region_skips_full(self, tmp_path, monkeypatch):
        # The expensive `--refine full` pass must NOT run under the default ceiling (D-035).
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="never")
        eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                   narrow=False, refine="none")
        assert "full" not in calls

    def test_ceiling_full_escalates_to_full(self, tmp_path, monkeypatch):
        # The opt-in Thorough effort / deep-refine pass raises the ceiling to full.
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="never")
        res, _ = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                            narrow=False, refine="none",
                                            escalate_ceiling="full")
        assert calls == ["none", "region", "full"]                        # tried everything
        assert res is not None

    def test_disabled_does_not_escalate(self, tmp_path, monkeypatch):
        calls, ref, asm, eu = self._setup(tmp_path, monkeypatch, clean_at="region")
        res, status = eu.extract_locus_exonerate(asm, ref, str(tmp_path / "o"), "S", "RPB2",
                                                 narrow=False, refine="none", escalate_refine=False)
        assert calls == ["none"] and res["n_internal_stops"] == 9         # stayed frameshifted


class TestScanFlaggedCds:
    """D-035: stateless disk scan that finds flagged CDS for the targeted deep-refine pass."""

    def _write_cds(self, root, strain, locus, seq):
        d = root / strain / locus
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{locus}_CDS.fasta").write_text(f">{strain}_{locus}\n{seq}\n")
        return d

    def test_flags_internal_stop_and_frameshift_omits_clean(self, tmp_path):
        from phylofetch.exonerate_utils import scan_flagged_cds
        self._write_cds(tmp_path, "S1", "RPB2", "ATGTAAAAA")       # internal stop (TAA), len%3=0
        self._write_cds(tmp_path, "S2", "TUB2", "ATGAAAC")         # len%3 = 1 (frameshift)
        self._write_cds(tmp_path, "S3", "TEF1", "ATGAAACCCGGG")    # clean: M K P G
        rows = scan_flagged_cds(str(tmp_path))
        flagged = {(r["strain"], r["locus"]) for r in rows}
        assert ("S1", "RPB2") in flagged
        assert ("S2", "TUB2") in flagged
        assert ("S3", "TEF1") not in flagged           # clean → omitted
        rpb2 = next(r for r in rows if r["locus"] == "RPB2")
        assert rpb2["n_internal_stops"] == 1 and rpb2["len_mod3"] == 0

    def test_clean_cds_ending_in_terminal_stop_is_omitted(self, tmp_path):
        # ATG AAA TAA = M K * — a clean CDS *with* a terminal stop. validate_cds must not count the
        # terminal stop as internal, so the scan omits it (guards the subtraction the scan relies on).
        from phylofetch.exonerate_utils import scan_flagged_cds
        self._write_cds(tmp_path, "S1", "RPB2", "ATGAAATAA")
        assert all(r["locus"] != "RPB2" for r in scan_flagged_cds(str(tmp_path)))

    def test_reads_refine_used_from_log_real_format(self, tmp_path):
        # Use the ACTUAL log format write_exonerate_log emits (D-035 review: the old test
        # fabricated a token the writer never wrote).
        from phylofetch.exonerate_utils import scan_flagged_cds
        d = self._write_cds(tmp_path, "S1", "RPB2", "ATGTAAAAA")
        (d / "RPB2_extraction.log").write_text(
            "[QC]\n  reading_frame len % 3 == 0 (OK)\n  refine_used    region\n"
            "  verdict        REVIEW\n")
        rows = scan_flagged_cds(str(tmp_path))
        assert rows[0]["refine_used"] == "region"

    def test_write_exonerate_log_actually_records_refine_used(self, tmp_path):
        # The producer side: write_exonerate_log must emit the refine_used line the scan reads.
        from phylofetch.exonerate_utils import write_exonerate_log
        result = {
            "strain_id": "S1", "locus_name": "RPB2", "n_introns": 0, "introns": [],
            "len_mod3": 0, "n_internal_stops": 0, "blast_type": "exonerate",
            "contig_id": "c1", "amplicon_start": 1, "amplicon_end": 99, "strand": "+",
            "n_exons": 1, "cds_length": 99, "amplicon_length": 99, "refine_used": "region",
        }
        log = write_exonerate_log(result, str(tmp_path))
        with open(log) as f:
            assert "refine_used    region" in f.read()

    def test_revalidates_under_logged_genetic_code_not_the_argument(self, tmp_path):
        # ATG TGA AAA = M * K under standard code 1 (TGA = internal stop → flagged), but M W K under
        # mold-mito code 4 (TGA = Trp → clean). The scan must re-validate under the code the locus was
        # EXTRACTED with (logged), not the caller's current widget, so a code change alone can't flip
        # the verdict (D-036 review). Here the log says 4, the arg says 1 → recorded code 4 wins → omit.
        from phylofetch.exonerate_utils import scan_flagged_cds
        d = self._write_cds(tmp_path, "S1", "RPB2", "ATGTGAAAA")
        (d / "RPB2_extraction.log").write_text("[QC]\n  geneticcode    4\n  verdict PASS\n")
        assert all(r["locus"] != "RPB2" for r in scan_flagged_cds(str(tmp_path), geneticcode=1))

    def test_falls_back_to_argument_code_when_log_has_none(self, tmp_path):
        # No geneticcode recorded (pre-fix extraction) → fall back to the caller's default; under
        # code 1 the same TGA-bearing CDS is flagged, and the row carries the code actually used.
        from phylofetch.exonerate_utils import scan_flagged_cds
        self._write_cds(tmp_path, "S1", "RPB2", "ATGTGAAAA")        # no log written
        rows = scan_flagged_cds(str(tmp_path), geneticcode=1)
        rpb2 = next(r for r in rows if r["locus"] == "RPB2")
        assert rpb2["n_internal_stops"] == 1 and rpb2["geneticcode"] == 1

    def test_write_exonerate_log_records_genetic_code(self, tmp_path):
        # Producer side: the code must be persisted for the scan above to read it back.
        from phylofetch.exonerate_utils import write_exonerate_log
        result = {
            "strain_id": "S1", "locus_name": "RPB2", "n_introns": 0, "introns": [],
            "len_mod3": 0, "n_internal_stops": 0, "blast_type": "exonerate",
            "contig_id": "c1", "amplicon_start": 1, "amplicon_end": 99, "strand": "+",
            "n_exons": 1, "cds_length": 99, "amplicon_length": 99, "refine_used": "none",
            "geneticcode": 4,
        }
        with open(write_exonerate_log(result, str(tmp_path))) as f:
            assert "geneticcode    4" in f.read()

    def test_empty_or_missing_dir_returns_empty(self, tmp_path):
        from phylofetch.exonerate_utils import scan_flagged_cds
        assert scan_flagged_cds(str(tmp_path / "nope")) == []
        assert scan_flagged_cds(str(tmp_path)) == []


class TestResolveGuidePath:
    """D-035 / D-036 review: guide lookup for the deep-refine re-run."""

    def _scratch(self, root, locus):
        d = root / "scratch" / "guides"
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{locus}_guide.fasta").write_text(">g\nMKL\n")

    def test_prefers_scratch_guide_for_catalogue_locus(self, tmp_path):
        from phylofetch.exonerate_utils import resolve_guide_path
        self._scratch(tmp_path, "RPB2")
        got = resolve_guide_path("RPB2", str(tmp_path), tmp_path / "goi", tmp_path / "refs")
        assert got.endswith("RPB2_guide.fasta")

    def test_goi_takes_precedence_over_same_named_catalogue_guide(self, tmp_path):
        from phylofetch.exonerate_utils import resolve_guide_path
        self._scratch(tmp_path, "RPB2")                       # catalogue guide exists
        goi = tmp_path / "goi"
        goi.mkdir()
        (goi / "RPB2.fasta").write_text(">user\nMKLV\n")       # user's GOI ortholog
        got = resolve_guide_path("RPB2", str(tmp_path), goi, tmp_path / "refs", is_goi=True)
        assert got == str(goi / "RPB2.fasta")                 # GOI wins

    def test_none_when_no_guide(self, tmp_path):
        from phylofetch.exonerate_utils import resolve_guide_path
        assert resolve_guide_path("XX", str(tmp_path), tmp_path / "goi", tmp_path / "refs") is None

    def test_no_directory_side_effect(self, tmp_path):
        from phylofetch.exonerate_utils import resolve_guide_path
        ref = tmp_path / "refs"
        resolve_guide_path("ZZ", str(tmp_path), tmp_path / "goi", ref)
        assert not (ref / "ZZ").exists()                      # must not create the dir


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
