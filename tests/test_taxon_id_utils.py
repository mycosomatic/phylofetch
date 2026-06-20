"""
Tests for src/phylofetch/taxon_id_utils.py — ITS-based provisional taxon ID (RM-007 step 3).

Network-free: covers remote-BLAST command construction, tabular-output parsing, organism
ranking/dedup, and longest-ITS-record selection. The live ITSx + remote-BLAST orchestration
(identify_taxon_from_assembly) needs binaries + internet and is not exercised here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.taxon_id_utils import (
    REMOTE_BLAST_FIELDS,
    _pick_longest_record,
    build_remote_blast_cmd,
    parse_blast_organism_hits,
    rank_unique_organisms,
)


class TestBuildRemoteBlastCmd:
    def test_core_remote_flags(self):
        cmd = build_remote_blast_cmd("q.fasta", "out.tsv")
        assert "-remote" in cmd
        assert cmd[cmd.index("-db") + 1] == "nt"
        assert "-num_threads" not in cmd               # rejected with -remote
        assert cmd[cmd.index("-query") + 1] == "q.fasta"
        assert cmd[cmd.index("-out") + 1] == "out.tsv"
        # outfmt carries the organism columns we parse
        outfmt = cmd[cmd.index("-outfmt") + 1]
        assert outfmt.startswith("6 ")
        assert "sscinames" in outfmt and "staxids" in outfmt

    def test_entrez_query_included_and_omittable(self):
        with_q = build_remote_blast_cmd("q.fasta", "o.tsv", entrez_query="fungi[ORGN]")
        assert with_q[with_q.index("-entrez_query") + 1] == "fungi[ORGN]"
        no_q = build_remote_blast_cmd("q.fasta", "o.tsv", entrez_query="")
        assert "-entrez_query" not in no_q
        blank_q = build_remote_blast_cmd("q.fasta", "o.tsv", entrez_query="   ")
        assert "-entrez_query" not in blank_q

    def test_params_threaded_through(self):
        cmd = build_remote_blast_cmd("q.fasta", "o.tsv", evalue=1e-5,
                                     max_target_seqs=5, task="blastn")
        assert cmd[cmd.index("-evalue") + 1] == "1e-05"
        assert cmd[cmd.index("-max_target_seqs") + 1] == "5"
        assert cmd[cmd.index("-task") + 1] == "blastn"


def _row(qseqid, acc, pident, length, evalue, bits, taxids, sciname, stitle):
    return "\t".join([qseqid, acc, str(pident), str(length), evalue, str(bits),
                      taxids, sciname, stitle])


class TestParseBlastOrganismHits:
    def test_parses_well_formed_rows(self):
        text = "\n".join([
            _row("q1", "MN123.1", 99.5, 600, "0.0", 1100, "5599",
                 "Alternaria alternata", "Alternaria alternata strain A 18S ... ITS"),
            _row("q1", "MN999.1", 98.0, 590, "0.0", 1050, "5599",
                 "Alternaria arborescens", "Alternaria arborescens isolate B ITS"),
        ])
        hits = parse_blast_organism_hits(text)
        assert len(hits) == 2
        assert hits[0]["accession"] == "MN123.1"
        assert hits[0]["pident"] == 99.5
        assert hits[0]["length"] == 600
        assert hits[0]["bitscore"] == 1100.0
        assert hits[0]["staxids"] == "5599"
        assert hits[0]["organism"] == "Alternaria alternata"      # from sscinames

    def test_skips_blank_comment_and_short_rows(self):
        text = "\n".join([
            "# a comment", "", "too\tfew\tcols",
            _row("q1", "X.1", 97.0, 500, "1e-50", 900, "1",
                 "Fusarium oxysporum", "Fusarium oxysporum strain Z ITS"),
        ])
        hits = parse_blast_organism_hits(text)
        assert len(hits) == 1
        assert hits[0]["organism"] == "Fusarium oxysporum"

    def test_non_numeric_pident_row_dropped(self):
        text = _row("q1", "X.1", "NA", 500, "1e-50", 900, "1", "Some fungus", "Some fungus ITS")
        assert parse_blast_organism_hits(text) == []

    def test_fields_constant_matches_parser(self):
        assert REMOTE_BLAST_FIELDS[-1] == "stitle"
        assert REMOTE_BLAST_FIELDS[-2] == "sscinames"
        assert len(REMOTE_BLAST_FIELDS) == 9


class TestOrganismResolution:
    """sscinames is N/A without a local taxdb (the common case) → derive from stitle."""

    def test_falls_back_to_title_when_sscinames_na(self):
        text = _row("q1", "KX463014.1", 100.0, 518, "0.0", 957, "5599",
                    "N/A", "Alternaria alternata strain ALT1246 18S ribosomal RNA gene, ITS")
        hit = parse_blast_organism_hits(text)[0]
        assert hit["organism"] == "Alternaria alternata"
        assert hit["sscinames"] == "N/A"

    def test_title_sp_and_uncultured_and_status_tag(self):
        rows = {
            "Alternaria sp. ALT1 internal transcribed spacer":          "Alternaria sp.",
            "Uncultured fungus clone 7 ITS region":                     "Uncultured fungus",
            "UNVERIFIED: Cladosporium herbarum isolate Q ITS":          "Cladosporium herbarum",
        }
        for title, expected in rows.items():
            text = _row("q", "A.1", 99.0, 500, "0.0", 900, "1", "N/A", title)
            assert parse_blast_organism_hits(text)[0]["organism"] == expected

    def test_sscinames_preferred_over_title_when_present(self):
        text = _row("q", "A.1", 99.0, 500, "0.0", 900, "1",
                    "Fusarium oxysporum", "Fungal sp. clone messy title")
        assert parse_blast_organism_hits(text)[0]["organism"] == "Fusarium oxysporum"

    def test_joined_na_sscinames_falls_back_to_title(self):
        # subject identical to several taxids → "N/A;N/A" without a local taxdb
        text = _row("q", "KX463014.1", 100.0, 518, "0.0", 957, "5599;5599",
                    "N/A;N/A", "Alternaria alternata strain ALT1246 18S ribosomal RNA gene")
        assert parse_blast_organism_hits(text)[0]["organism"] == "Alternaria alternata"


class TestRankUniqueOrganisms:
    def test_sorts_by_bitscore_and_dedupes_organism(self):
        hits = [
            {"organism": "Alternaria alternata", "bitscore": 900.0, "pident": 98.0},
            {"organism": "Alternaria arborescens", "bitscore": 1100.0, "pident": 99.0},
            {"organism": "Alternaria alternata", "bitscore": 1200.0, "pident": 99.9},
        ]
        ranked = rank_unique_organisms(hits)
        assert [h["organism"] for h in ranked] == [
            "Alternaria alternata", "Alternaria arborescens"]
        # best (highest-bitscore) row kept for the deduped organism
        assert ranked[0]["bitscore"] == 1200.0

    def test_case_insensitive_dedupe_and_blank_dropped(self):
        hits = [
            {"organism": "Fusarium oxysporum", "bitscore": 500.0},
            {"organism": "fusarium oxysporum", "bitscore": 400.0},
            {"organism": "", "bitscore": 999.0},
        ]
        ranked = rank_unique_organisms(hits)
        assert len(ranked) == 1
        assert ranked[0]["organism"] == "Fusarium oxysporum"


class TestPickLongestRecord:
    def test_picks_longest(self, tmp_path):
        fa = tmp_path / "its.fasta"
        fa.write_text(">a\nACGT\n>b\nACGTACGTAC\n>c\nACG\n")
        rec = _pick_longest_record(str(fa))
        assert rec.id == "b"
        assert len(rec.seq) == 10

    def test_empty_returns_none(self, tmp_path):
        fa = tmp_path / "empty.fasta"
        fa.write_text("")
        assert _pick_longest_record(str(fa)) is None
