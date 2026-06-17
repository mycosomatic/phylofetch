"""
Shared pytest fixtures for phylofetch tests.
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


# ── Synthetic FASTA helpers ────────────────────────────────────────────────

def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    """Write a list of (id, sequence) pairs to a FASTA file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(
        [SeqRecord(Seq(seq), id=rid, description="") for rid, seq in records],
        str(path),
        "fasta",
    )
    return path


@pytest.fixture
def tmp_fasta(tmp_path):
    """Factory: create a temporary FASTA file from (id, seq) pairs."""
    def _make(name: str, records: list[tuple[str, str]]) -> Path:
        return _write_fasta(tmp_path / name, records)
    return _make


@pytest.fixture
def simple_assembly(tmp_path) -> Path:
    """A small synthetic assembly with two contigs."""
    contig1 = "ATGCATGCATGC" * 100   # 1200 bp
    contig2 = "GTATGCTAGCTA" * 50    # 600 bp
    fp = tmp_path / "assembly.fasta"
    _write_fasta(fp, [("contig_1", contig1), ("contig_2", contig2)])
    return fp


# ── Mock BLAST hit dicts ───────────────────────────────────────────────────

def _make_hit(
    qseqid: str,
    sseqid: str,
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    bitscore: float = 200.0,
    pident: float = 95.0,
    evalue: float = 1e-50,
    strand: str = "+",
) -> dict:
    return {
        "qseqid":  qseqid,
        "sseqid":  sseqid,
        "qstart":  qstart,
        "qend":    qend,
        "sstart":  sstart if strand == "+" else send,
        "send":    send if strand == "+" else sstart,
        "bitscore": bitscore,
        "pident":  pident,
        "evalue":  evalue,
        "length":  abs(send - sstart) + 1,
        "qlen":    qend,
    }


@pytest.fixture
def single_ref_hsps():
    """Two HSPs from the same reference on the same contig — should be stitched."""
    return [
        _make_hit("ref1", "contig_1", 1,   250, 1001, 1250, bitscore=300.0),
        _make_hit("ref1", "contig_1", 251, 500, 1351, 1600, bitscore=280.0),
    ]


@pytest.fixture
def two_ref_hsps():
    """
    HSPs from TWO different references on the same contig.
    LXD-002: only one reference's HSPs must win; they must NOT be merged.
    """
    return [
        # ref1 — lower total bitscore
        _make_hit("ref1", "contig_1", 1,   250, 2001, 2250, bitscore=150.0),
        _make_hit("ref1", "contig_1", 251, 500, 2351, 2600, bitscore=140.0),
        # ref2 — higher total bitscore → should win
        _make_hit("ref2", "contig_1", 1,   250, 1001, 1250, bitscore=300.0),
        _make_hit("ref2", "contig_1", 251, 500, 1351, 1600, bitscore=280.0),
    ]
