"""
Tests for project persistence in src/phylofetch/project_manager.py
(assembly registry save/load, manifest, project listing).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from phylofetch.project_manager import (
    ASSEMBLY_MANIFEST_FIELDS,
    init_project,
    list_projects,
    load_assembly_registry,
    save_assembly_registry,
)


# Two records exercising both schemas the app produces.
_REGISTRY = {
    # Assembly Manager schema: nested stats{} with QUAST
    "NS26-3-C2": {
        "strain_id": "NS26-3-C2",
        "assembly_path": "/data/NS26-3-C2/NS26-3-C2_final_EGAP_assembly.fasta",
        "reads_r1": "/data/NS26-3-C2/r1.fastq",
        "reads_r2": "",
        "registered_at": "2026-06-18T00:00:00+00:00",
        "stats": {
            "assembler": "unknown",
            "num_contigs": 118,
            "n50": 1204567,
            "total_length_mb": 34.52,
            "mean_gc": 51.23,
            "quast": {"N50": 1204567},
            "quast_report": "/data/NS26-3-C2/..._quast/report.tsv",
        },
    },
    # Project Setup schema: flat keys, busco_dir linked
    "S9-1B-A2": {
        "assembly_path": "/data/S9-1B-A2/S9-1B-A2_final_EGAP_assembly.fasta",
        "assembler": "spades",
        "num_contigs": 90,
        "n50": 900000,
        "total_length_mb": 33.1,
        "gc_percent": 50.8,
        "busco_dir": "/data/S9-1B-A2/..._ascomycota_busco",
        "quast_report": "/data/S9-1B-A2/..._quast/report.tsv",
    },
}


class TestAssemblyRegistry:
    def test_roundtrip(self, tmp_path):
        proj = tmp_path / "proj1"
        save_assembly_registry(proj, _REGISTRY)
        loaded = load_assembly_registry(proj)
        assert set(loaded.keys()) == {"NS26-3-C2", "S9-1B-A2"}
        assert loaded["NS26-3-C2"]["stats"]["n50"] == 1204567

    def test_load_missing_is_empty(self, tmp_path):
        assert load_assembly_registry(tmp_path / "nope") == {}

    def test_manifest_tsv_written(self, tmp_path):
        proj = tmp_path / "proj2"
        save_assembly_registry(proj, _REGISTRY)
        manifest = proj / "metadata" / "assembly_manifest.tsv"
        assert manifest.exists()
        text = manifest.read_text()
        header = text.splitlines()[0].split("\t")
        assert header == ASSEMBLY_MANIFEST_FIELDS
        # Both schemas flatten correctly
        assert "NS26-3-C2" in text
        assert "1204567" in text          # n50 pulled from nested stats
        assert "spades" in text           # assembler pulled from flat record
        assert "report.tsv" in text       # QUAST report path logged
        assert "_ascomycota_busco" in text  # busco_dir logged from flat record

    def test_quast_presence_fallback(self, tmp_path):
        # Record with a parsed QUAST dict but no stored report path → "(found)"
        proj = tmp_path / "proj_q"
        save_assembly_registry(proj, {
            "X": {"assembly_path": "/d/x.fasta", "stats": {"quast": {"N50": 1}}},
        })
        text = (proj / "metadata" / "assembly_manifest.tsv").read_text()
        assert "(found)" in text

    def test_manifest_has_one_row_per_assembly(self, tmp_path):
        proj = tmp_path / "proj3"
        save_assembly_registry(proj, _REGISTRY)
        lines = (proj / "metadata" / "assembly_manifest.tsv").read_text().splitlines()
        assert len(lines) == 1 + len(_REGISTRY)  # header + 2 rows


class TestListProjects:
    def test_lists_created_projects(self, tmp_path):
        init_project(tmp_path / "alpha")
        save_assembly_registry(tmp_path / "beta", _REGISTRY)
        projects = list_projects(tmp_path)
        names = {p["name"] for p in projects}
        assert names == {"alpha", "beta"}
        beta = next(p for p in projects if p["name"] == "beta")
        assert beta["n_assemblies"] == 2

    def test_empty_root(self, tmp_path):
        assert list_projects(tmp_path / "does_not_exist") == []
