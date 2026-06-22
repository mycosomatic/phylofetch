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
    PROJECT_MANIFEST_SCHEMA_VERSION,
    WORKFLOW_STEPS,
    _migrate_assembly_record,
    clear_project_data,
    delete_project,
    effective_taxon,
    get_workflow,
    init_project,
    list_projects,
    load_assembly_registry,
    load_json,
    load_project_manifest,
    project_data_summary,
    project_output_dir,
    reset_workflow,
    set_output_dir,
    save_assembly_registry,
    save_json,
    set_assembly_taxon,
    set_default_taxon,
    set_workflow_loci,
    set_workflow_strategy,
    update_step,
)


class TestOutputDir:
    def test_default_is_project_results(self, tmp_path):
        p = tmp_path / "proj"; init_project(p)
        assert project_output_dir(p) == p / "results"

    def test_set_and_clear_override(self, tmp_path):
        p = tmp_path / "proj"; init_project(p)
        custom = tmp_path / "shared" / "out"
        set_output_dir(p, str(custom))
        assert project_output_dir(p) == custom
        set_output_dir(p, "")                       # blank reverts to default
        assert project_output_dir(p) == p / "results"

    def test_manifest_has_output_dir_default(self, tmp_path):
        assert load_project_manifest(tmp_path / "ghost")["output_dir"] == ""

    def test_summary_honours_override(self, tmp_path):
        p = tmp_path / "proj"; init_project(p)
        custom = tmp_path / "out"; set_output_dir(p, str(custom))
        (custom / "loci" / "combined").mkdir(parents=True)
        (custom / "loci" / "combined" / "ITS_combined.fasta").write_text(">x\nACGT\n")
        s = project_data_summary(p)
        assert s["n_combined"] == 1 and s["output_dir"] == str(custom)


class TestProjectDataManagement:
    def _proj(self, tmp_path):
        p = tmp_path / "proj"
        init_project(p)
        (p / "references" / "ITS").mkdir(parents=True)
        (p / "references" / "ITS" / "ITS_refs.fasta").write_text(">x\nACGT\n")
        (p / "results" / "loci" / "combined").mkdir(parents=True)
        (p / "results" / "loci" / "combined" / "ITS_combined.fasta").write_text(">x\nACGT\n")
        (p / "runs" / "r1").mkdir(parents=True)
        (p / "runs" / "r1" / "x.log").write_text("hi")
        return p

    def test_summary_counts(self, tmp_path):
        s = project_data_summary(self._proj(tmp_path))
        assert s["n_ref_loci"] == 1 and s["n_combined"] == 1 and s["n_runs"] == 1
        assert s["ref_bytes"] > 0

    def test_clear_references_recreates_empty(self, tmp_path):
        p = self._proj(tmp_path)
        assert clear_project_data(p, "references") is True
        assert project_data_summary(p)["n_ref_loci"] == 0
        assert (p / "references").exists()        # recreated empty, not just deleted

    def test_clear_results_and_runs(self, tmp_path):
        p = self._proj(tmp_path)
        clear_project_data(p, "results")
        clear_project_data(p, "runs")
        s = project_data_summary(p)
        assert s["n_combined"] == 0 and s["n_runs"] == 0

    def test_clear_rejects_bad_subdir(self, tmp_path):
        with pytest.raises(ValueError):
            clear_project_data(self._proj(tmp_path), "metadata")

    def test_clear_rejects_non_project(self, tmp_path):
        with pytest.raises(ValueError):
            clear_project_data(tmp_path / "nope", "references")

    def test_reset_workflow(self, tmp_path):
        p = self._proj(tmp_path)
        update_step(p, "coding", status="done")
        reset_workflow(p)
        assert get_workflow(p)["steps"]["coding"]["status"] == "pending"

    def test_delete_project(self, tmp_path):
        p = self._proj(tmp_path)
        assert delete_project(p) is True
        assert not p.exists()

    def test_delete_refuses_non_project(self, tmp_path):
        (tmp_path / "plain").mkdir()
        with pytest.raises(ValueError):
            delete_project(tmp_path / "plain")


# Two records exercising both schemas the app produces.
_REGISTRY = {
    # Assembly Manager schema: nested stats{} with QUAST
    "NS26-3-C2": {
        "strain_id": "NS26-3-C2",
        "assembly_path": "/data/NS26-3-C2/NS26-3-C2_final_EGAP_assembly.fasta",
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


class TestMigrateAssemblyRecord:
    def test_nested_record_preserved_with_taxon_defaults(self):
        # Canonical (nested-stats) records keep their content but now gain taxon defaults
        # (D-012); the original dict is not mutated.
        rec = {"assembly_path": "/d/x.fasta", "stats": {"n50": 100}, "reads_r1": ""}
        result = _migrate_assembly_record("X", rec)
        assert result["stats"] == {"n50": 100}
        assert result["assembly_path"] == "/d/x.fasta"
        assert result["reads_r1"] == ""
        assert result["taxon"] == "" and result["taxon_source"] == ""
        assert result["strain_id"] == "X"
        assert "taxon" not in rec          # caller's dict untouched

    def test_flat_record_gets_stats_nested(self):
        flat = {
            "assembly_path": "/d/x.fasta",
            "assembler": "spades",
            "n50": 900000,
            "num_contigs": 90,
            "total_length_mb": 33.1,
            "gc_percent": 50.8,
        }
        result = _migrate_assembly_record("X", flat)
        assert isinstance(result["stats"], dict)
        assert result["stats"]["n50"] == 900000
        assert result["stats"]["mean_gc"] == 50.8
        assert "gc_percent" not in result["stats"]
        assert result["assembly_path"] == "/d/x.fasta"

    def test_strain_id_backfilled(self):
        flat = {"assembly_path": "/d/x.fasta", "n50": 1}
        result = _migrate_assembly_record("MySample", flat)
        assert result["strain_id"] == "MySample"

    def test_roundtrip_migrates_on_load(self, tmp_path):
        """Save a flat-schema record; load_assembly_registry must normalize it."""
        flat_registry = {
            "OldSample": {
                "assembly_path": "/d/old.fasta",
                "assembler": "egap",
                "n50": 500000,
                "gc_percent": 51.0,
            }
        }
        # Save the raw (un-migrated) JSON directly so it's in the old format.
        import json
        meta = (tmp_path / "metadata")
        meta.mkdir(parents=True, exist_ok=True)
        (meta / "assemblies.json").write_text(json.dumps(flat_registry))

        loaded = load_assembly_registry(tmp_path)
        assert isinstance(loaded["OldSample"]["stats"], dict)
        assert loaded["OldSample"]["stats"]["mean_gc"] == 51.0
        assert "gc_percent" not in loaded["OldSample"]["stats"]


class TestProjectManifestSchema:
    def test_init_project_writes_v2_manifest(self, tmp_path):
        init_project(tmp_path / "p")
        m = load_json(tmp_path / "p" / "metadata" / "project_manifest.json")
        assert m["schema_version"] == PROJECT_MANIFEST_SCHEMA_VERSION
        assert m["default_taxon"] == ""
        assert set(m["workflow"]["steps"]) == set(WORKFLOW_STEPS)
        assert m["workflow"]["strategy"] is None
        assert m["workflow"]["loci"] == []
        assert all(s["status"] == "pending" for s in m["workflow"]["steps"].values())

    def test_load_fills_defaults_when_absent(self, tmp_path):
        m = load_project_manifest(tmp_path / "ghost")   # no project on disk yet
        assert m["default_taxon"] == ""
        assert set(m["workflow"]["steps"]) == set(WORKFLOW_STEPS)

    def test_load_upgrades_v1_manifest_without_data_loss(self, tmp_path):
        meta = tmp_path / "old" / "metadata"
        meta.mkdir(parents=True)
        save_json(meta / "project_manifest.json", {
            "name": "old", "created_at": "2026-01-01T00:00:00+00:00",
            "schema_version": 1, "notes": "legacy",
        })
        m = load_project_manifest(tmp_path / "old")
        assert m["name"] == "old"                                      # preserved
        assert m["notes"] == "legacy"                                  # preserved
        assert m["schema_version"] == PROJECT_MANIFEST_SCHEMA_VERSION  # upgraded
        assert set(m["workflow"]["steps"]) == set(WORKFLOW_STEPS)      # added

    def test_ensure_preserves_existing_values_and_extra_steps(self, tmp_path):
        meta = tmp_path / "x" / "metadata"
        meta.mkdir(parents=True)
        save_json(meta / "project_manifest.json", {
            "workflow": {"steps": {"references": {"status": "done"},
                                   "custom_step": {"status": "running"}}},
        })
        m = load_project_manifest(tmp_path / "x")
        assert m["workflow"]["steps"]["references"]["status"] == "done"  # value kept
        assert "custom_step" in m["workflow"]["steps"]                   # extra kept
        assert m["workflow"]["steps"]["coding"]["status"] == "pending"   # canonical added


class TestWorkflowStateHelpers:
    def test_set_default_taxon_strips(self, tmp_path):
        set_default_taxon(tmp_path / "p", "  Alternaria ")
        assert load_project_manifest(tmp_path / "p")["default_taxon"] == "Alternaria"

    def test_set_strategy_and_loci(self, tmp_path):
        set_workflow_strategy(tmp_path / "p", "fungal-barcodes")
        set_workflow_loci(tmp_path / "p", ["ITS", "RPB2"])
        wf = get_workflow(tmp_path / "p")
        assert wf["strategy"] == "fungal-barcodes"
        assert wf["loci"] == ["ITS", "RPB2"]

    def test_clearing_strategy(self, tmp_path):
        set_workflow_strategy(tmp_path / "p", "x")
        set_workflow_strategy(tmp_path / "p", None)
        assert get_workflow(tmp_path / "p")["strategy"] is None

    def test_update_step_merges_outputs_and_stamps(self, tmp_path):
        update_step(tmp_path / "p", "coding", status="running")
        update_step(tmp_path / "p", "coding", outputs={"cds": "RPB2_CDS.fasta"})
        update_step(tmp_path / "p", "coding", status="done",
                    outputs={"protein": "RPB2_protein.fasta"}, notes="2 loci")
        step = get_workflow(tmp_path / "p")["steps"]["coding"]
        assert step["status"] == "done"
        assert step["outputs"] == {"cds": "RPB2_CDS.fasta", "protein": "RPB2_protein.fasta"}
        assert step["notes"] == "2 loci"
        assert step["updated_at"]      # timestamp stamped

    def test_update_step_rejects_unknown_step(self, tmp_path):
        with pytest.raises(ValueError):
            update_step(tmp_path / "p", "nonsense", status="done")

    def test_update_step_rejects_bad_status(self, tmp_path):
        with pytest.raises(ValueError):
            update_step(tmp_path / "p", "coding", status="finished")


class TestAssemblyTaxonomy:
    def test_effective_taxon_prefers_override(self):
        assert effective_taxon({"taxon": "Alternaria alternata"}, "Fungi") == "Alternaria alternata"
        assert effective_taxon({"taxon": ""}, "Fungi") == "Fungi"
        assert effective_taxon({}, "Fungi") == "Fungi"
        assert effective_taxon({}, "") == ""

    def test_set_assembly_taxon_persists_and_logs_tsv(self, tmp_path):
        save_assembly_registry(tmp_path / "p", _REGISTRY)
        set_assembly_taxon(tmp_path / "p", "S9-1B-A2", "Alternaria alternata", source="manual")
        reg = load_assembly_registry(tmp_path / "p")
        assert reg["S9-1B-A2"]["taxon"] == "Alternaria alternata"
        assert reg["S9-1B-A2"]["taxon_source"] == "manual"
        tsv = (tmp_path / "p" / "metadata" / "assembly_manifest.tsv").read_text()
        assert "Alternaria alternata" in tsv

    def test_its_blast_source_allowed(self, tmp_path):
        save_assembly_registry(tmp_path / "p", _REGISTRY)
        set_assembly_taxon(tmp_path / "p", "S9-1B-A2", "Alternaria sp.", source="its_blast")
        assert load_assembly_registry(tmp_path / "p")["S9-1B-A2"]["taxon_source"] == "its_blast"

    def test_unknown_strain_raises(self, tmp_path):
        save_assembly_registry(tmp_path / "p", _REGISTRY)
        with pytest.raises(KeyError):
            set_assembly_taxon(tmp_path / "p", "NOPE", "X")

    def test_bad_source_raises(self, tmp_path):
        save_assembly_registry(tmp_path / "p", _REGISTRY)
        with pytest.raises(ValueError):
            set_assembly_taxon(tmp_path / "p", "S9-1B-A2", "X", source="guess")

    def test_migrated_records_gain_taxon_fields(self, tmp_path):
        save_assembly_registry(tmp_path / "p", _REGISTRY)
        for rec in load_assembly_registry(tmp_path / "p").values():
            assert rec.get("taxon") == ""
            assert rec.get("taxon_source") == ""


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
