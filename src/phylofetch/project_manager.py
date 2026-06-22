"""
project_manager.py
------------------
Reproducible run infrastructure: per-command directories, full logging, environment
and tool-version snapshots, append-only command history.

All subprocess calls in phylofetch go through RunManager so every analysis step
is traceable, re-runnable from the terminal, and accompanied by a log of exact
software versions used.
"""

from __future__ import annotations

import csv
import json
import os
import platform
import re
import shlex
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


DEFAULT_PROJECTS_ROOT = Path.home() / ".phylofetch" / "projects"
DEFAULT_PROJECT_DIR = DEFAULT_PROJECTS_ROOT / "default"

PROJECT_SUBDIRS = ("metadata", "runs", "results", "scratch", "logs")

# Project-manifest schema (metadata/project_manifest.json). Bumped to 2 for RM-007/D-012 to
# add a project-level default taxon and a workflow/step-state block. Older (v1) manifests are
# read tolerantly — load_project_manifest() fills the new keys in memory and the next save
# upgrades the file on disk. See DECISIONS.md D-012.
PROJECT_MANIFEST_SCHEMA_VERSION = 2

# Ordered extraction steps the Workflow page chains, and the status vocabulary each may hold.
WORKFLOW_STEPS = ("references", "rDNA", "coding", "primers", "combine")
STEP_STATUSES = ("pending", "running", "done", "error", "skipped")

# How an assembly's taxon was assigned (D-012): typed by the user, or picked from the
# ITS-extraction BLAST result. Empty string = unset.
TAXON_SOURCES = ("manual", "its_blast")

# Columns for the human-readable assembly manifest (metadata/assembly_manifest.tsv).
ASSEMBLY_MANIFEST_FIELDS = [
    "strain_id", "taxon", "taxon_source", "assembly_path", "assembler",
    "num_contigs", "n50", "total_length_mb", "gc_percent", "quast_report",
    "busco_dir", "busco_completeness", "registered_at",
]


@dataclass
class RunResult:
    """Structured return value from a logged external command."""

    run_id: str
    run_dir: str
    command: str
    returncode: int
    stdout_path: str
    stderr_path: str
    terminal_log_path: str
    command_json_path: str
    started_at: str
    finished_at: str


@dataclass
class ToolStatus:
    label: str
    executable: str
    available: bool
    resolved_path: str | None
    version: str | None = None


def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def timestamp_for_path() -> str:
    return datetime.now().strftime("%Y%m%d_%H%M%S")


def safe_slug(text: str, fallback: str = "item") -> str:
    cleaned = []
    for ch in text.strip():
        if ch.isalnum() or ch in ("-", "_"):
            cleaned.append(ch)
        elif ch in (" ", ".", "/", "\\", ":"):
            cleaned.append("_")
    slug = "".join(cleaned).strip("_")
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug or fallback


def _default_workflow() -> dict:
    """Fresh workflow/step-state block for a project manifest (D-012)."""
    return {
        "strategy": None,
        "loci": [],
        "steps": {
            step: {"status": "pending", "updated_at": "", "outputs": {}, "notes": ""}
            for step in WORKFLOW_STEPS
        },
    }


def _ensure_manifest_defaults(manifest: Mapping[str, Any] | None) -> dict:
    """
    Return a project-manifest dict with every current-schema key present, without
    discarding existing values. Tolerant of v1 manifests (no default_taxon / workflow):
    missing keys are filled in memory, and the canonical workflow steps are guaranteed to
    exist while any extra/unknown steps already present are preserved.
    """
    m = dict(manifest or {})
    m["schema_version"] = PROJECT_MANIFEST_SCHEMA_VERSION
    m.setdefault("default_taxon", "")
    m.setdefault("output_dir", "")        # "" → default <project>/results (see project_output_dir)
    wf = m.get("workflow")
    if not isinstance(wf, dict):
        wf = _default_workflow()
    else:
        wf.setdefault("strategy", None)
        wf.setdefault("loci", [])
        steps = wf.get("steps")
        steps = dict(steps) if isinstance(steps, dict) else {}
        for step in WORKFLOW_STEPS:
            s = steps.get(step)
            s = dict(s) if isinstance(s, dict) else {}
            s.setdefault("status", "pending")
            s.setdefault("updated_at", "")
            s.setdefault("outputs", {})
            s.setdefault("notes", "")
            steps[step] = s
        wf["steps"] = steps
    m["workflow"] = wf
    return m


def init_project(project_dir: str | Path) -> Path:
    root = Path(project_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    for subdir in PROJECT_SUBDIRS:
        (root / subdir).mkdir(parents=True, exist_ok=True)

    manifest = root / "metadata" / "project_manifest.json"
    if not manifest.exists():
        save_json(manifest, _ensure_manifest_defaults({
            "name": root.name,
            "created_at": now_iso(),
            "project_dir": str(root),
            "notes": "Created by phylofetch.",
        }))

    history = root / "metadata" / "command_history.tsv"
    if not history.exists():
        write_tsv_rows(history, [], fieldnames=[
            "run_id", "started_at", "finished_at", "module",
            "action", "returncode", "workdir", "command",
        ])

    return root


def load_json(path: str | Path, default: Any = None) -> Any:
    p = Path(path)
    if not p.exists():
        return default
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return default


def save_json(path: str | Path, data: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True, default=str)
        f.write("\n")


def append_jsonl(path: str | Path, item: Mapping[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(dict(item), sort_keys=True, default=str) + "\n")


def write_tsv_rows(path: str | Path, rows: Iterable[Mapping[str, Any]],
                   fieldnames: Sequence[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def append_tsv_row(path: str | Path, row: Mapping[str, Any],
                   fieldnames: Sequence[str]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    exists = p.exists()
    with open(p, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        if not exists:
            w.writeheader()
        w.writerow(row)


def command_to_display(command: str | Sequence[str]) -> str:
    if isinstance(command, str):
        return command
    return shlex.join([str(x) for x in command])


def _probe_version(executable: str) -> str:
    """Try common version flags; return first non-empty line of output."""
    for flag in ("--version", "-version", "version"):
        try:
            r = subprocess.run(
                [executable, flag],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, timeout=5,
            )
            first = (r.stdout or "").strip().splitlines()
            if first:
                return first[0][:120]
        except Exception:
            pass
    return "unknown"


def capture_tool_versions(tool_map: Mapping[str, str]) -> dict[str, str]:
    """Return {label: version_string} for each tool in tool_map."""
    return {label: _probe_version(exe) for label, exe in tool_map.items() if exe}


def snapshot_environment(tool_versions: dict[str, str] | None = None) -> dict[str, Any]:
    snap = {
        "captured_at": now_iso(),
        "platform": platform.platform(),
        "python_version": sys.version,
        "executable": sys.executable,
        "cwd": os.getcwd(),
        "conda_default_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "virtual_env": os.environ.get("VIRTUAL_ENV", ""),
        "path": os.environ.get("PATH", ""),
    }
    if tool_versions:
        snap["tool_versions"] = tool_versions
    return snap


def check_tools(tool_map: Mapping[str, str]) -> list[ToolStatus]:
    statuses = []
    for label, executable in tool_map.items():
        resolved = shutil.which(executable) if executable else None
        version = _probe_version(executable) if resolved else None
        statuses.append(ToolStatus(
            label=label, executable=executable,
            available=resolved is not None,
            resolved_path=resolved, version=version,
        ))
    return statuses


# ── Assembly registry persistence ───────────────────────────────────────────
#
# Each project stores its own re-openable assembly registry:
#   <project>/metadata/assemblies.json        machine-readable (full records)
#   <project>/metadata/assembly_manifest.tsv  human-readable log of file paths
#
# This is the durable record of *where every assembly and associated file lives*
# — it survives across sessions and can be re-opened as a named project.

def _assembly_record_to_row(strain_id: str, rec: Mapping[str, Any]) -> dict[str, Any]:
    """
    Flatten one assembly record into a manifest row. Tolerant of both the
    flat schema (Project Setup import) and the nested ``stats`` schema
    (Assembly Manager import).
    """
    stats = rec.get("stats") if isinstance(rec.get("stats"), Mapping) else {}
    stats = stats or {}

    def pick(*keys: str, default: Any = "") -> Any:
        for k in keys:
            if rec.get(k) not in (None, ""):
                return rec[k]
            if stats.get(k) not in (None, ""):
                return stats[k]
        return default

    quast = stats.get("quast") if isinstance(stats.get("quast"), Mapping) else None
    quast_report = (rec.get("quast_report") or stats.get("quast_report")
                    or ("(found)" if quast else ""))

    busco = rec.get("busco_dir") or rec.get("busco_dirs") or ""
    if isinstance(busco, (list, tuple)):
        busco = ";".join(str(b) for b in busco)

    busco_info = rec.get("busco")
    busco_pct = ""
    if isinstance(busco_info, Mapping):
        cp = busco_info.get("completeness_pct")
        busco_pct = cp if cp not in (None, "") else ""

    return {
        "strain_id":          strain_id,
        "taxon":              rec.get("taxon", ""),
        "taxon_source":       rec.get("taxon_source", ""),
        "assembly_path":      pick("assembly_path"),
        "assembler":          pick("assembler"),
        "num_contigs":        pick("num_contigs"),
        "n50":                pick("n50"),
        "total_length_mb":    pick("total_length_mb"),
        "gc_percent":         pick("gc_percent", "mean_gc"),
        "quast_report":       quast_report,
        "busco_dir":          busco,
        "busco_completeness": busco_pct,
        "registered_at":      rec.get("registered_at", ""),
    }


def save_assembly_registry(project_dir: str | Path,
                           assemblies: Mapping[str, Any]) -> Path:
    """
    Persist the assembly registry into a project as both a JSON record and a
    human-readable TSV manifest. Returns the project root.
    """
    root = init_project(project_dir)
    meta = root / "metadata"
    save_json(meta / "assemblies.json", dict(assemblies))
    rows = [_assembly_record_to_row(sid, rec) for sid, rec in assemblies.items()]
    write_tsv_rows(meta / "assembly_manifest.tsv", rows, ASSEMBLY_MANIFEST_FIELDS)
    return root


# Stats fields written by the old page-0 flat schema (before consolidation).
_FLAT_STAT_KEYS = frozenset({
    "assembler", "num_contigs", "total_length_bp", "total_length_mb",
    "n50", "l50", "largest_contig", "mean_gc", "gc_percent",
    "mean_coverage_header", "contigs", "quast", "quast_report",
})


def _migrate_assembly_record(sid: str, record: dict) -> dict:
    """
    Normalize legacy flat-format records to the nested-stats schema.

    Old page-0 schema stored stats at the top level; canonical schema wraps
    them under ``"stats"``. Records already in the canonical form pass through.
    Also renames ``gc_percent`` → ``mean_gc`` for consistency with
    ``assembly_utils.get_assembly_stats()`` output.
    """
    if isinstance(record.get("stats"), dict):
        rec = dict(record)
        rec.setdefault("strain_id", sid)
        rec.setdefault("taxon", "")
        rec.setdefault("taxon_source", "")
        return rec
    stats = {k: v for k, v in record.items() if k in _FLAT_STAT_KEYS}
    clean = {k: v for k, v in record.items() if k not in _FLAT_STAT_KEYS}
    if "gc_percent" in stats and "mean_gc" not in stats:
        stats["mean_gc"] = stats.pop("gc_percent")
    else:
        stats.pop("gc_percent", None)
    clean["stats"] = stats
    clean.setdefault("strain_id", sid)
    clean.setdefault("taxon", "")
    clean.setdefault("taxon_source", "")
    return clean


def load_assembly_registry(project_dir: str | Path) -> dict:
    """Load a project's saved assembly registry (empty dict if none)."""
    raw = load_json(Path(project_dir) / "metadata" / "assemblies.json", {})
    return {sid: _migrate_assembly_record(sid, rec) for sid, rec in raw.items()}


def list_projects(projects_root: str | Path | None = None) -> list[dict]:
    """
    Enumerate phylofetch projects under ``projects_root`` (default
    ~/.phylofetch/projects). Returns name, path, created_at, and assembly count.
    """
    root = Path(projects_root) if projects_root else DEFAULT_PROJECTS_ROOT
    projects: list[dict] = []
    if not root.exists():
        return projects
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        manifest = d / "metadata" / "project_manifest.json"
        if not manifest.exists():
            continue
        info = load_json(manifest, {})
        registry = load_json(d / "metadata" / "assemblies.json", {})
        projects.append({
            "name":         info.get("name", d.name),
            "path":         str(d),
            "created_at":   info.get("created_at", "?"),
            "n_assemblies": len(registry),
        })
    return projects


# ── Project manifest: taxonomy + workflow state (D-012) ──────────────────────
#
# The project manifest (metadata/project_manifest.json) is the chaining backbone for the
# component-page workflow: it carries the project-level default taxon and a workflow block
# (chosen strategy, selected loci, per-step status/outputs). Per-assembly taxon overrides
# live on the assembly records in the registry. Reads tolerate v1 manifests; writes upgrade
# the schema. See PLANNING.md RM-007.

def _project_manifest_path(project_dir: str | Path) -> Path:
    return Path(project_dir).expanduser().resolve() / "metadata" / "project_manifest.json"


def load_project_manifest(project_dir: str | Path) -> dict:
    """Load a project manifest with all current-schema keys filled (v1-tolerant)."""
    return _ensure_manifest_defaults(load_json(_project_manifest_path(project_dir), {}))


def save_project_manifest(project_dir: str | Path, manifest: Mapping[str, Any]) -> Path:
    """Persist a project manifest (schema kept current). Returns the manifest path."""
    init_project(project_dir)
    path = _project_manifest_path(project_dir)
    save_json(path, _ensure_manifest_defaults(manifest))
    return path


def set_default_taxon(project_dir: str | Path, taxon: str) -> dict:
    """Set the project-level default taxon; returns the updated manifest."""
    m = load_project_manifest(project_dir)
    m["default_taxon"] = (taxon or "").strip()
    save_project_manifest(project_dir, m)
    return m


def project_output_dir(project_dir: str | Path) -> Path:
    """
    Root for portable analysis artifacts (loci FASTAs, GFF3, partitions, alignments, …) that
    feed downstream tools. Defaults to ``<project>/results``; overridable via the manifest
    ``output_dir`` (e.g. a shared analysis folder). Pure path resolution — writers create their
    own subdirs. (D-021)
    """
    od = (load_project_manifest(project_dir).get("output_dir") or "").strip()
    return Path(od).expanduser() if od else (Path(project_dir) / "results")


def set_output_dir(project_dir: str | Path, path: str) -> dict:
    """Set (or clear with "") the project's output-directory override; returns the manifest."""
    m = load_project_manifest(project_dir)
    m["output_dir"] = (path or "").strip()
    save_project_manifest(project_dir, m)
    return m


def get_workflow(project_dir: str | Path) -> dict:
    """Return the workflow/step-state block of the project manifest."""
    return load_project_manifest(project_dir)["workflow"]


def set_workflow_strategy(project_dir: str | Path, strategy: str | None) -> dict:
    """Set the chosen strategy name (or None to clear); returns the updated manifest."""
    m = load_project_manifest(project_dir)
    m["workflow"]["strategy"] = strategy
    save_project_manifest(project_dir, m)
    return m


def set_workflow_loci(project_dir: str | Path, loci: Sequence[str]) -> dict:
    """Set the selected loci list; returns the updated manifest."""
    m = load_project_manifest(project_dir)
    m["workflow"]["loci"] = list(loci)
    save_project_manifest(project_dir, m)
    return m


def update_step(project_dir: str | Path, step: str, *,
                status: str | None = None,
                outputs: Mapping[str, Any] | None = None,
                notes: str | None = None) -> dict:
    """
    Update one workflow step's state. ``status`` (if given) must be in STEP_STATUSES;
    ``outputs`` is merged into the step's existing outputs; ``notes`` replaces. Always
    stamps ``updated_at``. Returns the updated manifest. Raises ValueError on an unknown
    step or invalid status.
    """
    if step not in WORKFLOW_STEPS:
        raise ValueError(f"unknown workflow step {step!r}; expected one of {WORKFLOW_STEPS}")
    if status is not None and status not in STEP_STATUSES:
        raise ValueError(f"invalid status {status!r}; expected one of {STEP_STATUSES}")
    m = load_project_manifest(project_dir)
    s = m["workflow"]["steps"][step]
    if status is not None:
        s["status"] = status
    if outputs:
        merged = dict(s.get("outputs") or {})
        merged.update(outputs)
        s["outputs"] = merged
    if notes is not None:
        s["notes"] = notes
    s["updated_at"] = now_iso()
    save_project_manifest(project_dir, m)
    return m


# ── Per-assembly taxonomy ────────────────────────────────────────────────────

def effective_taxon(record: Mapping[str, Any], default_taxon: str = "") -> str:
    """The assembly's own taxon if set, else the project default (both stripped)."""
    return (record.get("taxon") or "").strip() or (default_taxon or "").strip()


def set_assembly_taxon(project_dir: str | Path, strain_id: str, taxon: str,
                       source: str = "manual") -> dict:
    """
    Set ``taxon`` + ``taxon_source`` on a registered assembly and persist the registry.
    Returns the updated registry. Raises KeyError if the strain isn't registered and
    ValueError on an unrecognised source.
    """
    if source not in TAXON_SOURCES:
        raise ValueError(f"invalid taxon source {source!r}; expected one of {TAXON_SOURCES}")
    registry = load_assembly_registry(project_dir)
    if strain_id not in registry:
        raise KeyError(f"assembly {strain_id!r} is not registered in this project")
    registry[strain_id]["taxon"] = (taxon or "").strip()
    registry[strain_id]["taxon_source"] = source
    save_assembly_registry(project_dir, registry)
    return registry


# ── Project data management (inspect / clear caches, delete project) ──────────

_CLEARABLE_SUBDIRS = ("references", "results", "runs", "scratch", "logs")


def _is_phylofetch_project(project_dir: str | Path) -> bool:
    return (Path(project_dir) / "metadata" / "project_manifest.json").exists()


def project_data_summary(project_dir: str | Path) -> dict:
    """Counts + on-disk byte sizes of a project's caches, for the data-management UI."""
    root = Path(project_dir)

    def _bytes(p: Path) -> int:
        return sum(f.stat().st_size for f in p.rglob("*") if f.is_file()) if p.exists() else 0

    refs = root / "references"
    runs = root / "runs"
    out = project_output_dir(project_dir)          # honours the output_dir override
    combined = out / "loci" / "combined"
    n_ref_loci = (len([d for d in refs.iterdir()
                       if d.is_dir() and (d / f"{d.name}_refs.fasta").exists()])
                  if refs.exists() else 0)
    return {
        "n_assemblies": len(load_assembly_registry(project_dir)),
        "n_ref_loci": n_ref_loci,
        "ref_bytes": _bytes(refs),
        "n_combined": len(list(combined.glob("*.fasta"))) if combined.exists() else 0,
        "results_bytes": _bytes(out),
        "n_runs": len([d for d in runs.iterdir() if d.is_dir()]) if runs.exists() else 0,
        "runs_bytes": _bytes(runs),
        "output_dir": str(out),
    }


def clear_project_data(project_dir: str | Path, subdir: str) -> bool:
    """
    Remove and recreate one project cache subdir (references / results / runs / scratch / logs).
    Returns True if it existed. Refuses unknown subdirs and non-project paths.
    """
    if subdir not in _CLEARABLE_SUBDIRS:
        raise ValueError(f"not a clearable subdir: {subdir!r}")
    root = Path(project_dir)
    if not _is_phylofetch_project(root):
        raise ValueError("not a phylofetch project (no metadata/project_manifest.json)")
    target = root / subdir
    existed = target.exists()
    if existed:
        shutil.rmtree(target)
    target.mkdir(parents=True, exist_ok=True)
    return existed


def reset_workflow(project_dir: str | Path) -> dict:
    """Reset the manifest's workflow/step state to fresh (keeps identity + default_taxon)."""
    m = load_project_manifest(project_dir)
    m["workflow"] = _default_workflow()
    save_project_manifest(project_dir, m)
    return m


def delete_project(project_dir: str | Path) -> bool:
    """
    Delete an entire project directory. Guarded: must be a phylofetch project and not a
    protected directory (home, the projects root, or '/'). Returns True on deletion.
    """
    root = Path(project_dir).expanduser().resolve()
    if not _is_phylofetch_project(root):
        raise ValueError("refusing to delete: not a phylofetch project")
    protected = {Path.home().resolve(), DEFAULT_PROJECTS_ROOT.expanduser().resolve(), Path("/")}
    if root in protected:
        raise ValueError(f"refusing to delete protected directory: {root}")
    shutil.rmtree(root)
    return True


def clear_global_reference_cache() -> bool:
    """Remove the shared global reference library (~/.phylofetch/references). True if it existed."""
    cache = Path.home() / ".phylofetch" / "references"
    existed = cache.exists()
    if existed:
        shutil.rmtree(cache)
    return existed


class RunManager:
    """Execute and log external commands with full reproducibility artifacts."""

    _history_fields = [
        "run_id", "started_at", "finished_at", "module",
        "action", "returncode", "workdir", "command",
    ]

    def __init__(self, project_dir: str | Path):
        self.project_dir = init_project(project_dir)
        self.runs_dir = self.project_dir / "runs"
        self.command_history = self.project_dir / "metadata" / "command_history.tsv"

    def make_run_dir(self, module: str, action: str) -> tuple[str, Path]:
        run_id = (
            f"{timestamp_for_path()}_{safe_slug(module)}_{safe_slug(action)}"
            f"_{uuid.uuid4().hex[:8]}"
        )
        run_dir = self.runs_dir / run_id
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_id, run_dir

    def run(
        self,
        command: str | Sequence[str],
        *,
        module: str,
        action: str,
        workdir: str | Path | None = None,
        inputs: Mapping[str, Any] | None = None,
        outputs: Mapping[str, Any] | None = None,
        params: Mapping[str, Any] | None = None,
        tool_version_keys: Sequence[str] | None = None,
        shell: bool | None = None,
        timeout: int | None = None,
    ) -> RunResult:
        """
        Execute a command and write a reproducible run folder.

        tool_version_keys: executable names whose versions should be probed and
        stored in environment.json (e.g. ["blastn", "ITSx"]).
        """
        run_id, run_dir = self.make_run_dir(module, action)
        started_at = now_iso()
        display_cmd = command_to_display(command)
        cwd = str(Path(workdir).expanduser().resolve()) if workdir else None

        tool_versions = {}
        if tool_version_keys:
            tool_versions = {k: _probe_version(k) for k in tool_version_keys}

        command_json: dict[str, Any] = {
            "run_id": run_id,
            "module": module,
            "action": action,
            "started_at": started_at,
            "workdir": cwd,
            "command": display_cmd,
            "shell": bool(shell) if shell is not None else isinstance(command, str),
            "inputs": dict(inputs or {}),
            "outputs": dict(outputs or {}),
            "params": dict(params or {}),
        }
        save_json(run_dir / "command.json", command_json)
        save_json(run_dir / "environment.json", snapshot_environment(tool_versions))

        effective_shell = bool(shell) if shell is not None else isinstance(command, str)
        run_command: str | Sequence[str]
        if isinstance(command, str) and not effective_shell:
            run_command = shlex.split(command)
        else:
            run_command = command

        try:
            completed = subprocess.run(
                run_command,
                shell=effective_shell,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
                check=False,
            )
            stdout = completed.stdout or ""
            stderr = completed.stderr or ""
            returncode = int(completed.returncode)
        except Exception as exc:
            stdout = ""
            stderr = f"{type(exc).__name__}: {exc}\n"
            returncode = 999

        finished_at = now_iso()
        stdout_path = run_dir / "stdout.log"
        stderr_path = run_dir / "stderr.log"
        terminal_path = run_dir / "terminal.log"

        stdout_path.write_text(stdout, encoding="utf-8")
        stderr_path.write_text(stderr, encoding="utf-8")
        terminal_path.write_text(
            f"# Command\n{display_cmd}\n\n# STDOUT\n{stdout}\n\n# STDERR\n{stderr}",
            encoding="utf-8",
        )
        save_json(run_dir / "run_status.json", {
            "run_id": run_id, "started_at": started_at, "finished_at": finished_at,
            "returncode": returncode,
            "stdout_log": str(stdout_path), "stderr_log": str(stderr_path),
            "terminal_log": str(terminal_path),
        })

        append_tsv_row(self.command_history, {
            "run_id": run_id, "started_at": started_at, "finished_at": finished_at,
            "module": module, "action": action, "returncode": returncode,
            "workdir": cwd or "", "command": display_cmd,
        }, self._history_fields)
        append_jsonl(self.project_dir / "metadata" / "command_history.jsonl",
                     {**command_json, "finished_at": finished_at, "returncode": returncode})

        if returncode != 0:
            append_tsv_row(run_dir / "problems.tsv", {
                "run_id": run_id, "severity": "error",
                "problem": "Command returned nonzero exit status",
                "detail": f"returncode={returncode}", "created_at": finished_at,
            }, ["run_id", "severity", "problem", "detail", "created_at"])

        return RunResult(
            run_id=run_id, run_dir=str(run_dir), command=display_cmd,
            returncode=returncode,
            stdout_path=str(stdout_path), stderr_path=str(stderr_path),
            terminal_log_path=str(terminal_path),
            command_json_path=str(run_dir / "command.json"),
            started_at=started_at, finished_at=finished_at,
        )

    def dry_run(
        self,
        command: str | Sequence[str],
        *,
        module: str,
        action: str,
        workdir: str | Path | None = None,
        inputs: Mapping[str, Any] | None = None,
        outputs: Mapping[str, Any] | None = None,
        params: Mapping[str, Any] | None = None,
    ) -> tuple[str, Path]:
        run_id, run_dir = self.make_run_dir(module, action)
        save_json(run_dir / "command.json", {
            "run_id": run_id, "module": module, "action": action,
            "created_at": now_iso(), "workdir": str(workdir or ""),
            "command": command_to_display(command),
            "inputs": dict(inputs or {}), "outputs": dict(outputs or {}),
            "params": dict(params or {}), "dry_run": True,
        })
        return run_id, run_dir
