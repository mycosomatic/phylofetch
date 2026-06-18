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

# Columns for the human-readable assembly manifest (metadata/assembly_manifest.tsv).
ASSEMBLY_MANIFEST_FIELDS = [
    "strain_id", "assembly_path", "assembler", "num_contigs", "n50",
    "total_length_mb", "gc_percent", "quast_report", "busco_dir",
    "busco_completeness", "registered_at",
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


def init_project(project_dir: str | Path) -> Path:
    root = Path(project_dir).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    for subdir in PROJECT_SUBDIRS:
        (root / subdir).mkdir(parents=True, exist_ok=True)

    manifest = root / "metadata" / "project_manifest.json"
    if not manifest.exists():
        save_json(manifest, {
            "name": root.name,
            "created_at": now_iso(),
            "project_dir": str(root),
            "schema_version": 1,
            "notes": "Created by phylofetch.",
        })

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
        return record
    stats = {k: v for k, v in record.items() if k in _FLAT_STAT_KEYS}
    clean = {k: v for k, v in record.items() if k not in _FLAT_STAT_KEYS}
    if "gc_percent" in stats and "mean_gc" not in stats:
        stats["mean_gc"] = stats.pop("gc_percent")
    else:
        stats.pop("gc_percent", None)
    clean["stats"] = stats
    clean.setdefault("strain_id", sid)
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
