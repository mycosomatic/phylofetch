"""
Tests for the BLAST runner's RunManager provenance routing and timeout/launch guards
(D-033): run_blast must log through a manager when given, and must never hang or raise on
a timeout / missing binary — it surfaces a non-zero rc instead.
"""

import subprocess

from phylofetch import blast_loci_utils as bl


class _FakeRR:
    def __init__(self, returncode, stderr_path=""):
        self.returncode = returncode
        self.stdout_path = ""
        self.stderr_path = stderr_path


class _FakeManager:
    def __init__(self):
        self.calls = []

    def run(self, cmd, *, module, action, tool_version_keys=None, inputs=None,
            outputs=None, params=None, timeout=None, shell=None, workdir=None):
        self.calls.append({
            "cmd": list(cmd), "module": module, "action": action,
            "version_keys": tool_version_keys, "timeout": timeout, "inputs": inputs,
        })
        return _FakeRR(0)


def _fastas(tmp_path):
    q = tmp_path / "q.fa"; q.write_text(">q\nACGTACGT\n")
    s = tmp_path / "s.fa"; s.write_text(">s\nACGTACGT\n")
    return str(q), str(s)


def test_run_blast_routes_through_manager_with_version_key(tmp_path):
    mgr = _FakeManager()
    q, s = _fastas(tmp_path)
    rc, err, tsv = bl.run_blast(q, s, str(tmp_path), blast_bin="tblastn",
                                manager=mgr, module="blast", action="narrow_locus:TEF1")
    assert rc == 0
    assert len(mgr.calls) == 1
    call = mgr.calls[0]
    assert call["module"] == "blast" and call["action"] == "narrow_locus:TEF1"
    assert call["version_keys"] == ["tblastn"]   # tblastn binary → tblastn version key
    assert call["timeout"] == 600
    assert tsv.endswith("blast_hsps.tsv")


def test_run_blast_blastn_version_key(tmp_path):
    mgr = _FakeManager()
    q, s = _fastas(tmp_path)
    bl.run_blast(q, s, str(tmp_path), blast_bin="blastn", manager=mgr)
    assert mgr.calls[0]["version_keys"] == ["blastn"]


def test_run_blast_timeout_returns_124_not_raises(tmp_path, monkeypatch):
    def _boom(*a, **k):
        raise subprocess.TimeoutExpired(cmd="blastn", timeout=1)
    monkeypatch.setattr(bl.subprocess, "run", _boom)
    q, s = _fastas(tmp_path)
    rc, err, tsv = bl.run_blast(q, s, str(tmp_path), timeout=1)
    assert rc == 124
    assert "timed out" in err
    assert tsv.endswith("blast_hsps.tsv")


def test_run_blast_missing_binary_returns_127_not_raises(tmp_path, monkeypatch):
    def _boom(*a, **k):
        raise FileNotFoundError("blastn not found")
    monkeypatch.setattr(bl.subprocess, "run", _boom)
    q, s = _fastas(tmp_path)
    rc, err, _ = bl.run_blast(q, s, str(tmp_path))
    assert rc == 127
    assert "failed to launch" in err


def test_run_blast_no_manager_passes_timeout_to_subprocess(tmp_path, monkeypatch):
    seen = {}
    def _fake_run(cmd, **kw):
        seen.update(kw)
        class P:  # noqa: D401
            returncode = 0
            stderr = ""
        return P()
    monkeypatch.setattr(bl.subprocess, "run", _fake_run)
    q, s = _fastas(tmp_path)
    bl.run_blast(q, s, str(tmp_path), timeout=42)
    assert seen.get("timeout") == 42
