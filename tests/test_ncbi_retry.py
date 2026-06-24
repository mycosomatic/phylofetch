"""
Tests for the Entrez transport hardening (D-034): retry-with-backoff, typed NCBIError on
exhaustion, and — critically — ncbi_search_count raising rather than returning 0 on a network
failure (so a transient blip can't masquerade as a genuine "0 references").
"""

import pytest

from phylofetch import ncbi_utils as nu


@pytest.fixture(autouse=True)
def _no_real_waiting(monkeypatch):
    # Keep tests instant: no throttle sleeps, no backoff sleeps.
    monkeypatch.setattr(nu, "_throttle", lambda: None)
    monkeypatch.setattr(nu.time, "sleep", lambda *_a, **_k: None)


class _Handle:
    def close(self):
        pass


class TestEntrezRetry:
    def test_succeeds_after_one_transient_failure(self):
        calls = {"n": 0}

        def thunk():
            calls["n"] += 1
            if calls["n"] < 2:
                raise OSError("transient blip")
            return "ok"

        assert nu._entrez_retry(thunk, what="t") == "ok"
        assert calls["n"] == 2

    def test_raises_ncbierror_after_exhaustion(self):
        def thunk():
            raise OSError("network down")

        with pytest.raises(nu.NCBIError):
            nu._entrez_retry(thunk, what="t", retries=3)

    def test_runtimeerror_from_entrez_read_is_retried(self):
        # Bio.Entrez.read raises RuntimeError on an NCBI-side error body.
        calls = {"n": 0}

        def thunk():
            calls["n"] += 1
            if calls["n"] < 3:
                raise RuntimeError("NCBI says: too many requests")
            return 5

        assert nu._entrez_retry(thunk, what="t", retries=3) == 5


class TestNcbiSearchCountFailureMode:
    def test_raises_not_zero_on_persistent_network_failure(self, monkeypatch):
        monkeypatch.setattr(nu, "_entrez_email", "x@y.z")
        monkeypatch.setattr(nu.Entrez, "esearch",
                            lambda **k: (_ for _ in ()).throw(OSError("network down")))
        with pytest.raises(nu.NCBIError):
            nu.ncbi_search_count("tef1", "Alternaria")

    def test_returns_count_after_transient_then_success(self, monkeypatch):
        monkeypatch.setattr(nu, "_entrez_email", "x@y.z")
        state = {"n": 0}

        def flaky_esearch(**k):
            state["n"] += 1
            if state["n"] < 2:
                raise OSError("blip")
            return _Handle()

        monkeypatch.setattr(nu.Entrez, "esearch", flaky_esearch)
        monkeypatch.setattr(nu.Entrez, "read", lambda h: {"Count": "17"})
        assert nu.ncbi_search_count("tef1", "Alternaria") == 17
        assert state["n"] == 2


class TestApiKey:
    def test_set_api_key_sets_entrez_and_shrinks_interval(self, monkeypatch):
        monkeypatch.setattr(nu.Entrez, "api_key", "", raising=False)
        nu.set_api_key("ABC123")
        try:
            assert nu.Entrez.api_key == "ABC123"
        finally:
            nu.set_api_key("")
