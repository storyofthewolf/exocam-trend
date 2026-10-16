#!/usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test_cache.py
#
# Checks for the per-file global-mean cache in read_monthly_fields:
# cached results match a fresh read, fully cached files are never
# opened, new variables and rewritten files are recomputed, and a
# partial last line from a killed run is tolerated. Reuses the
# synthetic cam.h0 files from test_profiles.py.
#
# Run: `python3 -m pytest test_cache.py`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import json
import os
from types import SimpleNamespace

import numpy as np
import pytest

import trend_core as core
from test_profiles import case_dir, _weights, _write_month, NLAT, NLON  # noqa: F401  (fixture)


def _read(case_dir, cache, varnames=("TS",), profile_vars=()):
    return core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", list(varnames),
                                    1, 3, _weights(case_dir), list(profile_vars),
                                    cache_path=cache)


@pytest.fixture
def opened(monkeypatch):
    """Records every netCDF file the reader opens."""
    seen = []
    real = core.nc.Dataset

    def spy(path, *a, **k):
        seen.append(os.path.basename(str(path)))
        return real(path, *a, **k)

    # Patch trend_core's own reference only: the tests' _weights() opens files too.
    monkeypatch.setattr(core, "nc", SimpleNamespace(Dataset=spy))
    return seen


def test_cached_rerun_matches_and_opens_nothing(case_dir, tmp_path, opened):
    cache = str(tmp_path / "cache" / "c.cam.h0.jsonl")
    out1, prof1, _ = _read(case_dir, cache, profile_vars=["T"])
    assert len(opened) == 3
    opened.clear()
    out2, prof2, _ = _read(case_dir, cache, profile_vars=["T"])
    assert opened == []
    np.testing.assert_array_equal(out1, out2)
    for k in prof1:
        np.testing.assert_array_equal(prof1[k], prof2[k])
    out0, prof0, _ = _read(case_dir, None, profile_vars=["T"])
    np.testing.assert_array_equal(out0, out2)
    np.testing.assert_array_equal(prof0["PMID"], prof2["PMID"])


def test_new_variable_appends_only_that_variable(case_dir, tmp_path, opened):
    cache = str(tmp_path / "c.jsonl")
    _read(case_dir, cache, varnames=["TS"])
    opened.clear()
    out, _, _ = _read(case_dir, cache, varnames=["TS", "PS"])
    assert len(opened) == 3
    last = json.loads(open(cache).read().splitlines()[-1])
    assert set(last["vars"]) == {"PS"}
    np.testing.assert_allclose(out[:, 0], [301, 302, 303], rtol=1e-6)
    opened.clear()
    _read(case_dir, cache, varnames=["PS", "TS"])
    assert opened == []


def test_rewritten_file_is_read_again(case_dir, tmp_path, opened):
    cache = str(tmp_path / "c.jsonl")
    _read(case_dir, cache)
    f = case_dir / "c.cam.h0.0001-02.nc"
    st = os.stat(f)
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 1_000_000_000))
    opened.clear()
    _read(case_dir, cache)
    assert opened == ["c.cam.h0.0001-02.nc"]


def test_other_weights_ignore_the_cache(case_dir, tmp_path, opened):
    cache = str(tmp_path / "c.jsonl")
    _read(case_dir, cache)
    opened.clear()
    w = _weights(case_dir)
    w2 = np.roll(w, 1, axis=0)
    core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", ["TS"], 1, 3, w2,
                             cache_path=cache)
    assert len(opened) == 3


def test_missing_variable_is_cached_as_nan(case_dir, tmp_path, opened):
    cache = str(tmp_path / "c.jsonl")
    out1, _, _ = _read(case_dir, cache, varnames=["TS", "NOPE"])
    opened.clear()
    out2, _, _ = _read(case_dir, cache, varnames=["TS", "NOPE"])
    assert opened == []
    assert np.isnan(out1[:, 1]).all() and np.isnan(out2[:, 1]).all()


def test_partial_last_line_from_killed_run(case_dir, tmp_path, opened):
    cache = tmp_path / "c.jsonl"
    _read(case_dir, str(cache))
    lines = cache.read_text().splitlines()
    # a run killed while writing month 3: its record is cut off mid-line
    cache.write_text("\n".join(lines[:2]) + "\n" + lines[2][:20])
    opened.clear()
    out, _, _ = _read(case_dir, str(cache))
    assert opened == ["c.cam.h0.0001-03.nc"]
    np.testing.assert_allclose(out[:, 0], [301, 302, 303], rtol=1e-6)
    opened.clear()
    _read(case_dir, str(cache))
    assert opened == []


def test_extended_series_reads_only_new_months(case_dir, tmp_path, opened):
    # A finished run over 3 months; days later the GCM has written 2 more.
    cache = str(tmp_path / "c.jsonl")
    _read(case_dir, cache)
    for m in (4, 5):
        _write_month(case_dir / f"c.cam.h0.0001-{m:02d}.nc", m, 300.0 + m,
                     [200.0, 250.0, 290.0], np.full((NLAT, NLON), 1.0e5))
    opened.clear()
    out, _, files = core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", ["TS"],
                                             1, 6000, _weights(case_dir),
                                             cache_path=cache)
    assert opened == ["c.cam.h0.0001-04.nc", "c.cam.h0.0001-05.nc"]
    assert len(files) == 5
    np.testing.assert_allclose(out[:, 0], [301, 302, 303, 304, 305], rtol=1e-6)
