#!/usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test_profiles.py
#
# Checks for the per-level (3D) global means: read_monthly_fields with
# profile_vars, and print_profiles2text. Builds a few tiny synthetic
# cam.h0 files (hybrid levels, PS, T) in a temporary directory.
#
# Run: `python3 -m pytest test_profiles.py`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
import netCDF4 as nc
import pytest

import trend_core as core
import trend_utils as trend

NLEV, NLAT, NLON = 3, 4, 6
HYAM = np.array([0.01, 0.05, 0.0])
HYBM = np.array([0.0, 0.4, 0.95])
P0 = 1.0e5


def _write_month(path, month, ts, t_levels, ps):
    with nc.Dataset(path, "w") as d:
        for name, n in (("time", None), ("lev", NLEV), ("lat", NLAT), ("lon", NLON)):
            d.createDimension(name, n)
        d.createVariable("lat", "f8", ("lat",))[:] = np.linspace(-67.5, 67.5, NLAT)
        d.createVariable("lon", "f8", ("lon",))[:] = np.arange(NLON) * 360.0 / NLON
        d.createVariable("lev", "f8", ("lev",))[:] = np.arange(NLEV)
        d.createVariable("hyam", "f8", ("lev",))[:] = HYAM
        d.createVariable("hybm", "f8", ("lev",))[:] = HYBM
        d.createVariable("P0", "f8", ())[...] = P0
        d.createVariable("TS", "f4", ("time", "lat", "lon"))[0] = np.full((NLAT, NLON), ts)
        d.createVariable("PS", "f4", ("time", "lat", "lon"))[0] = ps
        T = np.broadcast_to(np.asarray(t_levels, dtype=float)[:, None, None],
                            (NLEV, NLAT, NLON))
        d.createVariable("T", "f4", ("time", "lev", "lat", "lon"))[0] = T


@pytest.fixture
def case_dir(tmp_path):
    for m in range(1, 4):
        ps = np.full((NLAT, NLON), 1.0e5 + 100.0 * m)
        ps[0, :] += 400.0  # a non-uniform PS, so <PS> must be area-weighted
        _write_month(tmp_path / f"c.cam.h0.0001-{m:02d}.nc", m, 300.0 + m,
                     [200.0 + m, 250.0 + 2 * m, 290.0 + 3 * m], ps)
    return tmp_path


def _weights(case_dir):
    with nc.Dataset(case_dir / "c.cam.h0.0001-01.nc") as d:
        return core.build_area_weights(d["lon"][:], d["lat"][:])


def test_profiles_per_level_and_pmid(case_dir):
    w = _weights(case_dir)
    out, prof, files = core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", ["TS"],
                                                1, 3, w, ["T"])
    assert len(files) == 3
    np.testing.assert_allclose(out[:, 0], [301, 302, 303], rtol=1e-6)
    assert set(prof) == {"T", "PMID"}
    assert prof["T"].shape == (3, NLEV)
    for m in range(1, 4):
        np.testing.assert_allclose(prof["T"][m - 1], [200 + m, 250 + 2 * m, 290 + 3 * m],
                                   rtol=1e-6)
        ps = np.full((NLAT, NLON), 1.0e5 + 100.0 * m)
        ps[0, :] += 400.0
        ps_mean = float((ps * w).sum())
        np.testing.assert_allclose(prof["PMID"][m - 1], HYAM * P0 + HYBM * ps_mean,
                                   rtol=1e-6)


def test_no_profiles_is_unchanged(case_dir):
    w = _weights(case_dir)
    out, files = core.read_monthly_files(str(case_dir), "c", ".cam.h0.", ["TS"], 1, 3, w)
    assert out.shape == (3, 1) and len(files) == 3
    _, prof, _ = core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", ["TS"], 1, 3, w)
    assert prof == {}


def test_profile_var_must_be_3d(case_dir):
    with pytest.raises(ValueError, match="not 3D"):
        core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", [], 1, 3,
                                 _weights(case_dir), ["TS"])


def test_missing_profile_var(case_dir):
    with pytest.raises(KeyError, match="Q"):
        core.read_monthly_fields(str(case_dir), "c", ".cam.h0.", [], 1, 3,
                                 _weights(case_dir), ["Q"])


def test_print_profiles2text(tmp_path):
    prof = {"T": np.array([[200.0, 250.5], [201.0, 251.25]]),
            "PMID": np.array([[1000.0, 9.5e4], [1000.0, 9.6e4]])}
    written = trend.print_profiles2text(prof, "0001-01", "0001-02", "c", outdir=str(tmp_path))
    assert [p.split("/")[-1] for p in written] == ["c_0001-01-0001-02_camlev_T.txt",
                                                   "c_0001-01-0001-02_camlev_PMID.txt"]
    lines = open(written[0]).read().splitlines()
    assert lines[0].split() == ["month", "L01", "L02"]
    assert lines[1].split() == ["1", "200", "250.5"]
    assert lines[2].split() == ["2", "201", "251.25"]
