#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trend_core.py
#
#  Author: Wolf, E.T.
#
#  Core computation functions: area weights, netCDF4 file I/O,
#  global means, and running statistics.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import netCDF4 as nc
import glob
import hashlib
import json
import re
import os
from tqdm import tqdm

# Module-level set to suppress duplicate missing-variable warnings
_warned_missing = set()

# Bump when the per-file reduction changes (masking, weighting, PMID formula),
# so every cached value computed the old way is ignored and recomputed.
CACHE_VERSION = 1


def build_area_weights(lon, lat):
    """
    Compute normalized 2D area weights for a regular lat/lon grid.
    Returns array of shape (nlat, nlon) with weights summing to 1.0.
    Uses staggered grid: cell edges at lat midpoints, poles at ±90°.
    """
    lat_rad = np.deg2rad(lat)

    # Cell edge latitudes: midpoints between adjacent points, poles at boundaries
    lat_edges = np.empty(len(lat) + 1)
    lat_edges[0]    = -np.pi / 2.0
    lat_edges[-1]   =  np.pi / 2.0
    lat_edges[1:-1] = 0.5 * (lat_rad[:-1] + lat_rad[1:])

    # Uniform longitude spacing (CAM uses regular lon grids)
    dlon_rad = np.deg2rad(lon[1] - lon[0]) if len(lon) > 1 else 2.0 * np.pi

    # Cell area proportional to dlon * delta(sin(lat)); broadcast to (nlat, nlon)
    d_sin_lat = np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1])
    weights   = d_sin_lat[:, np.newaxis] * np.full(len(lon), dlon_rad)

    return weights / weights.sum()


def global_mean_2d(var2d, weights):
    """
    Compute area-weighted global mean of a single 2D field.
    weights is the normalized (nlat, nlon) array from build_area_weights.
    Masked cells are excluded from both numerator and denominator.
    Returns np.nan if the entire field is masked.
    """
    ma = np.ma.asarray(var2d)
    if ma.mask is np.ma.nomask:
        return float(np.ma.average(ma, weights=weights))
    if ma.mask.all():
        return np.nan
    return float(np.ma.average(ma, weights=weights))


def global_mean_profile(var3d, weights):
    """
    Area-weighted global mean of each model level of a 3D field.
    var3d has shape (nlev, nlat, nlon); returns a (nlev,) array, one
    global_mean_2d per level, so masking behaves exactly as for 2D fields.
    No vertical interpolation: level k is model level k (top first).
    """
    return np.array([global_mean_2d(var3d[k], weights)
                     for k in range(var3d.shape[0])], dtype=float)


def select_monthly_files(root_path, case_id, prefix, start_year, n_months):
    """
    Monthly files for a case, selected by their YYYY-MM stamp: year >=
    start_year, keeping the first n_months.

    Select by the actual YYYY-MM stamp rather than by position in the sorted
    glob. Positional slicing is fragile: any stray file that sorts ahead of
    start_year (e.g. a year-0000 spin-up file, or another history stream the
    regex still matches) shifts the window and silently drops a real month
    off the tail. The file scan in trend.py selects by constructing expected
    names, so it never sees that error; matching by date here keeps the two
    phases in agreement.
    """
    date_pattern = re.compile(r'\.(\d{4})-(\d{2})\.nc$')
    all_files = sorted(glob.glob(f"{root_path}/{case_id}{prefix}*.nc"))
    files = []
    for f in all_files:
        m = date_pattern.search(f)
        if m is None:
            continue
        if int(m.group(1)) < start_year:
            continue
        files.append(f)
        if len(files) == n_months:
            break
    return files


def _weights_key(weights):
    """Short fingerprint of the area weights: a cached mean is only valid on the grid it was computed on."""
    return hashlib.sha1(np.ascontiguousarray(weights, dtype=float).tobytes()).hexdigest()[:16]


def _file_key(filepath):
    st = os.stat(filepath)
    return st.st_size, st.st_mtime_ns


def load_cache(cache_path, wkey):
    """
    Per-file global means saved by earlier runs, {basename: record}.

    The cache is append-only JSON lines, one record per file read. Records for
    the same file merge in order, later values winning, so a run that adds a
    variable appends only that variable. A record is ignored unless its
    CACHE_VERSION and weights fingerprint match; the caller also checks the
    file's size and mtime, so a rewritten history file is read again. A run
    killed mid-write leaves at most one partial last line, which is skipped.
    """
    entries = {}
    if cache_path is None or not os.path.exists(cache_path):
        return entries
    with open(cache_path) as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except ValueError:
                continue
            if rec.get("version") != CACHE_VERSION or rec.get("weights") != wkey:
                continue
            ent = entries.get(rec["file"])
            if ent is None or (ent["size"], ent["mtime_ns"]) != (rec["size"], rec["mtime_ns"]):
                ent = {"size": rec["size"], "mtime_ns": rec["mtime_ns"], "vars": {}, "profiles": {}}
                entries[rec["file"]] = ent
            ent["vars"].update(rec.get("vars", {}))
            ent["profiles"].update(rec.get("profiles", {}))
    return entries


def _open_cache_for_append(cache_path):
    d = os.path.dirname(cache_path)
    if d:
        os.makedirs(d, exist_ok=True)
    fh = open(cache_path, "a+")
    # A killed run can leave a partial last line; start the next record on a fresh one.
    if fh.tell() > 0:
        fh.seek(fh.tell() - 1)
        if fh.read(1) != "\n":
            fh.write("\n")
    return fh


def _warn_missing(prefix, vname, filepath):
    key = (prefix, vname)
    if key not in _warned_missing:
        print(f"  WARNING: variable '{vname}' not found in {os.path.basename(filepath)}, storing NaN")
        _warned_missing.add(key)


def _reduce_file(ncid, filepath, prefix, varnames, profile_vars, weights):
    """
    Global means for one open file: {var: float or None} for 2D fields
    (None = absent from the file) and {var: [per-level means]} for profile
    fields plus PMID when any are requested.
    """
    vals = {}
    for vname in varnames:
        if vname in ncid.variables:
            raw = ncid.variables[vname][0, :, :]
            vals[vname] = global_mean_2d(np.ma.masked_equal(raw, -999.0), weights)
        else:
            vals[vname] = None
    profs = {}
    if profile_vars:
        for name in list(profile_vars) + ['PS', 'hyam', 'hybm', 'P0']:
            if name not in ncid.variables:
                raise KeyError(f"profile variable '{name}' not found in "
                               f"{os.path.basename(filepath)}")
        for vname in profile_vars:
            var = ncid.variables[vname]
            if var.ndim != 4:
                raise ValueError(f"profile variable '{vname}' is not 3D "
                                 f"(time, lev, lat, lon): dims {var.dimensions}")
            raw = np.ma.masked_equal(var[0, :, :, :], -999.0)
            profs[vname] = global_mean_profile(raw, weights).tolist()
        ps = global_mean_2d(np.ma.masked_equal(ncid.variables['PS'][0, :, :], -999.0),
                            weights)
        hyam = np.asarray(ncid.variables['hyam'][:], dtype=float)
        hybm = np.asarray(ncid.variables['hybm'][:], dtype=float)
        p0 = float(np.asarray(ncid.variables['P0'][...]))
        profs['PMID'] = (hyam * p0 + hybm * ps).tolist()
    return vals, profs


def read_monthly_fields(root_path, case_id, prefix, varnames,
                        start_year, n_months, weights, profile_vars=(),
                        cache_path=None):
    """
    Read monthly netCDF files sequentially, in one pass, and compute
    area-weighted global means: of each 2D variable in varnames, and of each
    model level of each 3D variable in profile_vars.

    With cache_path, each file's means are appended to that JSON-lines file
    as soon as the file is read, and a file whose means are all already
    cached (same size, mtime, weights and CACHE_VERSION) is not opened. A
    killed run therefore loses at most the file it was on, and a rerun or a
    longer time series reads only the files not yet seen. See load_cache.

    Returns:
        out      -- numpy array, shape (n_months, len(varnames)), global means
        profiles -- dict {name: (n_months, nlev) array} for each profile
                    variable, plus 'PMID' (the mean pressure of each model
                    level, Pa: hyam*P0 + hybm*<PS>, with <PS> the area-weighted
                    global-mean surface pressure) when profile_vars is given;
                    {} otherwise
        files    -- list of file paths that were actually read

    Assumptions to verify against actual model output:
      - 2D variables are stored as (time, lat, lon) and 3D ones as
        (time, lev, lat, lon); index 0 is the single monthly snapshot.
      - weights shape (nlat, nlon) matches the spatial dims of each variable.
      - profile files carry hyam, hybm, P0 and PS (standard cam.h0 output).
    """
    files = select_monthly_files(root_path, case_id, prefix, start_year, n_months)

    out = np.zeros((len(files), len(varnames)), dtype=float)
    profiles = {}
    prof_names = list(profile_vars) + ['PMID'] if profile_vars else []

    wkey = _weights_key(weights)
    cache = load_cache(cache_path, wkey)
    cache_fh = _open_cache_for_append(cache_path) if cache_path else None
    n_hit = 0

    try:
        for i, filepath in tqdm(enumerate(files), total=len(files),
                                desc="reading files", unit="file"):
            base = os.path.basename(filepath)
            size, mtime_ns = _file_key(filepath)
            ent = cache.get(base)
            if ent is None or (ent["size"], ent["mtime_ns"]) != (size, mtime_ns):
                ent = {"vars": {}, "profiles": {}}
            need_vars = [v for v in varnames if v not in ent["vars"]]
            need_prof = [v for v in profile_vars if v not in ent["profiles"]]
            if prof_names and 'PMID' not in ent["profiles"]:
                need_prof = list(profile_vars)

            if need_vars or need_prof:
                with nc.Dataset(filepath, 'r') as ncid:
                    vals, profs = _reduce_file(ncid, filepath, prefix, need_vars,
                                               need_prof, weights)
                if cache_fh is not None:
                    cache_fh.write(json.dumps({"version": CACHE_VERSION, "weights": wkey,
                                               "file": base, "size": size,
                                               "mtime_ns": mtime_ns,
                                               "vars": vals, "profiles": profs}) + "\n")
                    cache_fh.flush()
                ent = {"vars": {**ent["vars"], **vals},
                       "profiles": {**ent["profiles"], **profs}}
            else:
                n_hit += 1

            for j, vname in enumerate(varnames):
                v = ent["vars"][vname]
                if v is None:
                    _warn_missing(prefix, vname, filepath)
                    out[i, j] = np.nan
                else:
                    out[i, j] = v
            for name in prof_names:
                prof = np.asarray(ent["profiles"][name], dtype=float)
                if name not in profiles:
                    profiles[name] = np.full((len(files), prof.size), np.nan)
                profiles[name][i] = prof
    finally:
        if cache_fh is not None:
            cache_fh.close()

    if cache_path is not None:
        print(f"  cache: {n_hit} of {len(files)} files reused ({cache_path})")

    return out, profiles, files


def read_monthly_files(root_path, case_id, prefix, varnames,
                       start_year, n_months, weights, cache_path=None):
    """
    Area-weighted global means of 2D variables only; see read_monthly_fields.

    Returns:
        out   -- numpy array, shape (n_months, len(varnames)), global means
        files -- list of file paths that were actually read
    """
    out, _, files = read_monthly_fields(root_path, case_id, prefix, varnames,
                                        start_year, n_months, weights,
                                        cache_path=cache_path)
    return out, files


def compute_running_means(vavg_vec, int1, int2):
    """
    Compute two causal rolling-window means and their per-year slopes.

    Window int1 (annual, 12 months) and int2 (decadal, 120 months).
    Before the window is full the mean is taken over all available data
    and the slope is computed relative to the first value, matching the
    original per-timestep logic.

    Slope units: [variable units] / year.

    Returns (intavg1, intavg2, slope1, slope2), each a 1D array of
    length len(vavg_vec).
    """
    N   = len(vavg_vec)
    cs  = np.cumsum(np.insert(vavg_vec, 0, 0.0))
    idx = np.arange(N)

    def _running_mean(window):
        # i < window: mean over vavg_vec[0 : i+1]  (cs[i+1] - cs[0]) / (i+1)
        # i >= window: mean over vavg_vec[i-window : i]  (cs[i] - cs[i-window]) / window
        end   = np.where(idx < window, idx + 1, idx)
        start = np.where(idx < window, 0,        idx - window)
        denom = np.where(idx < window, idx + 1,  window)
        return (cs[end] - cs[start]) / denom

    intavg1 = _running_mean(int1)
    intavg2 = _running_mean(int2)

    def _slope(intavg, window):
        # i < window: (intavg[i] - intavg[0]) / (window/12)
        # i >= window: (intavg[i] - intavg[i-window]) / (window/12)
        ref = np.where(idx < window, 0, idx - window)
        return (intavg - intavg[ref]) / (window / 12.0)

    slope1 = _slope(intavg1, int1)
    slope2 = _slope(intavg2, int2)

    return intavg1, intavg2, slope1, slope2
