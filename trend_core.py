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
import re
import os
from tqdm import tqdm

# Module-level set to suppress duplicate missing-variable warnings
_warned_missing = set()


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


def read_monthly_files(root_path, case_id, prefix, varnames,
                       start_year, n_months, weights):
    """
    Read monthly netCDF files sequentially and compute area-weighted global
    means for each variable at each timestep.

    Files are matched by glob and filtered to YYYY-MM date stamps, then
    selected by the date stamp itself (year >= start_year, keeping the first
    n_months) so the result aligns exactly with the file scan in trend.py
    that determined N_actual.

    Returns:
        out   -- numpy array, shape (n_months, len(varnames)), global means
        files -- list of file paths that were actually read

    Assumptions to verify against actual model output:
      - Variables are stored as (time, lat, lon); var[0, :, :] is the single
        monthly snapshot in each file.
      - weights shape (nlat, nlon) matches the spatial dims of each variable.
    """
    date_pattern = re.compile(r'\.(\d{4})-(\d{2})\.nc$')
    all_files = sorted(glob.glob(f"{root_path}/{case_id}{prefix}*.nc"))

    # Select by the actual YYYY-MM stamp rather than by position in the sorted
    # glob. Positional slicing is fragile: any stray file that sorts ahead of
    # start_year (e.g. a year-0000 spin-up file, or another history stream the
    # regex still matches) shifts the window and silently drops a real month
    # off the tail. The file scan in trend.py selects by constructing expected
    # names, so it never sees that error; matching by date here keeps the two
    # phases in agreement.
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

    out = np.zeros((len(files), len(varnames)), dtype=float)

    for i, filepath in tqdm(enumerate(files), total=len(files), 
                        desc="reading files", unit="file"):
    #for i, filepath in enumerate(files):
        ncid = nc.Dataset(filepath, 'r')
        for j, vname in enumerate(varnames):
            if vname in ncid.variables:
                raw = ncid.variables[vname][0, :, :]
                masked = np.ma.masked_equal(raw, -999.0)
                out[i, j] = global_mean_2d(masked, weights)
            else:
                key = (prefix, vname)
                if key not in _warned_missing:
                    print(f"  WARNING: variable '{vname}' not found in {os.path.basename(filepath)}, storing NaN")
                    _warned_missing.add(key)
                out[i, j] = np.nan
        ncid.close()

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
