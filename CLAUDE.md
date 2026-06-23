# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`trend` is a Python tool for producing time-series diagnostics from NCAR CESM (CAM/ExoCAM) climate model simulations. It reads monthly-mean NetCDF output files, computes area-weighted global means, calculates running averages and trend slopes over configurable short and long windows, then optionally prints results to screen, writes text files, and/or generates line plots. The primary use case is diagnosing model convergence and residual trends in energy balance, surface temperature, sea ice, and the hydrological cycle.

## Dependencies

- `netCDF4`, `numpy`, `matplotlib`, `tqdm`
- `xarray` and `exocampy_tools` are not used

## Running

```bash
python trend.py <case_id> --cam [--cice] [--clm] \
    [-y START_YEAR] [-n N_MONTHS] [-p PRINT_INTERVAL] [-a AVG_FREQ] \
    [--int1 YEARS] [--int2 YEARS] \
    [--plots] [--data] [--timing] \
    [--rundir] [--testdir PATH]
```

**Example:**
```bash
python trend.py $my_case_name --cam -p 100 -y 10 --plots
```

`trend.py` must be run from the repository root (it opens `vars.in` by relative path). Output text files go to `data/`, plots to `plots/`.

The hardcoded base directory `dir` near the top of `trend.py` must point to the CESM scratch space (currently `/gpfsm/dnb33/etwolf/cesm_scratch/`). Adjust this for your system before running.

### Command-line arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `case_id` | (required) | Simulation case name; must match the NetCDF filename prefix |
| `-y Y` | 1 | Start year for the time series |
| `-n N` | 6000 | Maximum number of months to process |
| `-p P` | (none) | Print to screen every P months; omit to suppress running output |
| `-a A` | 2 | Average shown in screen output: `0`=monthly, `1`=short-window, `2`=long-window |
| `--int1 N` | 1 | Short averaging window in years |
| `--int2 N` | 10 | Long averaging window in years |
| `--cam` | off | Read atmosphere model (`cam.h0`) files |
| `--cice` | off | Read sea ice model (`cice.h`) files |
| `--clm` | off | Read land model (`clm2.h0`) files |
| `--rundir` | off | Read from run directory instead of archive |
| `--testdir PATH` | off | Read all components from a flat directory (local testing) |
| `--plots` | off | Generate time-series line plots after processing |
| `--data` | off | Write global mean time series to text files in `data/` |
| `--timing` | off | Print wall-clock timing summary at end of run |

At least one of `--cam`, `--cice`, or `--clm` must be specified.

### File path modes

Three modes control where files are read from, checked in this order:

| Mode | Flag | Path |
|------|------|------|
| Test/local | `--testdir PATH` | `PATH/` flat directory, all components |
| Run directory | `--rundir` | `$dir/rundir/<case_id>/run/` |
| Archive (default) | _(none)_ | `$dir/archive/<case_id>/atm/hist/` etc. |

In `--testdir` mode, files must still follow CAM naming conventions (`case_id.cam.h0.YYYY-MM.nc`); no subdirectory structure is expected.

## File structure

| File | Description |
|------|-------------|
| `trend.py` | Main driver script |
| `trend_core.py` | Core computation: area weights, file I/O, global mean averaging, running statistics |
| `trend_utils.py` | `vars.in` parsing, screen/text output, time-series plotting, final summary |
| `vars.in` | Variable namelist: which fields to read, print, and plot for each component model |
| `data/` | Output directory for text files (must exist; not auto-created) |
| `plots/` | Output directory for plot files (must exist; not auto-created) |

## Architecture

**`trend.py`** — main driver. Execution proceeds in six phases:

1. **File scan** — iterates month-by-month checking file existence to determine `N_actual`, `firstDate`, `lastDate`, and `time_vec*` arrays. No data is read here. Exits cleanly when the next expected file is absent.
2. **Interval validation** — after the scan, `int1` and `int2` are checked against each other and against `N_actual`, with warnings and automatic capping if they exceed available data.
3. **Read and average** — calls `core.read_monthly_files()` once per active component to read all monthly files sequentially and compute area-weighted global means.
4. **Energy balance** — if `energy` is in the print list, computes TOA and surface energy balance for all `N_actual` timesteps from the averaged flux variables. Must happen before running means.
5. **Running statistics** — calls `core.compute_running_means()` once per variable column over the full filled series (O(N) via cumsum).
6. **Output loop** — post-processing pass over already-computed arrays: prints to screen at the `-p` interval, then optionally writes text files and generates plots. Always ends with `print_final_summary()`.

**`trend_core.py`** — core computation functions:

- `build_area_weights(lon, lat)` — staggered-grid normalized area weights, shape `(nlat, nlon)`, called once at startup. Cell edges are midpoints between adjacent lat values, with ±90° at the poles. Assumes a regular longitude grid (uniform `dlon`).
- `global_mean_2d(var2d, weights)` — area-weighted mean of a single 2D field via `np.ma.average`; masked cells (including -999.0 sentinels) are excluded from both numerator and denominator. Returns `np.nan` if the entire field is masked.
- `read_monthly_files(root_path, case_id, prefix, varnames, start_year, n_months, weights)` — globs and date-filters all matching files (regex `\.\d{4}-\d{2}\.nc$`), slices to `start_year` and `n_months`, reads each file with netCDF4, calls `global_mean_2d` per variable, and returns a `(n_months, len(varnames))` numpy array plus the file list. Progress shown via `tqdm`. Missing variables store `np.nan` and are warned once per `(prefix, varname)` pair via the module-level `_warned_missing` set.
- `compute_running_means(vavg_vec, int1, int2)` — O(N) causal rolling mean via cumsum; returns `(intavg1, intavg2, slope1, slope2)`. Before the window is full, the mean is taken over all available data and slope is computed relative to the first value. Slopes are in units per year.

**`trend_utils.py`** — parsing, output, and plotting functions:

- `read_request_var()` — parses `vars.in` into 9 lists (read/print/plot × atm/ice/lnd) and validates that print/plot variables are a subset of read variables. `energy` is exempt from this check (it is derived, not read).
- `print2screen()` — tabular running output; `-a` flag selects which average is shown (0=monthly, 1=short-window, 2=long-window). The `print_offset` variable maps `avgfreq` to the correct column index within the packed `atmout` array.
- `atm_energy_calc()` — derives `etop = FSNT - FLNT` and `ebot = FSNS - FLNS - LHFLX - SHFLX`. Requires those six variables in the atmosphere read list.
- `timeSeriesPlots()` — renders matplotlib line plots per requested variable, showing monthly, short-window, and long-window averages. Ice/land plotting branches are not yet implemented.
- `print2text()` — writes `data/<case_id>_<firstDate>-<lastDate>_cam.txt`; partially implemented — writes only the first atmosphere variable (hardcoded column index 4). The ice branch has a bug (writes `vavg_vecI[i,4]` which is wrong for most ice configurations).
- `print_final_summary()` — prints a formatted table of final-timestep values and running averages for all active components; always called at the end of `trend.py` regardless of `-p`.

**`vars.in`** — plain-text namelist with three blocks (read / print / plot), each with one line per component model (atm, ice, lnd). Comment lines beginning with `#` after the data blocks are ignored. The special token `energy` in print/plot blocks triggers `atm_energy_calc` rather than a direct variable lookup.

## Array indexing conventions

Variable arrays carry a leading offset of coordinate dimensions: atmosphere arrays have offset 4 (`time, lon, lat, lev`); ice and land arrays have offset 1 (`time`). `atm_vars_offset = 4` and `ice_vars_offset = 1` are module-level constants in `trend_utils.py`. When indexing `vavg_vecA`, add `atm_vars_offset` to the 0-based index into `atmvars_in`.

The derived `energy` variables are appended at the very end of the atmosphere array: `nvtotA-2` = `etop` (TOA), `nvtotA-1` = `ebot` (surface). `nvtotA` is set two larger than the normal variable count when `energy` is requested.

The `read_monthly_files` return array has shape `(n_months, len(varnames))` with no coordinate offset — index 0 is `varnames[0]`. `trend.py` writes these into `vavg_vec*` starting at `nv1dA`/`nv1dI`/`nv1dL` to maintain the offset convention expected by `trend_utils`.

## File naming conventions

CESM archive NetCDF files are expected at:
```
<dir>/archive/<case_id>/atm/hist/<case_id>.cam.h0.<YYYY>-<MM>.nc
<dir>/archive/<case_id>/ice/hist/<case_id>.cice.h.<YYYY>-<MM>.nc
<dir>/archive/<case_id>/lnd/hist/<case_id>.clm2.h0.<YYYY>-<MM>.nc
```
Years are zero-padded to 4 digits (e.g., `0010`).

`read_monthly_files` uses `(start_year - 1) * 12` as a slice offset into the sorted glob result to skip pre-`START_YEAR` files. **This assumes all years from 1 to `START_YEAR - 1` have exactly 12 complete monthly files.** Gaps or partial years before `START_YEAR` will silently misalign the slice.

## Constraints and assumptions

- **Monthly cadence only.** Only `nhtfrq=0` (monthly mean) output is supported. Sub-monthly or annual history files will not work.
- **2D variables only.** `read_monthly_files` reads `var[0, :, :]` — the first time index, all lats and lons. 3D variables (e.g., pressure-level fields) are not yet supported.
- **Regular lon grid.** `build_area_weights` assumes uniform longitude spacing; the first `dlon` value is used for all cells.
- **Grid consistency.** All files in a run are assumed to share the same `(nlat, nlon)` grid as the first file peeked at startup. No checking is done.
- **No partial years before START_YEAR.** The glob-slice offset assumes complete 12-file years prior to the start year.
- **`vars.in` must be present** in the working directory. `read_request_var` opens it by relative path with no fallback.
- **`data/` and `plots/` must exist.** These directories are not auto-created; running with `--data` or `--plots` will fail with a `FileNotFoundError` if they are absent.
- **GPFS / parallel filesystem latency.** The tool is designed for NCAR/NASA HPC GPFS scratch filesystems where individual file opens are expensive. The file scan phase (checking existence without reading) and the sequential-read design reflect this — random access patterns would be slower.

## Known limitations and bugs

- **`timeSeriesPlots` only handles atmosphere variables.** Ice and land branches exist but do nothing (noted `#!! routine is incomplete !!` in the source). The `do_ice`/`do_lnd` plotting paths remain intentionally stubbed for later.
- **`--int1` must be less than `--int2`.** If `int1 >= int2`, `trend.py` resets `int1` to 1 year with a warning but does not exit. Both are capped to `N_actual` if the run is shorter than the requested window.
