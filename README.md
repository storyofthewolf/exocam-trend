# trend

This family of scripts produces time-series diagnostics from CAM or ExoCAM
simulations using monthly mean netCDF output. The primary use case is evaluating
model convergence and diagnosing residual trends in energy balance, surface
temperature, sea ice, and the hydrological cycle.

The package assumes directory structure follows NCAR standards (`$RUNDIR` and
`$DOUT_S_ROOT` / archive directory), or a flat local test directory via `--testdir`.

**Requires:** `numpy`, `matplotlib`, `netCDF4`, `tqdm`

## Files

| File | Description |
|------|-------------|
| `trend.py` | Main driver script |
| `trend_core.py` | Core computation: area weights, file I/O, global mean averaging, running statistics |
| `trend_utils.py` | `vars.in` parsing, screen/text output, and time-series plotting |
| `vars.in` | Variable namelist: which fields to read, print, and plot for each component model |

## vars.in format

`vars.in` has three blocks separated by blank lines: **read**, **print**, and **plot**.
Each block has one line per component model (atmosphere, ice, land), with variable
names separated by spaces.

```
  // input fields //
TS FLNT FSNT FLNS FSNS LHFLX SHFLX ICEFRAC FSDTOA
Tsfc qi qs hi hs vicen005
TG
  // output fields //
TS FLNT FSNT FLNS FSNS LHFLX SHFLX ICEFRAC energy
Tsfc qi qs hi hs
TG
  // create plots //
TS energy
Tsfc
TG
```

The special token `energy` in the print or plot block triggers automatic
computation of TOA and surface energy balance from primary flux variables
(requires FSNT, FLNT, FSNS, FLNS, LHFLX, SHFLX in the read list).

## Usage

```
usage: trend.py [-h] [-y Y] [-n N] [-p P] [-a A]
                [--cam] [--cice] [--clm]
                [--rundir] [--testdir PATH]
                [--plots] [--save-data] [--timing]
                [--int1 INT1] [--int2 INT2] [--profile VARS]
                case_id
```

### Positional arguments

| Argument | Description |
|----------|-------------|
| `case_id` | Simulation case name (must match the filename prefix of the netCDF files) |

### Optional arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-y Y` | 1 | Start year for the time series |
| `-n N` | 6000 | Maximum number of months to process |
| `-p P` | (none) | Print to screen every P months; omit to suppress running output |
| `-a A` | 2 | Value shown in screen output: `0`=instantaneous monthly, `1`=1-year running mean, `2`=10-year running mean |
| `--cam` | off | Read atmosphere model (`cam.h0`) files |
| `--cice` | off | Read sea ice model (`cice.h`) files |
| `--clm` | off | Read land model (`clm2.h0`) files |
| `--rundir` | off | Read from run directory (`$dir/rundir/<case_id>/run/`) instead of archive |
| `--testdir PATH` | off | Read all component files from a single flat directory (local testing); files must still follow CAM naming conventions (`case_id.cam.h0.YYYY-MM.nc`) |
| `--plots` | off | Generate time-series line plots after processing |
| `--save-data` | off | Write global mean time series to text files in `data/` |
| `--timing` | off | Print wall-clock timing summary at end of run |
| `--int1 INT1` | 1 | Short averaging window in years |
| `--int2 INT2` | 10 | Long averaging window in years |
| `--profile VARS` | (none) | Comma-separated 3D `cam.h0` fields (e.g. `T,Q`): also compute per-level global means on model levels (see [Vertical profiles](#vertical-profiles)); written with `--save-data`. Needs `--cam`. |

At least one of `--cam`, `--cice`, or `--clm` must be specified.

### File path modes

| Mode | Path constructed |
|------|-----------------|
| Archive (default) | `$dir/archive/<case_id>/atm/hist/` etc. |
| `--rundir` | `$dir/rundir/<case_id>/run/` for all components |
| `--testdir PATH` | `PATH/` for all components (no subdirectory structure) |

The base directory `dir` is hardcoded near the top of `trend.py` and must be
updated for your system.

## Architecture

**`trend.py`** execution is structured in five sequential phases:

1. **File scan** — iterates month-by-month to check file existence, determine
   `N_actual` (number of available months), and record date strings. No data is read.
2. **Read and average** — calls `core.read_monthly_files()` once per active component
   to read all monthly files and compute area-weighted global means.
3. **Energy balance** — if `energy` is requested in `vars.in`, computes TOA and
   surface energy balance for all timesteps from the averaged flux variables.
4. **Running statistics** — calls `core.compute_running_means()` once per variable
   column to compute 1-year and 10-year rolling means and their slopes (O(N) via cumsum).
5. **Output loop** — iterates over timesteps printing to screen at the requested interval;
   optionally writes text files and generates plots.

**`trend_core.py`** functions:

| Function | Description |
|----------|-------------|
| `build_area_weights(lon, lat)` | Returns a normalized `(nlat, nlon)` weight array using a staggered lat grid (cell edges at midpoints between grid points, poles at ±90°). Called once at startup. |
| `global_mean_2d(var2d, weights)` | Area-weighted mean of a single 2D field. Uses `np.ma.average` so masked cells (including -999.0 sentinels) are excluded from both numerator and denominator. |
| `global_mean_profile(var3d, weights)` | Area-weighted mean of each model level of a `(nlev, nlat, nlon)` field: `global_mean_2d` per level, no vertical interpolation. |
| `read_monthly_fields(..., profile_vars=())` | Reads monthly netCDF files sequentially via netCDF4 in one pass: `global_mean_2d` for each 2D variable and `global_mean_profile` for each 3D profile variable, plus `PMID`, the mean pressure of each level. Returns `(out, profiles, files)`. Progress is shown via `tqdm`. |
| `read_monthly_files(root_path, case_id, prefix, varnames, start_year, n_months, weights)` | 2D-only wrapper of `read_monthly_fields`: returns a `(n_months, len(varnames))` array of global means plus the list of files read. |
| `compute_running_means(vavg_vec, int1, int2)` | Computes causal rolling-window means and slopes for annual (`int1=12`) and decadal (`int2=120`) windows using a cumsum trick. Returns `(intavg1, intavg2, slope1, slope2)`, each a 1D array of the same length as the input. |

**`trend_utils.py`** functions:

| Function | Description |
|----------|-------------|
| `read_request_var()` | Parses `vars.in` into 9 lists (read/print/plot × atm/ice/lnd) and validates that all print and plot variables are present in the read list. |
| `print2screen(...)` | Prints a formatted table row at the current timestep. The `-a` flag selects which average is displayed (monthly/annual/decadal). |
| `atm_energy_calc(atmvars, vavg_vecA)` | Derives `etop = FSNT − FLNT` and `ebot = FSNS − FLNS − LHFLX − SHFLX` from already-averaged values. |
| `timeSeriesPlots(...)` | Renders matplotlib line plots for each requested variable, showing monthly, 1-year, and 10-year averages. Ice and land plotting not yet implemented. |
| `print2text(...)` | Writes time-series data to `data/<case_id>_<firstDate>-<lastDate>_{cam,cice,clm}.txt`. |
| `print_profiles2text(...)` | Writes one per-level file per profile variable (and `PMID`) to `data/<case_id>_<firstDate>-<lastDate>_camlev_<VAR>.txt`. |

## Examples

Standard run from the archive directory, atmosphere only, starting at year 10:
```bash
python trend.py $my_case_name --cam -y 10 -p 100 --plots
```

Read both atmosphere and sea ice, print decadal averages every 12 months, write data files:
```bash
python trend.py $my_case_name --cam --cice -p 12 -a 2 --save-data
```

Also write per-level global-mean T and Q profiles:
```bash
python trend.py $my_case_name --cam -p 12 --save-data --profile T,Q
```

Local development with a flat test directory:
```bash
python trend.py my_case --cam -p 12 --testdir /path/to/test/data --plots
```

## Output

- **Screen** — formatted table with timestep index and global mean values at the
  interval set by `-p`. The column shown per variable is selected by `-a`.
- **`data/`** — text files named `<case_id>_<firstDate>-<lastDate>_cam.txt`
  (written with `--save-data`), and with `--profile` one
  `<case_id>_<firstDate>-<lastDate>_camlev_<VAR>.txt` per profile field plus
  `..._camlev_PMID.txt`.
- **`plots/`** — time-series line plots saved or displayed interactively
  (written with `--plots`).

## Vertical profiles

`--profile T,Q` reduces 3D `cam.h0` fields to area-weighted global means on
each **model level** (no interpolation to pressure surfaces), read in the same
pass over the files as the 2D variables, so the cost is the extra bytes of the
named fields (about 0.7 MB per field per month at 4x5 with 51 levels), not
another set of file opens. Each field is written to its own file, one row per
month, native monthly means only (no running means: the consumer picks its
own averaging window):

```
month  L01  L02  ...  L51
1      269.5  267.3  ...  378.0
```

`L01` is the model top. `..._camlev_PMID.txt` holds each level's mean pressure
in Pa, `hyam*P0 + hybm*<PS>` with `<PS>` the area-weighted global-mean surface
pressure of that month. exocam-accelerate reads these files to measure how the
atmosphere warms per level relative to the surface.
