#!/usr/bin/env bash
#
# run_trend_batch.sh — batch-generate exocam-trend time series for a list of
# cases, in the one consistent configuration that exocam-accelerate's Tier-0
# hindcast harness expects to consume.
#
# This script lives in exocam-trend (it orchestrates trend.py) and is the
# generation half of the exocam-trend -> exocam-accelerate text interface.
# exocam-accelerate never runs it and never imports from here; the two repos
# meet only at the data/*.txt output.
#
# It runs on the HPC, where the model archives live. For each case it invokes
# trend.py with:
#   --cam --cice --clm   all three components (sea-ice + snowpack vars live in
#                        the cice.h stream -> the *_cice.txt file)
#   --save-data          write data/<case>_<first>-<last>_{cam,cice,clm}.txt
#   -p 1                 emit an entry every month (monthly cadence: the
#                        hindcast needs the raw monthly series, not decadal)
#   -y / -n              start year and number of months per case
#   --int1 / --int2      short/long running-mean windows in YEARS
#
# The variables extracted are whatever exocam-trend/vars.in lists. For the
# cold-case sea-ice/snowpack spin-up study that must include (cice line):
#   Tsfc qi qs hi hs vicen005      and (cam line) ICEFRAC + the energy budget.
# This script does NOT edit vars.in; verify it before the run (see --check).
#
# Usage:
#   ./run_trend_batch.sh [options] CASE [CASE ...]
#   ./run_trend_batch.sh [options] --case-file cases.txt
#
# Options:
#   --start-year N     first model year of the series          (default 1)
#   --nmonths N        number of months to integrate           (default 1800 = 150 yr)
#   --int1 N           short running-mean window, years        (default 1)
#   --int2 N           long running-mean window, years         (default 10)
#   --cam              read atmosphere (cam.h0) — enables explicit component mode
#   --cice             read sea ice (cice.h)   — enables explicit component mode
#   --clm              read land (clm2.h0)     — enables explicit component mode
#   --case-file FILE   read case IDs from FILE (one per line, # comments ok)
#   --profile VARS     comma-separated 3D cam.h0 fields (e.g. T,Q): also write
#                      per-level global-mean series on model levels, one file
#                      per field plus PMID: data/<case>_<first>-<last>_camlev_<VAR>.txt
#   --outdir DIR       collect the .txt files THIS run produced into DIR (copied
#                      from data/; older series of the same case are left behind)
#   --trend-dir DIR    path to the exocam-trend checkout       (default: script dir)
#   --python BIN       python interpreter to use               (default: python3)
#   --check            print the resolved plan and vars.in, then exit
#   --dry-run          print each trend.py command without running it
#   -h, --help         this help
#
# Components: with no --cam/--cice/--clm given, all three are read (the old
# default). Passing ANY of them switches to explicit mode — only the named
# components are read. Use this for cases missing a component: a run with no
# land model takes "--cam --cice", so trend.py never looks for clm2.h0 files
# that do not exist (a missing component makes the file scan find 0 timesteps
# and trend.py then errors on the empty date range).
#
# Example (Tier-0 generation, aquaplanet cases with no land model):
#   ./run_trend_batch.sh --case-file cold15.txt --cam --cice --nmonths 1800 \
#       --int1 1 --int2 10 --outdir data/tier0_series
#
# Notes on --int1/--int2 (they matter for the hindcast, see below):
#   int1/int2 are running-mean windows exocam-trend writes as the _int1/_int2
#   columns alongside each variable's _native column. The Tier-0 harness fits
#   tendencies off whichever column you tell it to; a longer int2 suppresses
#   more interannual noise and can stabilise the long-term-drift slope, at the
#   cost of lagging turning points. Because the columns are additive (one run
#   emits native+int1+int2), generating with a couple of int2 values costs
#   little and lets the harness compare them. Keep native (-p 1) always: it is
#   the ground truth the hindcast scores against.

set -euo pipefail

# ---- defaults -------------------------------------------------------------
START_YEAR=1
NMONTHS=1800          # 150 years
INT1=1
INT2=10
CASE_FILE=""
OUTDIR=""
PROFILE=""
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TREND_DIR="$SCRIPT_DIR"
PYTHON="python3"
CHECK=0
DRY_RUN=0
CASES=()

# Component selection. Passing any of --cam/--cice/--clm switches from the
# default (all three) to an explicit set: only the components you name are
# read. Cases without a land model, for instance, take "--cam --cice" so
# trend.py never scans for the clm2.h0 files that do not exist (a missing
# component file makes the scan find 0 timesteps and trend.py then crashes on
# the empty date range).
COMP_CAM=""
COMP_CICE=""
COMP_CLM=""
COMP_EXPLICIT=0

# ---- arg parsing ----------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-year) START_YEAR="$2"; shift 2 ;;
    --nmonths)    NMONTHS="$2";    shift 2 ;;
    --int1)       INT1="$2";       shift 2 ;;
    --int2)       INT2="$2";       shift 2 ;;
    --case-file)  CASE_FILE="$2";  shift 2 ;;
    --outdir)     OUTDIR="$2";     shift 2 ;;
    --profile)    PROFILE="$2";    shift 2 ;;
    --trend-dir)  TREND_DIR="$2";  shift 2 ;;
    --python)     PYTHON="$2";     shift 2 ;;
    --cam)        COMP_CAM="--cam";   COMP_EXPLICIT=1; shift ;;
    --cice)       COMP_CICE="--cice"; COMP_EXPLICIT=1; shift ;;
    --clm)        COMP_CLM="--clm";   COMP_EXPLICIT=1; shift ;;
    --check)      CHECK=1;         shift ;;
    --dry-run)    DRY_RUN=1;       shift ;;
    -h|--help)    sed -n '2,71p' "${BASH_SOURCE[0]}"; exit 0 ;;
    --*)          echo "unknown option: $1" >&2; exit 2 ;;
    *)            CASES+=("$1");   shift ;;
  esac
done

# Default to all three components; any explicit --cam/--cice/--clm narrows it.
if [[ "$COMP_EXPLICIT" -eq 0 ]]; then
  COMP_CAM="--cam"; COMP_CICE="--cice"; COMP_CLM="--clm"
fi
COMPONENTS=()
SUFFIXES=()   # output-file suffixes trend.py writes per selected component
[[ -n "$COMP_CAM" ]]  && { COMPONENTS+=("$COMP_CAM");  SUFFIXES+=("cam"); }
[[ -n "$COMP_CICE" ]] && { COMPONENTS+=("$COMP_CICE"); SUFFIXES+=("cice"); }
[[ -n "$COMP_CLM" ]]  && { COMPONENTS+=("$COMP_CLM");  SUFFIXES+=("clm"); }
if [[ ${#COMPONENTS[@]} -eq 0 ]]; then
  echo "no components selected (need at least one of --cam/--cice/--clm)" >&2
  exit 2
fi

# case IDs from a file, if given (skip blank lines and # comments)
if [[ -n "$CASE_FILE" ]]; then
  if [[ ! -f "$CASE_FILE" ]]; then
    echo "case file not found: $CASE_FILE" >&2; exit 1
  fi
  while IFS= read -r line; do
    line="${line%%#*}"                       # strip trailing comment
    line="$(echo "$line" | xargs || true)"   # trim whitespace
    [[ -n "$line" ]] && CASES+=("$line")
  done < "$CASE_FILE"
fi

if [[ ${#CASES[@]} -eq 0 ]]; then
  echo "no cases given (positional args or --case-file)" >&2
  echo "run with --help for usage" >&2
  exit 2
fi

TREND_PY="$TREND_DIR/trend.py"
VARS_IN="$TREND_DIR/vars.in"
if [[ ! -f "$TREND_PY" ]]; then
  echo "trend.py not found at $TREND_PY (set --trend-dir)" >&2; exit 1
fi

# ---- plan summary ---------------------------------------------------------
echo "== exocam-trend batch plan =="
echo "  trend.py    : $TREND_PY"
echo "  python      : $PYTHON"
echo "  start year  : $START_YEAR"
echo "  nmonths     : $NMONTHS  ($(python3 -c "print(f'{$NMONTHS/12:.1f}')") yr)"
echo "  int1/int2   : $INT1 / $INT2  (years)"
echo "  cadence     : -p 1 (monthly)"
echo "  components  : ${COMPONENTS[*]}"
[[ -n "$PROFILE" ]] && echo "  profiles    : $PROFILE  (per model level)"
echo "  cases       : ${#CASES[@]}"
printf '                %s\n' "${CASES[@]}"
[[ -n "$OUTDIR" ]] && echo "  collect ->  : $OUTDIR"
echo

if [[ "$CHECK" -eq 1 ]]; then
  echo "== vars.in (verify the sea-ice/snowpack vars are present) =="
  if [[ -f "$VARS_IN" ]]; then
    cat "$VARS_IN"
  else
    echo "  (vars.in not found at $VARS_IN)"
  fi
  echo
  echo "check complete; no runs performed."
  exit 0
fi

[[ -n "$OUTDIR" ]] && mkdir -p "$OUTDIR"

# ---- run ------------------------------------------------------------------
# trend.py writes into <cwd>/data/, so run from the exocam-trend dir.
cd "$TREND_DIR"

failed=()
for case_id in "${CASES[@]}"; do
  echo "---- $case_id ----"
  cmd=("$PYTHON" "$TREND_PY" "$case_id"
       -y "$START_YEAR" -n "$NMONTHS" -p 1
       --int1 "$INT1" --int2 "$INT2"
       "${COMPONENTS[@]}" --save-data)
  [[ -n "$PROFILE" ]] && cmd+=(--profile "$PROFILE")

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '  '; printf '%q ' "${cmd[@]}"; echo
    continue
  fi

  # Marker for "written by this run": data/ keeps every earlier series of the
  # case (other spans), and consumers merge all of a case's files in a
  # directory, so only files newer than this marker are collected.
  marker="data/.collect-marker.$$"
  touch "$marker"
  if "${cmd[@]}"; then
    if [[ -n "$OUTDIR" ]]; then
      # collect the files this case just produced (one per selected component,
      # plus one per profile field)
      shopt -s nullglob
      produced=()
      for f in data/"${case_id}"_*.txt; do
        [[ "$f" -nt "$marker" ]] || continue
        for suffix in "${SUFFIXES[@]}"; do
          if [[ "$f" == *_"${suffix}".txt || ( "$suffix" == cam && "$f" == *_camlev_*.txt ) ]]; then
            produced+=("$f"); break
          fi
        done
      done
      shopt -u nullglob
      if [[ ${#produced[@]} -gt 0 ]]; then
        cp -v "${produced[@]}" "$OUTDIR"/
      else
        echo "  WARNING: no data/${case_id}_*.txt produced" >&2
      fi
    fi
  else
    echo "  FAILED: $case_id" >&2
    failed+=("$case_id")
  fi
  rm -f "$marker"
done

echo
if [[ ${#failed[@]} -gt 0 ]]; then
  echo "== done, with ${#failed[@]} failure(s): ${failed[*]} =="
  exit 1
fi
echo "== done: ${#CASES[@]} case(s) processed =="
