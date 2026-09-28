#!/usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test_trend_utils.py
#
# Self-contained regression check for print2text() / _write_component()
# in trend_utils.py. No netCDF data is needed: it builds synthetic
# time_vec / vavg_vec arrays the same shape and convention that trend.py
# produces, calls print2text(), and checks the written data/*.txt file.
#
# Regression covered: for a case with N_actual == n_months complete
# monthly files, print2text() used to write only N_actual - 1 rows, with
# every row's data shifted one month off from its printed label (the
# first month's data was silently dropped and the last labeled row
# actually held the second-to-last month's data). See CLAUDE.md commit
# history / bug report for exocam-accelerate integration.
#
# Run directly: `python3 test_trend_utils.py`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import numpy as np
import trend_utils as trend

atm_vars_offset = trend.atm_vars_offset


def _build_synthetic_atm(n_months, nvars=1):
    """
    Build time_vecA / vavg_vecA the way trend.py does for a fully complete
    run of n_months, with NT == N_actual == n_months (the exact condition
    that triggered the bug: -n N with a complete N-month archive).

    vavg_vecA[k, atm_vars_offset] is set to the 1-based month number so a
    written row's data value can be checked against its label directly.
    """
    NT = n_months
    time_vecA = np.zeros(NT, dtype=float)
    nvtotA = atm_vars_offset + nvars
    vavg_vecA = np.zeros((NT, nvtotA), dtype=float)
    for i in range(n_months):
        it = i + 1
        time_vecA[i] = it
        vavg_vecA[i, atm_vars_offset] = it  # 'TS' column carries month number
    intavg1_vecA = vavg_vecA.copy()
    intavg2_vecA = vavg_vecA.copy()
    slope1 = np.zeros_like(vavg_vecA)
    slope2 = np.zeros_like(vavg_vecA)
    return time_vecA, vavg_vecA, intavg1_vecA, intavg2_vecA, slope1, slope2


def _run_case(n_months):
    (time_vecA, vavg_vecA, intavg1_vecA, intavg2_vecA,
     slope1A, slope2A) = _build_synthetic_atm(n_months)

    atmvars_in  = ['TS']
    atmprint_in = ['TS']
    vnamesA = ['time', 'lon', 'lat', 'lev', 'TS']

    case_id   = '__test_trend_utils__'
    firstDate = '0001-01'
    lastDate  = f'{(n_months - 1) // 12 + 1:04d}-{(n_months - 1) % 12 + 1:02d}'

    outfile = f"data/{case_id}_{firstDate}-{lastDate}_cam.txt"
    if os.path.exists(outfile):
        os.remove(outfile)

    # dummy ice/land args -- disabled
    zeros = np.zeros((n_months, 1))
    trend.print2text(atmvars_in, [], [], atmprint_in, [], [],
                      True, vnamesA, time_vecA, vavg_vecA, intavg1_vecA, intavg2_vecA, slope1A, slope2A,
                      False, [], np.zeros(n_months), zeros, zeros, zeros, zeros, zeros,
                      False, [], np.zeros(n_months), zeros, zeros, zeros, zeros, zeros,
                      firstDate, lastDate, case_id)

    with open(outfile) as f:
        lines = f.readlines()
    os.remove(outfile)

    header = lines[0]
    data_lines = lines[1:]
    return header, data_lines


def test_row_count_matches_n_months():
    for n_months in (12, 359, 360, 1680):
        header, data_lines = _run_case(n_months)
        assert len(data_lines) == n_months, (
            f"n_months={n_months}: expected {n_months} data rows, got {len(data_lines)}")


def test_labels_and_data_are_not_shifted():
    n_months = 360
    header, data_lines = _run_case(n_months)

    assert len(data_lines) == n_months

    first_fields = data_lines[0].split()
    last_fields  = data_lines[-1].split()

    # month label is column 0; TS_native is column 1
    assert int(first_fields[0]) == 1, f"first row label should be month 1, got {first_fields[0]}"
    assert float(first_fields[1]) == 1.0, f"first row data should be month 1's value, got {first_fields[1]}"

    assert int(last_fields[0]) == n_months, f"last row label should be month {n_months}, got {last_fields[0]}"
    assert float(last_fields[1]) == float(n_months), (
        f"last row data should be month {n_months}'s value, got {last_fields[1]}")

    # spot-check every row: label i must carry value i (no one-month shift)
    for line in data_lines:
        fields = line.split()
        label = int(fields[0])
        value = float(fields[1])
        assert label == value, f"row label {label} does not match its data value {value}"


if __name__ == '__main__':
    test_row_count_matches_n_months()
    test_labels_and_data_are_not_shifted()
    print("All checks passed.")
