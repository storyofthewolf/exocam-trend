#!/usr/bin/env python  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trend.py
#
# Author Eric Wolf
# May 2023
#
# Purpose:  Analyze time-series data from NCAR CESM simulations.
#           Used for evalulating time-dependent trends, energy balance,
#           hydrological cycle balance, and equilibrium acceleration
#           approaches.
#
# Notes:    Presently it only works with monthly mean outputs (nhtfrq=0).
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import netCDF4 as nc
import numpy   as np
import os
import sys
import time
import trend_utils as trend
import trend_core  as core
import argparse

def print_timing(label, elapsed):
    print(f"  {label:<40} {elapsed:6.1f}s")

# input arguments and options
parser = argparse.ArgumentParser()
parser.add_argument('case_id'     , type=str,   nargs=1, default=' ',  help='Set simulation time series case name')
parser.add_argument('-y'          , type=int,   default='1', help='Start year over which to begin timeseries')
parser.add_argument('-n'          , type=int,   default='6000', help='Number of months to integrate over')
parser.add_argument('-p'          , type=int,   default=None, help='Interval for screen output, in number of months (omit to suppress running output)')
parser.add_argument('-a'          , type=int,   default='2',  help='Average frequency for screen output, 0:monthly, 1:yearly, 2:decadal')
parser.add_argument('--cam',        action='store_true', help='read atmosphere model data')
parser.add_argument('--cice',       action='store_true', help='read sea ice model data')
parser.add_argument('--clm',        action='store_true', help='read land model data')
parser.add_argument('--rundir',     action='store_true', help='read files from run directory instead of archive')
parser.add_argument('--testdir',    type=str, default=None, help='read files directly from this fixed directory path (for local testing)')
parser.add_argument('--plots',      action='store_true', help='do lineplots at end of sequence')
parser.add_argument('--show',       action='store_true', help='display plots interactively in addition to saving')
parser.add_argument('--save-data',  action='store_true', dest='save_data', help='write global mean time series to text files in data/')
parser.add_argument('--timing',     action='store_true', help='print wall-clock timing summary at end of run')
parser.add_argument('--int1',       type=int, default=1,  help='Short averaging window in years (default: 1)')
parser.add_argument('--int2',       type=int, default=10, help='Long averaging window in years (default: 10)')
parser.add_argument('--profile',    type=str, default='', help='comma-separated 3D cam.h0 fields (e.g. T,Q): write per-level global-mean monthly series on model levels, plus PMID (needs --cam; written with --save-data)')
args = parser.parse_args()

# define case
case_id     = str(args.case_id[0])
START_YEAR  = int(args.y)
NT          = int(args.n)

# averaging interval for output
avgfreq     = int(args.a)

# 3D fields reduced to per-level global means
profile_vars = [v.strip() for v in args.profile.split(',') if v.strip()]
if profile_vars and not args.cam:
    print("--profile needs --cam (the profile fields are read from cam.h0)")
    quit()

read_rundir   = False
if args.rundir:  read_rundir = args.rundir

# the root path to your working directory
dir = '/gpfsm/dnb33/etwolf/cesm_scratch/'
# dir = '/discover/nobackup/tfauchez/cesm_scratch/'

# model prefixs
prefixA = '.cam.h0.'
prefixL = '.clm2.h0.'
prefixI = '.cice.h.'

# component models activated
do_atm = False
do_ice = False
do_lnd = False
if args.cam:  do_atm = True
if args.cice: do_ice = True
if args.clm:  do_lnd = True


# define time averaging intervals
int1 = args.int1 * 12
int2 = args.int2 * 12
#  interval to print outputs to screen (None means suppress running output)
print_int = args.p


#
# read vars.in
#
atmvars_in, icevars_in, lndvars_in, atmprint_in, iceprint_in, lndprint_in, atmplot_in, iceplot_in, lndplot_in  = trend.read_request_var()

# data from atmosphere model
# include standard quartet of time and space indexes
group1d = ['time', 'lon','lat','lev']
group2d = atmvars_in
nv1dA   = len(group1d)
vnamesA = np.concatenate((group1d, group2d))
nv2dA   = len(vnamesA)  - nv1dA
nvtotA  = len(vnamesA)

# add two variables to atmosphere array if "energy" requested
# for print output
value_to_find = 'energy'
indexp = np.where(np.array(atmprint_in) == value_to_find)[0]
energy_requested = len(indexp) > 0
if energy_requested:
    nvtotA = nvtotA + 2

# data from ice model
group1d = ['time']
group2d = icevars_in
nv1dI   = len(group1d)
vnamesI = np.concatenate((group1d, group2d))
nv2dI   = len(vnamesI)  - nv1dI
nvtotI  = len(vnamesI)

# data from land model
group1d = ['time']
group2d = lndvars_in
nv1dL   = len(group1d)
vnamesL = np.concatenate((group1d, group2d))
nv2dL   = len(vnamesL)  - nv1dL
nvtotL  = len(vnamesL)

#-----------------------------------------------------------------------------------------------------------------
# End User Specification Section
#-----------------------------------------------------------------------------------------------------------------

# Determine the root path based on CESM standard for
# $RUNDIR and $ARCHIVEDIR, or a flat testdir path.
if args.testdir is not None:
    # testdir mode: all components read from the same flat directory.
    # Files must still follow CAM naming conventions: case_id.cam.h0.YYYY-MM.nc
    # No case_id subdirectory structure is assumed.
    root_atm = args.testdir
    root_ice = args.testdir
    root_lnd = args.testdir
elif args.rundir:
    root_atm = dir + "rundir/" + case_id + "/run"
    root_ice = dir + "rundir/" + case_id + "/run"
    root_lnd = dir + "rundir/" + case_id + "/run"
else:
    root_atm = dir + "archive/" + case_id + "/atm/hist"
    root_ice = dir + "archive/" + case_id + "/ice/hist"
    root_lnd = dir + "archive/" + case_id + "/lnd/hist"

root_atm = ' '.join(root_atm.split())
root_ice = ' '.join(root_ice.split())
root_lnd = ' '.join(root_lnd.split())

#------------------------------------------
# Peak in first file to get lon, lat, lev
#------------------------------------------
timing = {}

t0 = time.time()

file_atm = f"{root_atm}/{case_id}.cam.h0.{START_YEAR:04d}-01.nc"

ncid = nc.Dataset(file_atm, 'r')
lon = ncid.variables['lon'][:]
nlon = lon.size
lat = ncid.variables['lat'][:]
nlat = lat.size
lev = ncid.variables['lev'][:]
nlev = lev.size
ncid.close()

timing['file peek'] = time.time() - t0
#print_timing('file peek', timing['file peek'])


#------------------------------------------------------
# Define data arrays based on user inputs
#-------------------------------------------------------
# time vectors
time_vecA = np.zeros(NT, dtype=float)
time_vecI = np.zeros(NT, dtype=float)
time_vecL = np.zeros(NT, dtype=float)

# global mean variables in time series
vavg_vecA = np.zeros((NT, nvtotA), dtype=float)
vavg_vecI = np.zeros((NT, nvtotI), dtype=float)
vavg_vecL = np.zeros((NT, nvtotL), dtype=float)

# time-interval averages of global mean quantities
intavg1_vecA = np.zeros((NT, nvtotA), dtype=float)
intavg1_vecI = np.zeros((NT, nvtotI), dtype=float)
intavg1_vecL = np.zeros((NT, nvtotL), dtype=float)
intavg2_vecA = np.zeros((NT, nvtotA), dtype=float)
intavg2_vecI = np.zeros((NT, nvtotI), dtype=float)
intavg2_vecL = np.zeros((NT, nvtotL), dtype=float)

# slopes
slope_intavg1_vecA = np.zeros((NT, nvtotA), dtype=float)
slope_intavg1_vecI = np.zeros((NT, nvtotI), dtype=float)
slope_intavg1_vecL = np.zeros((NT, nvtotL), dtype=float)
slope_intavg2_vecA = np.zeros((NT, nvtotA), dtype=float)
slope_intavg2_vecI = np.zeros((NT, nvtotI), dtype=float)
slope_intavg2_vecL = np.zeros((NT, nvtotL), dtype=float)

#------------------------------------------------------
# Print setup information to screen
#-------------------------------------------------------
print("\n")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~ ExoCAM trend analysis ~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("=== File series descriptors ===")
print(case_id, " ", START_YEAR, " ", root_atm)
if do_atm == False and do_ice == False and do_lnd == False:
    print("must choose a source model; cam, cice, clm")
    quit()
if do_atm == True:
    print(prefixA)
if do_ice == True:
    print(prefixI)
if do_lnd == True:
    print(prefixL)

if do_atm == True:
    print("atmosphere model variables")
    print(atmvars_in)
    print("atmosphere print variables")
    print(atmprint_in)
    print("atmosphere plot variables")
    print(atmplot_in)

if do_ice == True:
    print("ice model variables")
    print(icevars_in)
    print("ice print variables")
    print(iceprint_in)
    print("ice plot variables")
    print(iceplot_in)

if do_lnd == True:
    print("land model variables")
    print(lndvars_in)
    print("land print variables")
    print(lndprint_in)
    print("land plot variables")
    print(lndplot_in)

print("=== Resolution ===");
print("nlon  ", nlon)
print("nlat  ", nlat)
print("nlev  ", nlev)

#-----------------------------------------------------------------------------------------------------------------
# File scan: determine N_actual, populate time vectors and date strings
# No data is read here; this pass only checks file existence.
#-----------------------------------------------------------------------------------------------------------------
print("========================================")
print("=========  scanning file series ========")
print("========================================")
t0 = time.time()

# Pre-scan sanity check on the first month of each enabled component. Without
# it, a component that is entirely absent (e.g. --cice or --clm requested for a
# boundary condition that has no sea-ice or land model) makes the scan below
# break at i=0, yielding N_actual=0 and a downstream crash on the empty date
# range. cam.h0 files are always present in every configuration, so a missing
# atm file is a genuine path/naming error; cice/clm may legitimately not exist,
# in which case we tell the user to drop that component and stop cleanly.
_first_month = f"{START_YEAR:04d}-01"
_missing = []
if do_atm and not os.path.isfile(f"{root_atm}/{case_id}{prefixA}{_first_month}.nc"):
    _missing.append(("atmosphere (--cam)", f"{root_atm}/{case_id}{prefixA}{_first_month}.nc"))
if do_ice and not os.path.isfile(f"{root_ice}/{case_id}{prefixI}{_first_month}.nc"):
    _missing.append(("sea ice (--cice)", f"{root_ice}/{case_id}{prefixI}{_first_month}.nc"))
if do_lnd and not os.path.isfile(f"{root_lnd}/{case_id}{prefixL}{_first_month}.nc"):
    _missing.append(("land (--clm)", f"{root_lnd}/{case_id}{prefixL}{_first_month}.nc"))
if _missing:
    print("ERROR: requested component(s) have no data at the start year:")
    for label, path in _missing:
        print(f"  {label}: not found\n    {path}")
    print("Nothing was written. If this boundary condition has no such")
    print("component, rerun without that flag (cam data always exists; cice")
    print("and clm may legitimately be absent).")
    quit()

firstDate = None
lastDate  = None
i = 0
while True:

    if i == NT:
        lastDate = f"{year:04d}-{month}"
        break

    it       = i + 1
    yr_count = (it - 1) // 12
    month_i  = it - (yr_count * 12)
    year     = START_YEAR + yr_count
    month    = f"{month_i:02d}"

    file_atm = f"{root_atm}/{case_id}{prefixA}{year:04d}-{month}.nc"
    file_ice = f"{root_ice}/{case_id}{prefixI}{year:04d}-{month}.nc"
    file_lnd = f"{root_lnd}/{case_id}{prefixL}{year:04d}-{month}.nc"

    if do_atm == True:
        if not os.path.isfile(file_atm):
            lastDate = f"{year:04d}-{month}"
            break
        time_vecA[i] = it
    if do_ice == True:
        if not os.path.isfile(file_ice):
            lastDate = f"{year:04d}-{month}"
            break
        time_vecI[i] = it
    if do_lnd == True:
        if not os.path.isfile(file_lnd):
            lastDate = f"{year:04d}-{month}"
            break
        time_vecL[i] = it

    if i == 0:
        firstDate = f"{year:04d}-{month}"
        print("Date of first data read =", firstDate)

    i += 1

N_actual = i
timing['file scan'] = time.time() - t0
print("Scan complete. Timesteps found:", N_actual)
#print_timing('file scan', timing['file scan'])

if int1 >= int2:
    print(f"  WARNING: --int1 ({args.int1} yr) must be less than --int2 ({args.int2} yr). Setting int1 = 1 year.")
    int1 = 12
if int1 > N_actual:
    print(f"  WARNING: --int1 ({args.int1} years, {int1} months) exceeds available timesteps ({N_actual}). "
          f"Capping int1 at {N_actual} months.")
    int1 = N_actual
if int2 > N_actual:
    print(f"  WARNING: --int2 ({args.int2} years, {int2} months) exceeds available timesteps ({N_actual}). "
          f"Capping int2 at {N_actual} months "
          f"({N_actual//12} years, {N_actual%12} months remainder).")
    int2 = N_actual

#-----------------------------------------------------------------------------------------------------------------
# Read all monthly files and compute area-weighted global means.
#
# Assumptions to verify against actual model output:
#   - All years from 1 to START_YEAR-1 have exactly 12 monthly files each, so
#     the (START_YEAR-1)*12 offset in read_monthly_files correctly aligns with
#     the file scan above that determined N_actual.
#   - Variables are stored as (time, lat, lon); var[0,:,:] is correct for
#     single-snapshot monthly files.
#   - weights shape (nlat, nlon) matches the spatial dims of each variable.
#   - vavg_vec columns match variable index offsets in trend_utils
#     (atm: cols 4+, ice/lnd: cols 1+).
#-----------------------------------------------------------------------------------------------------------------
print("========================================")
print("=========  reading and averaging  ======")
print("========================================")

weights = core.build_area_weights(lon, lat)

if do_atm == True:
    t0 = time.time()
    gm_atm, profiles_atm, files_atm = core.read_monthly_fields(
        root_atm, case_id, prefixA, list(atmvars_in), START_YEAR, N_actual,
        weights, profile_vars)
    timing['read and average (atm)'] = time.time() - t0
    print(f"  found {len(files_atm)} files")
    #print_timing('read and average (atm)', timing['read and average (atm)'])
    for n in range(len(atmvars_in)):
        vavg_vecA[:N_actual, nv1dA + n] = gm_atm[:, n]

if do_ice == True:
    t0 = time.time()
    gm_ice, files_ice = core.read_monthly_files(root_ice, case_id, prefixI,
                                                 list(icevars_in), START_YEAR,
                                                 N_actual, weights)
    timing['read and average (ice)'] = time.time() - t0
    print(f"  found {len(files_ice)} files")
    #print_timing('read and average (ice)', timing['read and average (ice)'])
    for n in range(len(icevars_in)):
        vavg_vecI[:N_actual, nv1dI + n] = gm_ice[:, n]

if do_lnd == True:
    t0 = time.time()
    gm_lnd, files_lnd = core.read_monthly_files(root_lnd, case_id, prefixL,
                                                 list(lndvars_in), START_YEAR,
                                                 N_actual, weights)
    timing['read and average (lnd)'] = time.time() - t0
    print(f"  found {len(files_lnd)} files")
    #print_timing('read and average (lnd)', timing['read and average (lnd)'])
    for n in range(len(lndvars_in)):
        vavg_vecL[:N_actual, nv1dL + n] = gm_lnd[:, n]

#------------------------------------------------------
# Energy balance for all timesteps (if requested).
# Must be computed before running means so energy columns
# are populated when the stats loop reads vavg_vecA.
#------------------------------------------------------
if energy_requested and do_atm == True:
    for i in range(N_actual):
        etop, ebot = trend.atm_energy_calc(atmvars_in, vavg_vecA[i, :])
        vavg_vecA[i, nvtotA - 2] = etop
        vavg_vecA[i, nvtotA - 1] = ebot

#------------------------------------------------------
# Running means and slopes -- O(N) via cumsum, computed
# once per variable column over the full filled series.
#------------------------------------------------------
print("========================================")
print("=========  running means          ======")
print("========================================")
t0 = time.time()

if do_atm == True:
    for n in range(nv1dA, nvtotA):
        r1, r2, s1, s2 = core.compute_running_means(vavg_vecA[:N_actual, n], int1, int2)
        intavg1_vecA[:N_actual, n]       = r1
        intavg2_vecA[:N_actual, n]       = r2
        slope_intavg1_vecA[:N_actual, n] = s1
        slope_intavg2_vecA[:N_actual, n] = s2

if do_ice == True:
    for n in range(nv1dI, nvtotI):
        r1, r2, s1, s2 = core.compute_running_means(vavg_vecI[:N_actual, n], int1, int2)
        intavg1_vecI[:N_actual, n]       = r1
        intavg2_vecI[:N_actual, n]       = r2
        slope_intavg1_vecI[:N_actual, n] = s1
        slope_intavg2_vecI[:N_actual, n] = s2

if do_lnd == True:
    for n in range(nv1dL, nvtotL):
        r1, r2, s1, s2 = core.compute_running_means(vavg_vecL[:N_actual, n], int1, int2)
        intavg1_vecL[:N_actual, n]       = r1
        intavg2_vecL[:N_actual, n]       = r2
        slope_intavg1_vecL[:N_actual, n] = s1
        slope_intavg2_vecL[:N_actual, n] = s2

timing['running means'] = time.time() - t0
#print_timing('running means', timing['running means'])

#-----------------------------------------------------------------------------------------------------------------
# Post-processing loop: print to screen at requested interval.
# All global means and statistics are already computed; this loop only handles output.
#-----------------------------------------------------------------------------------------------------------------
print("========================================")
print("=========  starting output loop  =======")
print("========================================")
t0 = time.time()
firstPrintCall = True
for i in range(N_actual):
    if print_int is not None and (i + 1) % print_int == 0:
        trend.print2screen(atmvars_in, icevars_in, lndvars_in,
                           atmprint_in, iceprint_in, lndprint_in,
                           firstPrintCall, avgfreq,
                           do_atm, time_vecA[i], vavg_vecA[i,:],
                           intavg1_vecA[i,:], intavg2_vecA[i,:],
                           slope_intavg1_vecA[i,:], slope_intavg2_vecA[i,:],
                           do_ice, time_vecI[i], vavg_vecI[i,:],
                           intavg1_vecI[i,:], intavg2_vecI[i,:],
                           slope_intavg1_vecI[i,:], slope_intavg2_vecI[i,:],
                           do_lnd, time_vecL[i], vavg_vecL[i,:],
                           intavg1_vecL[i,:], intavg2_vecL[i,:],
                           slope_intavg1_vecL[i,:], slope_intavg2_vecL[i,:],
                           i)
        firstPrintCall = False

timing['print output'] = time.time() - t0
#print_timing('print output', timing['print output'])

#-----------------------------------------------------------------------------------------------------------------
# end output loop
#-----------------------------------------------------------------------------------------------------------------
print("Concluding             ", case_id)
print("Number of files read:  ", N_actual)
print("End date (year-month): ", lastDate)


#------------------------------------------------------
# Print data to text file
#-------------------------------------------------------
if args.save_data == True:
  print('Printing to text file...')
  trend.print2text(atmvars_in, lndvars_in, icevars_in, atmprint_in, lndprint_in, iceprint_in, \
                   do_atm, vnamesA, time_vecA, vavg_vecA, intavg1_vecA, intavg2_vecA, slope_intavg1_vecA, slope_intavg2_vecA, \
                   do_ice, vnamesI, time_vecI, vavg_vecI, intavg1_vecI, intavg2_vecI, slope_intavg1_vecI, slope_intavg2_vecI, \
                   do_lnd, vnamesL, time_vecL, vavg_vecL, intavg1_vecL, intavg2_vecL, slope_intavg1_vecL, slope_intavg2_vecL, \
                   firstDate, lastDate, case_id)
  if do_atm == True and profile_vars:
    trend.print_profiles2text(profiles_atm, firstDate, lastDate, case_id)

#------------------------------------------------------
# Call line plotting script
#-------------------------------------------------------
if args.plots == True:
  print('Plotting...')
  trend.timeSeriesPlots(atmvars_in, lndvars_in, icevars_in, atmplot_in, lndplot_in, iceplot_in, \
                        do_atm, time_vecA, vavg_vecA, intavg1_vecA, intavg2_vecA, slope_intavg1_vecA, slope_intavg2_vecA, \
                        do_ice, time_vecI, vavg_vecI, intavg1_vecI, intavg2_vecI, slope_intavg1_vecI, slope_intavg2_vecI, \
                        do_lnd, time_vecL, vavg_vecL, intavg1_vecL, intavg2_vecL, slope_intavg1_vecL, slope_intavg2_vecL, \
                        firstDate, lastDate, case_id, show=args.show)

trend.print_final_summary(atmvars_in, icevars_in, lndvars_in,
                          do_atm, do_ice, do_lnd,
                          vavg_vecA, intavg1_vecA, intavg2_vecA,
                          vavg_vecI, intavg1_vecI, intavg2_vecI,
                          vavg_vecL, intavg1_vecL, intavg2_vecL,
                          N_actual, case_id, firstDate, lastDate,
                          args.int1, args.int2)

if args.timing:
    print("========================================")
    print("=========  timing summary         ======")
    print("========================================")
    total = sum(timing.values())
    for label, elapsed in timing.items():
        print_timing(label, elapsed)
    print_timing('total', total)

sys.exit()
