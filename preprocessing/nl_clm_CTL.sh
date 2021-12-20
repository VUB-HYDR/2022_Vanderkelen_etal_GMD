#! /bin/bash
#Author: Inne Vanderkelen

#use this file to change the settings for the run
#i.e. source (. ./) it in the simulation-creation file


cat > user_nl_clm << EOF

check_dynpft_consistency = .false.
do_transient_crops=.true.
do_transient_pfts =.true.
do_transient_lakes =.true.

hist_fincl2 = 'QRUNOFF', 'QRGWL', 'QRUNOFF_TO_COUPLER','QIRRIG_DEMAND','QIRRIG_FROM_SURFACE','QIRRIG_DRIP','EFLX_LH_TOT', 'RAIN', 'RAIN_FROM_ATM', 'SNOW'
hist_fincl3 = 'EFLX_LH_TOT'

! h0 with monthly output, h1 daily and h2 subgrid monthly output.
hist_nhtfrq = 0, -24, -24
hist_mfilt = 12, 365, 365
hist_dov2xy = .true., .true., .false.

! new surfdat and timeseries
flanduse_timeseries='/glade/work/ivanderk/inputdata/clm_input/CLM_for_mizuRoute/landuse.timeseries_360x720cru_hist_16pfts_dynlakes_Irrig_simyr1900-2015_c200217.nc'

fsurdat = '/glade/work/ivanderk/inputdata/clm_input/CLM_for_mizuRoute/surfdata_360x720cru_hist_16pfts_Irrig_CMIP6_simyr1900_c210331.nc'
!fatmlndfrc = '/glade/work/ivanderk/inputdata/clm_input/CLM_for_mizuRoute/domain.lnd.360x720_cruncep_largelakes.210426.nc'

EOF


cat > user_nl_mosart << EOF

! yearly output of monthly timesteps
rtmhist_nhtfrq = 0
rtmhist_mfilt = 12

EOF

# notes
