
#!/bin/bash


# bash script to set up a CLM5.0 run with transient land use
# to test  dynamic lakes



#=====================================================
# initialisation
#=====================================================


# define path where cases are stored
CASEDIR=~/cases


# define path where CESM scripts are stored
SCRIPTSDIR=~/clm5.0/cime/scripts


# define machine
MACH=cheyenne


# define compiler
COMPILER=intel;  MPILIB=openmpi      # io


# define compset - more information below
COMP=IHistClm50Sp             #HIST_DATM%GSWP3v1_CLM50%SP_SICE_SOCN_MOSART_CISM2%NOEVOLVE_SWAV

# define resolution
RES=hcru_hcru    #360x720cru       #  Exact half-degree CRUNCEP datm forcing grid with CRUNCEP land-mas

# define other information
# in this run the lake extent of 1900 is used (no reservoirs)
type=CTL

# define run settings
start_year=1966  # start year of the simulation (1966)
end_year=2015    # end year of GSIM data set (obs)
simul_length=1   # 1 year (total 50 years)
comp=$( echo ${COMP:0:1} | tr '[:upper:]' '[:lower:]' )

# define namelist script
nl_file=nl_clm_${type}.sh

# set whether you are in production mode or test mode
production=true  # true=final production runs; false=testing
nresubmit=49    # if simul_length=10, this is number of years -1; use for final production runs including spinup


#=====================================================
# 1. create new case
#=====================================================


# get path where setup scripts are
SETUPDIR=$(pwd)


# Change into the new case directory
cd $SCRIPTSDIR


# generate casename
CASE=$comp.$COMP.$RES.$type.2


# create a new case in the directory 'cases'
./create_newcase --case $CASEDIR/$CASE  --res $RES --compset $COMP --mach $MACH --project P93300041 --run-unsupported


#=====================================================
# 2. invoke cesm_setup
#=====================================================


# Change into the new case directory
cd $CASEDIR/$CASE


# copy this script to the case directory for reference
cp $SETUPDIR/`basename "$0"` .


# modify env_run.xml
./xmlchange PROJECT=P93300041
./xmlchange JOB_QUEUE="economy"

# modify env_run.xml - MAY BE CHANGED ANYTIME during a run
./xmlchange --file env_run.xml --id STOP_OPTION       --val nyears
./xmlchange --file env_run.xml --id STOP_N            --val ${simul_length}


# modify env_run.xml - don in run_clm_historical script
./xmlchange --file env_run.xml --id RUN_TYPE          --val startup                # Set to run type to startup (is the default)
./xmlchange --file env_run.xml --id RUN_STARTDATE     --val ${start_year}-01-01    # year in which CESM starts running (1st January)


# introduce changes to CLM namelist
cp $SETUPDIR/$nl_file .
. ./$nl_file

# Configure case
./case.setup


#=====================================================
# 3. Build and Submit run to the batch queue
#====================================================


if [ "$production" = true ]; then

   # copy ncar script and necessary datm namelist files into casedirectory
    # build the model
  ./case.build
  ./xmlchange -file env_run.xml -id RESUBMIT -val ${nresubmit}
  ./case.submit


else                               # single year run

  # build the model
  ./case.build

  # change options for test run
  ./xmlchange STOP_N=1
  ./xmlchange STOP_OPTION=nyears
  ./xmlchange --file env_batch.xml --id JOB_WALLCLOCK_TIME --val 04:00:00
  ./xmlchange -file env_run.xml -id DOUT_S -val False
  ./case.submit

fi


