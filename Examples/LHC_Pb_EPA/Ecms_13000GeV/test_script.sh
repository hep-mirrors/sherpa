#!/bin/bash

###------------------------###
### Sherpa EPA test script ###
###------------------------###

## Description:
################################################################################
# This script is used to ensure, that the EPA code keeps working while it is
# under development. Furthermore it is used to produce plots which can be 
# compared to to the review "Photon and gluon induced processes in relativistic
# heavy ion collisions." by F. Krauss, M. Greiner, G. Soff published in
# Prog.Part.Nucl.Phys.39:503-564,1997 

## Dependencies: 
################################################################################
# SHERPA 1.2.3 + EPA
# gnuplot
# EPA.plt (gnuplot plot file)

## Variables
################################################################################
# where is SHERPA?
SHERPA_PATH=/home/felix/Research/EPA/Development/svn_branch

# where is the SHERPA binary
EXE=${SHERPA_PATH}/bin/Sherpa

# where are the examples
SAMPLES=${SHERPA_PATH}/Examples/LHC_Pb_EPA/Ecms_13000GeV

# where should the output files be placed
OUTPUT_DIR=${SHERPA_PATH}/Examples/LHC_Pb_EPA/Ecms_13000GeV/Logs

# after how many seconds a SHERPA execution should be stopped
KILLSECS=90

# svn revision (including indication if modified)
SVN_REV=`svnversion`

# date/time
DATE=`date +%Y-%m-%d_%H-%M-%S`

# working directory
WRK_DIR=${OUTPUT_DIR}/${DATE}_${SVN_REV}

# logfile
LOGFILE=${WRK_DIR}/log_${DATE}_${SVN_REV}.log

# svn diff
SVN_DIFF=${WRK_DIR}/svn_diff_${DATE}_${SVN_REV}.patch

# svn status
SVN_STATUS=${WRK_DIR}/svn_status_${DATE}_${SVN_REV}.txt

## Functions
################################################################################ 
controlRuntime() { 
  progID=$1   # first parameter
  killtime=$2 # second parameter
  
  mypid=`eval ps ax|grep "$progID"|grep -iv "grep"| awk '{print $1}'`
  echo "mypid $mypid"
  mins=0
  secs=0
  startsecs=`eval date +%-S `
  starthours=`eval date +%-H`
  startmins=`eval date +%-M`

  #echo "startsecs=$startsecs starthours=$starthours startmins=$startmins"

  while [[ $secs -lt $killtime &&  -d /proc/$mypid ]]; do
    cursecs=`eval date +%-S `
    curhours=`eval date +%-H`
    curmins=`eval date +%-M`
    secs=$(( (curhours-starthours)*3600 + (curmins-startmins)*60 - startsecs + cursecs ))
    echo "Secs: $secs"
    sleep 1
  done
  if [ -d /proc/$mypid ]; then
    kill -9 "$mypid" 
    if [ $? -eq 0 ]; then
      echo "Process $mypid was killed"
    fi
  fi
}

run() {
  echo "Running SHERPA with $1" >> $LOGFILE

  # copy Run.dat to the working directory
  cp ${SAMPLES}/$1 ${WRK_DIR}

  # run sherpa
  ${EXE} -f ${WRK_DIR}/$1 -m Internal 2>&1 >> $LOGFILE &

  # control the runtime
  controlRuntime Sherpa $KILLSECS
}

 
## Execution
################################################################################
# create and enter working folder
mkdir ${WRK_DIR}
cd ${WRK_DIR}

# create logfile
touch $LOGFILE

# save code status
cd ${SHERPA_PATH}
svn diff > $SVN_DIFF
svn status > $SVN_STATUS
cd ${WRK_DIR}

# run SHERPA for EPA -> ...
###########################
# Pions
run Run.dat.Pions
run Run.dat.PionsHCSGauss 
run Run.dat.PionsSmoothHCS

# Kaons
#run Run.dat.Kaons

# D
#run Run.dat.D

# Muons
#run Run.dat.Muons

# Tauons
#run Run.dat.Tauons

# gnuplot 
#########
# EPA spectra
echo "gnuplot EPA spectra" >> $LOGFILE
cp ${SAMPLES}/EPA.plt ${WRK_DIR}
gnuplot EPA.plt 2>&1 >> $LOGFILE
ps2pdf EPA_spectrum.ps EPA_spectrum.pdf
ps2pdf EPA_spectrum.ps EPA_spectrum2.pdf

# matrix element for pi^+ pi^-

# matrix element for K^+ K^-

# matrix element for mu^+ mu^-

# archive
#########
cd ..
tar cjfv ${WRK_DIR}.tar.bz2 ${WRK_DIR}/*

# EOF
