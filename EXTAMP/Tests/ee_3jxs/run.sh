#!/bin/bash

# Expect this to converge not too well due to
# missing alpha parameter implementation

../../../bin/Sherpa -f Run_amegic.dat
./makelibs -j 14
cp Process/Amegic/lib/libProc_fsrchannels4.so ./

mpiexec -n 4 ../../../bin/Sherpa -f Run_ol.dat -g
