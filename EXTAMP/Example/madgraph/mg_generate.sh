#!/bin/bash

set -e

$BUILDS/MG5_aMC_v2_5_1/bin/mg5_aMC mg_proc_card.dat

cp PROC*/SubProcesses/P1_*/CPPProcess.* ./
cp PROC*/src/* ./

#rm -r PROC*/
