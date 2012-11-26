#!/bin/bash

if test -z "$2"; then
  SHERPAPATH=$PWD/../../../bin
else
  if test "$2" = "clean"; then
    rm -rf "${1}/Process" "${1}/Results" "${1}/*.log" "${1}/errors" "${1}/Analysis*"
    rm -rf "${1}_hard/Process" "${1}_hard/Results" "${1}_hard/*.log" "${1}_hard/errors" "${1}_hard/Analysis*"
    rm -rf "${1}_ME/Process" "${1}_ME/Results" "${1}_ME/*.log" "${1}_ME/errors" "${1}_ME/Analysis*"
    exit 0;
  else
    SHERPAPATH=$2
  fi
fi

if test -z "$1"; then
  PROCESSES="z_tautau w_taunu wp_taunu h_tautau hm_tautau z_toptop ww z_llvvbb"
else
  PROCESSES=$1
fi

EVENTS=50000

rivetplotstr=""
for process in $PROCESSES; do
    echo "Testing $process";
    if test -d $process; then cd $process && for gen in Amegic Comix; do
        echo "  with production in $gen and tau decays in HADRONS";
        for sc in 1 0; do
            echo "    for SPINCORRTAG:=$sc";
            if ! test -d Analysis${gen}SC${sc}/; then
                rm -rf Process/ Results/Comix Results/Amegic;
                if ! ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc LOG_FILE=out.$gen.$sc.log >> errors 2>&1; then
                    ./makelibs -j3 > /dev/null 2>&1;
                    ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc LOG_FILE=out.$gen.$sc.log >> errors 2>&1 || exit 1;
                fi
            fi
            ../../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis${gen}SC${sc} > /dev/null || exit 1
            rivetplotstr="$rivetplotstr ${process}/Analysis${gen}SC${sc}.aida:Title=${gen}SSC${sc}:LineWidth=0.1pt"
        done
    done && cd -; fi

    if test -d ${process}_hard; then cd ${process}_hard && for gen in Amegic Comix; do
        echo "  with production in $gen and hard tau decays";
        for sc in 1 0; do
            echo "    for SPINCORRTAG:=$sc";
            if ! test -d Analysis${gen}SC${sc}/; then
                rm -rf Process/ Results/Comix Results/Amegic;
                if ! ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc LOG_FILE=out.$gen.$sc.log >> errors 2>&1; then
                    ./makelibs -j3 > /dev/null 2>&1;
                    ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc LOG_FILE=out.$gen.$sc.log >> errors 2>&1 || exit 1;
                fi
            fi
            ../../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis${gen}SC${sc} > /dev/null || exit 1
            rivetplotstr="$rivetplotstr ${process}_hard/Analysis${gen}SC${sc}.aida:Title=${gen}HSC${sc}:LineWidth=0.1pt"
        done
    done && cd -; fi

    cd ${process}_ME;
    echo "  with production and decay in ME";
    if ! test -d Analysis/; then
        rm -rf Process/ Results/Comix Results/Amegic;
        ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} LOG_FILE=out.log >> errors 2>&1;
    fi;
    ../../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis > /dev/null || exit 1
    rivetplotstr="$rivetplotstr ${process}_ME/Analysis.aida:Title=ComixME:LineWidth=0.1pt"
    cd -;
done

rivet-mkhtml --ignore-missing -s -c make-plots.conf $rivetplotstr
