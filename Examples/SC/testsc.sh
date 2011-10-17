#!/bin/bash -x

SHERPAPATH=$PWD/../../bin
EVENTS=100000

rivetplotstr=""
for process in z_tautau w_taunu wp_taunu h_tautau hm_tautau; do
    cd $process;
    for gen in Amegic Comix; do
        for sc in 1 0; do
            if ! test -d Analysis${gen}SC${sc}/; then
                if ! ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc; then
                    ./makelibs -j3;
                    ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc || exit 1;
                fi
            fi
            ../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis${gen}SC${sc} || exit 1
            rivetplotstr="$rivetplotstr ${process}/Analysis${gen}SC${sc}.aida:Title=${gen}SSC${sc}:LineWidth=0.1pt"
        done
    done
    cd -;

    cd ${process}_hard;
    for gen in Amegic Comix; do
        for sc in 1 0; do
            if ! test -d Analysis${gen}SC${sc}/; then
                if ! ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc; then
                    ./makelibs -j3;
                    ${SHERPAPATH}/Sherpa EVENTS=${EVENTS} GENTAG:=$gen SPINCORRTAG:=$sc || exit 1;
                fi
            fi
            ../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis${gen}SC${sc} || exit 1
            rivetplotstr="$rivetplotstr ${process}_hard/Analysis${gen}SC${sc}.aida:Title=${gen}HSC${sc}:LineWidth=0.1pt"
        done
    done
    cd -;

    cd ${process}_ME;
    if ! test -d Analysis/; then
        ${SHERPAPATH}/Sherpa EVENTS=${EVENTS};
    fi;
    ../../../AddOns/Rivet/sherpa2aida -a MC_$process Analysis || exit 1
    rivetplotstr="$rivetplotstr ${process}_ME/Analysis.aida:Title=ComixME:LineWidth=0.1pt"
    cd -;
done

rivet-mkhtml -s -c make-plots.conf $rivetplotstr
