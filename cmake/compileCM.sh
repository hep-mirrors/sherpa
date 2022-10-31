#!/bin/bash
cmake -S . -B BUILD -DCMAKE_INSTALL_PREFIX=$(pwd)/CM -DSHERPA-MC_PYTHON_VERSION=3.10 -DCMAKE_C_FLAGS="-O2" -DCMAKE_Fortran_FLAGS="-O2" -DCMAKE_CXX_FLAGS="-O2 -Wl,--no-as-needed -fcx-fortran-rules"
bear --output CM.json -- cmake --build BUILD -j 10
cmake --install BUILD
rm -rf CMldd.txt
for a in $(ls -1 CM/lib64/SHERPA-MC/*); do
echo $(basename $a) >> CMldd.txt;
ldd $a | cut -f 1 -d\( | sort >> CMldd.txt;
done
