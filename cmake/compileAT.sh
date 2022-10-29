#!/bin/bash
export PYTHON=/usr/bin/python
export PYTHON_VERSION=3.10
autoreconf -fi

#export PYTHON=/usr/bin/python3
#export PYTHON_VERSION=3.10
./configure  --prefix=$(pwd)/AT  --disable-rpath  --enable-rivet=/usr  --enable-cernlib=/usr/lib64/cernlib/2006  \
--enable-blackhat=$(blackhat-config --prefix)     --enable-gzip --enable-recola=/usr  --enable-hztool=/usr   \
--enable-hepevtsize=4000    --enable-hepmc3=/usr --disable-hepmc3root   --libdir=$(pwd)/AT/lib64  \
  --enable-openloops=/usr/lib64/openloops --enable-hepmc2=/usr  --enable-root=/usr  --enable-binreloc  \
 --enable-pythia --enable-lhole --enable-lhapdf=/usr --enable-manual --enable-ewsud --enable-gosam=/usr #--enable-dihiggs 
 #--enable-analysis
 #--enable-analysis
bear --output AT.json -- make -j 10
make install
rm -f  AT/lib64/SHERPA-MC/*la
rm -f  AT/lib64/SHERPA-MC/*so.0.0
rm -f  AT/lib64/SHERPA-MC/*so.0
rm -rf ATldd.txt
for a in $(ls -1 AT/lib64/SHERPA-MC/*); do
echo $(basename $a) >> ATldd.txt;
ldd $a | cut -f 1 -d\( | sort >> ATldd.txt;
done
