#!/bin/bash

version=5.8

tarname=mcfm-${version}.tar.gz
dirname=MCFM-${version}

echo "installing MCFM from ${tarname} to ${dirname}"

if ! test -f $tarname; then
  wget http://mcfm.fnal.gov/$tarname
fi

if ! test -d $dirname; then
  tar -xzf $tarname
  mv MCFM $dirname
  cd $dirname
  mkdir obj
  sed -e's/\/Users\/johnmc\/MCFM/'$(pwd | sed -e's/\//\\\//g')'/g' \
      -e's/\(FFLAGS.*=.*\)-fno-f2c/\1-fPIC -DPIC/g' -i makefile
  sed -e's/\(.*call pdfwrap\)/c\1/g' -e's/\(.*nlooprun=0\)/c\1/g' -i src/Need/*.f
  sed -e's/\(.*[ \t]stop\)/c\1/g' -i src/Need/chooser.f
  cd -
fi

cd $dirname

make

if ! test -d ../lib; then mkdir ../lib; fi
ar cr ../lib/libMCFM.a */*.o
