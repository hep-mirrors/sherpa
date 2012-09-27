#!/bin/bash

version=6
subversion=3

tarname=mcfm-${version}.${subversion}.tar.gz
dirname=MCFM-${version}.${subversion}

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
  sed -e's/\/scratch\/ellis\/play\/MCFMdevel/'$(pwd | sed -e's/\//\\\//g')'/g' \
      -e's/\(FFLAGS.*=.*\)-fno-f2c/\1-fPIC -DPIC/g' -i makefile
  sed -e's/\(.*call pdfwrap\)/c\1/g' -e's/\(.*nlooprun=0\)/c\1/g' -i src/Procdep/*.f
  sed -e's/\(.*[ \t]stop\)/c\1/g' -i src/Procdep/chooser.f
  sed -e's/epinv\*\*2/epinv2/g' -i src/*/*.f
  sed -e"/      include 'epinv2.f'/d" -i src/*/*.f
  sed -e"/      INCLUDE 'epinv2.f'/d" -i src/*/*.f
  sed -e"/^      include 'epinv.f'/a\      include \'epinv2.f\'" -i src/*/*.f
  sed -e"/^      INCLUDE 'epinv.f'/a\      include \'epinv2.f\'" -i src/*/*.f
  if [ $version -ge 6 ]; then
    if [ $subversion -ge 3 ]; then
      echo "patching specifics for MCFM-${version}.${subversion}"
      sed -e"/^      common\/spira\/spira/a\c---SHERPA: short-cut, we don't care about MCFM's BRs\n      return\nc---SHERPA: end short-cut" -i src/Need/sethparams.f
      sed -e"/^C---end statement functions/a\c---SHERPA: short-cut, h->gamgam ratio will be reweighted\n      msqgamgam=1d0\n      return\nc---SHERPA: end short-cut" -i src/ggHgaga/msqgamgam.f
    fi
  fi
  cd -
fi

cd $dirname

if test -d QCDLoop; then
  cd QCDLoop
  make
  cd -
fi

make

if ! test -d ../lib; then mkdir ../lib; fi
ar cr ../lib/libMCFM.a */*.o
