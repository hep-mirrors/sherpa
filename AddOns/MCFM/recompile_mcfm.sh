#!/bin/bash

cd MCFM-6.1

make

if ! test -d ../lib; then mkdir ../lib; fi
ar cr ../lib/libMCFM.a */*.o
