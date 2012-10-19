#!/bin/bash

cd HAWK

make

if ! test -d ../lib; then mkdir ../lib; fi
ar cr ../lib/libhawk.a *.o
