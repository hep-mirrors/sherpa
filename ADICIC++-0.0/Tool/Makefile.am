MAKE    = make
SUBDIRS = savemaketool

SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0
SHERPAINC  = $(SHERPASYS)/SHERPA-1.0



if CVS_INCLUDE
CVSINCLUDE = CVS
endif



lib_LTLIBRARIES = libAdicicTool.la
libAdicicTool_la_SOURCES = Xio.C \
                           Sudakov_Calculator.C Recoil_Calculator.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA
AM_CXXFLAGS = -pedantic -Wall -Winline



INCLUDES = -I../Handler \
           -I../Structure \
           -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org \
           -I$(MODELINC)/Main


  
noinst_HEADERS = Xenum.H Xio.H Xmath.H defextra enumextra ioextra mathextra \
                 Sudakov_Calculator.H Recoil_Calculator.H



EXTRA_DIST = $(CVSINCLUDE) \
             Sudakov_Calculator.dat.hh Sudakov_Calculator.inl.hh \
             Sudakov_Calculator.dat.cc Sudakov_Calculator.tpt.cc \
             Recoil_Calculator.inl.hh Recoil_Calculator.tpt.cc









