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



lib_LTLIBRARIES = libAdicicHdl.la
libAdicicHdl_la_SOURCES = Dipole_Handler.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA



INCLUDES = -I../Tool \
           -I../Structure \
           -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org


  
noinst_HEADERS = Dipole_Handler.H



EXTRA_DIST = $(CVSINCLUDE) \
             Dipole_Handler.inl.hh Dipole_Handler.dat.cc








