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



lib_LTLIBRARIES = libAdicicStruct.la
libAdicicStruct_la_SOURCES = Dipole_Flavour.C Information.cpp \
                             Dipole.C Chain.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA



INCLUDES = -I../Tool \
           -I../Handler \
           -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org


  
noinst_HEADERS = Dipole_Flavour.H Information.hpp \
                 Dipole.H Chain.H



EXTRA_DIST = $(CVSINCLUDE) \
             Dipole.inl.hh Dipole.exa.cc Chain.inl.hh






