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



lib_LTLIBRARIES = libAdicicInput.la
libAdicicInput_la_SOURCES = Dipole_Flavour.C Information.cpp \
                            Recoil_Strategy.cpp Dipole_Parameter.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA
AM_CXXFLAGS = -pedantic -Wall -Winline



INCLUDES = -I../Tool \
           -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org \
           -I$(MODELINC)/Main


  
noinst_HEADERS = Dipole_Flavour.H Information.hpp \
                 Sudakov_Strategy.hpp Recoil_Strategy.hpp \
                 Evolution_Strategy.hpp \
                 Dipole_Parameter.H



EXTRA_DIST = $(CVSINCLUDE) \
             Dipole_Parameter.inl.hh Dipole_Parameter.prm.cc



