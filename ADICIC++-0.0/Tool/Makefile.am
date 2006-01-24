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
libAdicicTool_la_SOURCES = Xio.C Xmath.C Xhisto.C \
                           Sudakov_Calculator.C \
                           Sudakov_Group.C IISudakov_Group.C \
                           Recoil_Calculator.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA
AM_CXXFLAGS = -pedantic -Wall -Winline -ftemplate-depth-100



INCLUDES = -I../Input \
           -I../Handler \
           -I../Structure \
           -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org \
           -I$(MODELINC)/Main \
           -I$(PDFINC)/Main


  
include_HEADERS = Xenum.H Xio.H Xmath.H Xhisto.H defextra enumextra ioextra mathextra histoextra \
                  Sudakov_Calculator.H \
                  Sudakov_Group.H IISudakov_Group.H \
                  Recoil_Calculator.H



EXTRA_DIST = $(CVSINCLUDE) \
             Sudakov_Calculator.inl.hh \
             Sudakov_Group.inl.hh Sudakov_Group.tpt.cc \
             IISudakov_Group.inl.hh IISudakov_Group.tpt.cc \
             Recoil_Calculator.inl.hh Recoil_Calculator.tpt.cc


