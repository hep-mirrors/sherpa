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
libAdicicTool_la_SOURCES = Xio.C



#AM_CXXFLAGS = -pedantic -Wall -Winline -DUSING_SHERPA



INCLUDES = -I$(ATOOLSINC)/Phys \
           -I$(ATOOLSINC)/Math \
           -I$(ATOOLSINC)/Org


  
noinst_HEADERS = Xenum.H Xio.H defextra enumextra ioextra



EXTRA_DIST = $(CVSINCLUDE)








