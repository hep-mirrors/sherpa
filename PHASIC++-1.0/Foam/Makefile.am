SHERPASYS     =	../..

PHASICFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
#PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
#BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
#MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libFoam.la
libFoam_la_SOURCES    =	\
			TPSEMAR.cxx \
			TFOAM_INTEGRAND.cxx \
			TFMAXWT.cxx \
			TFHST.cxx \
			TFOAM.cxx
noinst_HEADERS        =	\
			TPSEMAR.h \
			TFOAM_INTEGRAND.h \
			TFMAXWT.h \
			TFHST.h \
			TFOAM.h
AM_CXXFLAGS   =	$(PHASICFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(ATOOLSINCS) $(ROOTINCS)
EXTRA_DIST    = $(CVSINCLUDE)
