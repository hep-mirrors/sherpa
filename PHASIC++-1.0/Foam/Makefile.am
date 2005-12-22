SHERPASYS     =	../..

PHASICFLAGS   =	# -pedantic -Wall #-Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
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
			ROOT_DEF.h \
			TPSEMAR.h \
			TFOAM_INTEGRAND.h \
			TFMAXWT.h \
			TFHST.h \
			TFOAM.h
AM_CXXFLAGS   =	$(PHASICFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(ATOOLSINCS) $(ROOTINCS)
EXTRA_DIST    = $(CVSINCLUDE)
