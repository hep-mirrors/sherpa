SHERPASYS     =	../..

PDFFLAGS      =	# -pedantic -Wall -Winline -Woverloaded-virtual

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if LHAPDF_SUPPORT
LHAPDFFLAGS      = -DLHAPDF_SUPPORT
lib_LTLIBRARIES = libLHAPDF.la 
libLHAPDF_la_SOURCES  =	\
			LHAPDF_Fortran_Interface.C \
			lhapdf_inter.f 
noinst_HEADERS        =	LHAPDF_Fortran_Interface.H
endif
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS) $(LHAPDFFLAGS)
INCLUDES      =	$(PDFINCS) $(MODELINCS) $(ATOOLSINCS)
