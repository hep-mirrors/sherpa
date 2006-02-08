SHERPASYS     =	../..

PDFFLAGS      =	# -pedantic -Wall -Winline
DEFS          =

PDFINCS       =	-I../Main -I../Sudakov
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
APACICINCS    =	-I$(APACICDIR)/Main -I$(APACICDIR)/Showers

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

if ROOT_SUPPORT
ROOTINCS      =	-I@ROOTINCDIR@
ROOTFLAGS     = -Wno-long-long
endif

pkglib_LTLIBRARIES       = libDUPDF.la
libDUPDF_la_SOURCES   =	\
			LL_Branching.C \
			LL_Single_Sudakov.C \
			LL_Sudakov.C \
			Doubly_Unintegrated_PDF.C
pkginclude_HEADERS        = \
			LL_Branching.H \
			LL_Single_Sudakov.H \
			LL_Sudakov.H \
			Doubly_Unintegrated_PDF.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(APACICINCS) $(PDFINCS) $(MODELINCS) $(ATOOLSINCS) \
		$(ROOTINCS)
