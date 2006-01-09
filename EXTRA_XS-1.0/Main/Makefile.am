SHERPASYS     =	../..

EXTRAXSFLAGS  =	# -pedantic -Wall -Winline
DEFS          =

EXTRAXSINCS   =	-I../Two2Two -I../BFKL
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES = libExtraXS.la
libExtraXS_la_SOURCES =	\
			XS_Base.C \
			Single_XS.C \
			XS_Group.C \
			Simple_XS.C \
			XS_Selector.C
pkginclude_HEADERS        =	\
			XS_Base.H \
			Single_XS.H \
			XS_Group.H \
			Simple_XS.H \
			XS_Selector.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(EXTRAXSFLAGS)
INCLUDES      =	$(SHERPAINCS) $(EXTRAXSINCS) $(PHASICINCS) $(PDFINCS) \
		$(BEAMINCS) $(MODELINCS) $(ATOOLSINCS)
