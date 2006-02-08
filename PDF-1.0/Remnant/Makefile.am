SHERPASYS     =	../..

PDFFLAGS      =	# -pedantic -Wall -Winline
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
BEAMINCS      = -I$(SHERPASYS)/BEAM-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS       =	-I../LHAPDF -I../MRST -I../CTEQ -I../GRV -I../Sudakov \
		-I../KMR -I../Main
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Showers

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

pkglib_LTLIBRARIES       =	libRemnant.la
libRemnant_la_SOURCES =	\
			Remnant_Base.C \
			No_Remnant.C \
			Color_Dipole.C \
			QCD_Remnant_Base.C \
			Hadron_Remnant.C \
			Electron_Remnant.C \
			Photon_Remnant.C 
pkginclude_HEADERS        =	\
			Remnant_Base.H \
			No_Remnant.H \
			Color_Dipole.H \
			QCD_Remnant_Base.H \
			Hadron_Remnant.H \
			Electron_Remnant.H \
			Photon_Remnant.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(PDFFLAGS) $(LHAPDFFLAGS)
INCLUDES      =	$(SHERPAINCS) $(APACICINCS) $(PDFINCS) \
		$(MODELINCS) $(BEAMINCS) $(ATOOLSINCS)

