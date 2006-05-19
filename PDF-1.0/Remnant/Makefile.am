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

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@ $(LHAPDFFLAGS)

INCLUDES      =	@SHERPAINCS@ @APACICINCS@ @PDFINCS@ @CTEQINCS@ \
		@MODELINCS@ @BEAMINCS@ @ATOOLSINCS@
