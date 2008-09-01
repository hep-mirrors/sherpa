pkglib_LTLIBRARIES       =	libLundTools.la

if PYTHIA_SUPPORT
PYTHIASOURCES = pythia6214.f pydata.f pyalps.f pypdfu.f pyplst.f pymult.f \
		Lund_Interface.C
PYTHIAHEADERS = Lund_Wrapper.H Lund_Interface.H
endif


libLundTools_la_SOURCES       =	$(PYTHIASOURCES)

pkginclude_HEADERS        =	$(PYTHIAHEADERS) 

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@MODELINCS@ @ATOOLSINCS@ @PDFINCS@ @SHERPAINCS@
