pkglib_LTLIBRARIES       =	libLundTools.la

libLundTools_la_SOURCES       =	\
				pythia6214.f \
				pydata.f \
				pyalps.f \
				pypdfu.f \
				pyplst.f \
				pymult.f \
				Lund_Interface.C 

pkginclude_HEADERS        =	\
		 	Lund_Wrapper.H \
			Lund_Interface.H 

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@MODELINCS@ @ATOOLSINCS@ @PDFINCS@ @SHERPAINCS@
