pkglib_LTLIBRARIES = libApacicMain.la

libApacicMain_la_SOURCES      = \
				Knot.C \
				Tree.C \
				Apacic.C

pkginclude_HEADERS       =	\
			Knot.H \
			Tree.H \
			Veto_Info.H \
			Apacic.H

AM_CXXFLAGS   = @SHERPAFLAGS@

INCLUDES = @APACICINCS@ @PDFINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
