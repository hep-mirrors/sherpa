pkglib_LTLIBRARIES = libExtraXS.la

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

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@EXTRAXSINCS@ @PHASICINCS@ @PDFINCS@ \
		@BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
