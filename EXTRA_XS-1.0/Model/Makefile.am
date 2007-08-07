pkglib_LTLIBRARIES = libExtraXSModel.la

libExtraXSModel_la_SOURCES    =	\
				XS_Model_Base.C \
				XS_Model_QCD.C \
				XS_Model_EW.C \
				XS_Model_SM.C \
				XS_Model_Handler.C

pkginclude_HEADERS    =	\
			XS_Model_Base.H \
			XS_Model_QCD.H \
			XS_Model_EW.H \
			XS_Model_SM.H \
			XS_Model_Handler.C

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
