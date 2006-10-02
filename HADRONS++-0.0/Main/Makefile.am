pkglib_LTLIBRARIES = libHadronsMain.la

libHadronsMain_la_SOURCES = \
		Hadrons.C \
		Decay_Table_Reader.C \
		Hadron_Decay_Channel.C 

pkginclude_HEADERS = \
		Hadrons.H \
		Decay_Table_Reader.H \
		Hadron_Decay_Channel.H 

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @HADRONSINCS@ @PHASICINCS@ @ATOOLSINCS@ @MODELINCS@

EXTRA_DIST    = @SVNINCLUDE@
