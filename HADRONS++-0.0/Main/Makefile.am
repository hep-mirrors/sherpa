pkglib_LTLIBRARIES = libHadronsMain.la

libHadronsMain_la_SOURCES = \
		Hadrons.C \
		Decay_Map.C \
		Hadron_Decay_Table.C \
		Hadron_Decay_Channel.C \
		Mixing_Handler.C \
		Tools.C

pkginclude_HEADERS = \
		Hadrons.H \
		Decay_Map.H \
		Hadron_Decay_Table.H \
		Hadron_Decay_Channel.H \
		Mixing_Handler.H \
		Tools.H

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @HADRONSINCS@ @PHASICINCS@ @ATOOLSINCS@ @MODELINCS@ @HELICITIESINCS@

EXTRA_DIST    = @SVNINCLUDE@
