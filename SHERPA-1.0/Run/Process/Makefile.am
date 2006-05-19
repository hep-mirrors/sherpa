SUBDIRS = Dummy

pkglib_LTLIBRARIES = libProcess.la

libProcess_la_SOURCES =	\
			Set_Channel.C \
			Set_Values.C 

EXTRA_DIST    =	@SVNINCLUDE@ Set_Channel.C.clean Set_Values.C.clean

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      = @AMEGICINCS@ @PHASICINCS@ @MODELINCS@ @PDFINCS@ \
		@BEAMINCS@ @ATOOLSINCS@
