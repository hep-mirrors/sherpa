
PGS_SOURCES   =	\
		pgs_finterface.f \
		PGS_Interface.C

AM_FFLAGS     =	-fno-automatic

if PGS_SUPPORT
pkglib_LTLIBRARIES    =	libPGSSherpa.la

libPGSSherpa_la_SOURCES       =	$(PGS_SOURCES)
else
PGS_EXTRA_DIST    = $(PGS_SOURCES)
endif

EXTRA_DIST = @SVNINCLUDE@ $(PGS_EXTRA_DIST)
