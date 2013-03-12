
PGS_SOURCES   =	\
		pgs_finterface.f \
		Output_PGS.C

AM_FFLAGS     =	-fno-automatic

if PGS_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaPGSOutput.la

libSherpaPGSOutput_la_SOURCES         =	$(PGS_SOURCES)
libSherpaPGSOutput_la_LIBADD          =	@CONDITIONAL_PGSLIBS@
else
PGS_EXTRA_DIST    = $(PGS_SOURCES)
endif

EXTRA_DIST = @SVNINCLUDE@ $(PGS_EXTRA_DIST)
