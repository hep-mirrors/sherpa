
GITTAG = AddOns/PGS
include ../../git.make

PGS_SOURCES   =	\
		pgs_finterface.f \
		Git_Info.C \
		Output_PGS.C

AM_FFLAGS     =	-fno-automatic

if PGS_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaPGSOutput.la

libSherpaPGSOutput_la_SOURCES         =	$(PGS_SOURCES)
libSherpaPGSOutput_la_LIBADD          =	@CONDITIONAL_PGSLIBS@
else
PGS_EXTRA_DIST    = $(PGS_SOURCES)
endif

EXTRA_DIST = $(PGS_EXTRA_DIST)
