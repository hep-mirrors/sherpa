
Phox_SOURCES  =	\
		Phox_Interface.C


if Phox_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaPhox.la

libSherpaPhox_la_SOURCES      =	$(Phox_SOURCES)

libSherpaPhox_la_LIBADD       =	@CONDITIONAL_PhoxLIBS@

libSherpaPhox_la_CPPFLAGS     =	$(AM_CPPFLAGS)

Phox_EXTRA_DIST = $(Phox_ADDS)
else
Phox_EXTRA_DIST = $(Phox_SOURCES) $(Phox_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(Phox_EXTRA_DIST)
