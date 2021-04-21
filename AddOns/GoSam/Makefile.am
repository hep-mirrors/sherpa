
GITTAG = AddOns/GoSam
include ../../git.make

GOSAM_SOURCES      =	\
			Git_Info.C \
                        GoSam_Virtual.C \
                        GoSam_Interface.C

GOSAM_ADDS         =	\
                        GoSam_Virtual.H \
                        GoSam_Interface.H

if GOSAM_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaGoSam.la

libSherpaGoSam_la_SOURCES  = $(GOSAM_SOURCES)

GOSAM_EXTRA_DIST = $(GOSAM_ADDS)
else
GOSAM_EXTRA_DIST = $(GOSAM_SOURCES) $(GOSAM_ADDS)
endif

EXTRA_DIST = $(GOSAM_EXTRA_DIST)
