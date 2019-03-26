GITTAG = AddOns/DiHiggsNLO
include ../../git.make

DIHIGGS_SOURCES = Git_Info.C DiHiggs_Virtual.C hhgrid.c
DIHIGGS_ADDS = DiHiggs_Virtual.H hhgrid.h

if ENABLE_DIHIGGS

pkglib_LTLIBRARIES    =	libSherpaDiHiggsNLO.la
libSherpaDiHiggsNLO_la_SOURCES = $(DIHIGGS_SOURCES)
DIHIGGS_EXTRA_DIST = $(DIHIGGS_ADDS)
libSherpaDiHiggsNLO_la_LIBADD = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LIBS)
libSherpaDiHiggsNLO_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)

else
DIHIGGS_EXTRA_DIST = $(DIHIGGS_SOURCES) $(DIHIGGS_ADDS)
endif

EXTRA_DIST = $(DIHIGGS_EXTRA_DIST)

dihiggsdir = $(pythondir)

dist_dihiggs_PYTHON = creategrid.py
