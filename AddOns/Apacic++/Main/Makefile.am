APAMAIN_SOURCES       =	\
			Knot.C \
			Tree.C \
			Apacic_Cluster_Definitions.C \
			Apacic.C

APAMAIN_ADDS          =	\
			Knot.H \
			Tree.H \
			Apacic_Cluster_Definitions.H

if USING__Apacic
pkglib_LTLIBRARIES = libApacicMain.la

libApacicMain_la_SOURCES      =	$(APAMAIN_SOURCES)

libApacicMain_la_LIBADD       =	-L../Showers -lApacicShowers

APAMAIN_EXTRA_DIST = $(APAMAIN_ADDS)
else
APAMAIN_EXTRA_DIST = $(APAMAIN_SOURCES) $(APAMAIN_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(APAMAIN_EXTRA_DIST)
