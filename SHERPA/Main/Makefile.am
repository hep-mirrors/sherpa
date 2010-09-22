pkglib_LTLIBRARIES       =	libSherpaMain.la

libSherpaMain_la_SOURCES      =	Sherpa.C 

localinc_HEADERS          =	Sherpa.H 

SYSLIBS	      =	$(FCLIBS) -lm -ldl $(CONDITIONAL_THREADLIBS)

libSherpaMain_la_LIBADD = -L$(SHERPABUILDDIR)/Initialization -L$(SHERPABUILDDIR)/Single_Events -L$(SHERPABUILDDIR)/Tools \
                -L$(SHERPABUILDDIR)/PerturbativePhysics -L$(SHERPABUILDDIR)/SoftPhysics -L$(SHERPABUILDDIR)/LundTools \
                -lSherpaInitialization -lSherpaSingleEvents \
                -lSherpaPerturbativePhysics -lSherpaSoftPhysics -lLundTools -lSherpaTools \
		@AMISICLIBS@ @CSSLIBS@ @ATOOLSLIBS@ @AHADICLIBS@ \
		@BEAMLIBS@ @COMIXLIBS@ @AMEGICLIBS@ @MODELLIBS@ @PDFLIBS@ @POWHEGLIBS@ @PHASICLIBS@  \
		@METOOLSLIBS@\
		@CONDITIONAL_HADRONSLIBS@ @CONDITIONAL_PHOTONSLIBS@ \
		@CONDITIONAL_GZIPLIBS@ @CONDITIONAL_ROOTLIBS@ \
		@EXTRAXSLIBS@ \
                $(SYSLIBS)

libSherpaMain_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaMain_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
