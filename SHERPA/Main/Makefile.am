pkglib_LTLIBRARIES       =	libSherpaMain.la

SVNTAG = SHERPA/Main
include ../../svn.make

libSherpaMain_la_SOURCES      =	SVN_Info.C Sherpa.C 

localinc_HEADERS          =    Sherpa.H 

SYSLIBS	      =	$(FCLIBS) -lm -ldl $(CONDITIONAL_THREADLIBS)

libSherpaMain_la_LIBADD = -L$(SHERPABUILDDIR)/Initialization -L$(SHERPABUILDDIR)/Single_Events -L$(SHERPABUILDDIR)/Tools \
                -L$(SHERPABUILDDIR)/PerturbativePhysics -L$(SHERPABUILDDIR)/SoftPhysics -L$(SHERPABUILDDIR)/LundTools \
                -lSherpaInitialization -lSherpaSingleEvents \
                -lSherpaPerturbativePhysics -lSherpaSoftPhysics -lLundTools -lSherpaTools \
		@AMISICLIBS@ @CSSLIBS@ @ATOOLSLIBS@ @AHADICLIBS@ @SINICLIBS@ \
		@BEAMLIBS@ @COMIXLIBS@ @AMEGICLIBS@ @MODELLIBS@ @PDFLIBS@ \
		@MCATNLOLIBS@ @POWHEGLIBS@ @PHASICLIBS@  \
		@METOOLSLIBS@ @CONDITIONAL_HEPMC2LIBS@ \
		@CONDITIONAL_HADRONSLIBS@ @CONDITIONAL_PHOTONSLIBS@ \
		@CONDITIONAL_GZIPLIBS@ \
		@EXTRAXSLIBS@ \
                $(SYSLIBS)
