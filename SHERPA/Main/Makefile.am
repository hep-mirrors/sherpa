pkglib_LTLIBRARIES       =	libSherpaMain.la

SVNTAG = SHERPA/Main
include ../../svn.make

libSherpaMain_la_SOURCES      =	SVN_Info.C Sherpa.C Filter.C 

localinc_HEADERS          =    Sherpa.H Filter.H

SYSLIBS	      =	$(FCLIBS) -lm -ldl $(CONDITIONAL_THREADLIBS)

libSherpaMain_la_LIBADD = @SHERPALIBS@ @SHRIMPSLIBS@ \
		@HADRONSLIBS@ @AHADICLIBS@ \
		@AMISICLIBS@ @PHOTONSLIBS@ \
		@MCATNLOLIBS@ @CSSLIBS@ @DIRELIBS@ \
		@COMIXLIBS@ @AMEGICLIBS@ @EXTRAXSLIBS@ \
		@PHASICLIBS@ @METOOLSLIBS@ @MODELLIBS@ \
		@PDFLIBS@ @BEAMLIBS@ @ATOOLSLIBS@ \
		@CONDITIONAL_HEPMC2LIBS@ \
                $(SYSLIBS)
