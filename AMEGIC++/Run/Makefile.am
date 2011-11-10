bin_PROGRAMS = Amegic

Amegic_SOURCES = Amplitude.C

Amegic_LDADD = $(pkglibdir)/libProcess.la \
           -L$(pkglibdir) @AMEGICLIBS@ @ATOOLSLIBS@ @MODELLIBS@ \
           @PDFLIBS@ @PHASICLIBS@ @BEAMLIBS@ \
	       -lg2c -lm @CONDITIONAL_ROOTLIBS@
