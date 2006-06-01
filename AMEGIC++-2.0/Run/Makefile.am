bin_PROGRAMS = Amegic

Amegic_SOURCES = Amplitude.C

Amegic_LDADD = $(pkglibdir)/libProcess.la \
               -L$(pkglibdir) @AMEGICLIBS@ @ATOOLSLIBS@ @MODELLIBS@ \
               @HDECAYLIBS@ @ISAJETLIBS@ @PDFLIBS@ @PHASICLIBS@ @BEAMLIBS@ \
               -lg2c -lm @CONDITIONAL_ROOTLIBS@

Amegic_LDFLAGS = @ldflags@

INCLUDES = @AMEGICINCS@ @PHASICINCS@ @MODELINCS@ @HDECAYINCS@ @ISAJETINCS@ @BEAMINCS@ \
           @PDFINCS@ @EXTRAXSINCS@ @ATOOLSINCS@
