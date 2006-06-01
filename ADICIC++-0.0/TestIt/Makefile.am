SUBDIRS = data

bin_PROGRAMS = Flavtest Emissiontest Imissiontest Chaintest Cascadetest ICascadetest

Flavtest_SOURCES = Flavtest.C
Flavtest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                 @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                 -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                 @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

Emissiontest_SOURCES = Emissiontest.C
Emissiontest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                     @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                     -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                     @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

Imissiontest_SOURCES = Imissiontest.C
Imissiontest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                     @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                     -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                     @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

Chaintest_SOURCES = Chaintest.C
Chaintest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                  @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                  -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                  @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

Cascadetest_SOURCES = Cascadetest.C
Cascadetest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                     @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                     -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                     @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

ICascadetest_SOURCES = ICascadetest.C
ICascadetest_LDADD = -L$(pkglibdir) @ADICICLIBS@ @ATOOLSLIBS@ @PDFLIBS@ \
                     @CTEQLIBS@ @CONDITIONAL_LHAPDFLIBS@ @APACICLIBS@ \
                     -lg2c -lm -ldl @CONDITIONAL_ROOTLIBS@ @MODELLIBS@ \
                     @CONDITIONAL_HDECAYLIBS@ @CONDITIONAL_ISAJETLIBS@

AM_CXXFLAGS = @SHERPAFLAGS@ -ftemplate-depth-100

INCLUDES = @ADICICINCS@ @ATOOLSINCS@ @MODELINCS@ @PDFINCS@

EXTRA_DIST = @SVNINCLUDE@
##EXTRA_DIST = Testrun/Run.dat Testrun/Integration.dat
