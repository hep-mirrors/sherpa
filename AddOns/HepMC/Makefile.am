
HEPMC_SOURCES =	Output_HepMC2_Genevent.C Output_HepMC2_Short.C
HEPMC_ADDS    =	Output_HepMC2_Genevent.H Output_HepMC2_Short.H

if HEPMC2_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaHepMCOutput.la

libSherpaHepMCOutput_la_SOURCES       =	$(HEPMC_SOURCES)

libSherpaHepMCOutput_la_LIBADD        =	@CONDITIONAL_HEPMC2LIBS@

libSherpaHepMCOutput_la_CPPFLAGS      =	$(AM_CPPFLAGS) @CONDITIONAL_GZIPINCS@ @CONDITIONAL_HEPMC2INCS@
else
HEPMC_EXTRA_DIST = $(HEPMC_SOURCES)
endif

EXTRA_DIST = $(HEPMC_EXTRA_DIST) $(HEPMC_ADDS)
