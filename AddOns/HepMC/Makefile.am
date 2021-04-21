GITTAG = AddOns/HepMC
include ../../git.make

HEPMC2_SOURCES =	Git_Info.C Output_HepMC2_Genevent.C Output_HepMC2_Short.C
HEPMC3_SOURCES =	Git_Info.C Output_HepMC3_Genevent.C
HEPMC2_ADDS    =	Output_HepMC2_Genevent.H Output_HepMC2_Short.H
HEPMC3_ADDS    =	Output_HepMC3_Genevent.H

pkglib_LTLIBRARIES  = 
if HEPMC2_SUPPORT
pkglib_LTLIBRARIES    +=	libSherpaHepMCOutput.la
endif

if HEPMC3_SUPPORT
pkglib_LTLIBRARIES    +=	libSherpaHepMC3Output.la
endif


libSherpaHepMCOutput_la_SOURCES       =	$(HEPMC2_SOURCES)
libSherpaHepMC3Output_la_SOURCES       =	$(HEPMC3_SOURCES)

libSherpaHepMCOutput_la_LIBADD        =	@CONDITIONAL_HEPMC2LIBS@
libSherpaHepMC3Output_la_LIBADD        =	@CONDITIONAL_HEPMC3LIBS@ @CONDITIONAL_ROOTLIBS@

libSherpaHepMCOutput_la_CPPFLAGS      =	$(AM_CPPFLAGS) @CONDITIONAL_HEPMC2INCS@ @CONDITIONAL_GZIPINCS@  
libSherpaHepMC3Output_la_CPPFLAGS      = -std=c++11 $(AM_CPPFLAGS) @CONDITIONAL_HEPMC3INCS@ @CONDITIONAL_ROOTINCS@ 

HEPMC2_EXTRA_DIST = $(HEPMC2_SOURCES)
HEPMC3_EXTRA_DIST = $(HEPMC3_SOURCES)


EXTRA_DIST = $(HEPM2_EXTRA_DIST) $(HEPMC2_ADDS) $(HEPMC3_EXTRA_DIST) $(HEPMC3_ADDS)
