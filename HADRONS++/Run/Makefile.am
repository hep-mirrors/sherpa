bin_PROGRAMS  = SingleDecay \
	SingleDecayPhotons \
	HadronFullDecayNoAnalysis \
	HadronFullDecayMixingAnalysis \
	HadronFullDecayInclusiveAnalysis \
	PhotonAnalysis

FULL_LDADD = $(FCLIBS) -lm -ldl -L$(pkglibdir) -lSherpaMain \
  @CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@

SingleDecay_SOURCES = Main_SingleDecay.C
SingleDecay_LDADD   = $(FULL_LDADD)
SingleDecay_LDFLAGS = @ldflags@

SingleDecayPhotons_SOURCES = Main_SingleDecayPhotons.C
SingleDecayPhotons_LDADD   = $(FULL_LDADD)
SingleDecayPhotons_LDFLAGS = @ldflags@

HadronFullDecayNoAnalysis_SOURCES = Main_HadronFullDecayNoAnalysis.C
HadronFullDecayNoAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayNoAnalysis_LDFLAGS = @ldflags@

HadronFullDecayMixingAnalysis_SOURCES = Main_HadronFullDecayMixingAnalysis.C
HadronFullDecayMixingAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayMixingAnalysis_LDFLAGS = @ldflags@

HadronFullDecayInclusiveAnalysis_SOURCES = Main_HadronFullDecayInclusiveAnalysis.C
HadronFullDecayInclusiveAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayInclusiveAnalysis_LDFLAGS = @ldflags@

PhotonAnalysis_SOURCES = Main_PhotonFullDistributionAnalysis.C
PhotonAnalysis_LDADD   = $(FULL_LDADD)
PhotonAnalysis_LDFLAGS = @ldflags@


AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @HADRONSINCS@ @APACICINCS@ \
@PHASICINCS@ @MODELINCS@ @PDFINCS@ @CTEQINCS@ \
@BEAMINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@ \
@SHERPAINCS@ @AMEGICINCS@ @AHADICINCS@ @ANALYSISINCS@ \
@CONDITIONAL_PHOTONSINCS@

EXTRA_DIST    = @SVNINCLUDE@ Main_FullDecay.C Main.H
