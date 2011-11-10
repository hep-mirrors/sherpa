bin_PROGRAMS  = SingleDecay \
	SingleDecayPhotons \
	HadronFullDecayNoAnalysis \
	PhotonAnalysis \
	PhotonAnalysisBR \
	HadronFullDecayMixingAnalysis \
	HadronFullDecayInclusiveAnalysis

FULL_LDADD = $(FCLIBS) -lm -ldl -L$(pkglibdir) -lSherpaMain \
  @CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@

SingleDecay_SOURCES = Main_SingleDecay.C
SingleDecay_LDADD   = $(FULL_LDADD)
SingleDecay_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
SingleDecay_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

SingleDecayPhotons_SOURCES = Main_SingleDecayPhotons.C
SingleDecayPhotons_LDADD   = $(FULL_LDADD)
SingleDecayPhotons_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
SingleDecayPhotons_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

HadronFullDecayNoAnalysis_SOURCES = Main_HadronFullDecayNoAnalysis.C
HadronFullDecayNoAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayNoAnalysis_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
HadronFullDecayNoAnalysis_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

HadronFullDecayMixingAnalysis_SOURCES = Main_HadronFullDecayMixingAnalysis.C
HadronFullDecayMixingAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayMixingAnalysis_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
HadronFullDecayMixingAnalysis_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

HadronFullDecayInclusiveAnalysis_SOURCES = Main_HadronFullDecayInclusiveAnalysis.C
HadronFullDecayInclusiveAnalysis_LDADD   = $(FULL_LDADD)
HadronFullDecayInclusiveAnalysis_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
HadronFullDecayInclusiveAnalysis_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

PhotonAnalysis_SOURCES = Main_PhotonFullDistributionAnalysis.C
PhotonAnalysis_LDADD   = $(FULL_LDADD)
PhotonAnalysis_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
PhotonAnalysis_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

PhotonAnalysisBR_SOURCES = Main_PhotonFullDistributionBRAnalysis.C
PhotonAnalysisBR_LDADD   = $(FULL_LDADD)
PhotonAnalysisBR_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@
PhotonAnalysisBR_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@ Main_FullDecay.C Main.H Run.dat Analysis.dat
