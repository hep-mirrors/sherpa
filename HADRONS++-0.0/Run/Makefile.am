SYSLIBS	      =	$(FLIBS) -lm -ldl

bin_PROGRAMS  = HadronFullDecayNoAnalysis \
	HadronFullDecayMixingAnalysis \
	HadronFullDecayInclusiveAnalysis \
  PhotonAnalysis

HadronFullDecayNoAnalysis_SOURCES        =	Main_HadronFullDecayNoAnalysis.C
HadronFullDecayNoAnalysis_LDADD  =	$(SYSLIBS) -L$(pkglibdir) -lSherpaMain \
	@CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@
HadronFullDecayNoAnalysis_LDFLAGS        =	@ldflags@


HadronFullDecayMixingAnalysis_SOURCES        =	Main_HadronFullDecayMixingAnalysis.C
HadronFullDecayMixingAnalysis_LDADD  =	$(SYSLIBS) -L$(pkglibdir) -lSherpaMain \
	@CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@
HadronFullDecayMixingAnalysis_LDFLAGS        =	@ldflags@


HadronFullDecayInclusiveAnalysis_SOURCES        =	Main_HadronFullDecayInclusiveAnalysis.C
HadronFullDecayInclusiveAnalysis_LDADD  =	$(SYSLIBS) -L$(pkglibdir) -lSherpaMain \
	@CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@
HadronFullDecayInclusiveAnalysis_LDFLAGS        =	@ldflags@

PhotonAnalysis_SOURCES        =	Main_PhotonFullDistributionAnalysis.C
PhotonAnalysis_LDADD  =	$(SYSLIBS) -L$(pkglibdir) -lSherpaMain \
	@CONDITIONAL_ROOTLIBS@ @CONDITIONAL_LHAPDFLIBS@ @CONDITIONAL_HEPMC2LIBS@
PhotonAnalysis_LDFLAGS        =	@ldflags@


AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @HADRONSINCS@ @APACICINCS@ \
@PHASICINCS@ @MODELINCS@ @PDFINCS@ @CTEQINCS@ \
@BEAMINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@ \
@SHERPAINCS@ @AMEGICINCS@ @AHADICINCS@ @ANALYSISINCS@ \
@CONDITIONAL_PHOTONSINCS@

EXTRA_DIST    = @SVNINCLUDE@
