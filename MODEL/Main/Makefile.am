pkglib_LTLIBRARIES = libModelMain.la

libModelMain_la_SOURCES = Model_Base.C \
			  Standard_Model.C SM_AGC.C MSSM.C SM_MSSM_EHC.C THDM.C \
			  ADD.C SM_Phantom_U1.C \
			  Running_AlphaS.C Running_AlphaQED.C Running_Fermion_Mass.C \
			  Effective_Higgs_Coupling.C LesHouches_Interface.C \
			  Spectrum_Generator_Base.C	

pkginclude_HEADERS = Model_Base.H \
		     Standard_Model.H SM_AGC.H MSSM.H SM_MSSM_EHC.H THDM.H \
		     ADD.H SM_Phantom_U1.H \
		     Running_AlphaS.H Running_AlphaQED.H  Running_Fermion_Mass.H \
	             Effective_Higgs_Coupling.H LesHouches_Interface.H \
		     Spectrum_Generator_Base.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES = @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
