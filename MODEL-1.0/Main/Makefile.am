pkglib_LTLIBRARIES = libModelMain.la

libModelMain_la_SOURCES = Model_Handler.C Model_Base.C \
			  Standard_Model.C MSSM.C ADD.C\
			  Running_AlphaS.C Running_AlphaQED.C Running_Fermion_Mass.C \
			  Effective_Higgs_Coupling.C LesHouches_Interface.C \
			  Spectrum_Generator_Base.C	

pkginclude_HEADERS = Model_Handler.H Model_Base.H \
		 Standard_Model.H MSSM.H ADD.H \
		 Running_AlphaS.H Running_AlphaQED.H  Running_Fermion_Mass.H \
	         Effective_Higgs_Coupling.H LesHouches_Interface.H \
		 Spectrum_Generator_Base.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES = @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
