pkglib_LTLIBRARIES = libModelMain.la

SVNTAG = MODEL/Main
include ../../svn.make

libModelMain_la_SOURCES = SVN_Info.C \
			  Coupling_Data.C Model_Base.C Model_Base_GF.C \
			  Standard_Model.C SM_AGC.C MSSM.C SM_MSSM_EHC.C \
			  THDM.C ADD.C SM_Phantom_U1.C SM_AxiGluon.C \
			  SM_U1_B.C Hidden_Valley.C \
			  Running_AlphaS.C Running_AlphaQED.C \
			  Running_Alpha_HV.C Running_Fermion_Mass.C \
			  Effective_Higgs_Coupling.C LesHouches_Interface.C \
			  Spectrum_Generator_Base.C \
			  FeynRules_Model.C FeynRules_Spectrum.C

localinc_HEADERS =   Coupling_Data.H Model_Base.H \
		     Standard_Model.H SM_AGC.H MSSM.H SM_MSSM_EHC.H \
		     THDM.H ADD.H SM_Phantom_U1.H SM_AxiGluon.H \
	 	     SM_U1_B.H Hidden_Valley.H \
		     Running_AlphaS.H Running_AlphaQED.H \
		     Running_Alpha_HV.H Running_Fermion_Mass.H \
	             Effective_Higgs_Coupling.H LesHouches_Interface.H \
		     Spectrum_Generator_Base.H \
		     FeynRules_Model.H FeynRules_Spectrum.H
