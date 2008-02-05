pkglib_LTLIBRARIES = libModelInteractions.la

libModelInteractions_la_SOURCES = Lorentz_Function.C Color_Function.C \
		      Single_Vertex.C Vertex.C \
                      Interaction_Model_Base.C \
		      Interaction_Model_QCD.C Interaction_Model_EW.C \
		      Interaction_Model_SM.C \
		      Interaction_Model_Higgs.C Interaction_Model_THDM.C \
		      Interaction_Model_EHC_S.C\
		      Interaction_Model_QCD_Grav.C Interaction_Model_EW_Grav.C \
		      Interaction_Model_ADD.C \
		      Interaction_Model_EHC.C Interaction_Model_SMEHC.C \
		      Interaction_Model_SM_ZPrime.C \
		      Interaction_Model_AEW.C Interaction_Model_SM_AGC.C \
		      Interaction_Model_MUED_QCD.C \
		      Interaction_Model_SM_Phantom_U1.C \
		      Interaction_Model_4GenLep.C \
		      Interaction_Model_Handler.C
#		      Interaction_Model_Inos.C Interaction_Model_sLepton_EW.C \
#		      Interaction_Model_sQCD.C Interaction_Model_sQuark_EW.C \
#		      Interaction_Model_sLepton_sQuark.C \
#		      Interaction_Model_MSSM.C \
#		      Interaction_Model_MSSM_EHC.C 

pkginclude_HEADERS = Lorentz_Function.H Color_Function.H \
		 Single_Vertex.H Vertex.H \
		 Interaction_Model_Base.H Interaction_Model_Handler.H \
		 Interaction_Model_QCD.H Interaction_Model_EW.H \
		 Interaction_Model_SM.H \
		 Interaction_Model_Higgs.H Interaction_Model_THDM.H \
		 Interaction_Model_Inos.H Interaction_Model_sLepton_EW.H \
		 Interaction_Model_sQCD.H Interaction_Model_sQuark_EW.H \
		 Interaction_Model_sLepton_sQuark.H \
	         Interaction_Model_MSSM.H Interaction_Model_MSSM_EHC.H \
	         Interaction_Model_QCD_Grav.H Interaction_Model_EW_Grav.H \
		 Interaction_Model_ADD.H Interaction_Model_EHC_S.H \
	         Interaction_Model_EHC.H Interaction_Model_SMEHC.H \
		 Interaction_Model_AEW.H Interaction_Model_SM_AGC.H \
	 	 Interaction_Model_MUED_QCD.H \
		 Interaction_Model_SM_Phantom_U1.H \
		 Interaction_Model_4GenLep.H \
	         Interaction_Model_SM_ZPrime.H

CXXFLAGS   =-g @SHERPAFLAGS@

INCLUDES = @AMEGICINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
