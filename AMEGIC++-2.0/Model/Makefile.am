lib_LTLIBRARIES = libModel.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libModel_la_SOURCES = Interaction_Model_Base.C Interaction_Model_Handler.C \
		      Interaction_Model_QCD.C Interaction_Model_EW.C \
		      Interaction_Model_SM.C \
		      Interaction_Model_Higgs.C Interaction_Model_THDM.C \
		      Interaction_Model_Inos.C Interaction_Model_sLepton_EW.C \
		      Interaction_Model_sQCD.C Interaction_Model_sQuark_EW.C \
		      Interaction_Model_sLepton_sQuark.C \
		      Interaction_Model_MSSM.C \
		      Interaction_Model_QCD_Grav.C Interaction_Model_EW_Grav.C \
		      Interaction_Model_ADD.C \
		      Interaction_Model_EHC.C Interaction_Model_SMEHC.C \
		      Vertex.C Lorentz_Function.C
INCLUDES = -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Model \
	   -I../String \
           -I../../MODEL-1.0/Main \
           -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = Interaction_Model_Base.H Interaction_Model_Handler.H \
		 Interaction_Model_QCD.H Interaction_Model_EW.H \
		 Interaction_Model_SM.H \
		 Interaction_Model_Higgs.H Interaction_Model_THDM.H \
		 Interaction_Model_Inos.H Interaction_Model_sLepton_EW.H \
		 Interaction_Model_sQCD.H Interaction_Model_sQuark_EW.H \
		 Interaction_Model_sLepton_sQuark.H \
	         Interaction_Model_MSSM.H \
	         Interaction_Model_QCD_Grav.H Interaction_Model_EW_Grav.H \
		 Interaction_Model_ADD.H \
	         Interaction_Model_EHC.H Interaction_Model_SMEHC.H \
	         Vertex.H Lorentz_Function.H Color_Function.H 
CXXFLAGS=-g
EXTRA_DIST    = $(CVSINCLUDE)