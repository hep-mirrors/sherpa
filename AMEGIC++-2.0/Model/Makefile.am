SUBDIRS = Isajet
lib_LTLIBRARIES = libModel.la
libModel_la_SOURCES = Model.C Model_Handler.C \
                 Model_EE_QCD.C Model_EW.C Model_QCD.C \
		 Model_SM.C Model_NC_QED.C \
                 Model_Chi.C Model_sLep_Chi.C Model_THDM.C Model_SQCD.C \
		 Model_sLeptons.C Model_SQCD_Chi.C Model_Higgs.C \
		 Model_MSSM.C \
                 Spectrum_EW.C \
		 Spectrum_Higgs.C Spectrum_sDownquarks.C \
		 Spectrum_sUpquarks.C Spectrum_sLeptons.C \
		 Spectrum_sNeutrinos.C\
		 Spectrum_Charginos.C Spectrum_Neutralinos.C \
                 Couplings_EW.C \
	         Vertex.C Lorentz_Function.C
INCLUDES = -I./Isajet/ \
	   -I../Amplitude \
	   -I../Amplitude/AmplTools \
	   -I../Amplitude/Zfunctions \
	   -I../Model \
	   -I../String \
           -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org 
noinst_HEADERS = Model.H Model_Handler.H \
                 Model_EE_QCD.H Model_EW.H Model_Higgs.H Model_QCD.H \
                 Model_Chi.H Model_sLep_Chi.H Model_THDM.H Model_SQCD.H \
		 Model_sLeptons.H Model_SQCD_Chi.H \
		 Model_SM.H Model_NC_QED.H \
		 Model_MSSM.H \
	         Spectrum_EW.H Spectrum_Higgs.H \
		 Spectrum_Charginos.H Spectrum_Neutralinos.H \
		 Spectrum_sDownquarks.H Spectrum_sUpquarks.H \
		 Spectrum_sLeptons.H Spectrum_sNeutrinos.H \
                 Couplings_EW.H \
	         Vertex.H Lorentz_Function.H Color_Function.H 
CXXFLAGS=-g
