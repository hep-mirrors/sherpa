lib_LTLIBRARIES = libSherpaTools.la
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
libSherpaTools_la_SOURCES = Output_Handler.C \
			    HepMC_Interface.C HepEvt_Interface.C fhepevt.f \
			    Kinematics.C Sample_Analysis.C 
INCLUDES = -I../Single_Events \
	   -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main \
	   -I$(HEPMCINCDIR) \
	   -I$(CLHEPincdir)
noinst_HEADERS = Output_Handler.H \
		 HepMC_Interface.H HepEvt_Interface.H \
		 Kinematics.H Sample_Analysis.H 




