lib_LTLIBRARIES = libSherpaTools.la
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
libSherpaTools_la_SOURCES = Output_Handler.C \
			    HepMC_Interface.C HepEvt_Interface.C \
			    Sample_Analysis.C \
			    fhepevt.f
INCLUDES = -I../Single_Events \
	   -I../../HepMC_Version1.23 \
	   -I../../HepMC_Version1.23/HepMC \
	   -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main \
	   -I$(CLHEPincdir)
noinst_HEADERS = Output_Handler.H \
		 HepMC_Interface.H HepEvt_Interface.H \
		 Sample_Analysis.H




