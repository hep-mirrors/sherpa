lib_LTLIBRARIES = libSherpaTools.la
CLHEPlibdir   = $(CLHEPDIR)/lib
CLHEPincdir   = $(CLHEPDIR)/include
libSherpaTools_la_SOURCES = Output_Handler.C \
			    HepMC_Interface.C IO_HepEvt.C \
			    Sample_Analysis.C
INCLUDES = -I../Single_Events \
	   -I../../HepMC_Version1.23 \
	   -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org \
	   -I../Single_Events \
	   -I$(CLHEPincdir)
noinst_HEADERS = Output_Handler.H \
		 HepMC_Interface.H IO_HepEvt.H \
		 Sample_Analysis.H




