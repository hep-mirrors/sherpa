SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

SHERPAINCS    =	-I../Single_Events 
ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Analysis
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS	      = -I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libSherpaTools.la
libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				IO_Handler.C \
			    	HepMC_Interface.C \
				ISR_Statistics.C \
				Shower_Statistics.C \
			    	Analysis_Handler.C \
                                Event_Reader.C \
				fhepevt.f 
noinst_HEADERS        =	\
			Kinematics.H \
			IO_Handler.H \
		 	HepMC_Interface.H \
			ISR_Statistics.H \
			Shower_Statistics.H \
			Analysis_Handler.H \
			Event_Reader.H

EXTRA_DIST    =	$(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(CLHEPFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(APACICINCS) $(MODELINCS) \
		$(PDFINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ANALYSISINCS) \
		$(ROOTINCS)




