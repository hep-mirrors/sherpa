SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ANALYSISDIR   =	$(SHERPASYS)/ANALYSIS-1.0
SHERPAINCS    =	-I../Single_Events 
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Analysis
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 
ANALYSISINCS  = -I$(ANALYSISDIR)/Main -I$(ANALYSISDIR)/Observables

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libSherpaTools.la
libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				IO_Handler.C \
			    	HepMC_Interface.C \
			    	Analysis_Handler.C \
                                Event_Reader.C \
				fhepevt.f 
noinst_HEADERS        =	\
			Kinematics.H \
			IO_Handler.H \
		 	HepMC_Interface.H \
			Analysis_Handler.H \
			Event_Reader.H

EXTRA_DIST    =	$(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(CLHEPFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(MODELINCS) $(ATOOLSINCS) \
		$(CLHEPINCS) $(ANALYSISINCS)




