SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
SHERPAINCS    =	-I../Single_Events 
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Analysis
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math \
		-I$(ATOOLSDIR)/Org 

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
				Output_Handler.C \
			    	HepMC_Interface.C \
				HepEvt_Interface.C \
			    	Sample_Analysis.C \
			    	Analysis_Phase.C \
			    	EvtReadin_Phase.C \
                                Event_Reader.C \
				fhepevt.f 
noinst_HEADERS        =	\
			Kinematics.H \
			Output_Handler.H \
		 	HepMC_Interface.H \
			HepEvt_Interface.H \
			Sample_Analysis.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			Event_Reader.H

EXTRA_DIST    =	$(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(CLHEPFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(MODELINCS) $(ATOOLSINCS) \
		$(CLHEPINCS)




