SHERPASYS     =	../..

ANALYSISFLAGS =	-pedantic -Wall -Winline

ANALYSISINCS  =	-I../Main
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libObservables.la
libObservables_la_SOURCES     =	\
			    	Statistics_Observable.C \
			    	MI_Observables.C \
			    	ISR_Statistics.C \
			    	Shower_Statistics.C \
			    	One_Particle_Observables.C \
			    	Two_Particle_Observables.C \
			    	Four_Particle_Observables.C \
			    	Calorimeter_Cone.C \
			    	Jet_Cone_Distribution.C \
			    	Jet_Observables.C \
			    	Event_Shapes_EE.C \
			    	Shape_Observables_EE.C \
			    	Scaled_Observables.C
noinst_HEADERS        =	\
			Statistics_Observable.H \
			MI_Observables.H \
			ISR_Statistics.H \
			Shower_Statistics.H \
			One_Particle_Observables.H \
			Two_Particle_Observables.H \
			Four_Particle_Observables.H \
			Calorimeter_Cone.H \
			Jet_Cone_Distribution.H \
			Jet_Observables.H \
			Event_Shapes_EE.H \
			Shape_Observables_EE.H \
			Scaled_Observables.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ANALYSISFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(ANALYSISINCS) $(PDFINCS) $(ATOOLSINCS) $(ROOTINCS)
