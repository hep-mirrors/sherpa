SHERPASYS     =	../..

ANALYSISFLAGS =	-pedantic -Wall -Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
ANALYSISINCS  =	-I../Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES           = libObservables.la
libObservables_la_SOURCES = \
			    Statistics_Observable.C \
			    MI_Statistics.C \
			    One_Particle_Observables.C \
			    Two_Particle_Observables.C \
			    Four_Particle_Observables.C \
			    Jet_Observables.C \
			    Event_Shapes_EE.C \
			    Shape_Observables_EE.C 
noinst_HEADERS            = \
			    Statistics_Observable.H \
			    MI_Statistics.H \
			    One_Particle_Observables.H \
			    Two_Particle_Observables.H \
			    Four_Particle_Observables.H \
			    Jet_Observables.H \
			    Event_Shapes_EE.H \
			    Shape_Observables_EE.H

EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ANALYSISFLAGS)
INCLUDES      =	$(ANALYSISINCS) $(ATOOLSINCS)
