SHERPASYS     =	../..

ANALYSISFLAGS =	-pedantic -Wall -Winline

ANALYSISINCS  =	-I../Observables
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES        = libAnalysis.la
libAnalysis_la_SOURCES = \
			 Final_Selector.C \
			 Primitive_Observable_Base.C \
			 Primitive_Observable.C \
			 Primitive_Detector.C \
			 Primitive_Calorimeter.C \
			 Primitive_Analysis.C \
	 		 Particle_Selector.C
noinst_HEADERS         = \
			 Final_Selector.H \
			 Primitive_Observable_Base.H \
			 Primitive_Observable.H \
			 Primitive_Detector.H \
			 Primitive_Calorimeter.H \
			 Primitive_Analysis.H \
	 		 Particle_Selector.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ANALYSISFLAGS)
INCLUDES      =	$(ANALYSISINCS) $(MODELINCS) $(ATOOLSINCS)

