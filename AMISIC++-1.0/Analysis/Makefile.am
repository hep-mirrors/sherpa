SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
MODELINC   = $(SHERPASYS)/MODEL-1.0

lib_LTLIBRARIES               =	libAmisicAnalysis.la
libAmisicAnalysis_la_SOURCES  = \
				Primordial_KT.C \
				Jet_PT.C \
				Jet_Y.C \
				Jet_Multiplicity.C \
				Jet_Observables.C
AM_CXXFLAGS   =	-pedantic -Wall -Winline -g
INCLUDES      = \
		-I../Tools \
		-I$(MODELINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org
noinst_HEADERS        =	\
			Primordial_KT.H \
			Jet_PT.H \
			Jet_Y.H \
			Jet_Multiplicity.H \
			Jet_Observables.H
