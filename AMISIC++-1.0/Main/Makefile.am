SHERPASYS  = ../..

ATOOLSINC  = $(SHERPASYS)/ATOOLS-2.0
PDFINC     = $(SHERPASYS)/PDF-1.0
BEAMINC    = $(SHERPASYS)/BEAM-1.0
PHASICINC  = $(SHERPASYS)/PHASIC++-1.0
MODELINC   = $(SHERPASYS)/MODEL-1.0
EXTRAXSINC = $(SHERPASYS)/EXTRA_XS-1.0
SHERPAINC  = $(SHERPASYS)/SHERPA-1.0

lib_LTLIBRARIES       = libAmisic.la
libAmisic_la_SOURCES  = \
			MI_Base.C \
			Environment.C \
			Simple_Chain.C \
			Amisic.C 
AM_CXXFLAGS   =	-pedantic -Wall -Winline -DUSING__SHERPA
INCLUDES      =	\
		-I../Tools \
		-I../Other \
		-I$(SHERPAINC)/LundTools \
		-I$(SHERPAINC)/PerturbativePhysics \
		-I$(EXTRAXSINC)/Two2Two \
		-I$(EXTRAXSINC)/Main \
		-I$(PHASICINC)/Main \
		-I$(MODELINC)/Main \
		-I$(PDFINC)/Main \
		-I$(BEAMINC)/Main \
		-I$(ATOOLSINC)/Phys \
		-I$(ATOOLSINC)/Math \
		-I$(ATOOLSINC)/Org 
noinst_HEADERS        = \
			Grid_Handler.H \
			Grid_Handler.IC \
			Grid_Handler.C \
			Grid_Creator_Base.H \
			Grid_Creator_Base.IC \
		       	Grid_Creator_Base.C \
			Grid_Creator.H \
			Grid_Creator.IC \
		       	Grid_Creator.C \
			Fixed_Variable_Channel.H \
			Fixed_Variable_Channel.IC \
			Fixed_Variable_Channel.C \
			MI_Base.H \
			MI_Base.IC \
			Environment.H \
			Simple_Chain.H \
			Simple_Chain.IC \
			Amisic.H \
			Amisic.IC 
