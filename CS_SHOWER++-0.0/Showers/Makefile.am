SHERPASYS  = ../..
ATOOLSDIR  = $(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS  = -I$(SHERPASYS)/MODEL-1.0/Main
CSSDIR     = $(SHERPASYS)/CS_SHOWER++-0.0
CSSINCS    = -I$(CSSDIR)/Tools

CSSFLAGS  =	-pedantic -Wall -Winline 
AM_CXXFLAGS   = $(CSSFLAGS)

if CVS_INCLUDE
CVSINCLUDE    = CVS
endif

lib_LTLIBRARIES = libCSShowers.la
libCSShowers_la_SOURCES = \
	Splitting_Function_Base.C Splitting_Function_Group.C \
	QCD_Splitting_Functions.C Sudakov.C Shower.C \
	Kinematics_Base.C

INCLUDES = \
	$(ATOOLSINCS) $(MODELINCS) $(CSSINCS)

noinst_HEADERS = \
	Splitting_Function_Base.H Splitting_Function_Group.H \
	QCD_Splitting_Functions.H Sudakov.H Shower.H \
	Kinematics_Base.H

EXTRA_DIST    = $(CVSINCLUDE)
