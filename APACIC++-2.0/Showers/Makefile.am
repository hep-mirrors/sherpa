if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

APACICFLAGS   =	-pedantic -Wall -Winline 

SHERPASYS     =	../..

if SHERPA_SUPPORT
SHERPAFLAGS   = -DSHERPA_SUPPORT
SHERPAINCS    =	-I$(SHERPASYS)/SHERPA-1.0/SoftPhysics
endif

lib_LTLIBRARIES = libApacicShowers.la
libApacicShowers_la_SOURCES = Final_State_Shower.C Timelike_Sudakov.C \
			      Timelike_Kinematics.C  Spacelike_Kinematics.C \
			      Initial_State_Shower.C Spacelike_Sudakov.C \
			      Sudakov_Tools.C Splitting_Function.C \
			      Splitting_Group.C Backward_Splitting_Group.C \
			      QED_Splitting_Functions.C QCD_Splitting_Functions.C \
			      SUSY_QCD_Splitting_Functions.C
AM_CXXFLAGS   = $(APACICFLAGS) $(SHERPAFLAGS)
INCLUDES = $(SHERPAINCS) -I../../APACIC++-2.0/Main \
	   -I../../MODEL-1.0/Main \
	   -I../../PDF-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org
noinst_HEADERS = Final_State_Shower.H Timelike_Sudakov.H \
	         Timelike_Kinematics.H Spacelike_Kinematics.H \
                 Initial_State_Shower.H Spacelike_Sudakov.H \
                 Sudakov_Tools.H Splitting_Function.H \
		 Splitting_Group.H Backward_Splitting_Group.H \
		 QED_Splitting_Functions.H QCD_Splitting_Functions.H \
		 SUSY_QCD_Splitting_Functions.H
EXTRA_DIST    = $(CVSINCLUDE)
