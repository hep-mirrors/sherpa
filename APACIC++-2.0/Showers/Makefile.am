lib_LTLIBRARIES = libShowers.la
libShowers_la_SOURCES = Final_State_Shower.C Timelike_Sudakov.C \
                        Timelike_Kinematics.C  Spacelike_Kinematics.C \
                        Initial_State_Shower.C Spacelike_Sudakov.C \
			Sudakov_Tools.C \
		        QED_Splitting_Functions.C \
                        QCD_Splitting_Functions.C \
                        Splitting_Function.C Splitting_Group.C Backward_Splitting_Group.C
INCLUDES = -I../Main \
	   -I../Hard_Processes \
	   -I../Hard_Processes/Process_Library \
	   -I../Hard_Processes/Integration_Tools \
	   -I../Hard_Interface \
           -I../../PDF-1.0/Include \
           -I../../PHASIC++-1.0/ISR \
           -I../../ATOOLS-2.0/Phys \
           -I../../ATOOLS-2.0/Math \
           -I../../ATOOLS-2.0/Org
# last 8 lines only for test in Final_State_Shower !!!
noinst_HEADERS = Final_State_Shower.H Timelike_Sudakov.H \
	         Timelike_Kinematics.H Spacelike_Kinematics.H \
                 Initial_State_Shower.H Spacelike_Sudakov.H \
                 Sudakov_Tools.H \
		 QED_Splitting_Functions.H QCD_Splitting_Functions.H \
                 Splitting_Group.H Splitting_Function.H \
		 Backward_Splitting_Group.H
