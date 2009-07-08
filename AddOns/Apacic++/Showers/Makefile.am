APASHOWERS_SOURCES    = \
			Sudakov_Tools.C \
			Splitting_Function.C \
			Splitting_Group.C \
			Backward_Splitting_Group.C \
			QED_Splitting_Functions.C \
			QCD_Splitting_Functions.C \
			SUSY_QCD_Splitting_Functions.C \
			Timelike_Kinematics.C \
			Timelike_Sudakov.C \
			Spacelike_Kinematics.C \
			Spacelike_Sudakov.C \
			Jet_Veto.C \
			Final_State_Shower.C \
			Initial_State_Shower.C

APASHOWERS_ADDS       = \
			Sudakov_Tools.H \
			Splitting_Function.H \
			Splitting_Group.H \
			Backward_Splitting_Group.H \
		 	QED_Splitting_Functions.H \
			QCD_Splitting_Functions.H \
		 	SUSY_QCD_Splitting_Functions.H \
	         	Timelike_Kinematics.H \
			Timelike_Sudakov.H \
			Spacelike_Kinematics.H \
			Spacelike_Sudakov.H \
			Jet_Veto.H \
			Final_State_Shower.H \
			Initial_State_Shower.H 

if USING__Apacic
pkglib_LTLIBRARIES = libApacicShowers.la

libApacicShowers_la_SOURCES   =	$(APASHOWERS_SOURCES)

APASHOWERS_EXTRA_DIST = $(APASHOWERS_ADDS)
else
APASHOWERS_EXTRA_DIST = $(APASHOWERS_SOURCES) $(APASHOWERS_ADDS)
endif

EXTRA_DIST = @SVNINCLUDE@ $(APASHOWERS_EXTRA_DIST)
