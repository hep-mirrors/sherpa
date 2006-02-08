SHERPASYS     =	../..

APACICFLAGS   =	# -pedantic -Wall -Winline 
DEFS          =

APACICINCS    = -I../Main 
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

pkglib_LTLIBRARIES = libApacicShowers.la
libApacicShowers_la_SOURCES   = \
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
pkginclude_HEADERS        = \
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
AM_CXXFLAGS   = $(APACICFLAGS)
INCLUDES      = $(APACICINCS) $(MODELINCS) $(PDFINCS) $(BEAMINCS) \
		$(ATOOLSINCS)
EXTRA_DIST    = $(CVSINCLUDE)
