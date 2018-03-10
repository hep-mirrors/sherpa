
OPENLOOPS_SOURCES      =	\
			OpenLoops_Virtual.C \
n			OpenLoops_Born.C \
			OpenLoops_Interface.C \
			Color_Correlated_ME2.C \
			Spin_Color_Correlated_ME2.C \
			GGH_Process_Manager.C \
			GGH_KFactor_Setter.C \
			GGH_Decay_KFactor_Setter.C \
			EWVirtKFactor_Setter.C

OPENLOOPS_ADDS         =	\
			OpenLoops_Virtual.H \
			OpenLoops_Born.H \
			OpenLoops_Interface.H \
			Color_Correlated_ME2.H \
			Spin_Color_Correlated_ME2.H \
			GGH_Process_Manager.H \
			GGH_KFactor_Setter.H \
			EWVirtKFactor_Setter.H 

if OPENLOOPS_SUPPORT
pkglib_LTLIBRARIES    =	libSherpaOpenLoops.la

libSherpaOpenLoops_la_SOURCES  = $(OPENLOOPS_SOURCES)

OPENLOOPS_EXTRA_DIST = $(OPENLOOPS_ADDS)
else
OPENLOOPS_EXTRA_DIST = $(OPENLOOPS_SOURCES) $(OPENLOOPS_ADDS)
endif

EXTRA_DIST = $(OPENLOOPS_EXTRA_DIST)
