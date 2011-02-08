
pkglib_LTLIBRARIES    =	libPhasicScales.la

SVNTAG = PHASIC++/Scales
include ../../svn.make

libPhasicScales_la_SOURCES    =	\
				SVN_Info.C \
				Scale_Setter_Base.C \
				KFactor_Setter_Base.C \
				Tag_Setter.C \
				Variable_Scale_Setter.C \
				Fastjet_Scale_Setter.C \
				MPI_Scale_Setter.C \
				QCD_Scale_Setter.C \
				METS_Scale_Setter.C \
				Variable_KFactor_Setter.C \
				No_KFactor_Setter.C

localinc_HEADERS      =	\
			Scale_Setter_Base.H \
			KFactor_Setter_Base.H \
			Tag_Setter.H

libPhasicScales_la_LIBADD     =	@CONDITIONAL_FASTJETLIBS@
        
INCLUDES      =	@CONDITIONAL_FASTJETINCS@
