lib_LTLIBRARIES = libZfunctions.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libZfunctions_la_SOURCES = Basic_Func.C Basic_Sfuncs.C Basic_Pfunc.C \
			   Basic_Mfunc.C Basic_MassTermfunc.C \
			   Basic_Vfunc.C Basic_Xfunc.C Basic_Yfunc.C \
			   Basic_Zfunc.C \
			   G4_Calc.C V4_Calc.C \
			   VVV_Calc.C ZLED_Calc.C\
			   Zfunc_Calc.C \
			   Calculator.C
INCLUDES = -I../../Amplitude \
	   -I../../Amplitude/AmplTools \
	   -I../../Amplitude/Zfunctions \
	   -I../../Main \
	   -I../../Model \
	   -I../../../MODEL-1.0/Main \
	   -I../../String \
	   -I../../../ATOOLS-2.0/Phys \
	   -I../../../ATOOLS-2.0/Math \
	   -I../../../ATOOLS-2.0/Org 
noinst_HEADERS = Basic_Func.H  Basic_Sfuncs.H Basic_Func.icc Zfunc_Calc.H Calculator.H

EXTRA_DIST    = $(CVSINCLUDE)
