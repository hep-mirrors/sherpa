SHERPASYS = ../..

MODELINC = $(SHERPASYS)/MODEL-1.0

lib_LTLIBRARIES       =	libToolsMath.la

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

libToolsMath_la_SOURCES       =	\
				MathTools.C \
				Random.C \
				Function_Base.C \
				Poincare.C \
				Vector.C \
				Vegas.C \
				Matrix.C \
				Fast_Function.C \
				Gauss_Integrator.C \
				Histogram.C
noinst_HEADERS        =	\
			MathTools.H \
			Random.H \
			Function_Base.H \
			Poincare.H \
			Vector.H \
			MyComplex.H \
			Matrix.H \
			Vegas.H \
			Fast_Function.H \
			Gauss_Integrator.H \
			Histogram.H
INCLUDES      =	-I../Phys \
	   	-I../Org \
	   	-I$(MODELINC)/Main

EXTRA_DIST    = $(CVSINCLUDE)
