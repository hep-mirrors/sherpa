SHERPASYS = ../..

ATOOLSFLAGS   =	-pedantic -Wall -Winline

ATOOLSINCS    =	-I../Org -I../Phys
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if MODEL_INCLUDE
MODELFLAGS    =
else
MODELFLAGS    = -DUSING__ATOOLS_only
endif

lib_LTLIBRARIES       =	libToolsMath.la
libToolsMath_la_SOURCES       =	\
				MathTools.C \
				Random.C \
				Function_Base.C \
				Poincare.C \
				Vector.C \
				Matrix.C \
				Fast_Function.C \
				Gauss_Integrator.C \
				Scaling.C \
				Axis.C \
				Histogram.C \
				Variable.C \
				Permutation.C \
				Grid.C \
				Primitive_Integrator.C \
				Algebra_Interpreter.C
noinst_HEADERS        =	\
			MathTools.H \
			Random.H \
			Function_Base.H \
			Poincare.H \
			Vector.H \
			MyComplex.H \
			Matrix.H \
			Fast_Function.H \
			Gauss_Integrator.H \
			Histogram.H \
			Variable.H \
			Scaling.H \
			Axis.H \
			Histogram.H \
			Permutation.H \
			Grid.H \
			Primitive_Integrator.H \
			Algebra_Interpreter.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(MODELFLAGS) $(ATOOLSFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS)
