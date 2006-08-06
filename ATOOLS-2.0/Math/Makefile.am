pkglib_LTLIBRARIES       =	libToolsMath.la

if ROOT_SUPPORT
Root_Histo_Source  = Root_Histogram.C 
Root_Histo_Include = Root_Histogram.H 
endif

libToolsMath_la_SOURCES       =	\
				MathTools.C \
				Random.C \
				Function_Base.C \
				Poincare.C \
				Vector.C \
				ComplexVector.C \
				Matrix.C \
				Fast_Function.C \
				Gauss_Integrator.C \
				Scaling.C \
				Axis.C \
				Variable.C \
				Permutation.C \
				PermStore.C \
				Grid.C \
				Primitive_Integrator.C \
				Algebra_Interpreter.C \
				Histogram.C \
				$(Root_Histo_Source) 

pkginclude_HEADERS        =	\
			MathTools.H \
			Random.H \
			Function_Base.H \
			Poincare.H \
			Vector.H \
			ComplexVector.H \
			MyComplex.H \
			Matrix.H \
			Fast_Function.H \
			Gauss_Integrator.H \
			Variable.H \
			Scaling.H \
			Axis.H \
			Permutation.H \
			PermStore.H \
			Grid.H \
			Primitive_Integrator.H \
			Algebra_Interpreter.H \
			Cluster_Algorithm.H \
			Cluster_Algorithm.C \
			Histogram_Base.H \
			Histogram.H \
			$(Root_Histo_Include) 

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @MODELINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
