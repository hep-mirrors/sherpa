pkglib_LTLIBRARIES       =	libSherpaAnalysis.la

libSherpaAnalysis_la_SOURCES  =	\
				Analysis_Object.C \
			 	Primitive_Analysis.C \
				Analysis_Handler.C

noinst_HEADERS        =	\
			Analysis_Object.H \
			Primitive_Analysis.H \
			Analysis_Handler.H

libSherpaAnalysis_la_LIBADD   =	-L../Tools \
	-L../Triggers -L../Observables \
	-lSherpaAnalysisTools -lSherpaAnalysisTrigger \
	-lSherpaObservables
