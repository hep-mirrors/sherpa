pkglib_LTLIBRARIES       =	libAnalysis.la

libAnalysis_la_SOURCES        =	\
			 	Final_Selector.C \
	 		 	Particle_Selector.C \
			 	Primitive_Observable_Base.C \
			 	Primitive_Calorimeter.C \
			 	Primitive_Detector.C \
			 	Primitive_Analysis.C \
				Analysis_Handler.C

pkginclude_HEADERS        =	\
			Final_Selector.H \
	 		Particle_Selector.H \
			Primitive_Observable_Base.H \
			Primitive_Calorimeter.H \
			Primitive_Detector.H \
			Primitive_Analysis.H \
			Analysis_Handler.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
