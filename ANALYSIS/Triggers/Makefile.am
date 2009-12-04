pkglib_LTLIBRARIES	=	libAnalysisTriggers.la

libAnalysisTriggers_la_SOURCES        =	\
					Trigger_Base.C \
					One_Variable_Selector.C \
					One_Particle_Selector.C \
					Two_Particle_Selector.C \
					All_Particle_Selector.C \
					Four_Particle_Calculator.C \
					Fastjet_Interface.C \
					Final_Selector.C \
					Midpoint_Cone.C \
				    	Calorimeter_Cone.C \
					Field_Stuart_Cone_Finder.C \
	 				D0_BFKL_Cone_Finder.C \
					One_Particle_Extractor.C \
					Primitive_Detector.C \
					Primitive_Calorimeter.C \
					Tau_Reconstruction.C 

pkginclude_HEADERS    =	\
			Trigger_Base.H \
			Final_Selector.H \
			Midpoint_Cone.H \
		    	Calorimeter_Cone.H \
			Primitive_Detector.H \
			Primitive_Calorimeter.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@ @CONDITIONAL_FASTJETINCS@

EXTRA_DIST    = @SVNINCLUDE@

libAnalysisTriggers_la_LIBADD = \
	@CONDITIONAL_FASTJETLIBS@
