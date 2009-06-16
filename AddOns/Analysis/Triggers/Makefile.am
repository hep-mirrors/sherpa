pkglib_LTLIBRARIES	=	libSherpaAnalysisTrigger.la

libSherpaAnalysisTrigger_la_SOURCES   =	\
					Trigger_Base.C \
					One_Variable_Selector.C \
					One_Particle_Selector.C \
					Two_Particle_Selector.C \
					All_Particle_Selector.C \
					Four_Particle_Calculator.C \
					Durham_Algorithm.C \
					Cambridge_Algorithm.C \
					Kt_Algorithm.C \
					DIS_Algorithm.C \
					Final_Selector.C \
					Midpoint_Cone.C \
				    	Calorimeter_Cone.C \
					Field_Stuart_Cone_Finder.C \
	 				D0_BFKL_Cone_Finder.C \
					Fastjet_Interface.C \
					One_Particle_Extractor.C \
					Primitive_Detector.C \
					Primitive_Calorimeter.C \
					Tau_Reconstruction.C 

noinst_HEADERS        =	\
			Trigger_Base.H \
			Durham_Algorithm.H \
			Cambridge_Algorithm.H \
			Kt_Algorithm.H \
			DIS_Algorithm.H \
			Final_Selector.H \
			Midpoint_Cone.H \
		    	Calorimeter_Cone.H \
			Primitive_Detector.H \
			Primitive_Calorimeter.H

libSherpaAnalysisTrigger_la_LIBADD = \
	@CONDITIONAL_FASTJETLIBS@

INCLUDES = @CONDITIONAL_FASTJETINCS@