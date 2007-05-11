pkglib_LTLIBRARIES	=	libAnalysisTriggers.la

libAnalysisTriggers_la_SOURCES        =	\
					Trigger_Base.C \
					One_Particle_Selector.C \
					Two_Particle_Selector.C \
					All_Particle_Selector.C \
					Four_Particle_Calculator.C \
					Final_Selector.C \
					Midpoint_Cone.C \
					Field_Stuart_Cone_Finder.C \
					One_Particle_Extractor.C 

pkginclude_HEADERS    =	\
			Trigger_Base.H \
			Final_Selector.H \
			Midpoint_Cone.H 

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
