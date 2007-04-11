pkglib_LTLIBRARIES	=	libAnalysisDetector.la

libAnalysisDetector_la_SOURCES	=	\
				Primitive_Detector.C \
				Primitive_Detector_Element.C \
				Primitive_Tracker.C \
				Primitive_ElMag_Calorimeter.C \
				Primitive_Hadron_Calorimeter.C \
				Primitive_Muon_Chambers.C \
				Object_Definition_Base.C \
				Photon_Definition.C

pkginclude_HEADERS        	=	\
				Primitive_Detector.H \
				Primitive_Detector_Element.H \
				Primitive_Tracker.H \
				Primitive_ElMag_Calorimeter.H \
				Primitive_Hadron_Calorimeter.H \
				Primitive_Muon_Chambers.H \
				Object_Definition_Base.H \
				Photon_Definition.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
