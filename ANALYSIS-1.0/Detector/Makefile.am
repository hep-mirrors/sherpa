pkglib_LTLIBRARIES	=	libAnalysisDetector.la

libAnalysisDetector_la_SOURCES	=	\
					Cell.C \
					Detector_Segment.C \
					Detector_Element.C \
					Tracker.C \
					ElMag_Calorimeter.C \
					Had_Calorimeter.C \
					Muon_Chambers.C \
					Detector.C \
					Particle_Smearer_Base.C \
					Electron_Handler.C \
					Photon_Handler.C \
					Muon_Handler.C \
					Hadron_Handler.C \
					Object_Definition_Base.C \
					Reconstructed_Object.C \
					Electron_Maker.C \
					Muon_Maker.C \
					Photon_Maker.C \
					Jet_Maker.C \
					Simple_Cone.C \
					Primitive_Hadron_Calorimeter.C

pkginclude_HEADERS        	=	\
					Cell.H \
					Detector_Segment.H \
					Detector_Element.H \
					Tracker.H \
					ElMag_Calorimeter.H \
					Had_Calorimeter.H \
					Muon_Chambers.H \
					Detector.H \
					Particle_Smearer_Base.H \
					Electron_Handler.H \
					Photon_Handler.H \
					Muon_Handler.H \
					Hadron_Handler.H \
					Object_Definition_Base.H \
					Reconstructed_Object.H \
					Electron_Maker.H \
					Muon_Maker.H \
					Photon_Maker.H \
					Jet_Maker.H \
					Simple_Cone.H \
					Primitive_Hadron_Calorimeter.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@ANALYSISINCS@ @MODELINCS@ @ATOOLSINCS@ @PHASICINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
