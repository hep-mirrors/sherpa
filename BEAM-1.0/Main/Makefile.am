MAKE = make 
if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

BEAMFLAGS  =	-pedantic -Wall -Winline

lib_LTLIBRARIES = libBeam.la
libBeam_la_SOURCES    =	Beam_Base.C \
			Monochromatic.C \
			Laser_Backscattering.C \
			Spectrum_Reader.C \
			Beam_Spectra_Handler.C
noinst_HEADERS        =	Beam_Base.H \
		 	Monochromatic.H \
			Laser_Backscattering.H \
		 	Spectrum_Reader.H \
			Beam_Spectra_Handler.H
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(BEAMFLAGS)
