MAKE = make 
if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif
lib_LTLIBRARIES = libBeam.la
libBeam_la_SOURCES = Beam_Spectra_Handler.C Monochromatic.C Laser_Backscattering.C \
	             Spectrum_Reader.C
noinst_HEADERS = Beam_Spectra_Handler.H Beam_Base.H \
		 Monochromatic.H Laser_Backscattering.H \
		 Spectrum_Reader.H
INCLUDES = -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Math \
	   -I../../ATOOLS-2.0/Org \
	   -I../../MODEL-1.0/Main
EXTRA_DIST    = $(CVSINCLUDE)
