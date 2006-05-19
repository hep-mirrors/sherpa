MAKE = make

pkglib_LTLIBRARIES = libBeam.la

libBeam_la_SOURCES    =	Beam_Base.C \
			Monochromatic.C \
			Laser_Backscattering.C \
			Spectrum_Reader.C \
			Beam_Spectra_Handler.C

pkginclude_HEADERS        =	Beam_Base.H \
		 	Monochromatic.H \
			Laser_Backscattering.H \
		 	Spectrum_Reader.H \
			Beam_Spectra_Handler.H

INCLUDES = @ATOOLSINCS@ @MODELINCS@

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@
