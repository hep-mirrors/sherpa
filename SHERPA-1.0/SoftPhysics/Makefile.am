pkglib_LTLIBRARIES       =	libSherpaSoftPhysics.la
libSherpaSoftPhysics_la_SOURCES	      =	\
					Primordial_KPerp.C \
					Beam_Remnant_Handler.C \
				   	Fragmentation_Handler.C \
					Hadron_Decay_Handler.C
pkginclude_HEADERS    =	\
			Primordial_KPerp.H \
			Beam_Remnant_Handler.H \
		 	Fragmentation_Handler.H \
			Hadron_Decay_Handler.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @SHERPAINCS@ @CONDITIONAL_HADRONSINCS@ @APACICINCS@ @AMEGICINCS@ \
		@EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ @PDFINCS@ \
		@BEAMINCS@ @ATOOLSINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
