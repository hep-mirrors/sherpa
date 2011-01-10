pkglib_LTLIBRARIES       =	libSherpaSoftPhysics.la

SVNTAG = SHERPA__SoftPhysics
include ../../svn.make

libSherpaSoftPhysics_la_SOURCES	      =	\
					SVN_Info.C \
					Primordial_KPerp.C \
					Beam_Remnant_Handler.C \
				   	Fragmentation_Handler.C \
					Hadron_Decay_Handler.C \
					Soft_Photon_Handler.C

localinc_HEADERS      =	\
			SVN_Info.C.in \
			Primordial_KPerp.H \
			Beam_Remnant_Handler.H \
		 	Fragmentation_Handler.H \
			Hadron_Decay_Handler.H \
			Soft_Photon_Handler.H

libSherpaSoftPhysics_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaSoftPhysics_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
