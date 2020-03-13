pkglib_LTLIBRARIES       =	libSherpaPerturbativePhysics.la

GITTAG = SHERPA/PerturbativePhysics
include ../../git.make

libSherpaPerturbativePhysics_la_SOURCES       =	\
						Git_Info.C \
						Matrix_Element_Handler.C \
						MI_Handler.C \
						Shower_Handler.C \
						Hard_Decay_Handler.C

localinc_HEADERS          =	\
			Matrix_Element_Handler.H \
			MI_Handler.H \
			Shower_Handler.H \
			Hard_Decay_Handler.H

libSherpaPerturbativePhysics_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_GZIPINCS@
