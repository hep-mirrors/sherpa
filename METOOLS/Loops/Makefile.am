pkglib_LTLIBRARIES       =	libHelicitiesLoops.la

SVNTAG = METOOLS__Loops
include ../../svn.make

if GOLEM95_SUPPORT
EXTRASOURCES=golem95_wrap.f90
endif

libHelicitiesLoops_la_SOURCES       =	\
                        SVN_Info.C \
			Master_Tadpole.C \
                        Master_Bubble.C \
                        Master_Triangle.C \
                        Master_Box.C \
                        PV_Tadpoles.C \
                        PV_Bubbles.C \
                        PV_Triangles.C \
                        Triangle_2.C \
                        Box_4.C \
                        $(EXTRASOURCES) \
                        Golem95_Wrapper.C

localincdir = $(pkgincludedir)/METOOLS/Loops

localinc_HEADERS      =	\
			Divergence_Array.H \
                        Master_Integrals.H \
                        PV_Integrals.H \
                        Triangle_2.H \
                        Box_4.H \
                        Golem95_Wrapper.H

AM_FCFLAGS = -I@CONDITIONAL_GOLEM95DIR@
