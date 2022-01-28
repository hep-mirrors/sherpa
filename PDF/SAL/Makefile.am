pkglib_LTLIBRARIES    =	libSALSherpa.la

AM_CXXFLAGS="-Wno-write-strings"

GITTAG = PDF/SAL
include ../../git.make

libSALSherpa_la_SOURCES       =	\
				Git_Info.C \
				sal.cpp \
				SALph_CPP_Interface.C 

localinc_HEADERS      =	\
			SALph_CPP_Interface.H \
			sal.h \
			math_ut.h

nobase_dist_pkgdata_DATA = \
            SALGrid/SAL.dat