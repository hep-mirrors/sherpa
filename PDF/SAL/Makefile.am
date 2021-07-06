pkglib_LTLIBRARIES    =	libSALSherpa.la

GITTAG = PDF/SAL
include ../../git.make

libSALSherpa_la_SOURCES       =	\
				Git_Info.C \
				sal.cpp \
				SALph_Fortran_Interface.C 

localinc_HEADERS      =	\
			SALph_Fortran_Interface.H \
			sal.h \
			math_ut.h

libSALSherpa_la_LIBADD = @PDFLIBS@ @ATOOLSLIBS@