pkglib_LTLIBRARIES = libHdecay.la

if CVS_INCLUDE
CVSINCLUDE    =	.svn
endif

MODELFLAGS    =	# -pedantic -Wall #-Winline
DEFS          =

libHdecay_la_SOURCES = Hdecay_Fortran_Interface.C hdecay.f \
		       haber.f feynhiggs.f hdecay_inter.f 
INCLUDES = -I../../MODEL-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Org \
	   -I../../ATOOLS-2.0/Math 
noinst_HEADERS = Hdecay_Fortran_Interface.H
EXTRA_DIST    = $(CVSINCLUDE) hdecay.in
AM_CXXFLAGS   =	$(MODELFLAGS)
