pkglib_LTLIBRARIES = libHdecay.la

libHdecay_la_SOURCES = Hdecay_Fortran_Interface.C hdecay.f \
		       haber.f feynhiggs.f hdecay_inter.f 

INCLUDES = @MODELINCS@ @ATOOLSINCS@

noinst_HEADERS = Hdecay_Fortran_Interface.H

EXTRA_DIST    = @SVNINCLUDE@ hdecay.in

AM_CXXFLAGS   =	@SHERPAFLAGS@
