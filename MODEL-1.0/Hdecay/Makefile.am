pkglib_LTLIBRARIES = libHdecay.la

libHdecay_la_SOURCES = Hdecay_Fortran_Interface.C hdecay.f \
		       haber.f feynhiggs.f hdecay_inter.f 

noinst_HEADERS = Hdecay_Fortran_Interface.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES = @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@ hdecay.in
