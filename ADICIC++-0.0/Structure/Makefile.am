pkglib_LTLIBRARIES = libAdicicStruct.la

libAdicicStruct_la_SOURCES = Dipole.C Chain.C Cascade.C

pkginclude_HEADERS = Dipole.H Chain.H Cascade.H

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @ADICICINCS@ @ATOOLSINCS@ @MODELINCS@ @PDFINCS@

EXTRA_DIST = @SVNINCLUDE@ \
             Dipole.inl.hh Dipole.exa.cc \
             Chain.inl.hh Cascade.inl.hh
