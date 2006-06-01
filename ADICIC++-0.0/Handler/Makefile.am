pkglib_LTLIBRARIES = libAdicicHdl.la

libAdicicHdl_la_SOURCES = Dipole_Handler.C \
                          Chain_Handler.C \
                          Cascade_Handler.C

pkginclude_HEADERS = Dipole_Handler.H \
                  Chain_Handler.H \
                  Cascade_Handler.H

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @ADICICINCS@ @ATOOLSINCS@ @MODELINCS@ @PDFINCS@

EXTRA_DIST = @SVNINCLUDE@ \
             Dipole_Handler.inl.hh \
             Chain_Handler.inl.hh \
             Cascade_Handler.inl.hh
