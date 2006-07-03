MAKE    = make

SUBDIRS = savemaketool

pkglib_LTLIBRARIES = libAdicicTool.la

libAdicicTool_la_SOURCES = Xio.C Xmath.C Xhisto.C \
                           Sudakov_Calculator.C \
                           Sudakov_Group.C IISudakov_Group.C \
                           Recoil_Calculator.C

pkginclude_HEADERS = Xenum.H Xio.H Xmath.H Xhisto.H defextra enumextra ioextra mathextra histoextra \
                     Sudakov_Calculator.H \
                     Sudakov_Group.H IISudakov_Group.H \
                     Recoil_Calculator.H

AM_CXXFLAGS = @SHERPAFLAGS@ -ftemplate-depth-100

INCLUDES = @ADICICINCS@ @ATOOLSINCS@ @MODELINCS@ @PDFINCS@

EXTRA_DIST = @SVNINCLUDE@ \
             Sudakov_Calculator.inl.hh \
             Sudakov_Group.inl.hh Sudakov_Group.tpt.cc \
             IISudakov_Group.inl.hh IISudakov_Group.tpt.cc \
             Recoil_Calculator.inl.hh Recoil_Calculator.tpt.cc
