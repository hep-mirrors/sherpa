pkglib_LTLIBRARIES = libAdicicInput.la

libAdicicInput_la_SOURCES = Dipole_Flavour.C Information.cpp \
                            Sudakov_Utils.cpp Recoil_Strategy.cpp \
                            Evolution_Strategy.cpp \
                            Dipole_Parameter.C

pkginclude_HEADERS = Dipole_Flavour.H Information.hpp \
                     Sudakov_Utils.hpp Recoil_Strategy.hpp \
                     Evolution_Strategy.hpp \
                     Dipole_Parameter.H

AM_CXXFLAGS = @SHERPAFLAGS@

INCLUDES = @ADICICINCS@ @ATOOLSINCS@ @MODELINCS@

EXTRA_DIST = @SVNINCLUDE@ \
             Recoil_Strategy.inl.hh \
             Dipole_Parameter.inl.hh Dipole_Parameter.prm.cc
