pkglib_LTLIBRARIES       =	libSherpaInitialization.la

SVNTAG = SHERPA__Initialization
include ../../svn.make

libSherpaInitialization_la_SOURCES    =	SVN_Info.C Initialization_Handler.C 

localinc_HEADERS      =        Initialization_Handler.H

libSherpaInitialization_la_CXXFLAGS = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaInitialization_la_CPPFLAGS = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
