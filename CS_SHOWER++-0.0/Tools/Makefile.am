SHERPASYS  = ../..
ATOOLSDIR  = $(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS = -I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org

CSSFLAGS  =	-pedantic -Wall -Winline 
AM_CXXFLAGS   = $(CSSFLAGS)

if CVS_INCLUDE
CVSINCLUDE    = CVS
endif

lib_LTLIBRARIES = libCSSTools.la
libCSSTools_la_SOURCES = \
	Parton.C Singlet.C

INCLUDES = \
	$(ATOOLSINCS)

noinst_HEADERS = \
	Parton.H Singlet.H

EXTRA_DIST    = $(CVSINCLUDE)
