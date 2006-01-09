SHERPASYS = ../..

ATOOLSFLAGS   =	# -pedantic -Wall -Winline
DEFS          =

ATOOLSINCS    =	-I../Org -I../Phys -I../Math
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main 

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libToolsMHV.la
libToolsMHV_la_SOURCES        =	\
				Mom.C \
				MHV_Calculator.C \
				PermStore.C
pkginclude_HEADERS       =	\
			Mom.H \
			MHV_Calculator.H \
			PermStore.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(ATOOLSFLAGS)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS)
