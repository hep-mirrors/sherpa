SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
SHERPAINCS    =	-I../Single_Events 	
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libLundTools.la
libLundTools_la_SOURCES       =	\
				pythia6214.f \
				pydata.f \
				pyplst.f \
				Lund_Interface.C 
noinst_HEADERS        =	\
		 	Lund_Wrapper.H \
			Lund_Interface.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      =	$(MODELINCS) $(ATOOLSINCS)




