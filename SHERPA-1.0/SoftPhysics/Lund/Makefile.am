SHERPASYS     = ../../..

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libLundInterface.la
libLundInterface_la_SOURCES   =	\
				Lund_Fortran_Interface.C \
			      	pythia6214.f \
				pydata.f \
				fhepinter.f \
				finterf.f   
noinst_HEADERS        =	Lund_Fortran_Interface.H 
EXTRA_DIST    =	$(CVSINCLUDE)
INCLUDES      = $(MODELINCS) $(ATOOLSINCS)
