SHERPASYS     =	../..

SHERPAFLAGS   =	# -pedantic -Wall #-Winline
DEFS          =

SHERPAINCS    =	-I../Single_Events 
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS	      = -I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -Wno-long-long
endif

lib_LTLIBRARIES       =	libSherpaTools.la
libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				Input_Output_Handler.C \
			    	HepMC_Interface.C \
                                Event_Reader.C
include_HEADERS        =	\
			Kinematics.H \
			Input_Output_Handler.H \
		 	HepMC_Interface.H \
			Event_Reader.H
EXTRA_DIST    =	$(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(APACICINCS) $(MODELINCS) \
		$(PDFINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)




