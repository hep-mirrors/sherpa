SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline

SHERPAINCS    =	-I../Single_Events 
APACICINCS    =	-I$(SHERPASYS)/APACIC++-2.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
PDFINCS	      = -I$(SHERPASYS)/PDF-1.0/Main
ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org 

if CLHEP_SUPPORT
CLHEPINCS     =	-I$(CLHEPDIR)/include
CLHEPFLAGS    = -DCLHEP_SUPPORT
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libSherpaTools.la
libSherpaTools_la_SOURCES     =	\
				Kinematics.C \
				IO_Handler.C \
			    	HepMC_Interface.C \
                                Event_Reader.C
noinst_HEADERS        =	\
			Kinematics.H \
			IO_Handler.H \
		 	HepMC_Interface.H \
			Event_Reader.H
EXTRA_DIST    =	$(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(CLHEPFLAGS) $(ROOTFLAGS)
INCLUDES      =	$(SHERPAINCS) $(AMISICINCS) $(APACICINCS) $(MODELINCS) \
		$(PDFINCS) $(ATOOLSINCS) $(CLHEPINCS) $(ROOTINCS)




