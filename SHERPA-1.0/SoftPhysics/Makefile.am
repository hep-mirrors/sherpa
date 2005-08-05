SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall #-Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../Single_Events -I../PerturbativePhysics \
		-I../LundTools -I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFDIR        =	$(SHERPASYS)/PDF-1.0
PDFINCS       =	-I$(PDFDIR)/Main -I$(PDFDIR)/Remnant
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if HADRONS_SUPPORT
HADRONSDIR    =	$(SHERPASYS)/HADRONS++-0.0
HADRONSINCS   = -I$(HADRONSDIR)/Main -I$(HADRONSDIR)/ME_Library -I$(HADRONSDIR)/PS_Library
HADRONFLAGS   = -DUSING__Hadrons
endif
if AHADIC_SUPPORT
AHADICDIR     =	$(SHERPASYS)/AHADIC++-0.0
AHADICINCS    = -I$(AHADICDIR)/Main -I$(AHADICDIR)/Formation -I$(AHADICDIR)/Decays -I$(AHADICDIR)/Tools
AHADICFLAGS   = -DUSING__Ahadic
endif
if AMISIC_SUPPORT
AMISICDIR     =	$(SHERPASYS)/AMISIC++-1.0
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMISICFLAGS   =	-DUSING__Amisic
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

if ROOT_SUPPORT
ROOTINCS      =	-I$(ROOTSYS)/include
ROOTFLAGS     = -DROOT_SUPPORT -Wno-long-long
endif

lib_LTLIBRARIES       =	libSherpaSoftPhysics.la
libSherpaSoftPhysics_la_SOURCES	      =	\
					Primordial_KPerp.C \
					Beam_Remnant_Handler.C \
				   	Fragmentation_Handler.C \
					Hadron_Decay_Handler.C
noinst_HEADERS        =	\
			Primordial_KPerp.H \
			Beam_Remnant_Handler.H \
		 	Fragmentation_Handler.H \
			Hadron_Decay_Handler.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(HADRONFLAGS) $(ROOTFLAGS) $(AHADICFLAGS)
INCLUDES      = $(SHERPAINCS) $(HADRONSINCS) $(APACICINCS) $(AMEGICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS) $(AHADICINCS) $(ROOTINCS)

