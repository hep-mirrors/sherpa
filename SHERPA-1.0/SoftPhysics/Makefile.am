SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall #-Winline

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../Single_Events -I../PerturbativePhysics -I../LundTools -I../Tools 
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFINCS       =	-I$(SHERPASYS)/PDF-1.0/Main
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libSherpaSoftPhysics.la
libSherpaSoftPhysics_la_SOURCES	      =	\
					Primordial_KPerp.C \
					Remnant_Base.C \
					No_Remnant.C \
					Color_Dipole.C \
					QCD_Remnant_Base.C \
					Hadron_Remnant.C \
					Electron_Remnant.C \
					Photon_Remnant.C \
					Beam_Remnant_Handler.C \
				   	Fragmentation_Handler.C \
					Hadron_Decay_Handler.C
noinst_HEADERS        =	\
			Primordial_KPerp.H \
			Primordial_KPerp.IC \
			Remnant_Base.H \
			No_Remnant.H \
			Color_Dipole.H \
			QCD_Remnant_Base.H \
			Hadron_Remnant.H \
			Electron_Remnant.H \
			Photon_Remnant.H \
			Beam_Remnant_Handler.H \
		 	Fragmentation_Handler.H \
			Hadron_Decay_Handler.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS)
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) $(EXTRAXSINCS) \
		$(PHASICINCS) $(MODELINCS) $(PDFINCS) $(BEAMINCS) $(ATOOLSINCS)
