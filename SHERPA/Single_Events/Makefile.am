pkglib_LTLIBRARIES	      =	libSherpaSingleEvents.la

if PYTHIA_SUPPORT
MCINTERFACESOURCES= MC_Interface.C
MCINTERFACEHEADERS= MC_Interface.H
endif

libSherpaSingleEvents_la_SOURCES      =	\
					Event_Handler.C \
					Event_Phase_Handler.C \
				   	Signal_Processes.C \
					Multiple_Interactions.C \
					Jet_Evolution.C \
					Hard_Decays.C \
					Beam_Remnants.C \
					Hadronization.C \
					Hadron_Decays.C \
					Analysis_Phase.C \
					EvtReadin_Phase.C \
					$(MCINTERFACESOURCES)

pkginclude_HEADERS    =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Hard_Decays.H \
			Beam_Remnants.H \
		 	Hadronization.H \
			Hadron_Decays.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			$(MCINTERFACEHEADERS)

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @SHERPAINCS@ @ANALYSISINCS@ @APACICINCS@ \
		@AMEGICINCS@ @CONDITIONAL_AMISICINCS@ \
		@EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ @AHADICINCS@ \
		@PDFINCS@ @BEAMINCS@ @CONDITIONAL_HADRONSINCS@ \
		@CONDITIONAL_PHOTONSINCS@ @ATOOLSINCS@ @HELICITIESINCS@ \
		@CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
