pkglib_LTLIBRARIES	      =	libSherpaSingleEvents.la

libSherpaSingleEvents_la_SOURCES      =	\
					Event_Handler.C \
					Event_Phase_Handler.C \
				   	Signal_Processes.C \
					Multiple_Interactions.C \
					Jet_Evolution.C \
					Hard_Decays.C \
					Hadronization.C \
					Hadron_Decays.C \
					Analysis_Phase.C \
					EvtReadin_Phase.C \
					MC_Interface.C

pkginclude_HEADERS        =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Hard_Decays.H \
		 	Hadronization.H \
			Hadron_Decays.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			MC_Interface.H

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      = @SHERPAINCS@ @ANALYSISINCS@ @APACICINCS@ \
		@CONDITIONAL_ADICICINCS@ @CONDITIONAL_CSSINCS@ @AMEGICINCS@ @CONDITIONAL_AMISICINCS@ \
		@CONDITIONAL_AHADICINCS@ @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ \
		@PDFINCS@ @BEAMINCS@ @CONDITIONAL_HADRONSINCS@ @ATOOLSINCS@ \
		@CONDITIONAL_ROOTINCS@ @CONDITIONAL_CLHEPINCS@

EXTRA_DIST    = @SVNINCLUDE@ 
