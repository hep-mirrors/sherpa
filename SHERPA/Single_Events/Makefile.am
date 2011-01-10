pkglib_LTLIBRARIES	      =	libSherpaSingleEvents.la

if PYTHIA_SUPPORT
MCINTERFACESOURCES= MC_Interface.C
MCINTERFACEHEADERS= MC_Interface.H
endif

SVNTAG = SHERPA__Single_Events
include ../../svn.make

libSherpaSingleEvents_la_SOURCES      =	\
					SVN_Info.C \
					Event_Handler.C \
					Event_Phase_Handler.C \
				   	Signal_Processes.C \
					Multiple_Interactions.C \
					Jet_Evolution.C \
					Signal_Process_FS_QED_Correction.C \
					Hard_Decays.C \
					Beam_Remnants.C \
					Hadronization.C \
					Hadron_Decays.C \
					Analysis_Phase.C \
					EvtReadin_Phase.C \
					$(MCINTERFACESOURCES)

localinc_HEADERS      =	\
			Event_Handler.H \
			Event_Phase_Handler.H \
		 	Signal_Processes.H \
			Multiple_Interactions.H \
			Jet_Evolution.H \
			Signal_Process_FS_QED_Correction.H \
			Hard_Decays.H \
			Beam_Remnants.H \
		 	Hadronization.H \
			Hadron_Decays.H \
			Analysis_Phase.H \
			EvtReadin_Phase.H \
			$(MCINTERFACEHEADERS)

libSherpaSingleEvents_la_CXXFLAGS   = $(AM_CXXFLAGS) @CONDITIONAL_ROOTFLAGS@

libSherpaSingleEvents_la_CPPFLAGS   = $(AM_CPPFLAGS) @CONDITIONAL_ROOTINCS@
