pkglib_LTLIBRARIES       =	libPhasespace.la
libPhasespace_la_SOURCES      =	\
				Regulator_Base.C \
				Regulator.C \
				Integrable_Base.C \
			   	Channel_Basics.C \
				Channel_Elements.C \
				Channel_Elements_KK.C \
				Single_Channel.C \
				Simple_Pole.C \
			   	Resonance.C \
				Threshold.C \
				Leading_Log.C \
				LBS_Compton_Peak.C \
			   	Leading_Log_Z.C \
			   	LL_KPerp.C \
				FSR_Channel.C \
				Rambo.C \
				RamboKK.C \
				Sarge.C \
				VHAAG.C \
				Vegas.C ISR_Vegas.C \
				Channel_Interface.C \
			   	Multi_Channel.C \
				PI_Interface.C \
				Phase_Space_Handler.C \
				Phase_Space_Integrator.C 
pkginclude_HEADERS        =	\
			Regulator_Base.H \
			Regulator.H \
		 	Channel_Basics.H \
			Channel_Elements.H \
			Channel_Elements_KK.H \
			Single_Channel.H \
		 	Simple_Pole.H \
		 	Resonance.H \
			Threshold.H \
			Leading_Log.H \
			LBS_Compton_Peak.H \
		 	Leading_Log_Z.H \
		 	LL_KPerp.H \
			FSR_Channel.H \
			Channel_Interface.H \
			Rambo.H \
			RamboKK.H \
			Sarge.H \
			VHAAG.H \
			Vegas.H ISR_Vegas.H \
			Multi_Channel.H \
		 	Integrable_Base.H \
			PI_Interface.H \
			Phase_Space_Handler.H \
			Phase_Space_Integrator.H 

AM_CXXFLAGS   =	@SHERPAFLAGS@ @CONDITIONAL_ROOTFLAGS@

INCLUDES      =	@PHASICINCS@ @MODELINCS@ @PDFINCS@ @BEAMINCS@ \
		@ATOOLSINCS@ @CONDITIONAL_ROOTINCS@

EXTRA_DIST    = @SVNINCLUDE@
