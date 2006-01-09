SHERPASYS     =	../..

SHERPAFLAGS   =	# -pedantic -Wall -Winline 
DEFS          =

ATOOLSDIR     =	$(SHERPASYS)/ATOOLS-2.0
EXTRAXSDIR    =	$(SHERPASYS)/EXTRA_XS-1.0
AMEGICDIR     =	$(SHERPASYS)/AMEGIC++-2.0
APACICDIR     =	$(SHERPASYS)/APACIC++-2.0
SHERPAINCS    =	-I../SingleEvents -I../SoftPhysics -I../Tools
ATOOLSINCS    =	-I$(ATOOLSDIR)/Phys -I$(ATOOLSDIR)/Math -I$(ATOOLSDIR)/Org
PDFDIR        =	-I$(SHERPASYS)/PDF-1.0
PDFINCS       = $(PDFDIR)/Main $(PDFDIR)/Sudakov
BEAMINCS      =	-I$(SHERPASYS)/BEAM-1.0/Main
PHASICINCS    =	-I$(SHERPASYS)/PHASIC++-1.0/Main
MODELINCS     =	-I$(SHERPASYS)/MODEL-1.0/Main
EXTRAXSINCS   =	-I$(EXTRAXSDIR)/Two2Two -I$(EXTRAXSDIR)/Main
AMEGICINCS    =	-I$(AMEGICDIR)/Main -I$(AMEGICDIR)/Amplitude/Zfunctions \
		-I$(AMEGICDIR)/Amplitude/AmplTools -I$(AMEGICDIR)/Amplitude \
		-I$(AMEGICDIR)/Phasespace -I$(AMEGICDIR)/Model \
		-I$(AMEGICDIR)/String 
APACICINCS    = -I$(APACICDIR)/Main -I$(APACICDIR)/Showers

if AMISIC_SUPPORT
AMISICDIR     =	$(SHERPASYS)/AMISIC++-1.0
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
endif

if CVS_INCLUDE
CVSINCLUDE    =	CVS
endif

lib_LTLIBRARIES       =	libSherpaPerturbativePhysics.la
libSherpaPerturbativePhysics_la_SOURCES       =	\
						Matrix_Element_Handler.C \
						MI_Handler.C \
					   	Shower_Handler.C \
						Hard_Decay_Handler.C \
						Interface_Tools.C \
					 	Perturbative_Interface.C \
					   	SimpleXS_Apacic_Interface.C \
						Tree_Filler.C \
						Combine_Table_Base.C \
						Combine_Table_CKKW.C \
					   	Cluster_Partons_Base.C \
					   	Cluster_Partons_CKKW.C \
						Amegic_Apacic_Interface.C \
					   	$(ADICICCFILES) \
					   	$(CSSCFILES)
pkginclude_HEADERS        =	\
			Matrix_Element_Handler.H \
			MI_Handler.H \
		 	Shower_Handler.H \
			Hard_Decay_Handler.H \
		 	Perturbative_Interface.H \
		 	Interface_Tools.H \
			SimpleXS_Apacic_Interface.H \
			SimpleXS_Adicic_Interface.H \
			SimpleXS_CSS_Interface.H \
			Tree_Filler.H \
			Combine_Table_Base.H \
			Combine_Table_CKKW.H \
		   	Cluster_Partons_Base.H \
		   	Cluster_Partons_CKKW.H \
			Amegic_Apacic_Interface.H
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) 
INCLUDES      = $(SHERPAINCS) $(APACICINCS) $(AMEGICINCS) \
		$(AMISICINCS) $(EXTRAXSINCS) $(PHASICINCS) \
		$(MODELINCS) $(PDFINCS) \
		$(BEAMINCS) $(ATOOLSINCS)


