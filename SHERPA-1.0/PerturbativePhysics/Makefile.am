SHERPASYS     =	../..

SHERPAFLAGS   =	-pedantic -Wall -Winline 

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

if ADICIC_SUPPORT
ADICICDIR     =	$(SHERPASYS)/ADICIC++-0.0
ADICICCFILES  = SimpleXS_Adicic_Interface.C 
ADICICFLAGS   = -DUSING__Adicic
ADICICINCS    = -I$(ADICICDIR)/Main -I$(ADICICDIR)/Structure \
                -I$(ADICICDIR)/Handler -I$(ADICICDIR)/Tool -I$(ADICICDIR)/Input
endif

if AMISIC_SUPPORT
AMISICDIR     =	$(SHERPASYS)/AMISIC++-1.0
AMISICINCS    =	-I$(SHERPASYS)/AMISIC++-1.0/Main
AMISICFLAGS   =	-DUSING__Amisic
endif
if CSS_SUPPORT
CSSDIR        = $(SHERPASYS)/CS_SHOWER++-0.0
CSSCFILES     = SimpleXS_CSS_Interface.C 
CSSINCS       = -I$(CSSDIR)/Tools -I$(CSSDIR)/Showers -I$(CSSDIR)/Main
CSSFLAGS      =	-DUSING__CSS 
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
					   	Cluster_Partons.C \
					 	Perturbative_Interface.C \
					   	SimpleXS_Apacic_Interface.C \
					   	$(ADICICCFILES) \
					   	$(CSSCFILES) \
						Amegic_Apacic_Interface.C \
						Combine_Table.C 
noinst_HEADERS        =	\
			Matrix_Element_Handler.H \
			MI_Handler.H \
		 	Shower_Handler.H \
			Hard_Decay_Handler.H \
		 	Perturbative_Interface.H \
		 	Interface_Tools.H \
			Cluster_Partons.H \
			SimpleXS_Apacic_Interface.H \
			SimpleXS_Adicic_Interface.H \
			SimpleXS_CSS_Interface.H \
			Amegic_Apacic_Interface.H \
			Combine_Table.H 
EXTRA_DIST    = $(CVSINCLUDE)
AM_CXXFLAGS   =	$(SHERPAFLAGS) $(AMISICFLAGS) $(ADICICFLAGS) $(CSSFLAGS)
INCLUDES      = $(SHERPAINCS) $(ADICICINCS) $(APACICINCS) $(AMEGICINCS) $(AMISICINCS) \
		$(EXTRAXSINCS) $(PHASICINCS) $(MODELINCS) $(PDFINCS) $(CSSINCS) \
		$(BEAMINCS) $(ATOOLSINCS)


