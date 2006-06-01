if ADICIC_SUPPORT
ADICICCFILES  = SimpleXS_Adicic_Interface.C 
endif
if CSS_SUPPORT
CSSCFILES     = SimpleXS_CSS_Interface.C 
endif

pkglib_LTLIBRARIES       =	libSherpaPerturbativePhysics.la

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

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      = @SHERPAINCS@ @CONDITIONAL_ADICICINCS@ @APACICINCS@ @AMEGICINCS@ \
		@CONDITIONAL_AMISICINCS@ @EXTRAXSINCS@ @PHASICINCS@ \
		@MODELINCS@ @PDFINCS@ @CONDITIONAL_CSSINCS@ \
		@BEAMINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
