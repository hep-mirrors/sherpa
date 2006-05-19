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
						Amegic_Apacic_Interface.C

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

EXTRA_DIST    = @SVNINCLUDE@

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      = @SHERPAINCS@ @APACICINCS@ @AMEGICINCS@ \
		@CONDITIONAL_AMISICINCS@ @EXTRAXSINCS@ @PHASICINCS@ \
		@MODELINCS@ @PDFINCS@ \
		@BEAMINCS@ @ATOOLSINCS@
