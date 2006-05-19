pkglib_LTLIBRARIES       = libAmisicModel.la

libAmisicModel_la_SOURCES     = \
			       	Grid_Creator.C \
				Profile_Function_Base.C \
				Profile_Function.C \
				Simple_Chain.C \
				Reggeon_Trajectory.C \
				Simple_String.C 
pkginclude_HEADERS        = \
			Grid_Creator.H \
			Profile_Function_Base.H \
			Profile_Function.H \
			Simple_Chain.H \
			Reggeon_Trajectory.H \
			Simple_String.H 

INCLUDES      =	@CONDITIONAL_SHERPAINCS@ @AMISICINCS@ @EXTRAXSINCS@ @PHASICINCS@ \
		@PDFINCS@ @BEAMINCS@ @MODELINCS@ @ATOOLSINCS@

AM_CXXFLAGS   =	@SHERPAFLAGS@

EXTRA_DIST    = @SVNINCLUDE@
