CDBGINCS      =	-I../Tools

pkglib_LTLIBRARIES    =	libExtraXSCDBG.la
libExtraXSCDBG_la_SOURCES    = \
				Vertex.C \
				Current.C \
				Vector_Current.C \
				Tensor_Current.C \
				CDBG_Amplitude.C \
				N_Gluon_CDBG.C \
				XS_QCD_CDBG.C
pkginclude_HEADERS    =	\
			Vertex.H \
			Current.H \
			Vector_Current.H \
			Tensor_Current.H \
			CDBG_Amplitude.H \
			N_Gluon_CDBG.H \
			XS_QCD_CDBG.H

EXTRA_DIST    = @SVNINCLUDE@
AM_CXXFLAGS   =	@SHERPAFLAGS@
INCLUDES      = @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ \
		@PDFINCS@ @BEAMINCS@ @ATOOLSINCS@
