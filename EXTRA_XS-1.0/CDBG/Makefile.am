CDBGINCS      =	-I../Tools

pkglib_LTLIBRARIES    =	libExtraXSCDBG.la
libExtraXSCDBG_la_SOURCES     =	\
				Current_Base.C \
				Current.C \
				Vertex.C \
				CDBG_Amplitude.C \
				N_Parton_CDBG.C \
				N_Parton_CDBG_T.C \
				XS_QCD_CDBG.C \
				XS_QCD_CDBG_T.C \
				QCD_Currents.C \
				QCD_Vertices.C
pkginclude_HEADERS    =	\
			Current_Base.H \
			Current.H \
			Vertex.H \
			CDBG_Amplitude.H \
			N_Parton_CDBG.H \
			N_Parton_CDBG_T.H \
			XS_QCD_CDBG.H \
			XS_QCD_CDBG_T.H \
			QCD_Currents.H \
			QCD_Vertices.H

EXTRA_DIST    = @SVNINCLUDE@
AM_CXXFLAGS   =	@SHERPAFLAGS@
INCLUDES      = @EXTRAXSINCS@ @PHASICINCS@ @MODELINCS@ \
		@PDFINCS@ @BEAMINCS@ @ATOOLSINCS@
