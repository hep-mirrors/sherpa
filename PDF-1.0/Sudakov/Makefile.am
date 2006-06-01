pkglib_LTLIBRARIES       =	libSudakov.la

libSudakov_la_SOURCES =	\
			NLL_Branching_Probability_Base.C \
			NLL_Branching_Probabilities.C \
			NLL_Sudakov_Base.C \
			NLL_Single_Sudakov.C \
			NLL_Combined_Sudakov.C \
			NLL_Sudakov.C \
			NLL_JetRate.C

pkginclude_HEADERS        =	\
			NLL_Branching_Probability_Base.H \
			NLL_Branching_Probabilities.H \
			NLL_Sudakov_Base.H \
			NLL_Single_Sudakov.H \
			NLL_Combined_Sudakov.H \
			NLL_Sudakov.H \
			NLL_JetRate.H

AM_CXXFLAGS   =	@SHERPAFLAGS@

INCLUDES      =	@APACICINCS@ @PDFINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
