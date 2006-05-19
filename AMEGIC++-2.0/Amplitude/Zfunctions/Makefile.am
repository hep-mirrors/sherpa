pkglib_LTLIBRARIES = libZfunctions.la

libZfunctions_la_SOURCES = Basic_Func.C Basic_Sfuncs.C Basic_Pfunc.C \
			   Basic_Mfunc.C Basic_MassTermfunc.C \
			   Basic_Vfunc.C Basic_Xfunc.C Basic_Yfunc.C \
			   Basic_Zfunc.C \
			   G4_Calc.C V4_Calc.C \
			   AnomalousV3_Calc.C AnomalousV4_Calc.C \
			   VVV_Calc.C ZLED_Calc.C Loops_Calc.C\
			   Zfunc_Calc.C \
			   Calculator.C \
                           MHV_Calculator.C

pkginclude_HEADERS = Basic_Func.H  Basic_Sfuncs.H Basic_Func.icc Zfunc_Calc.H Calculator.H MHV_Calculator.H

INCLUDES = @AMEGICINCS@ @MODELINCS@ @ATOOLSINCS@

EXTRA_DIST    = @SVNINCLUDE@
AM_CXXFLAGS   =	@SHERPAFLAGS@
