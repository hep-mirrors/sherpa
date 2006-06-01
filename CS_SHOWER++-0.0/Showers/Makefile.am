AM_CXXFLAGS   = @SHERPAFLAGS@

pkglib_LTLIBRARIES = libCSShowers.la

libCSShowers_la_SOURCES = \
	Splitting_Function_Base.C Splitting_Function_Group.C \
	QCD_Splitting_Functions.C Sudakov.C Shower.C \
	Kinematics_Base.C

noinst_HEADERS = \
	Splitting_Function_Base.H Splitting_Function_Group.H \
	QCD_Splitting_Functions.H Sudakov.H Shower.H \
	Kinematics_Base.H

INCLUDES = @ATOOLSINCS@ @MODELINCS@ @PDFINCS@ @CSSINCS@
