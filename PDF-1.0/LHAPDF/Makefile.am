lib_LTLIBRARIES = libLHAPDF.la
libLHAPDF_la_SOURCES = LHAPDF_Fortran_Interface.C lhapdf_inter.f \
		       LHpdflib.f description.f \
		       evolution.f alphas.f  \
		       parameter.f inputPDF.f \
		       wrapQCDNUM.f qcdnum.f Sqcdnum.f
INCLUDES = -I../../PDF-1.0/Main \
	   -I../../ATOOLS-2.0/Phys \
	   -I../../ATOOLS-2.0/Org \
	   -I../../ATOOLS-2.0/Math \
	   -I../../MODEL-1.0/Main
noinst_HEADERS = LHAPDF_Fortran_Interface.H
EXTRA_DIST = parmsetup.f