pkglib_LTLIBRARIES = libSherpaHEFT.la

GITTAG = MODEL/HEFT
include ../../git.make

libSherpaHEFT_la_SOURCES = \
	Git_Info.C \
	Effective_Higgs_Coupling.C \
	Model.C \
	SVV_LC.C \
	SVVV_LC.C \
	SVVP_LC.C \
	SVVVV_LC.C

localinc_HEADERS = \
	Effective_Higgs_Coupling.H

libSherpaHEFT_la_LIBADD = @MODELMAINLIB@ @ATOOLSLIBS@
