pkglib_LTLIBRARIES = libSherpaSMH.la

GITTAG = MODEL/SMH
include ../../git.make

libSherpaSMH_la_SOURCES = \
	Git_Info.C \
	gpl.cpp \
	gggH.cpp \
	qqgH.cpp \
	Model.C \
	SVV_LC.C \
	HGG_LC.C \
	HGGG_LC.C \
	HQQG_LC.C

localinc_HEADERS = \
	gpl.h \
	gggH.h \
	qqgH.h

libSherpaSMH_la_LIBADD = \
	$(CONDITIONAL_GINACLIBS) \
	$(CONDITIONAL_CLNLIBS)
