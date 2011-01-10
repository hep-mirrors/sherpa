pkglib_LTLIBRARIES    = \
			libMSTW08Sherpa.la

SVNTAG = PDF__MSTW
include ../../svn.make

libMSTW08Sherpa_la_SOURCES      =	\
				SVN_Info.C \
				mstwpdf.cc \
				PDF_MSTW.C

localinc_HEADERS      = \
			mstwpdf.h

nobase_dist_pkgdata_DATA         =	\
			MSTW08Grid/mstw2008lo.00.dat \
			MSTW08Grid/mstw2008nlo.00.dat \
			MSTW08Grid/mstw2008nnlo.00.dat
