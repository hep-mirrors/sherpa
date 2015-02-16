pkglib_LTLIBRARIES    = \
			libNNPDFSherpa.la

SVNTAG = PDF/NNPDF
include ../../svn.make

libNNPDFSherpa_la_SOURCES      =	\
				SVN_Info.C \
				NNPDFDriver.cc \
				PDF_NN.C

localinc_HEADERS      = \
			NNPDFDriver.h

NNPDF3archive = nnpdf30.tar.gz


nobase_dist_pkgdata_DATA         =	$(NNPDF3archive)

install-data-hook:
	tar xzf $(srcdir)/$(NNPDF3archive) -C $(DESTDIR)$(pkgdatadir)

