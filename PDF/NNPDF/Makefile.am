pkglib_LTLIBRARIES    = \
			libNNPDFSherpa.la

GITTAG = PDF/NNPDF
include ../../git.make

libNNPDFSherpa_la_SOURCES      =	\
				Git_Info.C \
				NNPDFDriver.cc \
				PDF_NN.C

localinc_HEADERS      = \
			NNPDFDriver.h

NNPDF3archive = nnpdf30.tar.bz2


nobase_dist_pkgdata_DATA         =	$(NNPDF3archive)

install-data-hook:
	tar xjf $(srcdir)/$(NNPDF3archive) -C $(DESTDIR)$(pkgdatadir)

