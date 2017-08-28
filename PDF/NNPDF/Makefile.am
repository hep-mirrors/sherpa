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

NNPDF3archive = nnpdf30.tar.bz2


nobase_dist_pkgdata_DATA         =	$(NNPDF3archive)

install-data-hook:
	tar xjf $(srcdir)/$(NNPDF3archive) -C $(DESTDIR)$(pkgdatadir)

uninstall-hook:
	-rm -rf $(DESTDIR)/$(pkgdatadir)/NNPDF30_nnlo_as_0118
	-rm -rf $(DESTDIR)/$(pkgdatadir)/NNPDF30_nlo_as_0118
