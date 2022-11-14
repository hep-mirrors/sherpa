##### Automatically or semiautomaticaly generated, do not edit.#########
#NNPDF/Makefile.am
 
 #pkglib_LTLIBRARIES = libNNPDFSherpa.la)
 set(GITTAG  PDF/NNPDF)
set(NNPDF_SOURCES  NNPDFDriver.cc PDF_NN.C)
set(NNPDF_HEADERS NNPDFDriver.h)
 set(NNPDF3archive   NNPDF31_nnlo_as_0118_mc.tar.bz2)
 #nobase_#dist_pkgdata_DATA = $(NNPDF3archive)
 #install-data-hook:)
 #tar xjf $(srcdir)/$(NNPDF3archive) -C $(DESTDIR)$(pkgdatadir)
 #uninstall-hook:)
 #-rm -rf $(DESTDIR)/$(pkgdatadir)/NNPDF31_nnlo_as_0118_mc)
sherpa_create_git_info(NNPDF ../NNPDF/)
list(TRANSFORM NNPDF_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../NNPDF/")
list(APPEND NNPDF_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(NNPDFSherpa SHARED ${NNPDF_SOURCES} )
sherpa_handle_shared_library(NNPDFSherpa NNPDFSherpa)
