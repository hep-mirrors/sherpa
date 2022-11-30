set(NNPDF_SOURCES NNPDFDriver.cc
                  PDF_NN.C)

set(NNPDF_HEADERS NNPDFDriver.h)

 set(NNPDF3archive   NNPDF31_nnlo_as_0118_mc.tar.bz2)

sherpa_create_git_info(NNPDF)
 
add_library(NNPDFSherpa SHARED ${NNPDF_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NNPDFSherpa NNPDFSherpa)
