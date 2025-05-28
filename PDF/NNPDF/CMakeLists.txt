set(NNPDF_SOURCES NNPDFDriver.cc
                  PDF_NN.C)

set(NNPDF_HEADERS NNPDFDriver.h)

 set(NNPDF3archive   NNPDF31_nnlo_as_0118_mc.tar.bz2)

sherpa_create_git_info(NNPDF ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(NNPDFSherpa SHARED ${NNPDF_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(NNPDFSherpa NNPDFSherpa)
add_dependencies(NNPDFSherpa NNPDF_Git_Info)
