##### Automatically or semiautomaticaly generated, do not edit.#########
#LHAPDF/Makefile.am
 
set(GITTAG  PDF/LHAPDF)
set(LHAPDF_SOURCES  LHAPDF_CPP_Interface.C)
sherpa_create_git_info(LHAPDF ../LHAPDF/)
list(TRANSFORM LHAPDF_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LHAPDF/")
list(APPEND LHAPDF_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(LHAPDFSherpa SHARED ${LHAPDF_SOURCES} )
sherpa_handle_shared_library(LHAPDFSherpa LHAPDFSherpa)

target_link_libraries(LHAPDFSherpa PRIVATE LHAPDF::LHAPDF)
