 
set(LHAPDF_SOURCES  LHAPDF_CPP_Interface.C)
sherpa_create_git_info(LHAPDF)
 
add_library(LHAPDFSherpa SHARED ${LHAPDF_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(LHAPDFSherpa LHAPDFSherpa)

target_link_libraries(LHAPDFSherpa PRIVATE LHAPDF::LHAPDF)
