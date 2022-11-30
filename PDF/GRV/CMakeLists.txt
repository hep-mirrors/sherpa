set(GRV_SOURCES grvphoton.f
                GRVph_Fortran_Interface.C)

set(GRV_HEADERS GRVph_Fortran_Interface.H)

sherpa_create_git_info(GRV)
 
add_library(GRVSherpa SHARED ${GRV_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(GRVSherpa GRVSherpa)
