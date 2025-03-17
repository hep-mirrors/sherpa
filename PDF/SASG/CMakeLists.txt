set(SASG_SOURCES sasgam2.f
                 SASGph_Fortran_Interface.C)

set(SASG_HEADERS SASGph_Fortran_Interface.H)

sherpa_create_git_info(SASG ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SASGSherpa SHARED ${SASG_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SASGSherpa SASGSherpa)
add_dependencies(SASGSherpa SASG_Git_Info)
