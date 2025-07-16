set(H1_SOURCES i_2006_fita.f
        i_2006_fitb.f
        qcd_2006.f
        H1_Fortran_Interface.C)

set(H1_HEADERS H1_Fortran_Interface.H)

sherpa_create_git_info(H1 ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(H1Sherpa SHARED ${H1_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(H1Sherpa H1Sherpa)
add_dependencies(H1Sherpa H1_Git_Info)
