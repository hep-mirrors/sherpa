set(SAL_SOURCES sal.cpp
                SALph_CPP_Interface.C)

set(SAL_HEADERS SALph_CPP_Interface.H
                sal.h
                math_ut.h)

sherpa_create_git_info(SAL ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SALSherpa SHARED ${SAL_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SALSherpa SALSherpa)
add_dependencies(SALSherpa SAL_Git_Info)
