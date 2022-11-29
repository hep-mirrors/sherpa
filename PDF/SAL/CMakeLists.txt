##### Automatically or semiautomaticaly generated, do not edit.#########
#SAL/Makefile.am
 
set(GITTAG  PDF/SAL)

set(SAL_SOURCES sal.cpp
                SALph_CPP_Interface.C)

set(SAL_HEADERS SALph_CPP_Interface.H
                sal.h
                math_ut.h)

sherpa_create_git_info(SAL ../SAL/)
list(TRANSFORM SAL_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SAL/")
list(APPEND SAL_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SALSherpa SHARED ${SAL_SOURCES} )
sherpa_handle_shared_library(SALSherpa SALSherpa)
