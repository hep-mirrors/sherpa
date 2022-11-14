##### Automatically or semiautomaticaly generated, do not edit.#########
#SASG/Makefile.am
 
set(GITTAG  PDF/SASG)
set(SASG_SOURCES  sasgam2.f SASGph_Fortran_Interface.C)
set(SASG_HEADERS SASGph_Fortran_Interface.H)
sherpa_create_git_info(SASG ../SASG/)
list(TRANSFORM SASG_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SASG/")
list(APPEND SASG_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SASGSherpa SHARED ${SASG_SOURCES} )
sherpa_handle_shared_library(SASGSherpa SASGSherpa)
