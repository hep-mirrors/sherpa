##### Automatically or semiautomaticaly generated, do not edit.#########
#GRV/Makefile.am
 
 set(GITTAG  PDF/GRV)
set(GRV_SOURCES  grvphoton.f GRVph_Fortran_Interface.C)
set(GRV_HEADERS GRVph_Fortran_Interface.H)
sherpa_create_git_info(GRV ../GRV/)
list(TRANSFORM GRV_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GRV/")
list(APPEND GRV_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(GRVSherpa SHARED ${GRV_SOURCES} )
sherpa_handle_shared_library(GRVSherpa GRVSherpa)
