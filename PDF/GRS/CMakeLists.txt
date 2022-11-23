##### Automatically or semiautomaticaly generated, do not edit.#########
#GRS/Makefile.am
 
set(GITTAG  PDF/GRS)
set(GRS_SOURCES  grsg99.f GRSph_Fortran_Interface.C)
set(GRS_HEADERS GRSph_Fortran_Interface.H)
sherpa_create_git_info(GRS ../GRS/)
list(TRANSFORM GRS_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../GRS/")
list(APPEND GRS_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(GRSSherpa SHARED ${GRS_SOURCES} )
sherpa_handle_shared_library(GRSSherpa GRSSherpa)
check_Fortran_compiler_flag(-std=legacy Fortran_COMPILER_SUPPORTS_legacy)
if (Fortran_COMPILER_SUPPORTS_legacy)
  target_compile_options(GRSSherpa PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)
endif()
