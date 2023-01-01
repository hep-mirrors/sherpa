set(GRS_SOURCES grsg99.f
                GRSph_Fortran_Interface.C)

set(GRS_HEADERS GRSph_Fortran_Interface.H)

sherpa_create_git_info(GRS)
 
add_library(GRSSherpa SHARED ${GRS_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(GRSSherpa GRSSherpa)
check_Fortran_compiler_flag(-std=legacy Fortran_COMPILER_SUPPORTS_legacy)
if (Fortran_COMPILER_SUPPORTS_legacy)
  target_compile_options(GRSSherpa PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)
endif()
