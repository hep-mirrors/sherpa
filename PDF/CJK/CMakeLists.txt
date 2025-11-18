set(CJK_SOURCES CJKL.f
        cjkhogrid.f
        cjk1grid1.f
        cjk2grid1.f
        CJKph_Fortran_Interface.C)

set(CJK_HEADERS CJKph_Fortran_Interface.H)

sherpa_create_git_info(CJK ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(CJKSherpa SHARED ${CJK_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CJKSherpa CJKSherpa)
check_Fortran_compiler_flag(-std=legacy Fortran_COMPILER_SUPPORTS_legacy)
if (Fortran_COMPILER_SUPPORTS_legacy)
    target_compile_options(CJKSherpa PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>)
endif()
add_dependencies(CJKSherpa CJK_Git_Info)
