set(CJK_SOURCES CJKL.f
        cjkhogrid.f
        cjk1grid1.f
        cjk2grid1.f
        CJKph_Fortran_Interface.C)

set(CJK_HEADERS CJKph_Fortran_Interface.H)

sherpa_create_git_info(CJK)

add_library(CJKSherpa SHARED ${CJK_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(CJKSherpa CJKSherpa)
