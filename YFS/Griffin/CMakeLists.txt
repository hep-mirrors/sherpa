set(Griffin_SOURCES B0.cc 
                C0.cc 
                classes.cc 
                D0.cc 
                delrho.cc 
                deltar.cc 
                EWPOZ2.cc 
                EWPOZ.cc 
                ff.cc 
                li.cc 
                linex.cc 
                xscnnlo.cc )

set(Griffin_HEADERS  classes.h
                    Cplx.h
                    delrho.h
                    deltar.h
                    EWPOZ2.h
                    EWPOZ.h
                    ff0.h
                    ff.h
                    li.h
                    linex.h
                    oneloop.h
                    SMvalG.h
                    SMval.h
                    tools.h
                    xscnnlo.h )


sherpa_create_git_info(Griffin)
 
add_library(YFSGriffin SHARED ${Griffin_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(YFSGriffin YFSGriffin)
