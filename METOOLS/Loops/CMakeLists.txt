##### Automatically or semiautomaticaly generated, do not edit.#########
#Loops/Makefile.am
 
set(GITTAG  METOOLS/Loops)

set(Loops_SOURCES Master_Tadpole.C
                  Master_Bubble.C
                  Master_Triangle.C
                  Master_Box.C
                  PV_Tadpoles.C
                  PV_Bubbles.C
                  PV_Triangles.C)

set(Loops_HEADERS Divergence_Array.H
                  Master_Integrals.H
                  PV_Integrals.H)

sherpa_create_git_info(Loops)
 
add_library(METoolsLoops SHARED ${Loops_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(METoolsLoops METoolsLoops)
