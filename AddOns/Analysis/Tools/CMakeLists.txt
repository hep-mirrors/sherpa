##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  AddOns/Analysis/Tools)

set(Tools_SOURCES Lorentz.C
                  Breit_Frame.C
                  MET.C
                  Blob_Selector.C
                  List_Creator.C
                  List_Merger.C
                  Particle_Ordering.C
                  Particle_Qualifier.C
                  Particle_Selector.C)

set(Tools_HEADERS Particle_Qualifier.H)

sherpa_create_git_info(Tools)
 
add_library(SherpaAnalysisTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaAnalysisTools SherpaAnalysisTools)
