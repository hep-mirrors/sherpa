##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  REMNANTS/Tools)

set(Tools_SOURCES Kinematics_Generator.C
                  Colour_Generator.C
                  Beam_Decorrelator.C
                  Primordial_KPerp.C
                  Form_Factor.C)

set(Tools_HEADERS Kinematics_Generator.H
                  Colour_Generator.H
                  Beam_Decorrelator.H
                  Primordial_KPerp.H
                  Form_Factor.H)

sherpa_create_git_info(Tools)
 
add_library(RemnantsTools SHARED ${Tools_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(RemnantsTools RemnantsTools)
