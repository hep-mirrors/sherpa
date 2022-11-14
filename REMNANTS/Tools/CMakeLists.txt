##### Automatically or semiautomaticaly generated, do not edit.#########
#Tools/Makefile.am
 
set(GITTAG  REMNANTS/Tools)
set(Tools_SOURCES  Kinematics_Generator.C Colour_Generator.C Beam_Decorrelator.C Primordial_KPerp.C Form_Factor.C)
set(Tools_HEADERS Kinematics_Generator.H Colour_Generator.H Beam_Decorrelator.H Primordial_KPerp.H Form_Factor.H)
sherpa_create_git_info(Tools ../Tools/)
list(TRANSFORM Tools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Tools/")
list(APPEND Tools_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(RemnantsTools SHARED ${Tools_SOURCES} )
sherpa_handle_shared_library(RemnantsTools RemnantsTools)
