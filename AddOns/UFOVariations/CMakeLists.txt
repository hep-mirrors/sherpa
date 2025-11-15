set(UFOVariations_SOURCES   Variation_Generator.C
                            Combinations.C
                            Variations.C)

set(UFOVariations_HEADERS   Variation_Generator.H
                            UFO_Variation_Key.H
                            Combinations.H
                            Variations.H)

sherpa_create_git_info(UFOVariations ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaUFOVariations SHARED ${UFOVariations_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaUFOVariations SherpaUFOVariations)
add_dependencies(SherpaUFOVariations UFOVariations_Git_Info)