set(MODEL_VARIATIONS_SOURCES    Variation_Generator.C
                                Variation_Key.C
                                Variations.C)

set(MODEL_VARIATIONS_HEADERS    Variation_Generator.H
                                Variation_Key.H
                                Variations.H)

sherpa_create_git_info(MODEL_VARIATIONS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ModelVariations SHARED ${MODEL_VARIATIONS_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ModelVariations ModelVariations)
add_dependencies(ModelVariations MODEL_VARIATIONS_Git_Info)