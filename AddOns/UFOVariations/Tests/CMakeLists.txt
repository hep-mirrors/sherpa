set(UFOVariations_SOURCES   Model.C)

set(UFOVariations_HEADERS   Model.H)

sherpa_create_git_info(UFOVariationsSMTest ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(SherpaUFOVariationsSMTest SHARED ${UFOVariations_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(SherpaUFOVariationsSMTest SherpaUFOVariationsSMTest)
add_dependencies(SherpaUFOVariationsSMTest UFOVariationsSMTest_Git_Info)