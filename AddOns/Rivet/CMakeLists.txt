set(Rivet_SOURCES Rivet_Interface.C)
sherpa_create_git_info(Rivet ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(SherpaRivetAnalysis SHARED ${Rivet_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaRivetAnalysis SherpaRivetAnalysis)
target_link_libraries(SherpaRivetAnalysis PRIVATE Rivet::Rivet HepMC3::HepMC3)
add_dependencies(SherpaRivetAnalysis Rivet_Git_Info)
