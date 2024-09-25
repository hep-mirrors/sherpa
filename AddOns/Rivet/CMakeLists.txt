set(Rivet_SOURCES Rivet_Interface.C)
sherpa_create_git_info(Rivet)

add_library(SherpaRivetAnalysis SHARED ${Rivet_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaRivetAnalysis SherpaRivetAnalysis)
target_link_libraries(SherpaRivetAnalysis PRIVATE Rivet::Rivet HepMC3::HepMC3)
