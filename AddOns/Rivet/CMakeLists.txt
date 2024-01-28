set(Rivet_SOURCES Rivet_Interface.C)
sherpa_create_git_info(Rivet)

add_library(SherpaRivetAnalysis SHARED ${Rivet_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
sherpa_handle_shared_library(SherpaRivetAnalysis SherpaRivetAnalysis)
target_compile_features(SherpaRivetAnalysis PRIVATE cxx_std_14)
target_link_libraries(SherpaRivetAnalysis PRIVATE ${RIVET_LIBRARIES})
target_include_directories(SherpaRivetAnalysis PRIVATE ${RIVET_INCLUDE_DIR} ${RIVET_CONFIG_CPPFLAGS_DIRS})
target_include_directories(SherpaRivetAnalysis PRIVATE ${HEPMC3_INCLUDE_DIR})
