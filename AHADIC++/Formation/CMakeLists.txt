set(Formation_SOURCES Beam_Particles_Shifter.C
                      Singlet_Former.C
                      Singlet_Checker.C
                      Gluon_Decayer.C
                      Gluon_Splitter.C
                      Trivial_Splitter.C)

set(Formation_HEADERS Beam_Particles_Shifter.H
                      Singlet_Former.H
                      Singlet_Checker.H
                      Gluon_Decayer.H
                      Gluon_Splitter.H
                      Trivial_Splitter.H)

sherpa_create_git_info(Formation ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(AhadicFormation SHARED ${Formation_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(AhadicFormation AhadicFormation)
add_dependencies(AhadicFormation Formation_Git_Info)
