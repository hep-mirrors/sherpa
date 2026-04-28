set(Ions_SOURCES Ion_Base.C
                 Ion_Propagation.C
                 Ion_Potentials.C
                 Ion_Parameters.C)

set(Ions_HEADERS Ion_Base.H
                 Ion_Propagation.H
                 Ion_Potentials.H
                 Ion_Parameters.H)

sherpa_create_git_info(Ions ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(BeamIons SHARED ${Ions_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(BeamIons BeamIons)
add_dependencies(BeamIons Ions_Git_Info)
