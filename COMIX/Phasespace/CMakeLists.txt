set(ComixPhasespace_SOURCES PS_Info.C
                       PS_Vertex.C
                       PS_Current.C
                       PS_Generator.C
                       PS_Channel.C)

set(ComixPhasespace_HEADERS PS_Info.H
                       PS_Vertex.H
                       PS_Current.H
                       PS_Generator.H
                       PS_Channel.H)

sherpa_create_git_info(ComixPhasespace ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
 
add_library(ComixPhasespace SHARED ${ComixPhasespace_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ComixPhasespace ComixPhasespace)
add_dependencies(ComixPhasespace ComixPhasespace_Git_Info)
if (SHERPA_ENABLE_THREADING)
  target_link_libraries(ComixPhasespace PRIVATE Threads::Threads)
endif()
