set(YAML_SOURCES binary.cpp
                 convert.cpp
                 depthguard.cpp
                 directives.cpp
                 emit.cpp
                 emitfromevents.cpp
                 emitter.cpp
                 emitterstate.cpp
                 emitterutils.cpp
                 exceptions.cpp
                 exp.cpp
                 memory.cpp
                 node.cpp
                 node_data.cpp
                 nodebuilder.cpp
                 nodeevents.cpp
                 null.cpp
                 ostream_wrapper.cpp
                 parse.cpp
                 parser.cpp
                 regex_yaml.cpp
                 scanner.cpp
                 scanscalar.cpp
                 scantag.cpp
                 scantoken.cpp
                 simplekey.cpp
                 singledocparser.cpp
                 stream.cpp
                 tag.cpp)

set(YAML_HEADERS collectionstack.h
                 directives.h
                 emitterstate.h
                 emitterutils.h
                 exp.h
                 indentation.h
                 nodebuilder.h
                 nodeevents.h
                 ptr_vector.h
                 regex_yaml.h
                 regeximpl.h
                 scanner.h
                 scanscalar.h
                 scantag.h
                 setting.h
                 singledocparser.h
                 stream.h
                 streamcharsource.h
                 stringsource.h
                 tag.h
                 token.h
                 yaml-cpp/anchor.h
                 yaml-cpp/binary.h
                 yaml-cpp/depthguard.h
                 yaml-cpp/dll.h
                 yaml-cpp/emitfromevents.h
                 yaml-cpp/emitter.h
                 yaml-cpp/emitterdef.h
                 yaml-cpp/emittermanip.h
                 yaml-cpp/emitterstyle.h
                 yaml-cpp/eventhandler.h
                 yaml-cpp/exceptions.h
                 yaml-cpp/mark.h
                 yaml-cpp/noexcept.h
                 yaml-cpp/null.h
                 yaml-cpp/ostream_wrapper.h
                 yaml-cpp/parser.h
                 yaml-cpp/stlemitter.h
                 yaml-cpp/traits.h
                 yaml-cpp/yaml.h
                 yaml-cpp/node/convert.h
                 yaml-cpp/node/emit.h
                 yaml-cpp/node/impl.h
                 yaml-cpp/node/iterator.h
                 yaml-cpp/node/node.h
                 yaml-cpp/node/parse.h
                 yaml-cpp/node/ptr.h
                 yaml-cpp/node/type.h
                 yaml-cpp/node/detail/impl.h
                 yaml-cpp/node/detail/iterator_fwd.h
                 yaml-cpp/node/detail/iterator.h
                 yaml-cpp/node/detail/memory.h
                 yaml-cpp/node/detail/node_data.h
                 yaml-cpp/node/detail/node_iterator.h
                 yaml-cpp/node/detail/node_ref.h
                 yaml-cpp/node/detail/node.h)

sherpa_create_git_info(YAML)
 
 
set(YAML_yaml-cpp_HEADERS anchor.h
                          binary.h
                          depthguard.h
                          dll.h
                          emitfromevents.h
                          emitter.h
                          emitterdef.h
                          emittermanip.h
                          emitterstyle.h
                          eventhandler.h
                          exceptions.h
                          mark.h
                          noexcept.h
                          null.h
                          ostream_wrapper.h
                          parser.h
                          stlemitter.h
                          traits.h
                          yaml.h)

 
 
set(YAML_yaml-cpp_node_HEADERS convert.h
                               emit.h
                               impl.h
                               iterator.h
                               node.h
                               parse.h
                               ptr.h
                               type.h)

 
 
set(YAML_yaml-cpp_node_detail_HEADERS impl.h
                                      iterator.h
                                      iterator_fwd.h
                                      memory.h
                                      node.h
                                      node_data.h
                                      node_iterator.h
                                      node_ref.h)

 
add_library(ToolsYaml SHARED ${YAML_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C )
sherpa_handle_shared_library(ToolsYaml ToolsYaml)
