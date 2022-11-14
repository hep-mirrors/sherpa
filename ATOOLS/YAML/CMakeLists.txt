##### Automatically or semiautomaticaly generated, do not edit.#########
#YAML/Makefile.am
 
 #SUBDIRS = yaml-cpp)
 #pkglib_LTLIBRARIES = libToolsYaml.la)
 set(GITTAG  ATOOLS/YAML)
 #EXTRA_DIST = yaml-cpp)
set(YAML_SOURCES  binary.cpp convert.cpp depthguard.cpp directives.cpp emit.cpp emitfromevents.cpp emitter.cpp emitterstate.cpp emitterutils.cpp exceptions.cpp exp.cpp memory.cpp node.cpp node_data.cpp
 nodebuilder.cpp nodeevents.cpp null.cpp ostream_wrapper.cpp parse.cpp parser.cpp regex_yaml.cpp scanner.cpp scanscalar.cpp scantag.cpp scantoken.cpp simplekey.cpp singledocparser.cpp stream.cpp tag.cpp)
set(YAML_HEADERS collectionstack.h directives.h emitterstate.h emitterutils.h exp.h indentation.h nodebuilder.h nodeevents.h ptr_vector.h regex_yaml.h regeximpl.h scanner.h scanscalar.h scantag.h setting.h
 singledocparser.h stream.h streamcharsource.h stringsource.h tag.h token.h)
sherpa_create_git_info(YAML ../YAML/)
list(TRANSFORM YAML_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../YAML/")
list(APPEND YAML_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
#YAML/yaml-cpp/Makefile.am
 
 #SUBDIRS = node)
set(YAML_yaml-cpp_HEADERS anchor.h binary.h depthguard.h dll.h emitfromevents.h emitter.h emitterdef.h emittermanip.h emitterstyle.h eventhandler.h exceptions.h mark.h noexcept.h null.h ostream_wrapper.h
 parser.h stlemitter.h traits.h yaml.h)
 
#YAML/yaml-cpp/node/Makefile.am
 
 #SUBDIRS = detail)
set(YAML_yaml-cpp_node_HEADERS convert.h emit.h impl.h iterator.h node.h parse.h ptr.h type.h)
 
#YAML/yaml-cpp/node/detail/Makefile.am
 
set(YAML_yaml-cpp_node_detail_HEADERS impl.h iterator.h iterator_fwd.h memory.h node.h node_data.h node_iterator.h node_ref.h)
 
add_library(ToolsYaml SHARED ${YAML_SOURCES} )
sherpa_handle_shared_library(ToolsYaml ToolsYaml)
