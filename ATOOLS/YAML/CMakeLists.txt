
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am
# yaml-cpp/Makefile.am
# yaml-cpp/node/Makefile.am
# yaml-cpp/node/detail/Makefile.am

########################################################################
#YAML/Makefile.am
 
 #SUBDIRS = yaml-cpp)
 #pkglib_LTLIBRARIES = libToolsYaml.la)
 #EXTRA_DIST = yaml-cpp)
set(YAML_SOURCES binary.cpp convert.cpp depthguard.cpp directives.cpp emit.cpp emitfromevents.cpp emitter.cpp emitterstate.cpp emitterutils.cpp exceptions.cpp exp.cpp memory.cpp node.cpp node_data.cpp
 nodebuilder.cpp nodeevents.cpp null.cpp ostream_wrapper.cpp parse.cpp parser.cpp regex_yaml.cpp scanner.cpp scanscalar.cpp scantag.cpp scantoken.cpp simplekey.cpp singledocparser.cpp stream.cpp tag.cpp)
set(YAML_HEADERS collectionstack.h directives.h emitterstate.h emitterutils.h exp.h indentation.h nodebuilder.h nodeevents.h ptr_vector.h regex_yaml.h regeximpl.h scanner.h scanscalar.h scantag.h setting.h
 singledocparser.h stream.h streamcharsource.h stringsource.h tag.h token.h)
creategitinfo(YAML ../YAML/)
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
 
set(YAML_esources )
if (SHERPA_BUILD_SHARED)
add_library(ToolsYaml SHARED ${ToolsYaml_esources}
                             ${YAML_SOURCES} 
)
target_include_directories(ToolsYaml PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ToolsYaml PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ToolsYaml DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ToolsYaml        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ToolsYaml SOVERSION 0.0.0)
set_target_properties(ToolsYaml        PROPERTIES DEFINE_SYMBOL "")
endif()
