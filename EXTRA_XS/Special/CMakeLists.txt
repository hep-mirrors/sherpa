##### Automatically or semiautomaticaly generated, do not edit.#########
#Special/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXSSpecial.la)
 set(GITTAG  EXTRA_XS/Special)
set(Special_SOURCES  Instanton.C)
creategitinfo(Special ../Special/)
list(TRANSFORM Special_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Special/")
list(APPEND Special_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ExtraXSSpecial SHARED ${Special_SOURCES} )
target_include_directories(ExtraXSSpecial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ExtraXSSpecial)
install(TARGETS ExtraXSSpecial DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXSSpecial        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXSSpecial SOVERSION 0.0.0)
set_target_properties(ExtraXSSpecial        PROPERTIES DEFINE_SYMBOL "")
