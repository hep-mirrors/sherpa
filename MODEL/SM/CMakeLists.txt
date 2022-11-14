##### Automatically or semiautomaticaly generated, do not edit.#########
#SM/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSM.la)
 set(GITTAG  MODEL/SM)
set(SM_SOURCES  Model.C)
creategitinfo(SM ../SM/)
list(TRANSFORM SM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SM/")
list(APPEND SM_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaSM SHARED ${SM_SOURCES} )
target_include_directories(SherpaSM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(SherpaSM)
install(TARGETS SherpaSM DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSM        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSM SOVERSION 0.0.0)
set_target_properties(SherpaSM        PROPERTIES DEFINE_SYMBOL "")
