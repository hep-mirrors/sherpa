
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#String/Makefile.am
 
 #pkglib_LTLIBRARIES = libString.la)
 set(GITTAG  AMEGIC++/String)
 #include ../../git.make)
set(String_SOURCES  String_Handler.C String_Generator.C String_Tree.C String_Output.C MyString.C String_Library.C)
set(String_HEADERS String_Handler.H String_Generator.H String_Tree.H String_Output.H MyString.H Values.H String_Library.H)
creategitinfo(String ../String/)
list(TRANSFORM String_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../String/")
list(APPEND String_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(String SHARED 
                             ${String_SOURCES} 
)
target_include_directories(String PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(String PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
sherpa_mpi_link_libraries(String)
install(TARGETS String DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(String        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME String SOVERSION 0.0.0)
set_target_properties(String        PROPERTIES DEFINE_SYMBOL "")
