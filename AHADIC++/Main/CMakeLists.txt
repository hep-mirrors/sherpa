
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libAhadicMain.la)
 set(GITTAG  AHADIC++/Main)
 #include ../../git.make)
set(Main_SOURCES  Ahadic.C)
 #localincdir = $(pkg#includedir)/AHADIC++/Main)
set(Main_HEADERS Ahadic.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(AhadicMain SHARED 
                             ${Main_SOURCES} 
)
target_include_directories(AhadicMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(AhadicMain)
install(TARGETS AhadicMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicMain SOVERSION 0.0.0)
set_target_properties(AhadicMain        PROPERTIES DEFINE_SYMBOL "")
