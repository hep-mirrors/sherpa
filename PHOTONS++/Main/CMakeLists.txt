
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhotonsMain.la)
 set(GITTAG  PHOTONS++/Main)
 #include ../../git.make)
set(Main_SOURCES  Photons.C Define_Dipole.C Dipole_Type.C)
set(Main_HEADERS Photons.H Dipole_Type.H Define_Dipole.H)
creategitinfo(Main ../Main/)
list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
list(APPEND Main_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Main_esources )
if (SHERPA_BUILD_SHARED)
add_library(PhotonsMain SHARED ${PhotonsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(PhotonsMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhotonsMain SOVERSION 0.0.0)
set_target_properties(PhotonsMain        PROPERTIES DEFINE_SYMBOL "")
endif()
