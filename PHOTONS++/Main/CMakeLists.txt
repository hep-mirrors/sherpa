
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
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(PhotonsMain_static STATIC ${PhotonsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(PhotonsMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME PhotonsMain)
endif()
if (SHERPA_BUILD_SHARED)
add_library(PhotonsMain SHARED ${PhotonsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(PhotonsMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(PhotonsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhotonsMain SOVERSION 0.0.0)
endif()
