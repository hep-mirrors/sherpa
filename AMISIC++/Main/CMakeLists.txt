
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libAmisic.la)
 set(GITTAG  AMISIC/Main)
 #include ../../git.make)
set(Main_SOURCES  Amisic.C)
set(Main_HEADERS Amisic.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(Amisic_static STATIC ${Amisic_esources}
                             ${Main_SOURCES} 
)
target_include_directories(Amisic_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amisic_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Amisic_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amisic_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME Amisic)
endif()
if (SHERPA_BUILD_SHARED)
add_library(Amisic SHARED ${Amisic_esources}
                             ${Main_SOURCES} 
)
target_include_directories(Amisic PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amisic PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS Amisic DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amisic        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Amisic SOVERSION 0.0.0)
endif()
