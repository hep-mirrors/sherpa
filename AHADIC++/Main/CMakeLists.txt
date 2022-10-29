
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
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
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(AhadicMain_static STATIC ${AhadicMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(AhadicMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME AhadicMain)
set_target_properties(AhadicMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(AhadicMain SHARED ${AhadicMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(AhadicMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(AhadicMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS AhadicMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(AhadicMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME AhadicMain SOVERSION 0.0.0)
set_target_properties(AhadicMain        PROPERTIES DEFINE_SYMBOL "")
endif()
