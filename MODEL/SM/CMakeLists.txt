
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#SM/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaSM.la)
 set(GITTAG  MODEL/SM)
 #include ../../git.make)
set(SM_SOURCES  Model.C)
creategitinfo(SM ../SM/)
  list(TRANSFORM SM_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../SM/")
 
set(SM_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaSM_static STATIC ${SherpaSM_esources}
                             ${SM_SOURCES} 
)
target_include_directories(SherpaSM_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSM_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSM_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSM_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaSM)
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaSM SHARED ${SherpaSM_esources}
                             ${SM_SOURCES} 
)
target_include_directories(SherpaSM PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaSM PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaSM DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaSM        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaSM SOVERSION 0.0.0)
endif()
