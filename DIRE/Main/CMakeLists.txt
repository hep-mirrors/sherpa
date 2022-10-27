
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libDireMain.la)
 set(GITTAG  DIRE/Main)
 #include ../../git.make)
set(Main_SOURCES  Jet_Criterion.C Dire.C)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(DireMain_static STATIC ${DireMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(DireMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DireMain)
endif()
if (SHERPA_BUILD_SHARED)
add_library(DireMain SHARED ${DireMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(DireMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DireMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS DireMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DireMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DireMain SOVERSION 0.0.0)
endif()
