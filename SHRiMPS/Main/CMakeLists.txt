
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsMain.la)
 set(GITTAG  SHRiMPS/Main)
 #include ../../git.make)
set(Main_SOURCES  Shrimps.C Cluster_Algorithm.C Soft_Jet_Criterion.C Hadron_Init.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Main)
set(Main_HEADERS Shrimps.H Cluster_Algorithm.H Soft_Jet_Criterion.H Hadron_Init.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsMain_static STATIC ${ShrimpsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ShrimpsMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsMain)
set_target_properties(ShrimpsMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsMain SHARED ${ShrimpsMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ShrimpsMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsMain SOVERSION 0.0.0)
set_target_properties(ShrimpsMain        PROPERTIES DEFINE_SYMBOL "")
endif()
