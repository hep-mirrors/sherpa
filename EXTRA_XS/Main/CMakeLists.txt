
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS.la)
 set(GITTAG  EXTRA_XS/Main)
 #include ../../git.make)
set(Main_SOURCES  Single_Process.C Process_Group.C Simple_XS.C ME_Tools.C ME2_Base.C)
set(Main_HEADERS Single_Process.H Process_Group.H ME_Tools.H ME2_Base.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(ExtraXS_static STATIC ${ExtraXS_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ExtraXS_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ExtraXS)
set_target_properties(ExtraXS_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ExtraXS SHARED ${ExtraXS_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ExtraXS PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS SOVERSION 0.0.0)
set_target_properties(ExtraXS        PROPERTIES DEFINE_SYMBOL "")
endif()
