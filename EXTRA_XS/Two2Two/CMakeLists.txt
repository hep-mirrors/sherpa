
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Two2Two/Makefile.am
 
 #pkglib_LTLIBRARIES = libExtraXS2_2.la)
 set(GITTAG  EXTRA_XS/Two2Two)
 #include ../../git.make)
set(Two2Two_SOURCES  XS_Drell_Yan.C XS_DiPhoton.C XS_QCD.C ee_Y4S_BB.C gg_yy.C yy_bobo.C DMDM_mumu.C)
creategitinfo(Two2Two ../Two2Two/)
  list(TRANSFORM Two2Two_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Two2Two/")
 
set(Two2Two_esources )
if (SHERPA_BUILD_STATIC)
add_library(ExtraXS2_2_static STATIC ${ExtraXS2_2_esources}
                             ${Two2Two_SOURCES} 
)
target_include_directories(ExtraXS2_2_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS2_2_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS2_2_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS2_2_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ExtraXS2_2)
set_target_properties(ExtraXS2_2_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ExtraXS2_2 SHARED ${ExtraXS2_2_esources}
                             ${Two2Two_SOURCES} 
)
target_include_directories(ExtraXS2_2 PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ExtraXS2_2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ExtraXS2_2 DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ExtraXS2_2        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ExtraXS2_2 SOVERSION 0.0.0)
set_target_properties(ExtraXS2_2        PROPERTIES DEFINE_SYMBOL "")
endif()
