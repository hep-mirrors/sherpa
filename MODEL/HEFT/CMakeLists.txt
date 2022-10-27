
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#HEFT/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaHEFT.la)
 set(GITTAG  MODEL/HEFT)
 #include ../../git.make)
set(HEFT_SOURCES  Effective_Higgs_Coupling.C Model.C SVV_LC.C SVVV_LC.C SVVP_LC.C SVVVV_LC.C)
set(HEFT_HEADERS Effective_Higgs_Coupling.H)
creategitinfo(HEFT ../HEFT/)
  list(TRANSFORM HEFT_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../HEFT/")
 
set(HEFT_esources )
if (SHERPA_BUILD_STATIC)
add_library(SherpaHEFT_static STATIC ${SherpaHEFT_esources}
                             ${HEFT_SOURCES} 
)
target_include_directories(SherpaHEFT_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaHEFT_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaHEFT_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHEFT_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME SherpaHEFT)
endif()
if (SHERPA_BUILD_SHARED)
add_library(SherpaHEFT SHARED ${SherpaHEFT_esources}
                             ${HEFT_SOURCES} 
)
target_include_directories(SherpaHEFT PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaHEFT PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaHEFT DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaHEFT        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaHEFT SOVERSION 0.0.0)
endif()
