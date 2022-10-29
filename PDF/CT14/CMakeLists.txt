
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#CT14/Makefile.am
 
 #pkglib_LTLIBRARIES = libCT14Sherpa.la)
 set(GITTAG  PDF/CT14)
 #include ../../git.make)
set(CT14_SOURCES  CT14Pdf.f CT14_Fortran_Interface.C)
#nobase_#dist_pkgdata_DATA = CT14Grid/ct14n.00.pds CT14Grid/ct14nn.00.pds
creategitinfo(CT14 ../CT14/)
  list(TRANSFORM CT14_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../CT14/")
 
set(CT14_esources )
if (SHERPA_BUILD_STATIC)
add_library(CT14Sherpa_static STATIC ${CT14Sherpa_esources}
                             ${CT14_SOURCES} 
)
target_include_directories(CT14Sherpa_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CT14Sherpa_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CT14Sherpa_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CT14Sherpa_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME CT14Sherpa)
set_target_properties(CT14Sherpa_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(CT14Sherpa SHARED ${CT14Sherpa_esources}
                             ${CT14_SOURCES} 
)
target_include_directories(CT14Sherpa PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(CT14Sherpa PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS CT14Sherpa DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(CT14Sherpa        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME CT14Sherpa SOVERSION 0.0.0)
set_target_properties(CT14Sherpa        PROPERTIES DEFINE_SYMBOL "")
endif()
