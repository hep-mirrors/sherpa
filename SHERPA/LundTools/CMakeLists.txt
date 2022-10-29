#Manually edited
########################################################################
#
#  Automatically or semiautomaticall generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#LundTools/Makefile.am
 
 #pkglib_LTLIBRARIES = libLundTools.la)
set(PYTHIASOURCES   pyplst.f Lund_Interface.C)
 set(PYTHIAHEADERS  Lund_Wrapper.H Lund_Interface.H)
 ##MD5_EXCLUDE = pythia-6.4.18.f)
 #GITTAG = SHERPA/LundTools)
 #include ../../git.make)
 if (SHERPA_ENABLE_PYTHIA6)
set(LundTools_SOURCES ${PYTHIASOURCES} )
 #localinc_HEADERS = $(PYTHIAHEADERS) Dummy.C)
  list(TRANSFORM LundTools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LundTools/")
  list(APPEND LundTools_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/../LundTools/pythia-6.4.18.f")
 else()
set(LundTools_SOURCES  Dummy.C)
 #localinc_HEADERS = ${PYTHIASOURCES} $(PYTHIAHEADERS)
  list(TRANSFORM LundTools_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../LundTools/")
 endif()
 
set(LundTools_esources )
if (SHERPA_BUILD_STATIC)
add_library(LundTools_static STATIC ${LundTools_esources}
                             ${LundTools_SOURCES} 
)
target_include_directories(LundTools_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(LundTools_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(LundTools_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS LundTools_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(LundTools_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME LundTools)
set_target_properties(LundTools_static       PROPERTIES DEFINE_SYMBOL "")

endif()
if (SHERPA_BUILD_SHARED)
add_library(LundTools SHARED ${LundTools_esources}
                             ${LundTools_SOURCES} 
)
target_include_directories(LundTools PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(LundTools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(LundTools PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS LundTools DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(LundTools        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME LundTools SOVERSION 0.0.0)
set_target_properties(LundTools       PROPERTIES DEFINE_SYMBOL "")

endif()
