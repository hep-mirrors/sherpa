
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#DipoleSubtraction/Makefile.am
 
 #pkglib_LTLIBRARIES = libDipoleSubtraction.la)
 set(GITTAG  AMEGIC++/DipoleSubtraction)
 #include ../../git.make)
set(DipoleSubtraction_SOURCES  DipoleSplitting_Base.C FF_DipoleSplitting.C FI_DipoleSplitting.C IF_DipoleSplitting.C II_DipoleSplitting.C Single_LOProcess.C Single_LOProcess_MHV.C Single_LOProcess_External.C
 Single_DipoleTerm.C Single_OSTerm.C Single_Virtual_Correction.C Single_Real_Correction.C)
set(DipoleSubtraction_HEADERS DipoleSplitting_Base.H FF_DipoleSplitting.H FI_DipoleSplitting.H IF_DipoleSplitting.H II_DipoleSplitting.H Single_LOProcess.H Single_LOProcess_MHV.H Single_LOProcess_External.H
 Single_DipoleTerm.H Single_OSTerm.H Single_Virtual_Correction.H Single_Real_Correction.H)
creategitinfo(DipoleSubtraction ../DipoleSubtraction/)
  list(TRANSFORM DipoleSubtraction_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../DipoleSubtraction/")
 
set(DipoleSubtraction_esources )
if (SHERPA_BUILD_STATIC)
add_library(DipoleSubtraction_static STATIC ${DipoleSubtraction_esources}
                             ${DipoleSubtraction_SOURCES} 
)
target_include_directories(DipoleSubtraction_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DipoleSubtraction_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DipoleSubtraction_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS DipoleSubtraction_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DipoleSubtraction_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME DipoleSubtraction)
endif()
if (SHERPA_BUILD_SHARED)
add_library(DipoleSubtraction SHARED ${DipoleSubtraction_esources}
                             ${DipoleSubtraction_SOURCES} 
)
target_include_directories(DipoleSubtraction PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(DipoleSubtraction PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(DipoleSubtraction PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS DipoleSubtraction DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(DipoleSubtraction        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME DipoleSubtraction SOVERSION 0.0.0)
endif()
