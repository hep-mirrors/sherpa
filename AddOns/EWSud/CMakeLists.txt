
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#EWSud/Makefile.am
 
 if (SHERPA_ENABLE_EWSUD)
 #SUBDIRS = Tests)
 endif()
 set(GITTAG  AddOns/EWSud)
 #include ../../git.make)
set(EWSud_SOURCES Variation_Generator.C Coefficient_Checker.C KFactor_Checker.C Comix_Interface.C HE_Comix_Interface.C EWGroupConstants.C EWSud.C Amplitudes.C  KFactor.C Calculator.C Clustered_Calculator.C)
set(EWSUD_ADDS   Variation_Generator.H Coefficient_Checker.H KFactor_Checker.H Comix_Interface.H HE_Comix_Interface.H EWGroupConstants.H EWSud.H Amplitudes.H KFactor.H Calculator.H Clustered_Calculator.H)
 if (SHERPA_ENABLE_EWSUD)
 #pkglib_LTLIBRARIES = libSherpaEWSud.la)
set(EWSud_SOURCES ${EWSud_SOURCES})
set(EWSud_HEADERS ${EWSUD_ADDS})
 #EWSUD_#EXTRA_DIST = ${EWSUD_ADDS})
 endif()
 #EXTRA_DIST = $(#EWSUD_#EXTRA_DIST)
creategitinfo(EWSud ../EWSud/)
list(TRANSFORM EWSud_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../EWSud/")
list(APPEND EWSud_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(EWSud_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaEWSud SHARED ${SherpaEWSud_esources}
                             ${EWSud_SOURCES} 
)
target_include_directories(SherpaEWSud PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaEWSud PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaEWSud DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaEWSud        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaEWSud SOVERSION 0.0.0)
set_target_properties(SherpaEWSud        PROPERTIES DEFINE_SYMBOL "")
endif()
