##### Automatically or semiautomaticaly generated, do not edit.#########
#EWSud/Makefile.am
 
 if (SHERPA_ENABLE_EWSUD)
 #SUBDIRS = Tests)
 endif()
 set(GITTAG  AddOns/EWSud)
set(EWSud_SOURCES Variation_Generator.C Coefficient_Checker.C KFactor_Checker.C Comix_Interface.C HE_Comix_Interface.C EWGroupConstants.C EWSud.C Amplitudes.C  KFactor.C Calculator.C Clustered_Calculator.C)
set(EWSUD_ADDS   Variation_Generator.H Coefficient_Checker.H KFactor_Checker.H Comix_Interface.H HE_Comix_Interface.H EWGroupConstants.H EWSud.H Amplitudes.H KFactor.H Calculator.H Clustered_Calculator.H)
 if (SHERPA_ENABLE_EWSUD)
 #pkglib_LTLIBRARIES = libSherpaEWSud.la)
set(EWSud_SOURCES ${EWSud_SOURCES})
set(EWSud_HEADERS ${EWSUD_ADDS})
 #EWSUD_#EXTRA_DIST = ${EWSUD_ADDS})
 endif()
 #EXTRA_DIST = $(#EWSUD_#EXTRA_DIST)
sherpa_create_git_info(EWSud ../EWSud/)
list(TRANSFORM EWSud_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../EWSud/")
list(APPEND EWSud_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(SherpaEWSud SHARED ${EWSud_SOURCES} )
sherpa_handle_shared_library(SherpaEWSud SherpaEWSud)
