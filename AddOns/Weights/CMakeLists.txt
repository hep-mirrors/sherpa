
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Weights/Makefile.am
 
 #pkglib_LTLIBRARIES = libSherpaWeightsOutput.la)
 set(GITTAG  AddOns/Weights)
 #include ../../git.make)
set(Weights_SOURCES  Output_Weights.C)
creategitinfo(Weights ../Weights/)
list(TRANSFORM Weights_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Weights/")
list(APPEND Weights_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Weights_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaWeights SHARED ${SherpaWeights_esources}
                             ${Weights_SOURCES} 
)
target_include_directories(SherpaWeights PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaWeights PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaWeights DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaWeights        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaWeights SOVERSION 0.0.0)
set_target_properties(SherpaWeights        PROPERTIES DEFINE_SYMBOL "")
endif()
