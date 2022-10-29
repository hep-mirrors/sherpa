
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Amplitude/Makefile.am
 
 #SUBDIRS = Zfunctions)
 #pkglib_LTLIBRARIES = libAmplitude.la)
 set(GITTAG  AMEGIC++/Amplitude)
 #include ../../git.make)
set(Amplitude_SOURCES  Single_Vertex.C Vertex.C CFColor.C Amplitude_Base.C Single_Amplitude.C Amplitude_Manipulator.C Amplitude_Generator.C Color_Generator.C Prop_Generator.C Zfunc_Generator.C Amplitude_Group.C
 Amplitude_Handler.C Single_Amplitude_Base.C Super_Amplitude.C Zfunc.C FullAmplitude_MHV_Base.C FullAmplitude_External.C Amplitude_Output.C Lorentz_Function.C Lorentz_Functions.C)
set(Amplitude_HEADERS Single_Vertex.H Vertex.H CFColor.H Single_Amplitude.H Amplitude_Generator.H Amplitude_Manipulator.H Color_Generator.H Prop_Generator.H Zfunc_Generator.H Pfunc.H Zfunc.H Amplitude_Handler.H
 Amplitude_Base.H Amplitude_Group.H Color_Group.H Super_Amplitude.H Single_Amplitude_Base.H FullAmplitude_MHV_Base.H FullAmplitude_External.H Amplitude_Output.H Lorentz_Function.H)
creategitinfo(Amplitude ../Amplitude/)
  list(TRANSFORM Amplitude_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Amplitude/")
 
set(Amplitude_esources )
if (SHERPA_BUILD_STATIC)
add_library(Amplitude_static STATIC ${Amplitude_esources}
                             ${Amplitude_SOURCES} 
)
target_include_directories(Amplitude_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amplitude_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Amplitude_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS Amplitude_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amplitude_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME Amplitude)
set_target_properties(Amplitude_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(Amplitude SHARED ${Amplitude_esources}
                             ${Amplitude_SOURCES} 
)
target_include_directories(Amplitude PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Amplitude PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Amplitude PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS Amplitude DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Amplitude        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Amplitude SOVERSION 0.0.0)
set_target_properties(Amplitude        PROPERTIES DEFINE_SYMBOL "")
endif()
