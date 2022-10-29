
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Zfunctions/Makefile.am
 
 #pkglib_LTLIBRARIES = libZfunctions.la)
 set(GITTAG  AMEGIC++/Amplitude/Zfunctions)
 #include ../../../git.make)
set(Zfunctions_SOURCES  Basic_Func.C Basic_Sfuncs.C Basic_Pfunc.C Basic_Mfunc.C Basic_MassTermfunc.C Basic_Vfunc.C Basic_Xfunc.C Basic_Yfunc.C Basic_Zfunc.C Basic_Epsilonfunc.C Unitarityfunc.C G4_Calc.C
 V4_Calc.C VVV_Calc.C Loops_Calc.C Zfunc_Calc.C Calculator.C Mom.C MHVCalculator.C)
set(Zfunctions_HEADERS Basic_Func.H Basic_Sfuncs.H Basic_Func.icc Zfunc_Calc.H Calculator.H Mom.H MHVCalculator.H)
creategitinfo(Zfunctions ../Zfunctions/)
  list(TRANSFORM Zfunctions_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Zfunctions/")
 
set(Zfunctions_esources )
if (SHERPA_BUILD_STATIC)
add_library(Zfunctions_static STATIC ${Zfunctions_esources}
                             ${Zfunctions_SOURCES} 
)
target_include_directories(Zfunctions_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Zfunctions_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Zfunctions_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(Zfunctions_static PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS Zfunctions_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Zfunctions_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME Zfunctions)
set_target_properties(Zfunctions_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(Zfunctions SHARED ${Zfunctions_esources}
                             ${Zfunctions_SOURCES} 
)
target_include_directories(Zfunctions PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(Zfunctions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Zfunctions PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
target_include_directories(Zfunctions PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/")
install(TARGETS Zfunctions DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(Zfunctions        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME Zfunctions SOVERSION 0.0.0)
set_target_properties(Zfunctions        PROPERTIES DEFINE_SYMBOL "")
endif()
