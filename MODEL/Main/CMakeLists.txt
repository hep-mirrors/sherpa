
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Main/Makefile.am
 
 #pkglib_LTLIBRARIES = libModelMain.la)
 set(GITTAG  MODEL/Main)
 #include ../../git.make)
set(Main_SOURCES  Color_Function.C Single_Vertex.C Coupling_Data.C Model_Base.C Running_AlphaS.C Strong_Coupling.C Running_AlphaQED.C Running_Fermion_Mass.C)
set(Main_HEADERS Color_Function.H Single_Vertex.H Coupling_Data.H Model_Base.H Running_AlphaS.H Strong_Coupling.H Running_AlphaQED.H Running_Fermion_Mass.H)
creategitinfo(Main ../Main/)
  list(TRANSFORM Main_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Main/")
 
set(Main_esources )
if (SHERPA_BUILD_STATIC)
add_library(ModelMain_static STATIC ${ModelMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ModelMain_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ModelMain_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ModelMain_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ModelMain_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ModelMain)
set_target_properties(ModelMain_static PROPERTIES DEFINE_SYMBOL "")
endif()
if (SHERPA_BUILD_SHARED)
add_library(ModelMain SHARED ${ModelMain_esources}
                             ${Main_SOURCES} 
)
target_include_directories(ModelMain PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ModelMain PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ModelMain DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ModelMain        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ModelMain SOVERSION 0.0.0)
set_target_properties(ModelMain        PROPERTIES DEFINE_SYMBOL "")
endif()
