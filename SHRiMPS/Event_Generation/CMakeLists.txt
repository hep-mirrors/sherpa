
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Event_Generation/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsEvents.la)
 set(GITTAG  SHRiMPS/Events)
 #include ../../git.make)
set(Event_Generation_SOURCES  Event_Generator.C Inelastic_Event_Generator.C Elastic_Event_Generator.C Soft_Diffractive_Event_Generator.C Collinear_Emission_Generator.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Event_Generation)
set(Event_Generation_HEADERS Event_Generator_Base.H Event_Generator.H Inelastic_Event_Generator.H Elastic_Event_Generator.H Soft_Diffractive_Event_Generator.H Collinear_Emission_Generator.H)
creategitinfo(Event_Generation ../Event_Generation/)
list(TRANSFORM Event_Generation_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Event_Generation/")
list(APPEND Event_Generation_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsEvents SHARED 
                             ${Event_Generation_SOURCES} 
)
target_include_directories(ShrimpsEvents PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsEvents DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsEvents        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsEvents SOVERSION 0.0.0)
set_target_properties(ShrimpsEvents        PROPERTIES DEFINE_SYMBOL "")
