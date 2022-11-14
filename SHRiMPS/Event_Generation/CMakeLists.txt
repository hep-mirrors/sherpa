##### Automatically or semiautomaticaly generated, do not edit.#########
#Event_Generation/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsEvents.la)
 set(GITTAG  SHRiMPS/Events)
set(Event_Generation_SOURCES  Event_Generator.C Inelastic_Event_Generator.C Elastic_Event_Generator.C Soft_Diffractive_Event_Generator.C Collinear_Emission_Generator.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Event_Generation)
set(Event_Generation_HEADERS Event_Generator_Base.H Event_Generator.H Inelastic_Event_Generator.H Elastic_Event_Generator.H Soft_Diffractive_Event_Generator.H Collinear_Emission_Generator.H)
sherpa_create_git_info(Event_Generation ../Event_Generation/)
list(TRANSFORM Event_Generation_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Event_Generation/")
list(APPEND Event_Generation_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsEvents SHARED ${Event_Generation_SOURCES} )
sherpa_handle_shared_library(ShrimpsEvents ShrimpsEvents)
