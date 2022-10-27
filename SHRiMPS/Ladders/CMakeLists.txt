
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Ladders/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsLadders.la)
set(Ladders_SOURCES Ladder_Particle.C T_Prop.C Ladder.C Ladder_Generator_Base.C Ladder_Generator_Eik.C Ladder_Generator_QT.C Ladder_Generator_KT.C Ladder_Generator_Seeded.C MEs.C Primary_Ladders.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Ladders)
set(Ladders_HEADERS Ladder_Particle.H T_Prop.H Ladder.H Ladder_Generator_Base.H Ladder_Generator_Eik.H Ladder_Generator_QT.H Ladder_Generator_KT.H Ladder_Generator_Seeded.H MEs.H Primary_Ladders.H)
creategitinfo(Ladders ../Ladders/)
  list(TRANSFORM Ladders_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Ladders/")
 
set(Ladders_esources )
if (SHERPA_BUILD_STATIC)
add_library(ShrimpsLadders_static STATIC ${ShrimpsLadders_esources}
                             ${Ladders_SOURCES} 
)
target_include_directories(ShrimpsLadders_static PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsLadders_static PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsLadders_static DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsLadders_static PROPERTIES POSITION_INDEPENDENT_CODE ${SHERPA_POSITION_INDEPENDENT_CODE} OUTPUT_NAME ShrimpsLadders)
endif()
if (SHERPA_BUILD_SHARED)
add_library(ShrimpsLadders SHARED ${ShrimpsLadders_esources}
                             ${Ladders_SOURCES} 
)
target_include_directories(ShrimpsLadders PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(ShrimpsLadders PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS ShrimpsLadders DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsLadders        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsLadders SOVERSION 0.0.0)
endif()
