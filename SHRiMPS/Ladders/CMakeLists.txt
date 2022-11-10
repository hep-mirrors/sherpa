
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Ladders/Makefile.am
 
 #pkglib_LTLIBRARIES = libShrimpsLadders.la)
 set(GITTAG  SHRiMPS/Ladders)
 #include ../../git.make)
set(Ladders_SOURCES  Ladder_Particle.C T_Prop.C Ladder.C Ladder_Generator_Base.C Ladder_Generator_Eik.C Ladder_Generator_QT.C Ladder_Generator_KT.C Ladder_Generator_Seeded.C MEs.C Primary_Ladders.C)
 #localincdir = $(pkg#includedir)/SHRiMPS/Ladders)
set(Ladders_HEADERS Ladder_Particle.H T_Prop.H Ladder.H Ladder_Generator_Base.H Ladder_Generator_Eik.H Ladder_Generator_QT.H Ladder_Generator_KT.H Ladder_Generator_Seeded.H MEs.H Primary_Ladders.H)
creategitinfo(Ladders ../Ladders/)
list(TRANSFORM Ladders_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Ladders/")
list(APPEND Ladders_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ShrimpsLadders SHARED 
                             ${Ladders_SOURCES} 
)
target_include_directories(ShrimpsLadders PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ShrimpsLadders)
install(TARGETS ShrimpsLadders DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ShrimpsLadders        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ShrimpsLadders SOVERSION 0.0.0)
set_target_properties(ShrimpsLadders        PROPERTIES DEFINE_SYMBOL "")
