
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Phys/Makefile.am
 
 set(GITTAG  ATOOLS/Phys)
 #include ../../git.make)
 #pkglib_LTLIBRARIES = libToolsPhys.la)
set(Phys_SOURCES  Fragmentation_Base.C fjcore.cc Flavour.C Ordering.C Particle.C Simple_Polarisation_Info.C Particle_List.C Particle_Dresser.C KF_Table.C NLO_Subevt.C NLO_Types.C Decay_Info.C Cluster_Leg.C
 Cluster_Amplitude.C Blob.C Blob_List.C Parton_Finder.C Color_Tester.C Momentum_Shifter.C Flow.C Color.C Spinor.C Selector_List.C Fastjet_Helpers.C Momenta_Stretcher.C Weight_Info.C Variations.C Weights.C
 Hard_Process_Variation_Generator.C)
set(Phys_HEADERS fjcore.hh Fragmentation_Base.H Ordering.H Flavour_Tags.H Flavour.H Particle.H Simple_Polarisation_Info.H Particle_List.H Particle_Dresser.H KF_Table.H NLO_Subevt.H NLO_Types.H Decay_Info.H
 Cluster_Leg.H Cluster_Amplitude.H Flow.H Blob.H Blob_List.H Parton_Finder.H Color_Tester.H Momentum_Shifter.H Color.H Spinor.H Selector_List.H Fastjet_Helpers.H Momenta_Stretcher.H Weight_Info.H Variations.H
 Weights.H Hard_Process_Variation_Generator.H)
creategitinfo(Phys ../Phys/)
list(TRANSFORM Phys_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Phys/")
list(APPEND Phys_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(ToolsPhys SHARED 
                             ${Phys_SOURCES} 
)
target_include_directories(ToolsPhys PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
sherpa_mpi_link_libraries(ToolsPhys)
install(TARGETS ToolsPhys DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(ToolsPhys        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME ToolsPhys SOVERSION 0.0.0)
set_target_properties(ToolsPhys        PROPERTIES DEFINE_SYMBOL "")

if (SHERPA_ENABLE_LHAPDF)
target_link_libraries(ToolsPhys PRIVATE ${LHAPDF_LIBRARIES})
target_include_directories(ToolsPhys PRIVATE ${LHAPDF_INCLUDE_DIRS})
endif()
