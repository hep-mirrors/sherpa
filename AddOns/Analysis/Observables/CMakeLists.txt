
########################################################################
#
#  Automatically or semiautomaticaly generated, do not edit.
#
########################################################################
# The following input was used
# Makefile.am

########################################################################
#Observables/Makefile.am
 
 set(GITTAG  AddOns/Analysis/Observables)
 #include ../../../git.make)
set(Observables_SOURCES  Event_Output.C Blobdata.C Global_Calculator.C Primitive_Observable_Base.C Normalized_Observable.C One_Variable_Observable.C HT.C XS.C TotalMass.C LMass.C EV_C_Observables.C EV_C2_Observables.C
 EV_C4_Observables.C EV_SC_Observables.C PSM_Observables.C Multiplicity.C One_Particle_Observables.C One_Particle_Observables_2.C Conditional_One_Particle_Observables.C Two_Particle_Observables.C Two_Particle_Observables_2.C
 Three_Particle_Observables.C Four_Particle_Observables.C Four_Particle_Observables_2.C Six_Particle_Observables.C Soft_Photon_Observables.C Jet_Cone_Distribution.C Jet_Observables.C Jet_NJet_Observables.C
 Jet_X_Observables.C Momentum_Tensor.C CD_Parameters.C Sphericity_Observables.C Event_Shapes_EE.C Shape_Observables_EE.C Jet_Mass_and_Broadening.C Scaled_Observables.C Four_Jet_Angles.C ChargedDEta.C)
set(ANAOBS_ADDS  Primitive_Observable_Base.H Normalized_Observable.H HT.H PSM_Observables.H Multiplicity.H MI_Observables.H One_Particle_Observables.H Two_Particle_Observables.H Three_Particle_Observables.H
 Four_Particle_Observables.H Six_Particle_Observables.H Soft_Photon_Observables.H Jet_Cone_Distribution.H Jet_Observables.H Momentum_Tensor.H CD_Parameters.H Sphericity_Observables.H Event_Shapes_EE.H
 Shape_Observables_EE.H Jet_Mass_and_Broadening.H Scaled_Observables.H Four_Jet_Angles.H ChargedDEta.H)
 if (SHERPA_ENABLE_ANALYSIS)
 #pkglib_LTLIBRARIES = libSherpaObservables.la)
set(Observables_SOURCES ${ANAOBS_SOURCES})
set(Observables_HEADERS ${ANAOBS_ADDS})
 #ANAOBS_#EXTRA_DIST = ${ANAOBS_ADDS})
 else()
 #ANAOBS_#EXTRA_DIST = ${ANAOBS_SOURCES} ${ANAOBS_ADDS})
 endif()
 #EXTRA_DIST = $(#ANAOBS_#EXTRA_DIST)
creategitinfo(Observables ../Observables/)
list(TRANSFORM Observables_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../Observables/")
list(APPEND Observables_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
set(Observables_esources )
if (SHERPA_BUILD_SHARED)
add_library(SherpaObservables SHARED ${SherpaObservables_esources}
                             ${Observables_SOURCES} 
)
target_include_directories(SherpaObservables PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(SherpaObservables PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS SherpaObservables DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(SherpaObservables        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME SherpaObservables SOVERSION 0.0.0)
set_target_properties(SherpaObservables        PROPERTIES DEFINE_SYMBOL "")
endif()