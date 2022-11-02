
########################################################################
#  Automatically or semiautomaticaly generated, do not edit.
########################################################################
# The following input was used
# Makefile.am

########################################################################
#MEs/Makefile.am
 
 #pkglib_LTLIBRARIES = libPhotonsMEs.la)
 set(GITTAG  PHOTONS++/MEs)
 #include ../../git.make)
set(MEs_SOURCES  PHOTONS_ME_Base.C Scalar_To_Scalar_Scalar.C Scalar_To_Fermion_Fermion.C Scalar_To_Scalar_Lepton_Neutrino.C Scalar_To_Vector_Lepton_Neutrino.C Tau_To_Lepton_Neutrinos.C W_To_Lepton_Neutrino.C
 Z_To_Fermion_Fermion.C Vector_To_Fermion_Fermion.C Vector_To_Scalar_Scalar.C)
set(MEs_HEADERS PHOTONS_ME_Base.H Scalar_To_Scalar_Scalar.H Scalar_To_Fermion_Fermion.H Scalar_To_Scalar_Lepton_Neutrino.H Scalar_To_Vector_Lepton_Neutrino.H Tau_To_Lepton_Neutrinos.H W_To_Lepton_Neutrino.H
 Z_To_Fermion_Fermion.H Vector_To_Fermion_Fermion.H Vector_To_Scalar_Scalar.H)
creategitinfo(MEs ../MEs/)
list(TRANSFORM MEs_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/../MEs/")
list(APPEND MEs_SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C)
 
add_library(PhotonsMEs SHARED 
                             ${MEs_SOURCES} 
)
target_include_directories(PhotonsMEs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS PhotonsMEs DESTINATION ${CMAKE_INSTALL_LIBDIR}/SHERPA-MC COMPONENT libs)
set_target_properties(PhotonsMEs        PROPERTIES POSITION_INDEPENDENT_CODE ON OUTPUT_NAME PhotonsMEs SOVERSION 0.0.0)
set_target_properties(PhotonsMEs        PROPERTIES DEFINE_SYMBOL "")
