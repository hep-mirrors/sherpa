### INSTRUCTIONS FOR ADDING AN EXAMPLE #########################################
#
# First consider if the example should be installed. If not, ignore these
# instructions. If yes, add the to-be-installed files in the lists below:
#
# - INSTALL_PATHS: If it should simply be installed as it is.
# - CONFIGURE_AND_INSTALL_PATHS: If it should be configured before installing.
#
# The latter list can therefore be used for *.in files. If the example relies
# on the python interface, look for the `if (SHERPA_ENABLE_PYTHON)` clause,
# where the lists above are extended conditionally.


### HELPER FUNCTIONS FOR INSTALLING EXAMPLE FILES ##############################

# Helper function to install an example file, preserving the subdirectories
function (INSTALL_WITHIN_TREE PATH)
  get_filename_component(target_dir "${PATH}" DIRECTORY)
  install(FILES "${PATH}"
          DESTINATION "${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING}/Examples/${target_dir}")
endfunction()

# Helper function to configure and install an example file, preserving the
# subdirectories
function (CONFIGURE_AND_INSTALL_WITHIN_TREE PATH)
  configure_file("${PATH}.in" "${CMAKE_CURRENT_BINARY_DIR}/${PATH}" @ONLY)
  get_filename_component(target_dir "${PATH}" DIRECTORY)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PATH}"
          DESTINATION "${CMAKE_INSTALL_DATADIR}/SHERPA-MC${VERSIONING}/Examples/${target_dir}")
endfunction()


### CONSTRUCT LISTS OF EXAMPLE FILES TO BE INSTALLED ###########################

# Example files that should be installed as they are
set (INSTALL_PATHS
     API/ME2-CPP/CMakeLists.txt
     API/ME2-CPP/Sherpa.yaml
     API/ME2-CPP/test.cpp
     Approximate_EW_Corrections/EWSudakov/Sherpa.DY.LOPS.yaml
     Approximate_EW_Corrections/EWSudakov/Sherpa.DY.MENLOPS.yaml
     BFactory/BDecay/Sherpa.yaml
     BFactory/Belle_QCDBackground/Sherpa.YFS.yaml
     BFactory/Belle_QCDBackground/Sherpa.yaml
     BFactory/Belle_Signal/Sherpa.YFS.yaml
     BFactory/Belle_Signal/Sherpa.yaml
     BSM/SM_ZPrime/Makefile
     BSM/SM_ZPrime/Model.C
     BSM/SM_ZPrime/Sherpa.yaml
     BSM/SM_ggV/Makefile
     BSM/SM_ggV/Model.C
     BSM/SM_ggV/Sherpa.WW.yaml
     BSM/SM_ggV/Sherpa.ZZ.yaml
     BSM/SM_ggV/Sherpa.h.yaml
     BSM/UFO_AQGC/AQGC/__init__.py
     BSM/UFO_AQGC/AQGC/coupling_orders.py
     BSM/UFO_AQGC/AQGC/couplings.py
     BSM/UFO_AQGC/AQGC/function_library.py
     BSM/UFO_AQGC/AQGC/lorentz.py
     BSM/UFO_AQGC/AQGC/object_library.py
     BSM/UFO_AQGC/AQGC/parameters.py
     BSM/UFO_AQGC/AQGC/particles.py
     BSM/UFO_AQGC/AQGC/propagators.py
     BSM/UFO_AQGC/AQGC/vertices.py
     BSM/UFO_AQGC/AQGC/write_param_card.py
     BSM/UFO_MSSM/MSSM/__init__.py
     BSM/UFO_MSSM/MSSM/build_restrict.py
     BSM/UFO_MSSM/MSSM/coupling_orders.py
     BSM/UFO_MSSM/MSSM/couplings.py
     BSM/UFO_MSSM/MSSM/function_library.py
     BSM/UFO_MSSM/MSSM/lorentz.py
     BSM/UFO_MSSM/MSSM/object_library.py
     BSM/UFO_MSSM/MSSM/parameters.py
     BSM/UFO_MSSM/MSSM/particles.py
     BSM/UFO_MSSM/MSSM/vertices.py
     BSM/UFO_MSSM/MSSM/write_param_card.py
     BSM/UFO_MUED/MUED/__init__.py
     BSM/UFO_MUED/MUED/coupling_orders.py
     BSM/UFO_MUED/MUED/couplings.py
     BSM/UFO_MUED/MUED/function_library.py
     BSM/UFO_MUED/MUED/lorentz.py
     BSM/UFO_MUED/MUED/model.pkl
     BSM/UFO_MUED/MUED/object_library.py
     BSM/UFO_MUED/MUED/parameters.py
     BSM/UFO_MUED/MUED/particles.py
     BSM/UFO_MUED/MUED/vertices.py
     BSM/UFO_MUED/MUED/write_param_card.py
     BSM/UFO_SMEFT/Sherpa.yaml
     EPA/LEP_EPA_BosonPairs.yaml
     EPA/LEP_EPA_JetProduction.yaml
     EPA/LHC_EPA_JetProduction.yaml
     EPA/MinBias.yaml
     FixedOrder/MINLO_ppenuj/Sherpa.yaml
     FixedOrder/NLO_EW_ppWW/Sherpa.LO.yaml
     FixedOrder/NLO_EW_ppWW/Sherpa.RS.EW.yaml
     FixedOrder/NLO_QCD_pptt/Sherpa.BI.QCD.yaml
     FixedOrder/NNLO_QCD_pph/Sherpa.yaml
     FixedOrder/NNLO_QCD_ppll/Sherpa.yaml
     FixedOrder/NNLO_QCD_pplv/Sherpa.yaml
     HPC/WJets/Sherpa.yaml
     HPC/WJets_Async/Sherpa.yaml
     H_in_AssociatedProduction/LHC_Hll_MCatNLO/Sherpa.yaml
     H_in_AssociatedProduction/LHC_WHJets/Sherpa.yaml
     H_in_AssociatedProduction/LHC_ZHJets/Sherpa.yaml
     H_in_GluonFusion/LHC_HInt_MCatNLO/Sherpa.yaml
     H_in_GluonFusion/LHC_HInterference/Sherpa.yaml
     H_in_GluonFusion/LHC_HJets/Sherpa.yaml
     H_in_GluonFusion/LHC_HJets_Finite_MTop/Sherpa.yaml
     H_in_GluonFusion/LHC_HJets_UN2LOPS/Sherpa.yaml
     H_in_TTBar/LHC_TTH_MCatNLO/Sherpa.yaml
     Instanton/LHC_13/Instanton.yaml
     Instanton/LHC_13/Instanton_InvRho_noMPI.yaml
     Instanton/LHC_13/Instanton_InvRho_wMPI.yaml
     Instanton/LHC_13/Instanton_Original.yaml
     Instanton/LHC_13/Instanton_shatoverN_wMPI.yaml
     Instanton/LHC_13/mc.Sh_Inst_m300_nf5.yaml
     Jets_at_HadronColliders/LHC_Jets_MCatNLO/Sherpa.yaml
     Jets_at_HadronColliders/LHC_Jets_MEPS/Sherpa.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.197.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.206.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.EPA.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.EPA_PDF.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.GammaPDF.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.NLO.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.PDF.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.Reconnect.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.YFS.yaml
     Jets_at_LeptonColliders/LEP_Jets/Sherpa.yaml
     Jets_at_LeptonColliders/PETRA_Jets/Sherpa.14.yaml
     Jets_at_LeptonColliders/PETRA_Jets/Sherpa.35.yaml
     Jets_at_LeptonColliders/PETRA_Jets/Sherpa.44.yaml
     Jets_at_LeptonColliders/PETRA_Jets/Sherpa.yaml
     Jets_in_DIS/HERA/Sherpa.yaml
     Photoproduction/HERA/Sherpa.Direct.yaml
     Photoproduction/HERA/Sherpa.Resolved.yaml
     Photoproduction/LEP/Sherpa.Direct.yaml
     Photoproduction/LEP/Sherpa.Double-Resolved.yaml
     Photoproduction/LEP/Sherpa.Single-Resolved.yaml
     RelicDensity/Sherpa.yaml
     RemnantTests/SimpleHera/Sherpa.yaml
     RemnantTests/SimpleLEP/Sherpa.EPA.yaml
     RemnantTests/SimpleLEP/Sherpa.ISR.yaml
     RemnantTests/SimpleLEP/Sherpa.yaml
     RemnantTests/SimpleLHC/DY.yaml
     RemnantTests/SimpleLHC/Sherpa.MPI.yaml
     RemnantTests/SimpleLHC/Sherpa.yaml
     SingleTop_Channels/Sherpa.tW.yaml
     SingleTop_Channels/Sherpa.tj-s_channel-nf4.yaml
     SingleTop_Channels/Sherpa.tj-s_channel.yaml
     SingleTop_Channels/Sherpa.tj-t_channel-nf4.yaml
     SingleTop_Channels/Sherpa.tj-t_channel.yaml
     Soft_QCD/LHC_7TeV_MinBias/Amisic.yaml
     Soft_QCD/LHC_7TeV_MinBias/Elastic.yaml
     Soft_QCD/LHC_7TeV_MinBias/Shrimps.yaml
     Tops_plus_Jets/LHC_Tops/Sherpa.EWapprox.yaml
     Tops_plus_Jets/LHC_Tops/Sherpa.Reconnect.yaml
     Tops_plus_Jets/LHC_Tops/Sherpa.yaml
     Tops_plus_V/LHC_TTW/Sherpa.yaml
     Userhook/Sherpa.yaml
     Userhook/Userhook_Example.C
     VBF_Interference/Loop2/Sherpa.yaml
     VV_plus_Jets/LHC_2l2nu2jJets_SameSign/Sherpa.yaml
     VV_plus_Jets/LHC_2l2nuJets/Sherpa.yaml
     VV_plus_Jets/LHC_2l2nuJets_GluonInitiated/Sherpa.tree.yaml
     VV_plus_Jets/LHC_2l2nuJets_GluonInitiated/Sherpa.yaml
     VV_plus_Jets/LHC_4lJets/Sherpa.yaml
     VV_plus_Jets/LHC_4lJets_GluonInitiated/Sherpa.tree.yaml
     VV_plus_Jets/LHC_4lJets_GluonInitiated/Sherpa.yaml
     VV_plus_Jets/LHC_WZ_polarization_at_nLO_PS/Sherpa.yaml
     VV_plus_Jets/LHC_WqqZnunuJets/Sherpa.yaml
     VV_plus_Jets/LHC_ssWW_polarization/Sherpa.yaml
     V_plus_Bs/LHC_Wbb/Sherpa.yaml
     V_plus_Bs/LHC_Wbb/Wbb_Interface.C
     V_plus_Bs/LHC_Wbb/Wbb_Virtual.f
     V_plus_Bs/LHC_Zbb/Sherpa.yaml
     V_plus_Jets/LHC_WJets/My_JetCriterion.C
     V_plus_Jets/LHC_WJets/Sherpa.yaml
     V_plus_Jets/LHC_WJets/cmake/Modules/FindFastJet.cmake
     V_plus_Jets/LHC_WJets_UN2LOPS/Sherpa.yaml
     V_plus_Jets/LHC_ZJets/Sherpa.LO.yaml
     V_plus_Jets/LHC_ZJets/Sherpa.yaml
     V_plus_Jets/LHC_ZJets_UN2LOPS/Sherpa.yaml
)

# Example files that should be first configured and then installed
set (CONFIGURE_AND_INSTALL_PATHS
     Userhook/CMakeLists.txt
     V_plus_Bs/LHC_Wbb/CMakeLists.txt
     V_plus_Jets/LHC_WJets/CMakeLists.txt
)

# Conditionally extend the above lists with files for examples that rely on the
# Python interface.
if (SHERPA_ENABLE_PYTHON)
  set (CONFIGURE_AND_INSTALL_PATHS ${CONFIGURE_AND_INSTALL_PATHS}
       API/Events/test.py
       API/ME2-Python/test.py
       API/MPIEvents/test.py
       VBF_Interference/BH_Interface/BlackHat.py
       VBF_Interference/BH_Interface/Comix.py
       VBF_Interference/Tests/Comix.py
  )
  set (INSTALL_PATHS ${INSTALL_PATHS}
       API/Events/Sherpa.yaml
       API/ME2-Python/Sherpa.yaml
       API/MPIEvents/Sherpa.yaml
  )
endif()


### INSTALL EXAMPLE FILES ######################################################

foreach (PATH ${INSTALL_PATHS})
  INSTALL_WITHIN_TREE(${PATH})
endforeach()

foreach (PATH ${CONFIGURE_AND_INSTALL_PATHS})
  CONFIGURE_AND_INSTALL_WITHIN_TREE(${PATH})
endforeach()
