(run){
  # general setting
  EVENTS 1M; INTEGRATION_ERROR 0.05;

  # tags for process setup
  NJET:=4; LJET:=2,3,4; QCUT:=20.;

  # me generator settings
  ME_SIGNAL_GENERATOR Comix Amegic LOOPGEN;
  LOOPGEN:=OpenLoops;

  # collider setup
  BEAM_1 2212; BEAM_ENERGY_1 = 6500.;
  BEAM_2 2212; BEAM_ENERGY_2 = 6500.;

  ### Fusing settings fragmentation component
  SHERPA_LDADD=SherpaFusing
  USERHOOK = Fusing_Fragmentation
  ## to set event weights instead of rejecting events use:
  # FUSING_FRAGMENTATION_STORE_AS_WEIGHT=1

  # for a consistent fusing, the parameters have to be chosen identical between
  # direct and fragmentation component
  MASS[5]=4.75 # consistent with PDF set
  COMIX_CLUSTER_ORDERED 1;
  COMIX_CLUSTER_RS_ORDERED 1;
  CSS_SCALE_SCHEME 2;
  CSS_EVOLUTION_SCHEME 3;
  PP_RS_SCALE VAR{H_TM2};
  PP_HPSMODE=0;
}(run)

(processes){
  Process 93 93 -> 11 -11 93{NJET};
  Order (*,2); CKKW sqr(QCUT/E_CMS);
  NLO_QCD_Mode MC@NLO {LJET};
  ME_Generator Amegic {LJET};
  RS_ME_Generator Comix {LJET};
  Loop_Generator LOOPGEN {LJET};
  End process;
}(processes)

(selector){
  Mass 11 -11 66 E_CMS
}(selector)
