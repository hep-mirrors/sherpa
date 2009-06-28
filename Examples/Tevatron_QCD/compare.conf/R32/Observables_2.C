(observables2){
  PIECE_SETUP rc.C (rc0){ }(rc0);
  PIECE_SETUP Observables_2.C (dr32){ }(dr32);
  PIECE_SETUP dc.C (dc0){ }(dc0);
  LEFT_MARGIN 0.15;
  Y_TITLE_OFFSET 1.4;
  X_TITLE_OFFSET 1.1;
  @@ PATHPIECE NJet_Cone/output_C07/;
  LEG_DELTAY 0.04;
  LEG_TEXT_SIZE 0.03;
}(observables2);

(dr32){
  PATH_PIECE DATAPATH/
  FILE_PIECE FinalState_R32_Data.dat;
  ALIAS_NAME FinalStateD32;
  LEGEND_TITLE D0 Data;
  DATA_TYPE ASCII;
  DRAW_OPTION P;
  X_VALUE 1;
  Y_VALUE 3;
  X_ERROR_PLUS 2;
  X_ERROR_MINUS 2;
  Y_ERROR_PLUS 5;
  Y_ERROR_MINUS 4;
}(dr32);

