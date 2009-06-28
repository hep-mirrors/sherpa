(rc0){
  PIECE_SETUP rc.C (djet){ }(djet);
  PIECE_SETUP rc.C (sjet){ }(sjet);
  DRAW YES;
  X_MIN 100;
  X_MAX 600;
  Y_MIN 0.001;
  Y_MAX 0.7;
  Y_AXIS_TITLE R_{32}(H_{T});
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  BOTTOM_MARGIN 0.3;
  X_AXIS_LABEL_SIZE 0;
  X_TITLE_SIZE 0;
  HISTOGRAM_NAME DR32_DHT;
  FIGURE_CAPTION $R_{32}(H_{T})$;
  WEBPAGE_CAPTION R<sub>32</sub>(H<sub>T</sub>);
  DRAW_LINE V 170 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 240 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 330 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 430 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LATEX SHERPA | LEFT 0.05 TOP 0.075 ALIGN 12 PRIORITY -20\;;
}(rc0);

(sjet){
  PIECE_SETUP rc.C (sjet20){ }(sjet20);
  PIECE_SETUP rc.C (sjet30){ }(sjet30);
  PIECE_SETUP rc.C (sjet50){ }(sjet50);
  PIECE_SETUP rc.C (sjet85){ }(sjet85);
  PIECE_SETUP rc.C (sjet115){ }(sjet115);
  DATA_TYPE ALGEBRA(y[0]/y[1]);
  LEG_LEFT 0.075;
  LEG_RIGHT 0.325;
  LEG_TOP 1.15;
}(sjet);
(sjet20){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_jet20_HT.dat c0_two_jet_inc_jet20_HT.dat;
  ALIAS_NAME jet20R32;
  ## ADOPT_BINS jet20D32[0];
}(sjet20);
(sjet30){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_jet30_HT.dat c0_two_jet_inc_jet30_HT.dat;
  ALIAS_NAME jet30R32;
  ## ADOPT_BINS jet30D32[0];
  DRAW_LEGEND NO;
}(sjet30);
(sjet50){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_jet50_HT.dat c0_two_jet_inc_jet50_HT.dat;
  ALIAS_NAME jet50R32;
  ## ADOPT_BINS jet50D32[0];
  DRAW_LEGEND NO;
}(sjet50);
(sjet85){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_jet85_HT.dat c0_two_jet_inc_jet85_HT.dat;
  ALIAS_NAME jet85R32;
  ## ADOPT_BINS jet85D32[0];
  DRAW_LEGEND NO;
}(sjet85);
(sjet115){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_jet115_HT.dat c0_two_jet_inc_jet115_HT.dat;
  ALIAS_NAME jet115R32;
  ## ADOPT_BINS jet115D32[0];
  DRAW_LEGEND NO;
}(sjet115);

(djet){
  PIECE_SETUP rc.C (djet20){ }(djet20);
  PIECE_SETUP rc.C (djet30){ }(djet30);
  PIECE_SETUP rc.C (djet50){ }(djet50);
  PIECE_SETUP rc.C (djet85){ }(djet85);
  PIECE_SETUP rc.C (djet115){ }(djet115);
  DATA_TYPE ASCII;
  DRAW_OPTION P;
  X_VALUE 1;
  Y_VALUE 3;
  X_ERROR_PLUS 2;
  X_ERROR_MINUS 2;
  Y_ERROR_PLUS 5;
  Y_ERROR_MINUS 4;
  LEG_LEFT 0.725;
  LEG_RIGHT 0.975;
  LEG_TOP 0.35;
}(djet);
(djet20){
  PATH_PIECE DATAPATH/
  FILE_PIECE jet20_R32_Data.dat;
  ALIAS_NAME jet20D32;
  MARKER_STYLE 20;
  LEGEND_TITLE Jet 20;
}(djet20);
(djet30){
  PATH_PIECE DATAPATH/
  FILE_PIECE jet30_R32_Data.dat;
  ALIAS_NAME jet30D32;
  MARKER_STYLE 24;
  LEGEND_TITLE Jet 30;
}(djet30);
(djet50){
  PATH_PIECE DATAPATH/
  FILE_PIECE jet50_R32_Data.dat;
  ALIAS_NAME jet50D32;
  LEGEND_TITLE Jet 50;
  MARKER_STYLE 21;
}(djet50);
(djet85){
  PATH_PIECE DATAPATH/
  FILE_PIECE jet85_R32_Data.dat;
  ALIAS_NAME jet85D32;
  LEGEND_TITLE Jet 85;
  MARKER_STYLE 25;
}(djet85);
(djet115){
  PATH_PIECE DATAPATH/
  FILE_PIECE jet115_R32_Data.dat;
  ALIAS_NAME jet115D32;
  LEGEND_TITLE Jet 115;
  MARKER_STYLE 22;
}(djet115);

