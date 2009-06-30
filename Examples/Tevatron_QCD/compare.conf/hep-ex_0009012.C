(initial){
  PIECE_SETUP hep-ex_0009012.C (observables3){ }(observables3);
  PIECE_SETUP hep-ex_0009012.C (observables2){ }(observables2);
  PIECE_SETUP hep-ex_0009012.C (observables1){ }(observables1);

  RIGHT_AXIS YES;
  TOP_AXIS YES;

  DEFINE_COLOURS VIOLET1 195 0 185;
  DEFINE_COLOURS BLUE1 25 25 205;
  DEFINE_COLOURS GREEN1 30 145 35;
  DEFINE_COLOURS RED1 195 15 15;
  DEFINE_COLOURS YELLOW1 240 215 10,;
  DEFINE_COLOURS BLUE2 15 20 180;
  DEFINE_COLOURS GREEN2 15 125 25;
  DEFINE_COLOURS RED2 175 5 5;
}(initial);

(observables1){
  PIECE_SETUP hep-ex_0009012.C (yrc0){ }(yrc0);
  PIECE_SETUP hep-ex_0009012.C (ydc0){ }(ydc0);
  LEFT_MARGIN 0.15;
  Y_TITLE_OFFSET 1.4;
  X_TITLE_OFFSET 1.1;
  @@ PATHPIECE hep-ex_0009012/;
  LEG_DELTAY 0.04;
  LEG_TEXT_SIZE 0.03;
}(observables1);

(yrc0){
  PIECE_SETUP hep-ex_0009012.C (ydr32){ }(ydr32);
  PIECE_SETUP hep-ex_0009012.C (ysr32){ }(ysr32);
  X_MIN 100;
  X_MAX 600;
  Y_MIN 0.001;
  Y_MAX 0.7;
  DRAW YES;
  LEG_LEFT 0.725; LEG_RIGHT 0.975; LEG_TOP 0.35;
  Y_AXIS_TITLE R_{32}(H_{T});
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  BOTTOM_MARGIN 0.3;
  X_AXIS_LABEL_SIZE 0;
  X_TITLE_SIZE 0;
  HISTOGRAM_NAME DR32_DHT;
  FIGURE_CAPTION $R_{32}(H_{T})$;
  WEBPAGE_CAPTION R<sub>32</sub>(H<sub>T</sub>) Untriggered;
  DRAW_LATEX SHERPA | LEFT 0.05 TOP 0.075 ALIGN 12 PRIORITY -20\;;
}(yrc0);
(ydc0){
  PIECE_SETUP hep-ex_0009012.C (yddr32){ }(yddr32);
  PIECE_SETUP hep-ex_0009012.C (ydsr32+){ }(ydsr32+);
  PIECE_SETUP hep-ex_0009012.C (ydsr32-){ }(ydsr32-);
  DRAW YES;
  DRAW_LEGEND NO;
  DIFF_PLOT YES;
  TOP_MARGIN 0.7;
  Y_TITLE_SIZE 0;
  Y_AXIS_TICK_LENGTH 0.08;
  Y_AXIS_NDIVISIONS 505;
  Y_AXIS_LABEL_DIVISIONS 1;
  DRAW_LINE H 0 | STYLE 1 COLOUR 1 PRIORITY -10;
  X_MIN 100;
  X_MAX 600;
  Y_MIN -0.25;
  Y_MAX 0.25;
  DRAW YES;
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  HISTOGRAM_NAME DDR32_DHT;
  FIGURE_CAPTION $R_{32 }(H_{T})$;
  WEBPAGE_CAPTION R<sub>32</sub>(H<sub>T</sub>);
}(ydc0);

(yddr32){ 
  FILE_PIECE R32[0] D32[0];
  LINE_COLOUR RED1;
  DATA_TYPE ALGEBRA(y[0]/y[1]-1);
  ## ALGEBRA_DATA 4; ## COLUMNS 2 3;
}(yddr32);
(ydsr32+){ 
  FILE_PIECE D32[0] D32[0];
  DATA_TYPE ALGEBRA(y[0]/y[1]);
  ## ALGEBRA_DATA 4; ## COLUMNS 5 3;
  LINE_COLOUR 5;
  FILL_COLOUR 5;
  FILL_STYLE 1000;
  DRAW_PRIORITY -20;
}(ydsr32+);
(ydsr32-){ 
  FILE_PIECE D32[0] D32[0];
  DATA_TYPE ALGEBRA(-y[0]/y[1]);
  ## ALGEBRA_DATA 4; ## COLUMNS 4 3;
  LINE_COLOUR 5;
  FILL_COLOUR 5;
  FILL_STYLE 1000;
  DRAW_PRIORITY -20;
}(ydsr32-);

(ysr32){
  PIECE_SETUP hep-ex_0009012.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_HT.dat c0_two_jet_inc_HT.dat;
  DATA_TYPE ALGEBRA(y[0]/y[1]);
  ALIAS_NAME R32;
  ## ADOPT_BINS D32[0];
}(ysr32);

(ydr32){
  PATH_PIECE data/hep-ex_0009012/
  FILE_PIECE R32_Data.dat;
  ALIAS_NAME D32; LEGEND_TITLE D0 Data;
  DATA_TYPE ASCII; DRAW_OPTION P;
  X_VALUE 1; Y_VALUE 3;
  X_ERROR_PLUS 2; X_ERROR_MINUS 2;
  Y_ERROR_PLUS 5; Y_ERROR_MINUS 4;
}(ydr32);


(observables3){
  PIECE_SETUP hep-ex_0009012.C (uc0_three){ }(uc0_three);
  PIECE_SETUP hep-ex_0009012.C (uc0_two){ }(uc0_two);
  @@ PATHPIECE hep-ex_0009012/;
}(observables3);

(uc0_three){
  PIECE_SETUP hep-ex_0009012.C (level1){ }(level1);
  FILE_PIECE c0_three_jet_inc_HT.dat;
}(uc0_three);

(uc0_two){
  PIECE_SETUP hep-ex_0009012.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_HT.dat;
}(uc0_two);

(level1){
  PIECE_SETUP hep-ex_0009012.C ([1[+1]2]){ }([1[+1]2]);
}(level1);

(1){
  PIECE_SETUP hep-ex_0009012.C (final){ }(final);
  PATH_PIECE BPATH1PATHPIECE;
  LEGEND_TITLE TITLE1;
  LINE_COLOUR RED1;
}(1);
(2){
  PIECE_SETUP hep-ex_0009012.C (final){ }(final);
  PATH_PIECE BPATH2PATHPIECE;
  LEGEND_TITLE TITLE2;
  LINE_COLOUR BLUE1; LINE_STYLE 2;
}(2);

(final){ 
  DATA_TYPE ATOOLS; 
}(final);
