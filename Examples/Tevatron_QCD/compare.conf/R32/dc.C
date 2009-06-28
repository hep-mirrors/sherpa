(dc0){
  PIECE_SETUP dc.C (ddr32_1){ }(ddr32_1);
  PIECE_SETUP dc.C (ddr32_2){ }(ddr32_2);
  PIECE_SETUP dc.C (ddr32_3){ }(ddr32_3);
  PIECE_SETUP dc.C (ddr32_4){ }(ddr32_4);
  PIECE_SETUP dc.C (dr32){ }(dr32);
  PIECE_SETUP dc.C (dsr32+){ }(dsr32+);
  PIECE_SETUP dc.C (dsr32-){ }(dsr32-);
  DRAW YES;
  DRAW_LEGEND NO;
  DIFF_PLOT YES;
  TOP_MARGIN 0.7;
  Y_TITLE_SIZE 0;
  Y_AXIS_TICK_LENGTH 0.08;
  Y_AXIS_NDIVISIONS 505;
  Y_AXIS_LABEL_DIVISIONS 1;
  X_MIN 100;
  X_MAX 600;
  Y_MIN -0.25;
  Y_MAX 0.25;
  DRAW YES;
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  HISTOGRAM_NAME DDR32_DHT;
  FIGURE_CAPTION $R_{32 }(H_{T})$;
  WEBPAGE_CAPTION R<sub>32</sub>(H<sub>T</sub>);
  DRAW_LINE H 0 | STYLE 1 COLOUR 1 PRIORITY 10\;;
  DRAW_LINE V 170 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 240 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 330 | STYLE 1 COLOUR 16 PRIORITY 10\;;
  DRAW_LINE V 430 | STYLE 1 COLOUR 16 PRIORITY 10\;;
}(dc0);

(ddr32_1){
  PIECE_SETUP dc.C (ddr32_1jet20){ }(ddr32_1jet20);
  PIECE_SETUP dc.C (ddr32_1jet30){ }(ddr32_1jet30);
  PIECE_SETUP dc.C (ddr32_1jet50){ }(ddr32_1jet50);
  PIECE_SETUP dc.C (ddr32_1jet85){ }(ddr32_1jet85);
  PIECE_SETUP dc.C (ddr32_1jet115){ }(ddr32_1jet115);
  LINE_COLOUR RED1;
  DATA_TYPE ALGEBRA(y[0]/y[1]-1);
  ## ALGEBRA_DATA 4; ## COLUMNS 2 3;
}(ddr32_1);
(ddr32_1jet20){ 
  FILE_PIECE jet20R32[0] jet20D32[0];
}(ddr32_1jet20);
(ddr32_1jet30){ 
  FILE_PIECE jet30R32[0] jet30D32[0];
}(ddr32_1jet30);
(ddr32_1jet50){ 
  FILE_PIECE jet50R32[0] jet50D32[0];
}(ddr32_1jet50);
(ddr32_1jet85){ 
  FILE_PIECE jet85R32[0] jet85D32[0];
}(ddr32_1jet85);
(ddr32_1jet115){ 
  FILE_PIECE jet115R32[0] jet115D32[0];
}(ddr32_1jet115);

(ddr32_2){
  PIECE_SETUP dc.C (ddr32_2jet20){ }(ddr32_2jet20);
  PIECE_SETUP dc.C (ddr32_2jet30){ }(ddr32_2jet30);
  PIECE_SETUP dc.C (ddr32_2jet50){ }(ddr32_2jet50);
  PIECE_SETUP dc.C (ddr32_2jet85){ }(ddr32_2jet85);
  PIECE_SETUP dc.C (ddr32_2jet115){ }(ddr32_2jet115);
  LINE_COLOUR GREEN1;
  DATA_TYPE ALGEBRA(y[0]/y[1]-1);
  ## ALGEBRA_DATA 4; ## COLUMNS 2 3;
}(ddr32_2);
(ddr32_2jet20){ 
  FILE_PIECE jet20R32[1] jet20D32[0];
}(ddr32_2jet20);
(ddr32_2jet30){ 
  FILE_PIECE jet30R32[1] jet30D32[0];
}(ddr32_2jet30);
(ddr32_2jet50){ 
  FILE_PIECE jet50R32[1] jet50D32[0];
}(ddr32_2jet50);
(ddr32_2jet85){ 
  FILE_PIECE jet85R32[1] jet85D32[0];
}(ddr32_2jet85);
(ddr32_2jet115){ 
  FILE_PIECE jet115R32[1] jet115D32[0];
}(ddr32_2jet115);

(ddr32_3){
  PIECE_SETUP dc.C (ddr32_3jet20){ }(ddr32_3jet20);
  PIECE_SETUP dc.C (ddr32_3jet30){ }(ddr32_3jet30);
  PIECE_SETUP dc.C (ddr32_3jet50){ }(ddr32_3jet50);
  PIECE_SETUP dc.C (ddr32_3jet85){ }(ddr32_3jet85);
  PIECE_SETUP dc.C (ddr32_3jet115){ }(ddr32_3jet115);
  LINE_COLOUR BLUE1;
  DATA_TYPE ALGEBRA(y[0]/y[1]-1);
  ## ALGEBRA_DATA 4; ## COLUMNS 2 3;
}(ddr32_3);
(ddr32_3jet20){ 
  FILE_PIECE jet20R32[2] jet20D32[0];
}(ddr32_3jet20);
(ddr32_3jet30){ 
  FILE_PIECE jet30R32[2] jet30D32[0];
}(ddr32_3jet30);
(ddr32_3jet50){ 
  FILE_PIECE jet50R32[2] jet50D32[0];
}(ddr32_3jet50);
(ddr32_3jet85){ 
  FILE_PIECE jet85R32[2] jet85D32[0];
}(ddr32_3jet85);
(ddr32_3jet115){ 
  FILE_PIECE jet115R32[2] jet115D32[0];
}(ddr32_3jet115);

(ddr32_4){
  PIECE_SETUP dc.C (ddr32_4jet20){ }(ddr32_4jet20);
  PIECE_SETUP dc.C (ddr32_4jet30){ }(ddr32_4jet30);
  PIECE_SETUP dc.C (ddr32_4jet50){ }(ddr32_4jet50);
  PIECE_SETUP dc.C (ddr32_4jet85){ }(ddr32_4jet85);
  PIECE_SETUP dc.C (ddr32_4jet115){ }(ddr32_4jet115);
  LINE_COLOUR YELLOW1;
  DATA_TYPE ALGEBRA(y[0]/y[1]-1);
  ## ALGEBRA_DATA 4; ## COLUMNS 2 3;
}(ddr32_4);
(ddr32_4jet20){ 
  FILE_PIECE jet20R32[3] jet20D32[0];
}(ddr32_4jet20);
(ddr32_4jet30){ 
  FILE_PIECE jet30R32[3] jet30D32[0];
}(ddr32_4jet30);
(ddr32_4jet50){ 
  FILE_PIECE jet50R32[3] jet50D32[0];
}(ddr32_4jet50);
(ddr32_4jet85){ 
  FILE_PIECE jet85R32[3] jet85D32[0];
}(ddr32_4jet85);
(ddr32_4jet115){ 
  FILE_PIECE jet115R32[3] jet115D32[0];
}(ddr32_4jet115);

(dsr32+){ 
  FILE_PIECE FinalStateD32[0] FinalStateD32[0];
  DATA_TYPE ALGEBRA(y[0]/y[1]);
  ## ALGEBRA_DATA 4; ## COLUMNS 5 3;
  LINE_COLOUR 5;
  FILL_COLOUR 5;
  FILL_STYLE 1000;
  DRAW_PRIORITY -20;
}(dsr32+);
(dsr32-){ 
  FILE_PIECE FinalStateD32[0] FinalStateD32[0];
  DATA_TYPE ALGEBRA(-y[0]/y[1]);
  ## ALGEBRA_DATA 4; ## COLUMNS 4 3;
  LINE_COLOUR 5;
  FILL_COLOUR 5;
  FILL_STYLE 1000;
  DRAW_PRIORITY -20;
}(dsr32-);

