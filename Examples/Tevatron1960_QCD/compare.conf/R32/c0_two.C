(c0_two){
  PIECE_SETUP c0_two.C (jet20){ }(jet20);
  PIECE_SETUP c0_two.C (jet30){ }(jet30);
  PIECE_SETUP c0_two.C (jet50){ }(jet50);
  PIECE_SETUP c0_two.C (jet85){ }(jet85);
  PIECE_SETUP c0_two.C (jet115){ }(jet115);
  DRAW YES;
  X_MIN 100;
  X_MAX 600;
  Y_MAX 1e5;
  Y_MIN 1e-2;
  LEG_LEFT 0.65;
  LEG_RIGHT 0.9;
  Y_SCALING Log_B_10;
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  Y_AXIS_TITLE d#sigma_{2}/dH_{T} #left[ pb/GeV #right];
  FIGURE_CAPTION $d\sigma_{2}/dH_{T}$;
  WEBPAGE_CAPTION d&sigma<sub>2</sub>/dH<sub>T</sub>;
  HISTOGRAM_NAME DR3_DHT;
  DRAW_LATEX Jet 20 Trigger | COLOUR 1 LEFT 0.07 TOP 0.05 ALIGN 12 ANGLE 90 PRIORITY 20 SIZE 0.03\;;
  DRAW_LATEX Jet 30 Trigger | COLOUR 1 LEFT 0.21 TOP 0.05 ALIGN 12 ANGLE 90 PRIORITY 20 SIZE 0.03\;;
  DRAW_LATEX Jet 50 Trigger | COLOUR 1 LEFT 0.37 TOP 0.05 ALIGN 12 ANGLE 90 PRIORITY 20 SIZE 0.03\;;
  DRAW_LATEX Jet 85 Trigger | COLOUR 1 LEFT 0.56 TOP 0.05 ALIGN 12 ANGLE 90 PRIORITY 20 SIZE 0.03\;;
  DRAW_LATEX Jet 115 Trigger | COLOUR 1 LEFT 0.83 TOP 0.05 ALIGN 12 ANGLE 90 PRIORITY 20 SIZE 0.03\;;
  DRAW_LINE V 170 | STYLE 1 PRIORITY 10\;;
  DRAW_LINE V 240 | STYLE 1 PRIORITY 10\;;
  DRAW_LINE V 330 | STYLE 1 PRIORITY 10\;;
  DRAW_LINE V 430 | STYLE 1 PRIORITY 10\;;
}(c0_two);
(jet20){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_jet20_HT.dat;
}(jet20);
(jet30){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_jet30_HT.dat;
  DRAW_LEGEND NO;
}(jet30);
(jet50){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_jet50_HT.dat;
  DRAW_LEGEND NO;
}(jet50);
(jet85){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_jet85_HT.dat;
  DRAW_LEGEND NO;
}(jet85);
(jet115){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_jet115_HT.dat;
  DRAW_LEGEND NO;
}(jet115);

(uc0_two){
  PIECE_SETUP Level_1.C (level1){ }(level1);
  FILE_PIECE c0_two_jet_inc_HT.dat;
  DRAW YES;
  X_MIN 100;
  X_MAX 600;
  Y_MAX 1e5;
  Y_MIN 1e-2;
  LEG_LEFT 0.65;
  LEG_RIGHT 0.9;
  Y_SCALING Log_B_10;
  X_AXIS_TITLE H_{T} #left[ GeV #right];
  Y_AXIS_TITLE d#sigma_{2}/dH_{T} #left[ pb/GeV #right];
  FIGURE_CAPTION $d\sigma_{2}/dH_{T}$ Untriggered;
  WEBPAGE_CAPTION d&sigma<sub>2</sub>/dH<sub>T</sub> Untriggered;
  HISTOGRAM_NAME DR2_DHT;
  DRAW_LATEX Untriggered | COLOUR 2 LEFT 0.5 TOP 0.15 ALIGN 22 PRIORITY 20\;;
}(uc0_two);

