(initial){
  PIECE_SETUP hep-ex_0701051.C (ptall){ }(ptall);
  PIECE_SETUP hep-ex_0701051.C (dptall){ }(dptall);
  PIECE_SETUP hep-ex_0701051.C (logptall){ }(logptall);
  PIECE_SETUP hep-ex_0701051.C (dlogptall){ }(dlogptall);
  TOP_AXIS YES; RIGHT_AXIS YES;
  DEFINE_COLOURS VIOLET1 195 0 185,;
  DEFINE_COLOURS BLUE1 25 25 205,;
  DEFINE_COLOURS GREEN1 30 145 35,;
  DEFINE_COLOURS YELLOW1 240 215 10,;
  DEFINE_COLOURS RED1 195 15 15,;
  DEFINE_COLOURS BLUE2 15 20 180,;
  DEFINE_COLOURS GREEN2 15 125 25,;
  DEFINE_COLOURS RED2 175 5 5,;
}(initial);

// standard plots

(ptall){
  PIECE_SETUP hep-ex_0701051.C (yranges){ }(yranges);
  DRAW YES;
  X_TITLE_OFFSET 1.1; BOTTOM_MARGIN 0.11;
  Y_TITLE_OFFSET 1.6; LEFT_MARGIN 0.14;
  RIGHT_MARGIN 0.035; TOP_MARGIN 0.025;
  X_MIN 0.001; X_MAX 799.999;
  X_AXIS_TITLE k_{#perp }^{jet} #left[ GeV #right]
  Y_MIN 1.001e-15; Y_MAX 0.999e9; Y_SCALING Log_B_10;
  @@ YSFAC 1000;
  Y_AXIS_TITLE d#sigma / dk_{#perp }^{jet}dy^{jet} #left[ nb/GeV #right]
  HISTOGRAM_NAME Jet_PT;
  FIGURE_CAPTION Jet $P_T$; WEBPAGE_CAPTION Jet P<sub>T</sub>;
  LEG_LEFT 0.65; LEG_RIGHT 0.85; LEG_TOP 0.925; LEG_TEXT_SIZE 0.03;
  WATERMARK K_{T }\, D = 0.7 | COLOUR 1 SIZE 0.03 LEFT 0.25 TOP 0.9 ALIGN 13\;;
  WATERMARK \\\|y^{jet }\\\|<0.1 (x10^{6 }) | COLOUR 1 SIZE 0.025 LEFT 0.775 TOP 0.525 ALIGN 22\;;
  WATERMARK 0.1<\\\|y^{jet }\\\|<0.7 (x10^{3 }) | COLOUR 1 SIZE 0.025 LEFT 0.775 TOP 0.4 ALIGN 22\;;
  WATERMARK 0.7<\\\|y^{jet }\\\|<1.1 | COLOUR 1 SIZE 0.025 LEFT 0.625 TOP 0.285 ALIGN 22\;;
  WATERMARK 1.1<\\\|y^{jet }\\\|<1.6 (x10^{-3 }) | COLOUR 1 SIZE 0.025 LEFT 0.475 TOP 0.2 ALIGN 22\;;
  WATERMARK 1.6<\\\|y^{jet }\\\|<2.1 (x10^{-6 }) | COLOUR 1 SIZE 0.025 LEFT 0.35 TOP 0.05 ALIGN 22\;;
  WATERMARK parton level | SIZE .035 COLOUR RED1 LEFT .925 TOP .1 ALIGN 32 PRIORITY 10\;;
  WATERMARK SHERPA | SIZE .04 COLOUR 19 LEFT .95 TOP .05 ALIGN 32 PRIORITY -10\;;
}(ptall);

(logptall){
  PIECE_SETUP hep-ex_0701051.C (yranges){ }(yranges);
  DRAW YES;
  X_TITLE_OFFSET 1.1; BOTTOM_MARGIN 0.11;
  Y_TITLE_OFFSET 1.6; LEFT_MARGIN 0.14;
  RIGHT_MARGIN 0.035; TOP_MARGIN 0.025;
  X_MIN 45; X_MAX 1500; X_SCALING Log_B_10;
  X_AXIS_TITLE k_{#perp }^{jet} #left[ GeV #right]
  Y_MIN 1.001e-15; Y_MAX 0.999e9; Y_SCALING Log_B_10;
  @@ YSFAC 1000;
  Y_AXIS_TITLE d#sigma / dk_{#perp }^{jet}dy^{jet} #left[ nb/GeV #right]
  HISTOGRAM_NAME Jet_PT;
  FIGURE_CAPTION Jet $P_T$; WEBPAGE_CAPTION Jet P<sub>T</sub>;
  LEG_LEFT 0.65; LEG_RIGHT 0.85; LEG_TOP 0.95; LEG_TEXT_SIZE 0.03;
  WATERMARK K_{T }\, D = 0.7 | COLOUR 1 SIZE 0.03 LEFT 0.35 TOP 0.925 ALIGN 13\;;
  WATERMARK \\\|y^{jet }\\\|<0.1 (x10^{6 }) | COLOUR 1 SIZE 0.025 LEFT 0.825 TOP 0.525 ALIGN 22\;;
  WATERMARK 0.1<\\\|y^{jet }\\\|<0.7 (x10^{3 }) | COLOUR 1 SIZE 0.025 LEFT 0.825 TOP 0.4 ALIGN 22\;;
  WATERMARK 0.7<\\\|y^{jet }\\\|<1.1 | COLOUR 1 SIZE 0.025 LEFT 0.79 TOP 0.285 ALIGN 22\;;
  WATERMARK 1.1<\\\|y^{jet }\\\|<1.6 (x10^{-3 }) | COLOUR 1 SIZE 0.025 LEFT 0.76 TOP 0.2 ALIGN 22\;;
  WATERMARK 1.6<\\\|y^{jet }\\\|<2.1 (x10^{-6 }) | COLOUR 1 SIZE 0.025 LEFT 0.675 TOP 0.05 ALIGN 22\;;
  WATERMARK parton level | SIZE .035 COLOUR RED1 LEFT .075 TOP .1 ALIGN 12 PRIORITY 10\;;
  WATERMARK SHERPA | SIZE .04 COLOUR 19 LEFT .05 TOP .05 ALIGN 12 PRIORITY -10\;;
}(logptall);

(yranges){ 
  PIECE_SETUP hep-ex_0701051.C (0.0-0.1){ }(0.0-0.1);
  PIECE_SETUP hep-ex_0701051.C (0.1-0.7){ }(0.1-0.7);
  PIECE_SETUP hep-ex_0701051.C (0.7-1.1){ }(0.7-1.1);
  PIECE_SETUP hep-ex_0701051.C (1.1-1.6){ }(1.1-1.6);
  PIECE_SETUP hep-ex_0701051.C (1.6-2.1){ }(1.6-2.1);
}(yranges);

(0.0-0.1){ 
  PIECE_SETUP hep-ex_0701051.C (samples){ }(samples);
  FILE_PIECE KTJets_0.0-0.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e6; @@ YRANGE (0.1-0.0)
  LEGEND_ENABLED YES;
}(0.0-0.1);
(0.1-0.7){ 
  PIECE_SETUP hep-ex_0701051.C (samples){ }(samples);
  FILE_PIECE KTJets_0.1-0.7_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e3; @@ YRANGE (0.7-0.1)
}(0.1-0.7);
(0.7-1.1){ 
  PIECE_SETUP hep-ex_0701051.C (samples){ }(samples);
  FILE_PIECE KTJets_0.7-1.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1; @@ YRANGE (1.1-0.7)
}(0.7-1.1);
(1.1-1.6){ 
  PIECE_SETUP hep-ex_0701051.C (samples){ }(samples);
  FILE_PIECE KTJets_1.1-1.6_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e-3; @@ YRANGE (1.6-1.1)
}(1.1-1.6);
(1.6-2.1){ 
  PIECE_SETUP hep-ex_0701051.C (samples){ }(samples);
  FILE_PIECE KTJets_1.6-2.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e-6; @@ YRANGE (2.1-1.6)
}(1.6-2.1);

(samples){
  PIECE_SETUP hep-ex_0701051.C (data){ }(data);
  PIECE_SETUP hep-ex_0701051.C (bfkl){ }(bfkl);
  LEGEND_ENABLED NO;
}(samples);

(data){ 
  DATA_TYPE ASCII;
  PATH_PIECE data/hep-ex_0701051/;
  Y_FUNCTION y*DATSFAC;
  DRAW_OPTION P; DRAW_PRIORITY 20;
  X_VALUE 1; Y_VALUE 2;
  X_ERROR_MINUS 3; X_ERROR_PLUS 4;
  Y_ERROR_MINUS 5; Y_ERROR_PLUS 6;
  LEGEND_TITLE CDF Data;
}(data);
(bfkl){
  PIECE_SETUP hep-ex_0701051.C (paths){ }(paths);
  READER_PARAMETERS ADOPT_BINS ./data/hep-ex_0701051/CURRENT_FILE_PIECE;
}(bfkl);

(paths){
  PIECE_SETUP hep-ex_0701051.C (path1){ }(path1);
  PIECE_SETUP hep-ex_0701051.C (path2){ }(path2);
}(paths);

(path1){
  PIECE_SETUP hep-ex_0701051.C (jets){ }(jets);
  PATH_PIECE BPATH1/hep-ex_0701051/; LINE_STYLE 1; @@ SUBJ 1;
  @@ LTITLE Comix #otimes CSS; @@ KFAC 1;
}(path1);
(path2){
  PIECE_SETUP hep-ex_0701051.C (jets){ }(jets);
  PATH_PIECE BPATH2/hep-ex_0701051/; LINE_STYLE 3; @@ SUBJ 0;
  @@ LTITLE CSS only; @@ KFAC 1;
}(path2);

(jets){ 
  if (SUBJ) PIECE_SETUP hep-ex_0701051.C (j[2[+1]5]){ }(j[2[+1]5]);
  PIECE_SETUP hep-ex_0701051.C (sum){ }(sum);
  Y_FUNCTION y*DATSFAC/YSFAC/YRANGE*KFAC;
  DRAW_OPTION L;
}(jets);
(sum){ 
  LINE_COLOUR 1; LINE_WIDTH 2;
  DRAW_PRIORITY 10;
  LEGEND_TITLE LTITLE;
}(sum);
(j2){ 
  PATH_PIECE j2/;
  DATA_TYPE ATOOLS;
  LINE_COLOUR RED1;
  LINE_STYLE 2;
  LEGEND_TITLE 2-jet;
}(j2);
(j3){ 
  PATH_PIECE j3/;
  DATA_TYPE ATOOLS;
  LINE_COLOUR GREEN1;
  LINE_STYLE 2;
  LEGEND_TITLE 3-jet;
}(j3);
(j4){ 
  PATH_PIECE j4/;
  DATA_TYPE ATOOLS;
  LINE_COLOUR BLUE1;
  LINE_STYLE 2;
  LEGEND_TITLE 4-jet;
}(j4);
(j5){ 
  PATH_PIECE j5/;
  DATA_TYPE ATOOLS;
  LINE_COLOUR YELLOW1;
  LINE_STYLE 2;
  LEGEND_TITLE 5-jet;
}(j5);

// diff plots

(dptall){
  PIECE_SETUP hep-ex_0701051.C (dyranges){ }(dyranges);
  X_MIN 0.001; X_MAX 799.999;
  X_AXIS_TITLE k_{#perp }^{jet} #left[ GeV #right]
  Y_MIN -0.999; Y_MAX 0.999; Y_SCALING Id;
  @@ YSFAC 1000;
  Y_AXIS_TITLE d#sigma_{th} / d#sigma_{exp} - 1;
  HISTOGRAM_NAME Jet_PT_D;
  FIGURE_CAPTION Jet $P_T$; WEBPAGE_CAPTION Jet P<sub>T</sub>;
  LEG_LEFT 0.65; LEG_RIGHT 0.85; LEG_TOP 0.95; LEG_TEXT_SIZE 0.03;
  Y_AXIS_NDIVISIONS 305; Y_AXIS_LABEL_SIZE 0.04;
  X_TITLE_OFFSET 1.1; Y_TITLE_OFFSET 1.6; 
  @@ GBM 0.11; @@ GTM 0.025;
  @@ JMAX 5; @@ SDIV (1-GBM-GTM)/JMAX;
  LEFT_MARGIN 0.14; RIGHT_MARGIN 0.035;
}(dptall);

(dlogptall){
  PIECE_SETUP hep-ex_0701051.C (dyranges){ }(dyranges);
  X_MIN 45; X_MAX 1500; X_SCALING Log_B_10;
  X_AXIS_TITLE k_{#perp }^{jet} #left[ GeV #right]
  Y_MIN -0.999; Y_MAX 0.999; Y_SCALING Id;
  @@ YSFAC 1000;
  Y_AXIS_TITLE d#sigma_{th} / d#sigma_{exp} - 1;
  HISTOGRAM_NAME Jet_PT_D;
  FIGURE_CAPTION Jet $P_T$; WEBPAGE_CAPTION Jet P<sub>T</sub>;
  LEG_LEFT 0.65; LEG_RIGHT 0.85; LEG_TOP 0.95; LEG_TEXT_SIZE 0.03;
  Y_AXIS_NDIVISIONS 305; Y_AXIS_LABEL_SIZE 0.04;
  X_TITLE_OFFSET 1.1; Y_TITLE_OFFSET 1.6; 
  @@ GBM 0.11; @@ GTM 0.025;
  @@ JMAX 5; @@ SDIV (1-GBM-GTM)/JMAX;
  LEFT_MARGIN 0.14; RIGHT_MARGIN 0.035;
}(dlogptall);

(dyranges){
  PIECE_SETUP hep-ex_0701051.C (d0.0-0.1){ }(d0.0-0.1);
  PIECE_SETUP hep-ex_0701051.C (d0.1-0.7){ }(d0.1-0.7);
  PIECE_SETUP hep-ex_0701051.C (d0.7-1.1){ }(d0.7-1.1);
  PIECE_SETUP hep-ex_0701051.C (d1.1-1.6){ }(d1.1-1.6);
  PIECE_SETUP hep-ex_0701051.C (d1.6-2.1){ }(d1.6-2.1);
}(dyranges);

(d0.0-0.1){ 
  PIECE_SETUP hep-ex_0701051.C (dsamples){ }(dsamples);
  FILE_PIECE KTJets_0.0-0.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e6; @@ YRANGE (0.1-0.0)
  DRAW YES; TOP_MARGIN GTM;
  BOTTOM_MARGIN GBM+(JMAX-1)*SDIV;
  Y_TITLE_SIZE 0.04; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK \\\|y^{jet }\\\|<0.1 | COLOUR 1 SIZE 0.025 LEFT 0.85 TOP 0.65 ALIGN 22 PRIORITY 20\;;
}(d0.0-0.1);
(d0.1-0.7){ 
  PIECE_SETUP hep-ex_0701051.C (dsamples){ }(dsamples);
  FILE_PIECE KTJets_0.1-0.7_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e3; @@ YRANGE (0.7-0.1)
  DRAW YES; DIFF_PLOT YES; TOP_MARGIN GTM+SDIV;
  BOTTOM_MARGIN GBM+(JMAX-2)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK 0.1<\\\|y^{jet }\\\|<0.7 | COLOUR 1 SIZE 0.025 LEFT 0.85 TOP 0.65 ALIGN 22 PRIORITY 20\;;
}(d0.1-0.7);
(d0.7-1.1){ 
  PIECE_SETUP hep-ex_0701051.C (dsamples){ }(dsamples);
  FILE_PIECE KTJets_0.7-1.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1; @@ YRANGE (1.1-0.7)
  DRAW YES; DIFF_PLOT YES; TOP_MARGIN GTM+2*SDIV;
  BOTTOM_MARGIN GBM+(JMAX-3)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK 0.7<\\\|y^{jet }\\\|<1.1 | COLOUR 1 SIZE 0.025 LEFT 0.85 TOP 0.65 ALIGN 22 PRIORITY 20\;;
}(d0.7-1.1);
(d1.1-1.6){ 
  PIECE_SETUP hep-ex_0701051.C (dsamples){ }(dsamples);
  FILE_PIECE KTJets_1.1-1.6_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e-3; @@ YRANGE (1.6-1.1)
  DRAW YES; DIFF_PLOT YES;
  TOP_MARGIN GTM+3*SDIV;
  BOTTOM_MARGIN GBM+(JMAX-4)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK 1.1<\\\|y^{jet }\\\|<1.6 | COLOUR 1 SIZE 0.025 LEFT 0.85 TOP 0.65 ALIGN 22 PRIORITY 20\;;
}(d1.1-1.6);
(d1.6-2.1){ 
  PIECE_SETUP hep-ex_0701051.C (dsamples){ }(dsamples);
  FILE_PIECE KTJets_1.6-2.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 1e-6; @@ YRANGE (2.1-1.6)
  DRAW YES; DIFF_PLOT YES;
  TOP_MARGIN GTM+4*SDIV;
  BOTTOM_MARGIN GBM+(JMAX-5)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0.04; X_AXIS_LABEL_SIZE 0.04;
  WATERMARK 1.6<\\\|y^{jet }\\\|<2.1 | COLOUR 1 SIZE 0.025 LEFT 0.85 TOP 0.65 ALIGN 22 PRIORITY 20\;;
  WATERMARK parton level | SIZE .035 COLOUR RED1 LEFT .975 TOP .25 ALIGN 32 PRIORITY 10\;;
  WATERMARK SHERPA | SIZE .04 COLOUR 19 LEFT .05 TOP .25 ALIGN 12 PRIORITY -20\;;
}(d1.6-2.1);


(dsamples){
  PIECE_SETUP hep-ex_0701051.C (ddata){ }(ddata);
  PIECE_SETUP hep-ex_0701051.C (ddatap){ }(ddatap);
  PIECE_SETUP hep-ex_0701051.C (dbfkl){ }(dbfkl);
  LEGEND_ENABLED NO;
}(dsamples);

(ddata){
  PATH_PIECE data/hep-ex_0701051/;
  Y_FUNCTION 0|dy/y;
  X_VALUE 1; Y_VALUE 2;
  X_ERROR_MINUS 3; X_ERROR_PLUS 4;
  Y_ERROR_MINUS 5; Y_ERROR_PLUS 6;
  DRAW_OPTION P2;
  FILL_COLOUR 5;
  MARKER_COLOUR 5;
  DRAW_PRIORITY -5;
  DRAW_LINES H 0 | STYLE 1 COLOUR 1\;;
}(ddata);
(ddatap){
  PATH_PIECE data/hep-ex_0701051/;
  Y_FUNCTION 0|dy/y;
  X_VALUE 1; Y_VALUE 2;
  X_ERROR_MINUS 3; X_ERROR_PLUS 4;
  Y_ERROR_MINUS 5; Y_ERROR_PLUS 6;
  DRAW_OPTION P;
  MARKER_COLOUR 12;
  LINE_COLOUR 12;
  DRAW_PRIORITY 5;
  DRAW_LINES H 0 | STYLE 1 COLOUR 1\;;
}(ddatap);
(dbfkl){
  PIECE_SETUP hep-ex_0701051.C (dpaths){ }(dpaths);
  READER_PARAMETERS ADOPT_BINS ./data/hep-ex_0701051/CURRENT_FILE_PIECE;
  DATA_TYPE ALGEBRA(y[1]/y[0]-1);
  DRAW_PRIORITY 10;
}(dbfkl);

(dpaths){
  PIECE_SETUP hep-ex_0701051.C (dpath1){ }(dpath1);
  PIECE_SETUP hep-ex_0701051.C (dpath2){ }(dpath2);
}(dpaths);

(dpath1){
  PIECE_SETUP hep-ex_0701051.C (djets){ }(djets);
  @@ BPPIECE BPATH1/hep-ex_0701051/; LINE_STYLE 1;
}(dpath1);
(dpath2){
  PIECE_SETUP hep-ex_0701051.C (djets){ }(djets);
  @@ BPPIECE BPATH2/hep-ex_0701051/; LINE_STYLE 3;
}(dpath2);

(djets){ 
  PIECE_SETUP hep-ex_0701051.C (dsum){ }(dsum);
}(djets);
(dsum){ 
  PATH_PIECE data/hep-ex_0701051/ BPPIECE;
  LINE_COLOUR 1; LINE_WIDTH 2;
  DRAW_PRIORITY 10;
  LEGEND_TITLE LTITLE;
}(dsum);
