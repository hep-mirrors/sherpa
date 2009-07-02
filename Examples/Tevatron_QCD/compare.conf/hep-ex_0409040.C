(initial){
  PIECE_SETUP hep-ex_0409040.C (ptall){ }(ptall);
  PIECE_SETUP hep-ex_0409040.C (dptall){ }(dptall);
  TOP_AXIS YES; RIGHT_AXIS YES;
  DEFINE_COLOURS VIOLET1 195 0 185,;
  DEFINE_COLOURS BLUE1 25 25 205,;
  DEFINE_COLOURS GREEN1 30 145 35,;
  DEFINE_COLOURS RED1 195 15 15,;
  DEFINE_COLOURS BLUE2 15 20 180,;
  DEFINE_COLOURS GREEN2 15 125 25,;
  DEFINE_COLOURS RED2 175 5 5,;
}(initial);

// standard plots

(ptall){
  PIECE_SETUP hep-ex_0409040.C (yranges){ }(yranges);
  DRAW YES;
  X_TITLE_OFFSET 1.1; BOTTOM_MARGIN 0.11;
  Y_TITLE_OFFSET 1.6; LEFT_MARGIN 0.14;
  RIGHT_MARGIN 0.035; TOP_MARGIN 0.025;
  X_MIN M_PI/2; X_MAX M_PI;
  X_AXIS_TITLE #Delta#phi_{dijet};
  Y_MIN 1.001e-3; Y_MAX 1.999e5; Y_SCALING Log_B_10;
  Y_AXIS_TITLE 1/#sigma_{dijet} d#sigma_{dijet}/d#Delta#phi_{dijet};
  HISTOGRAM_NAME Jet_DPhi;
  FIGURE_CAPTION Jet $\Delta\phi$; WEBPAGE_CAPTION Jet &Delta\;&phi\;;
  LEG_LEFT 0.075; LEG_RIGHT 0.35; LEG_TOP 0.95; LEG_TEXT_SIZE 0.0275; LEG_DELTA_Y 0.05;
  WATERMARK parton level | SIZE .035 COLOUR RED1 LEFT .925 TOP .1 ALIGN 32 PRIORITY 10\;;
  WATERMARK SHERPA | SIZE .04 COLOUR 19 LEFT .95 TOP .05 ALIGN 32 PRIORITY -10\;;
}(ptall);

(yranges){ 
  PIECE_SETUP hep-ex_0409040.C (0.0-0.1){ }(0.0-0.1);
  PIECE_SETUP hep-ex_0409040.C (0.1-0.7){ }(0.1-0.7);
  PIECE_SETUP hep-ex_0409040.C (0.7-1.1){ }(0.7-1.1);
  PIECE_SETUP hep-ex_0409040.C (1.1-1.6){ }(1.1-1.6);
}(yranges);

(0.0-0.1){ 
  PIECE_SETUP hep-ex_0409040.C (samples){ }(samples);
  FILE_PIECE TwoDPhi_j-0_j-1_C.dat; MARKER_STYLE 20;
  @@ DATSFAC 8000;
  LEGEND_ENABLED YES; @@ LEGENDTITLE p_{T}^{max} > 180 GeV (x DATSFAC)
}(0.0-0.1);
(0.1-0.7){ 
  PIECE_SETUP hep-ex_0409040.C (samples){ }(samples);
  FILE_PIECE TwoDPhi_j-0_j-1_B.dat; MARKER_STYLE 24;
  FILE_PIECE KTJets_0.7-1.1_jet_1_1_pt_0.dat;
  @@ DATSFAC 400;
  @@ LEGENDTITLE 130 GeV < p_{T}^{max} < 180 GeV (x DATSFAC)
}(0.1-0.7);
(0.7-1.1){ 
  PIECE_SETUP hep-ex_0409040.C (samples){ }(samples);
  FILE_PIECE TwoDPhi_j-0_j-1_A.dat; MARKER_STYLE 21;
  FILE_PIECE KTJets_0.1-0.7_jet_1_1_pt_0.dat;
  @@ DATSFAC 20;
  @@ LEGENDTITLE 100 GeV < p_{T}^{max} < 130 GeV (x DATSFAC)
}(0.7-1.1);
(1.1-1.6){ 
  PIECE_SETUP hep-ex_0409040.C (samples){ }(samples);
  FILE_PIECE TwoDPhi_j-0_j-1.dat; MARKER_STYLE 25;
  @@ DATSFAC 1;
  @@ LEGENDTITLE 75 GeV < p_{T}^{max} < 100 GeV
}(1.1-1.6);

(samples){
  PIECE_SETUP hep-ex_0409040.C (data){ }(data);
  PIECE_SETUP hep-ex_0409040.C (bfkl){ }(bfkl);
}(samples);

(data){ 
  DATA_TYPE ASCII;
  PATH_PIECE data/hep-ex_0409040/;
  Y_FUNCTION y*DATSFAC;
  LEGEND_TITLE LEGENDTITLE; LEGEND_ENABLED YES;
  DRAW_OPTION P; DRAW_PRIORITY 20;
  X_VALUE 1; Y_VALUE 2;
  X_ERROR_MINUS 3; X_ERROR_PLUS 4;
  Y_ERROR_MINUS 5; Y_ERROR_PLUS 6;
}(data);
(bfkl){
  PIECE_SETUP hep-ex_0409040.C (paths){ }(paths);
  READER_PARAMETERS ADOPT_BINS ./data/hep-ex_0409040/CURRENT_FILE_PIECE;
  LEGEND_ENABLED NO;
}(bfkl);

(paths){
  PIECE_SETUP hep-ex_0409040.C (path1){ }(path1);
  PIECE_SETUP hep-ex_0409040.C (path2){ }(path2);
}(paths);

(path1){
  PIECE_SETUP hep-ex_0409040.C (jets){ }(jets);
  PATH_PIECE BPATH1/hep-ex_0409040/; LINE_STYLE 1; @@ SUBJ 0;
  @@ LTITLE BTITLE1; LINE_COLOUR RED1;
}(path1);
(path2){
  PIECE_SETUP hep-ex_0409040.C (jets){ }(jets);
  PATH_PIECE BPATH2/hep-ex_0409040/; LINE_STYLE 3; @@ SUBJ 0;
  @@ LTITLE BTITLE2; LINE_COLOUR BLUE1;
}(path2);

(jets){ 
  if (SUBJ) PIECE_SETUP hep-ex_0409040.C (j[2[+1]5]){ }(j[2[+1]5]);
  PIECE_SETUP hep-ex_0409040.C (sum){ }(sum);
  Y_FUNCTION y*DATSFAC/HISTOGRAM_NORM;
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
  LINE_STYLE 2;
  LEGEND_TITLE 2-jet;
}(j2);
(j3){ 
  PATH_PIECE j3/;
  DATA_TYPE ATOOLS;
  LINE_STYLE 3;
  LEGEND_TITLE 3-jet;
}(j3);
(j4){ 
  PATH_PIECE j4/;
  DATA_TYPE ATOOLS;
  LINE_STYLE 4;
  LEGEND_TITLE 4-jet;
}(j4);
(j5){ 
  PATH_PIECE j5/;
  DATA_TYPE ATOOLS;
  LINE_STYLE 5;
  LEGEND_TITLE 5-jet;
}(j5);

// diff plots

(dptall){
  PIECE_SETUP hep-ex_0409040.C (dyranges){ }(dyranges);
  X_MIN M_PI/2; X_MAX M_PI;
  X_AXIS_TITLE #Delta#phi_{dijet};
  Y_MIN -0.999; Y_MAX 0.999; Y_SCALING Id;
  Y_AXIS_TITLE Theory / Data - 1;
  HISTOGRAM_NAME Jet_DPhi_D;
  FIGURE_CAPTION Jet $\Delta\phi$; WEBPAGE_CAPTION Jet &Delta\;&phi\;;
  LEG_LEFT 0.65; LEG_RIGHT 0.85; LEG_TOP 0.95; LEG_TEXT_SIZE 0.03;
  Y_AXIS_NDIVISIONS 305; Y_AXIS_LABEL_SIZE 0.04;
  X_TITLE_OFFSET 1.1; Y_TITLE_OFFSET 1.6; 
  @@ GBM 0.11; @@ GTM 0.025;
  @@ JMAX 4; @@ SDIV (1-GBM-GTM)/JMAX;
  LEFT_MARGIN 0.14; RIGHT_MARGIN 0.035;
}(dptall);

(dyranges){
  PIECE_SETUP hep-ex_0409040.C (d0.0-0.1){ }(d0.0-0.1);
  PIECE_SETUP hep-ex_0409040.C (d0.1-0.7){ }(d0.1-0.7);
  PIECE_SETUP hep-ex_0409040.C (d0.7-1.1){ }(d0.7-1.1);
  PIECE_SETUP hep-ex_0409040.C (d1.1-1.6){ }(d1.1-1.6);
}(dyranges);

(d0.0-0.1){ 
  PIECE_SETUP hep-ex_0409040.C (dsamples){ }(dsamples);
  FILE_PIECE TwoDPhi_j-0_j-1_C.dat;
  @@ DATSFAC 1e6;
  DRAW YES; TOP_MARGIN GTM;
  BOTTOM_MARGIN GBM+(JMAX-1)*SDIV;
  Y_TITLE_SIZE 0.04; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK p_{T}^{max} > 180 GeV | COLOUR 1 SIZE 0.025 LEFT 0.8 TOP 0.75 ALIGN 22 PRIORITY 20\;;
}(d0.0-0.1);
(d0.1-0.7){ 
  PIECE_SETUP hep-ex_0409040.C (dsamples){ }(dsamples);
  FILE_PIECE TwoDPhi_j-0_j-1_B.dat;
  @@ DATSFAC 1e3;
  DRAW YES; DIFF_PLOT YES; TOP_MARGIN GTM+SDIV;
  BOTTOM_MARGIN GBM+(JMAX-2)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK 130 GeV < p_{T}^{max} < 180 GeV | COLOUR 1 SIZE 0.025 LEFT 0.8 TOP 0.75 ALIGN 22 PRIORITY 20\;;
}(d0.1-0.7);
(d0.7-1.1){ 
  PIECE_SETUP hep-ex_0409040.C (dsamples){ }(dsamples);
  FILE_PIECE TwoDPhi_j-0_j-1_A.dat;
  @@ DATSFAC 1;
  DRAW YES; DIFF_PLOT YES; TOP_MARGIN GTM+2*SDIV;
  BOTTOM_MARGIN GBM+(JMAX-3)*SDIV;
  Y_TITLE_SIZE 0; X_TITLE_SIZE 0; X_AXIS_LABEL_SIZE 0;
  WATERMARK 100 GeV < p_{T}^{max} < 130 GeV | COLOUR 1 SIZE 0.025 LEFT 0.8 TOP 0.75 ALIGN 22 PRIORITY 20\;;
}(d0.7-1.1);
(d1.1-1.6){ 
  PIECE_SETUP hep-ex_0409040.C (dsamples){ }(dsamples);
  FILE_PIECE TwoDPhi_j-0_j-1.dat;
  @@ DATSFAC 1e-3;
  DRAW YES; DIFF_PLOT YES;
  TOP_MARGIN GTM+3*SDIV;
  BOTTOM_MARGIN GBM+(JMAX-4)*SDIV;
  Y_TITLE_SIZE 0;  X_TITLE_SIZE 0.04; X_AXIS_LABEL_SIZE 0.04;
  WATERMARK 75 GeV < p_{T}^{max} < 100 GeV | COLOUR 1 SIZE 0.025 LEFT 0.8 TOP 0.75 ALIGN 22 PRIORITY 20\;;
  WATERMARK parton level | SIZE .035 COLOUR RED1 LEFT .975 TOP .25 ALIGN 32 PRIORITY 10\;;
  WATERMARK SHERPA | SIZE .04 COLOUR 19 LEFT .05 TOP .25 ALIGN 12 PRIORITY -20\;;
}(d1.1-1.6);


(dsamples){
  PIECE_SETUP hep-ex_0409040.C (ddata){ }(ddata);
  PIECE_SETUP hep-ex_0409040.C (ddatap){ }(ddatap);
  PIECE_SETUP hep-ex_0409040.C (dbfkl){ }(dbfkl);
  LEGEND_ENABLED NO;
}(dsamples);

(ddata){
  PATH_PIECE data/hep-ex_0409040/;
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
  PATH_PIECE data/hep-ex_0409040/;
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
  PIECE_SETUP hep-ex_0409040.C (dpaths){ }(dpaths);
  READER_PARAMETERS ADOPT_BINS ./data/hep-ex_0409040/CURRENT_FILE_PIECE;
  DATA_TYPE ALGEBRA(y[1]/y[0]-1);
  DRAW_PRIORITY 10;
}(dbfkl);

(dpaths){
  PIECE_SETUP hep-ex_0409040.C (dpath1){ }(dpath1);
  PIECE_SETUP hep-ex_0409040.C (dpath2){ }(dpath2);
}(dpaths);

(dpath1){
  PIECE_SETUP hep-ex_0409040.C (djets){ }(djets);
  @@ BPPIECE BPATH1/hep-ex_0409040/;
  LINE_COLOUR RED1;
}(dpath1);
(dpath2){
  PIECE_SETUP hep-ex_0409040.C (djets){ }(djets);
  @@ BPPIECE BPATH2/hep-ex_0409040/; LINE_STYLE 2;
  LINE_COLOUR BLUE1;
}(dpath2);

(djets){ 
  PIECE_SETUP hep-ex_0409040.C (dsum){ }(dsum);
}(djets);
(dsum){ 
  PATH_PIECE data/hep-ex_0409040/ BPPIECE;
  LINE_COLOUR 1; LINE_WIDTH 2;
  DRAW_PRIORITY 10;
  LEGEND_TITLE LTITLE;
}(dsum);
