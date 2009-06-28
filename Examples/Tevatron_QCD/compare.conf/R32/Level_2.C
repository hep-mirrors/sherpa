(level2){
  PIECE_SETUP Level_2.C (.){ }(.);
  if (DIFFPLOT!=true) {
    PIECE_SETUP Level_2.C (j2){ }(j2);
    PIECE_SETUP Level_2.C (j3){ }(j3);
    PIECE_SETUP Level_2.C (j4){ }(j4);
  }
}(level2);

(j2){
  PIECE_SETUP Final.C (final){ }(final);
  PATH_PIECE j2/;
  LINE_STYLE 2;
  DRAW_LEGEND NO;
}(j2);

(j3){
  PIECE_SETUP Final.C (final){ }(final);
  PATH_PIECE j3/;
  LINE_STYLE 3;
  DRAW_LEGEND NO;
}(j3);

(j4){
  PIECE_SETUP Final.C (final){ }(final);
  PATH_PIECE j4/;
  LINE_STYLE 4;
  DRAW_LEGEND NO;
}(j4);

(.){
  PIECE_SETUP Final.C (final){ }(final);
  LINE_STYLE 1;
}(.);
