/* =============================================================
   The resonant gg -> H g -> gamma gamma g amplitudes,
   written in spinor products with phase given by specific definition
   of angle bracket spa(...) and square bracket spb(...)
   ============================================================= */

Complex ggH_pp(int i1, int i2){ return -spb(i1,i2)*spb(i1,i2); } 

Complex ggH_mp(int i1, int i2){ return 0; }

Complex ggH_mm(int i1, int i2){ return -spa(i1,i2)*spa(i1,i2); } 

// gggH:

Complex gggH_ppp(int i1, int i2, int i3){
 return pow( sij(i1,i2)+sij(i2,i3)+sij(i3,i1), 2) /spa(i1,i2)/spa(i2,i3)/spa(i3,i1);
} 

Complex gggH_mpp(int i1, int i2, int i3){
 return -pow( spb(i2,i3), 4) /spb(i1,i2)/spb(i2,i3)/spb(i3,i1);
}

Complex gggH_pmm(int i1, int i2, int i3){
 return  pow( spa(i2,i3), 4) /spa(i1,i2)/spa(i2,i3)/spa(i3,i1);
}

Complex gggH_mmm(int i1, int i2, int i3){
 return -pow( sij(i1,i2)+sij(i2,i3)+sij(i3,i1), 2) /spb(i1,i2)/spb(i2,i3)/spb(i3,i1);
} 

// Get ggH for all helicities using symmetries: 

Complex ggH_gen(int i1, int h1, int i2, int h2){
  int hsum = h1+h2;
  if (hsum == 2){ return ggH_pp(i1,i2); }
  else if (hsum == 0){
    if (h1 == -1){ return ggH_mp(i1,i2); }
    else { return ggH_mp(i2,i1); } }
  else { return ggH_mm(i1,i2); }
}

// Get gggH for all helicities using symmetries: 

Complex gggH_gen(int i1, int h1, int i2, int h2, int i3, int h3){
  int hsum = h1+h2+h3;
  if (hsum == 3){ return gggH_ppp(i1,i2,i3); }
  else if (hsum == 1){
    if (h1 == -1){ return gggH_mpp(i1,i2,i3); }
    else if (h2 == -1){ return gggH_mpp(i2,i3,i1); }
    else { return gggH_mpp(i3,i1,i2); } }
  else if (hsum == -1){
    if (h1 == 1){ return gggH_pmm(i1,i2,i3); }
    else if (h2 == 1){ return gggH_pmm(i2,i3,i1); }
    else { return gggH_pmm(i3,i1,i2); } }
  else { return gggH_mmm(i1,i2,i3); }
}

// q(1,i1,h1) qbar(2,i2,-h1) -> g(3,i3,h3) H(4+5) for all helicities:

Complex qqgH_gen(int i1, int h1, int i2, int i3, int h3){
  if (h3 == 1) {
    if (h1 == 1) { return -spb(i1,i3)*spb(i1,i3)/spb(i1,i2); }
    else { return spb(i2,i3)*spb(i2,i3)/spb(i1,i2); } }
  else if (h3 == -1) {
    if (h1 == 1) { return -spa(i2,i3)*spa(i2,i3)/spa(i1,i2); }
    else { return spa(i1,i3)*spa(i1,i3)/spa(i1,i2); } }
  else return 0.;
}

//=================================================================
// Special versions with fewer labels:

// g(1) g(2) -> H(3+4):
Complex ggH(int h1, int h2){ return ggH_gen(1,h1,2,h2); }

// g(1) g(2) -> H(3+4)  g(5):
Complex ggHg(int h1, int h2, int h5){ return gggH_gen(1,h1,2,h2,5,h5); }

// q(1) qbar(2) -> H(3+4)  g(5):
Complex qqbHg(int h1, int h5){ return qqgH_gen(1,h1,2,5,h5); }

// qbar(1) q(2) -> H(3+4)  g(5):
Complex qbqHg(int h2, int h5){ return qqgH_gen(2,h2,1,5,h5); }

// q(1) g(2) -> H(3+4)  q(5):
Complex qgHq(int h1, int h2){ return qqgH_gen(1,h1,5,2,h2); }

// g(1) q(2) -> H(3+4)  q(5):
Complex gqHq(int h1, int h2){ return qqgH_gen(2,h2,5,1,h1); }

//============ DECAY TO gamma gamma ======================================

// H(3+4) -> gamma(3) + gamma(4) :
Complex Hgamgam(int h3, int h4){ return ggH_gen(3,h3,4,h4); }
