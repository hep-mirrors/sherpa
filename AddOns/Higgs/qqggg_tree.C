Complex qqggamgam_tree_mpmpp(int i1, int i2, int i3, int i4, int i5) { return spa(i1,i2)*spa(i1,i3)*spa(i1,i3)/spa(i1,i4)/spa(i1,i5)/spa(i2,i4)/spa(i2,i5); }
Complex qqggamgam_tree_pmpmm(int i1, int i2, int i3, int i4, int i5) { return -spb(i1,i2)*spb(i1,i3)*spb(i1,i3)/spb(i1,i4)/spb(i1,i5)/spb(i2,i4)/spb(i2,i5); }
Complex qqggamgam_tree_pmmpp(int i1, int i2, int i3, int i4, int i5) { return -spa(i1,i2)*spa(i2,i3)*spa(i2,i3)/spa(i1,i4)/spa(i1,i5)/spa(i2,i4)/spa(i2,i5); }
Complex qqggamgam_tree_mppmm(int i1, int i2, int i3, int i4, int i5) { return spb(i1,i2)*spb(i2,i3)*spb(i2,i3)/spb(i1,i4)/spb(i1,i5)/spb(i2,i4)/spb(i2,i5); }

// h2 == -h1 otherwise the amplitude vanishes
Complex qqggamgam_tree_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4, int i5, int h5) {
  if ( h3 == h4 && h4 == h5 ) { return 0.; }
  else if ( h4 == h5 ) {
    if ( h3 == 1) {
      if ( h1 == 1 ) { return qqggamgam_tree_pmpmm(i1,i2,i3,i4,i5); }
      else { return qqggamgam_tree_mppmm(i1,i2,i3,i4,i5); }}
    else { // h3 == -1
      if ( h1 == 1 ){ return qqggamgam_tree_pmmpp(i1,i2,i3,i4,i5); }
      else { return qqggamgam_tree_mpmpp(i1,i2,i3,i4,i5); } }
  }
  else if ( h3 == h4 ) { return qqggamgam_tree_gen(i1,h1,i2,i5,h5,i3,h3,i4,h4); }
  else if ( h3 == h5 ) { return qqggamgam_tree_gen(i1,h1,i2,i4,h4,i3,h3,i5,h5); }
  else return 0.;
}

// q(1) qbar(2) -> gam(3) gam(4) g(5)
Complex qqbgamgamg_tree(int h1, int h3, int h4, int h5) { return qqggamgam_tree_gen(1,h1,2,5,h5,3,h3,4,h4); }
// qbar(1) q(2) -> gam(3) gam(4) g(5)
Complex qbqgamgamg_tree(int h2, int h3, int h4, int h5) { return qqggamgam_tree_gen(2,h2,1,5,h5,3,h3,4,h4); }
// q(1) g(2) -> gam(3) gam(4) q(5)
Complex qggamgamq_tree(int h1, int h2, int h3, int h4) { return qqggamgam_tree_gen(1,h1,5,2,h2,3,h3,4,h4); }
// g(1) q(2) -> gam(3) gam(4) q(5)
Complex gqgamgamq_tree(int h1, int h2, int h3, int h4) { return qqggamgam_tree_gen(2,h2,5,1,h1,3,h3,4,h4); }

