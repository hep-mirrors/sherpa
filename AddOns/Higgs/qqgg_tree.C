Complex qqgamgam_tree_mpmp(int i1, int i2, int i3, int i4) { return spa(i1,i2)*spa(i1,i3)*spa(i1,i3)/spa(i1,i4)/spa(i1,i2)/spa(i2,i4); }
Complex qqgamgam_tree_pmpm(int i1, int i2, int i3, int i4) { return -spb(i1,i2)*spb(i1,i3)*spb(i1,i3)/spb(i1,i4)/spb(i1,i2)/spb(i2,i4); }
Complex qqgamgam_tree_pmmp(int i1, int i2, int i3, int i4) { return -spa(i1,i2)*spa(i2,i3)*spa(i2,i3)/spa(i1,i4)/spa(i1,i2)/spa(i2,i4); }
Complex qqgamgam_tree_mppm(int i1, int i2, int i3, int i4) { return spb(i1,i2)*spb(i2,i3)*spb(i2,i3)/spb(i1,i4)/spb(i1,i2)/spb(i2,i4); }

// h2 == -h1 otherwise the amplitude vanishes
Complex qqgamgam_tree_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4) {
  if ( h3 == h4 ) { return 0.; }
  if ( h3 == 1) {
    if ( h1 == 1 ) { return qqgamgam_tree_pmpm(i1,i2,i3,i4); }
    else { return qqgamgam_tree_mppm(i1,i2,i3,i4); }
  }
  else { // h3 == -1
    if ( h1 == 1 ){ return qqgamgam_tree_pmmp(i1,i2,i3,i4); }
    else { return qqgamgam_tree_mpmp(i1,i2,i3,i4); }
  }
  return 0.;
}

// q(1) qbar(2) -> gam(3) gam(4)
Complex qqbgamgam_tree(int h1, int h3, int h4) { return qqgamgam_tree_gen(1,h1,2,3,h3,4,h4); }
// qbar(1) q(2) -> gam(3) gam(4)
Complex qbqgamgam_tree(int h2, int h3, int h4) { return qqgamgam_tree_gen(2,h2,1,3,h3,4,h4); }

