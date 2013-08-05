//=== FIVE-POINT AMPLITUDES ===

//the fermion loop for (-,+,+,+,+):
Complex qqggamgam_mpppp(int i1, int i2, int i3, int i4, int i5){
 return -1./sij(i1,i2)/spa(i4,i5)
 * ( spa(i1,i5)*spb(i2,i4)*spb(i3,i5)/spa(i3,i5)
   - spa(i1,i4)*spb(i2,i5)*spb(i3,i4)/spa(i3,i4) ) ;
}
//the fermion loop for (-,+,-,+,+):
Complex qqggamgam_mpmpp(int i1, int i2, int i3, int i4, int i5){
 return -2. * spa(i1,i4)*spa(i1,i5)*spa(i3,i4)*spa(i3,i5)*spb(i4,i5)*spb(i4,i5)*spb(i4,i5)*spb(i4,i5)/spa(i1,i2)/pow(sij(i1,i2),4)
   * ( Ls3(sij(i3,i4),sij(i1,i2),sij(i3,i5),sij(i1,i2))
     - 1./3. * sij(i1,i2)/sij(i4,i5) * ( L2(sij(i3,i4),sij(i1,i2)) + L2(sij(i3,i5),sij(i1,i2)) ) )
- spa(i1,i3)*spa(i1,i3)*spb(i4,i5)*spb(i4,i5)/spa(i1,i2) * Ls1(sij(i3,i4),sij(i1,i2),sij(i3,i5),sij(i1,i2))/pow(sij(i1,i2),2)
+ spa(i1,i3)*spa(i1,i5)*spa(i3,i4)*spb(i4,i5)*spb(i4,i5)/spa(i4,i5)/spa(i1,i2)
   * L1(sij(i3,i4),sij(i1,i2))/pow(sij(i1,i2),2)
- spa(i1,i3)*spa(i1,i4)*spa(i3,i5)*spb(i4,i5)*spb(i4,i5)/spa(i4,i5)/spa(i1,i2)
   * L1(sij(i3,i5),sij(i1,i2))/pow(sij(i1,i2),2)
- spa(i1,i3)*spa(i1,i3)*spb(i4,i5)/spa(i4,i5)/spa(i1,i2)/sij(i1,i2)
   * ( L0(sij(i3,i4),sij(i1,i2)) + L0(sij(i3,i5),sij(i1,i2)) )
+ 1./3. * spa(i1,i4)*spa(i1,i5)*spb(i4,i5)*spb(i4,i5)*spb(i4,i5)*spb(i4,i5) /pow(sij(i1,i2),2)/spa(i1,i2)/spb(i3,i4)/spb(i3,i5)
- 1./3. * spa(i1,i3)/sij(i1,i2) * spb(i4,i5)/spa(i4,i5)
     * ( spb(i2,i4)/spb(i3,i4) + spb(i2,i5)/spb(i3,i5) + 4. * spa(i1,i3)/spa(i1,i2) )
+ 1./3. * spb(i2,i4)*spb(i2,i5)*spb(i4,i5)/spa(i4,i5)/spb(i1,i2)/spb(i3,i4)/spb(i3,i5);
}
Complex qqggamgam_pmppp(int i1, int i2, int i3, int i4, int i5){ return qqggamgam_mpppp(i2,i1,i3,i4,i5); }
Complex qqggamgam_pmmpp(int i1, int i2, int i3, int i4, int i5){ return qqggamgam_mpmpp(i2,i1,i3,i4,i5); }
//the fermion loop for (+,-,-,-,-):
Complex qqggamgam_pmmmm(int i1, int i2, int i3, int i4, int i5){
 return 1./sij(i1,i2)/spb(i4,i5)
 * ( spb(i1,i5)*spa(i2,i4)*spa(i3,i5)/spb(i3,i5)
   - spb(i1,i4)*spa(i2,i5)*spa(i3,i4)/spb(i3,i4) ) ;
}
//the fermion loop for (+,-,+,-,-):
Complex qqggamgam_pmpmm(int i1, int i2, int i3, int i4, int i5){
 return 2. * spb(i1,i4)*spb(i1,i5)*spb(i3,i4)*spb(i3,i5)*spa(i4,i5)*spa(i4,i5)*spa(i4,i5)*spa(i4,i5)/spb(i1,i2)/pow(sij(i1,i2),4)
   * ( Ls3(sij(i3,i4),sij(i1,i2),sij(i3,i5),sij(i1,i2))
     - 1./3. * sij(i1,i2)/sij(i4,i5) * ( L2(sij(i3,i4),sij(i1,i2)) + L2(sij(i3,i5),sij(i1,i2)) ) )
+ spb(i1,i3)*spb(i1,i3)*spa(i4,i5)*spa(i4,i5)/spb(i1,i2) * Ls1(sij(i3,i4),sij(i1,i2),sij(i3,i5),sij(i1,i2))/pow(sij(i1,i2),2)
- spb(i1,i3)*spb(i1,i5)*spb(i3,i4)*spa(i4,i5)*spa(i4,i5)/spb(i4,i5)/spb(i1,i2)
   * L1(sij(i3,i4),sij(i1,i2))/pow(sij(i1,i2),2)
+ spb(i1,i3)*spb(i1,i4)*spb(i3,i5)*spa(i4,i5)*spa(i4,i5)/spb(i4,i5)/spb(i1,i2)
   * L1(sij(i3,i5),sij(i1,i2))/pow(sij(i1,i2),2)
+ spb(i1,i3)*spb(i1,i3)*spa(i4,i5)/spb(i4,i5)/spb(i1,i2)/sij(i1,i2)
   * ( L0(sij(i3,i4),sij(i1,i2)) + L0(sij(i3,i5),sij(i1,i2)) )
- 1./3. * spb(i1,i4)*spb(i1,i5)*spa(i4,i5)*spa(i4,i5)*spa(i4,i5)*spa(i4,i5) /pow(sij(i1,i2),2)/spb(i1,i2)/spa(i3,i4)/spa(i3,i5)
+ 1./3. * spb(i1,i3)/sij(i1,i2) * spa(i4,i5)/spb(i4,i5)
     * ( spa(i2,i4)/spa(i3,i4) + spa(i2,i5)/spa(i3,i5) + 4. * spb(i1,i3)/spb(i1,i2) )
- 1./3. * spa(i2,i4)*spa(i2,i5)*spa(i4,i5)/spb(i4,i5)/spa(i1,i2)/spa(i3,i4)/spa(i3,i5);
}
Complex qqggamgam_mpmmm(int i1, int i2, int i3, int i4, int i5){ return qqggamgam_pmmmm(i2,i1,i3,i4,i5); }
Complex qqggamgam_mppmm(int i1, int i2, int i3, int i4, int i5){ return qqggamgam_pmpmm(i2,i1,i3,i4,i5); }

// h2 == -h1 otherwise the amplitude vanishes
Complex qqggamgam_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4, int i5, int h5) {
  if ( h3 == h4 && h4 == h5 ) {
    if ( h3 == 1) {
      if ( h1 == 1 ) { return qqggamgam_pmppp(i1,i2,i3,i4,i5); }
      else { return qqggamgam_mpppp(i1,i2,i3,i4,i5); }}
    else { // h3 == -1
      if ( h1 == 1 ){ return qqggamgam_pmmmm(i1,i2,i3,i4,i5); }
      else { return qqggamgam_mpmmm(i1,i2,i3,i4,i5); } }
  }
  else if ( h4 == h5 ) {
    if ( h3 == 1) {
      if ( h1 == 1 ) { return qqggamgam_pmpmm(i1,i2,i3,i4,i5); }
      else { return qqggamgam_mppmm(i1,i2,i3,i4,i5); }}
    else { // h3 == -1
      if ( h1 == 1 ){ return qqggamgam_pmmpp(i1,i2,i3,i4,i5); }
      else { return qqggamgam_mpmpp(i1,i2,i3,i4,i5); } }
  }
  else if ( h3 == h4 ) { return qqggamgam_gen(i1,h1,i2,i5,h5,i3,h3,i4,h4); }
  else if ( h3 == h5 ) { return qqggamgam_gen(i1,h1,i2,i4,h4,i3,h3,i5,h5); }
  else return 0.;
}

//=============================================================================

// q(1) qbar(2) -> gam(3) gam(4) g(5)
Complex qqbgamgamg(int h1, int h3, int h4, int h5) { return qqggamgam_gen(1,h1,2,5,h5,3,h3,4,h4); }
// qbar(1) q(2) -> gam(3) gam(4) g(5)
Complex qbqgamgamg(int h2, int h3, int h4, int h5) { return qqggamgam_gen(2,h2,1,5,h5,3,h3,4,h4); }
// q(1) g(2) -> gam(3) gam(4) q(5)
Complex qggamgamq(int h1, int h2, int h3, int h4) { return qqggamgam_gen(1,h1,5,2,h2,3,h3,4,h4); }
// g(1) q(2) -> gam(3) gam(4) q(5)
Complex gqgamgamq(int h1, int h2, int h3, int h4) { return qqggamgam_gen(2,h2,5,1,h1,3,h3,4,h4); }

