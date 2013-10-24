/* =============================================================
   The resonant gg -> X -> gamma gamma  amplitudes, where X is spin-2 graviton,
   written in spinor products with phase given by specific definition
   of angle bracket spa(...) and square bracket spb(...)
   ============================================================= */

Complex ggXgamgam_mpmp(int i1, int i2, int i3, int i4) {
  return pow(spa(i1,i3)*spb(i2,i4),2);
}

Complex ggXgamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3, int i4, int h4) {
  if ( h1==h2 || h3==h4 ) return 0;
  else if ( h1==h3 ) { 
     if ( h1==1 ) return ggXgamgam_mpmp(i2,i1,i4,i3);
     else return ggXgamgam_mpmp(i1,i2,i3,i4);
  }
  else { // h1==h4
     if ( h1==1 ) return ggXgamgam_mpmp(i2,i1,i3,i4);
     else return ggXgamgam_mpmp(i1,i2,i4,i3);
  }
}

//======

Complex qqXgamgam_mpmp(int i1, int i2, int i3, int i4) {
    return pow(spa(i1,i3),2)*spb(i2,i4)*spb(i1,i4);
}
Complex qqXgamgam_pmpm(int i1, int i2, int i3, int i4) {
    return pow(spb(i1,i3),2)*spa(i2,i4)*spa(i1,i4);
}

Complex qqXgamgam_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4) {
    if ( h3==h4 ) return 0;
    else if ( h1==h3 ) {
        if ( h1==1 ) return qqXgamgam_pmpm(i1,i2,i3,i4);
        else return qqXgamgam_mpmp(i1,i2,i3,i4);
    }
    else { // h1==h4
        if ( h1==1 ) return qqXgamgam_pmpm(i1,i2,i4,i3);
        else return qqXgamgam_mpmp(i1,i2,i4,i3);
    }
}

//======

Complex qqgXgamgam_mppmp(int i1, int i2, int i3, int i4, int i5) {
  return -pow(spa(i1,i4),3)*pow(spb(i4,i5),2)*spa(i4,i2)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1);
}
Complex qqgXgamgam_pmmpm(int i1, int i2, int i3, int i4, int i5) {
  return pow(spb(i1,i4),3)*pow(spa(i4,i5),2)*spb(i4,i2)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1);
}
Complex qqgXgamgam_pmpmp(int i1, int i2, int i3, int i4, int i5) {
  return qqgXgamgam_mppmp(i2,i1,i3,i4,i5); }
Complex qqgXgamgam_mpmpm(int i1, int i2, int i3, int i4, int i5) {
  return qqgXgamgam_pmmpm(i2,i1,i3,i4,i5); }

Complex qqgXgamgam_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4, int i5, int h5) {
  if ( h4==h5 ) return 0;
  if ( h4==1 ) {
    if ( h1==1 ) {
      if ( h3==1 ) return qqgXgamgam_pmpmp(i1,i2,i3,i5,i4);
      else return qqgXgamgam_pmmpm(i1,i2,i3,i4,i5);
    }
    else { // h1==-1
      if ( h3==-1 ) return qqgXgamgam_mpmpm(i1,i2,i3,i4,i5);
      else return qqgXgamgam_mppmp(i1,i2,i3,i5,i4);
    }
  }
  else { // h4==-1
    if ( h1==1 ) {
      if ( h3==1 ) return qqgXgamgam_pmpmp(i1,i2,i3,i4,i5);
      else return qqgXgamgam_pmmpm(i1,i2,i3,i5,i4);
    }
    else { // h1==-1
      if ( h3==-1 ) return qqgXgamgam_mpmpm(i1,i2,i3,i5,i4);
      else return qqgXgamgam_mppmp(i1,i2,i3,i4,i5);
    }
  }
}

//======

Complex gggXgamgam_mppmp(int i1, int i2, int i3, int i4, int i5) {
  return -pow(spa(i1,i4),4)*pow(spb(i4,i5),2)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1);
}
Complex gggXgamgam_pmmpm(int i1, int i2, int i3, int i4, int i5) {
  return pow(spb(i1,i4),4)*pow(spa(i4,i5),2)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1);
}

Complex gggXgamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3, int i4, int h4, int i5, int h5) {
  if ( h4==h5 ) return 0;
  if ( h1==h2 && h2==h3 ) return 0;
  else if ( h4==1 ) {
    if ( h1==1 ) {
      if ( h2==1 ) return gggXgamgam_mppmp(i3,i1,i2,i5,i4);
      else if ( h3==1 ) return gggXgamgam_mppmp(i2,i3,i1,i5,i4);
      else return gggXgamgam_pmmpm(i1,i2,i3,i4,i5);
    }
    else { // h1==-1
      if ( h2==-1 ) return gggXgamgam_pmmpm(i3,i1,i2,i4,i5);
      else if ( h3==-1 ) return gggXgamgam_pmmpm(i2,i3,i1,i4,i5);
      else return gggXgamgam_mppmp(i1,i2,i3,i5,i4);
    }
  }
  else { // h4==-1
    if ( h1==1 ) {
      if ( h2==1 ) return gggXgamgam_mppmp(i3,i1,i2,i4,i5);
      else if ( h3==1 ) return gggXgamgam_mppmp(i2,i3,i1,i4,i5);
      else return gggXgamgam_pmmpm(i1,i2,i3,i5,i4);
    }
    else { // h1==-1
      if ( h2==-1 ) return gggXgamgam_pmmpm(i3,i1,i2,i5,i4);
      else if ( h3==-1 ) return gggXgamgam_pmmpm(i2,i3,i1,i5,i4);
      else return gggXgamgam_mppmp(i1,i2,i3,i4,i5);
    }
  }
}

//=================================================================
// Special versions with fewer labels:

// g(1) g(2) -> gam(3) gam(4):
Complex ggXgamgam(int h1, int h2, int h3, int h4){ return -ggXgamgam_gen(1,h1,2,h2,3,h3,4,h4); }

// q(1) qbar(2) -> gam(3) gam(4):
Complex qqbXgamgam(int h1, int h3, int h4){ return -qqXgamgam_gen(1,h1,2,3,h3,4,h4); }
// qbar(1) q(2) -> gam(3) gam(4):
Complex qbqXgamgam(int h2, int h3, int h4){ return -qqXgamgam_gen(2,h2,1,3,h3,4,h4); }

// g(1) g(2) -> gam(3) gam(4) g(5):
Complex ggXgamgamg(int h1, int h2, int h3, int h4, int h5){ return -gggXgamgam_gen(1,h1,2,h2,5,h5,3,h3,4,h4); }

// q(1) qbar(2) -> gam(3) gam(4) g(5)
Complex qqbXgamgamg(int h1, int h3, int h4, int h5) { return -qqgXgamgam_gen(1,h1,2,5,h5,3,h3,4,h4); }
// qbar(1) q(2) -> gam(3) gam(4) g(5)
Complex qbqXgamgamg(int h2, int h3, int h4, int h5) { return -qqgXgamgam_gen(2,h2,1,5,h5,3,h3,4,h4); }
// q(1) g(2) -> gam(3) gam(4) q(5)
Complex qgXgamgamq(int h1, int h2, int h3, int h4) { return -qqgXgamgam_gen(1,h1,5,2,h2,3,h3,4,h4); }
// g(1) q(2) -> gam(3) gam(4) q(5)
Complex gqXgamgamq(int h1, int h2, int h3, int h4) { return -qqgXgamgam_gen(2,h2,5,1,h1,3,h3,4,h4); }
