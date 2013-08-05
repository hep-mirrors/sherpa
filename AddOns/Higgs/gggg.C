//=== FOUR-POINT AMPLITUDES ===

Complex gamgamgamgam_pppp(int i1,int i2,int i3,int i4) {
  return spb(i3,i4)*spb(i3,i4)/spa(i1,i2)/spa(i1,i2) ; }

Complex gamgamgamgam_mppp(int i1,int i2,int i3,int i4) {
  return spa(i1,i3)*spb(i2,i3)*spb(i3,i4)
        /spb(i1,i3)/spa(i2,i3)/spa(i3,i4) ; }

Complex gamgamgamgam_mmpp(int i1,int i2,int i3,int i4) {
  return spa(i1,i2)*spa(i1,i2)*spa(i1,i2)/spa(i2,i3)/spa(i3,i4)/spa(i4,i1) * (
 0.5 * sij(i2,i3)/pow(sij(i1,i2),3) * ( pow(sij(i1,i3),2) + pow(sij(i2,i3),2) )
     * ( pow(Clog(sij(i1,i3),sij(i2,i3)),2) + PISQ )
+ sij(i2,i3)/sij(i1,i2)/sij(i1,i2) * (sij(i1,i3)-sij(i2,i3)) 
     * ( Clog(sij(i1,i3),sij(i2,i3)) )
+ sij(i2,i3)/sij(i1,i2) ) ; }

// mpmp is the same as mmpp but with 2 <-> 3
inline Complex gamgamgamgam_mpmp(int i1,int i2,int i3,int i4) { return gamgamgamgam_mmpp(i1,i3,i2,i4) ; }

Complex gamgamgamgam_mmmm(int i1,int i2,int i3,int i4) {
  return spa(i3,i4)*spa(i3,i4)/spb(i1,i2)/spb(i1,i2) ; }

Complex gamgamgamgam_pmmm(int i1,int i2,int i3,int i4) {
  return spb(i1,i3)*spa(i2,i3)*spa(i3,i4)
        /spa(i1,i3)/spb(i2,i3)/spb(i3,i4) ; }


Complex gamgamgamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3, int i4, int h4){
  int hsum = h1+h2+h3+h4;
//  more (+)'s than (-)'s
  if (hsum == 4){ return gamgamgamgam_pppp(i1,i2,i3,i4); }
  else if (hsum == 2){
    if (h1 == -1){ return gamgamgamgam_mppp(i1,i2,i3,i4); }
    else if (h2 == -1){ return gamgamgamgam_mppp(i2,i3,i4,i1); }
    else if (h3 == -1){ return gamgamgamgam_mppp(i3,i4,i1,i2); }
    else { return gamgamgamgam_mppp(i4,i1,i2,i3); } }
  else if (hsum == 0){ 
    if ((h1 == -1) && (h2 == -1)){ return gamgamgamgam_mmpp(i1,i2,i3,i4); }
    else if ((h2 == -1) && (h3 == -1)){ return gamgamgamgam_mmpp(i2,i3,i4,i1); }
    else if ((h3 == -1) && (h4 == -1)){ return gamgamgamgam_mmpp(i3,i4,i1,i2); }
    else if ((h4 == -1) && (h1 == -1)){ return gamgamgamgam_mmpp(i4,i1,i2,i3); }
    else if ((h1 == -1) && (h3 == -1)){ return gamgamgamgam_mpmp(i1,i2,i3,i4); }
    else { return gamgamgamgam_mpmp(i4,i1,i2,i3); } }
//  more (-)'s than (+)'s
  else if (hsum == -2){
    if (h1 == 1){ return gamgamgamgam_pmmm(i1,i2,i3,i4); }
    else if (h2 == 1){ return gamgamgamgam_pmmm(i2,i3,i4,i1); }
    else if (h3 == 1){ return gamgamgamgam_pmmm(i3,i4,i1,i2); }
    else { return gamgamgamgam_pmmm(i4,i1,i2,i3); } }
  else { return gamgamgamgam_mmmm(i1,i2,i3,i4); }
}

/* gggamgam is the same as gamgamgamgam,
   the light by light primitive scattering
   amplitudes, i.e. up to coupling factors. */

inline Complex gggamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3,
                     int i4, int h4) {
 return gamgamgamgam_gen(i1,h1,i2,h2,i3,h3,i4,h4); }

//================================================================

// Version without i1, ..., i4 labels, for g(1) g(2) -> gamma(3) gamma(4),
// since there is only one independent color structure, delta^{a_1 a_2}.

Complex gggamgam(int h1, int h2, int h3, int h4){
 return gggamgam_gen(1,h1,2,h2,3,h3,4,h4);
}
