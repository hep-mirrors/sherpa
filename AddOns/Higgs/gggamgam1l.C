// Amplitudes for leading-order (1-loop) contribution to 
// g g -> gamma gamma.

// hel:
// pppp = 1, mppp = 2, pmpp = 3, ppmp = 4,
// mmpp = 5, mpmp = 6, pmmp = 7, mmmp = 8 

// gggamgam1l represents the contribution of one quark flavor.
// The numerical prefactor omitted from gggamgam1l is
// (note "Catani" normalization conventions):
// 2 * (sqrt(2)*e*Q)^2 * g^2/(4*Pi)^2 = 4 * Q^2 * alpha * alpha_s 

Complex gggamgam1l_pppp(){ return 1.; }
Complex gggamgam1l_mppp(){ return 1.; }
Complex gggamgam1l_ppmp(){ return 1.; }

Complex gggamgam1l_mmpp(){ 
  double s = sij(1,2);   double t = sij(2,3);  double u = sij(1,3);
  return 
    - 0.5 * (t*t+u*u)/s/s * (log(t/u)*log(t/u) + PISQ) 
	     - (t-u)/s * log(t/u) - 1. ;
}

Complex gggamgam1l_mpmp(){ 
  double s = sij(1,2);   double t = sij(2,3);  double u = sij(1,3);
  return 
    - 0.5 * (t*t+s*s)/u/u * log(-t/s) * ( log(-t/s) + 2. * PI_DEF * I )
	     - (t-s)/u * ( log(-t/s) + PI_DEF * I ) - 1. ;
}

Complex gggamgam1l_mppm(){ 
  double s = sij(1,2);   double u = sij(2,3);  double t = sij(1,3);
  return 
    - 0.5 * (t*t+s*s)/u/u * log(-t/s) * ( log(-t/s) + 2. * PI_DEF * I )
	     - (t-s)/u * ( log(-t/s) + PI_DEF * I ) - 1. ;
}

/*=====================================================================
   Get g(1,h1) g(2,h2) -> gamma(3,h3) gamma(4,h4) 
   1-loop amplitudes for all helicities, with phases removed,
   using symmetries
======================================================================*/

Complex gggamgam1l(int h1, int h2, int h3, int h4){
  int hsum = h1+h2+h3+h4;
//  more (+)'s than (-)'s
  if (hsum == 4){ return gggamgam1l_pppp(); }
  else if (hsum == 2){
    if ((h1 == -1) || (h2 == -1)){ return gggamgam1l_mppp(); }
    else { return gggamgam1l_ppmp(); } }
  else if (hsum == 0){
    if (((h1 == -1) && (h2 == -1)) || ((h3 == -1) && (h4 == -1))){
      return gggamgam1l_mmpp(); } 
    else if (((h2 == -1) && (h3 == -1)) || ((h4 == -1) && (h1 == -1))){
      return gggamgam1l_mppm(); }
    else { return gggamgam1l_mpmp(); } }
//  more (-)'s than (+)'s
  else if (hsum == -2){
    if ((h1 == 1) || (h2 == 1)){ return gggamgam1l_mppp(); }
    else { return gggamgam1l_ppmp(); } }
  else { return gggamgam1l_pppp(); }
}
