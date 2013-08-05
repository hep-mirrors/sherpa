/* =============================================================
   The 1-loop and 2-loop continuum gg -> gamma gamma amplitudes.
   with light quark masses and top quark contribution neglected
   ============================================================= */

#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;

#include "gggamgam1l.C"  // gggamgam 1-loop
                         // amplitudes, with and without phases.
                         // as another contributing process, at NNLO
#include "gggamgam2l.C"  // gggamgam 2-loop amplitudes (finite parts)
                         // with phases removed.  Also reads in 1-loop
                         // amplitudes, with and without phases.
                         // as another contributing process, at NNLO

/* =============================================================
   The 1-loop continuum gg -> gamma gamma amplitudes.
   The light quark (u,d,c,s,b) masses included
   ============================================================= */

// NEW helicity LABELLING:
// hel = 1 means ++ -> ++ (or equivalently -- -> --) or "mmpp"
// hel = 2 means ++ -> -- (or equivalently -- -> ++) or "pppp" 
// Other helicity configurations are irrelevant for Higgs interference

// OLD helicity LABELLING:   pppp = 1, mppp = 2, mmpp = 3, mpmp = 4, mppm = 5
// Unless specified defaulting to new helicity labelling

/* The one-loop light-by-light amplitudes */

// Basic scalar integrals, with a uniform mass m in the loop:

Complex B0(double m, double s){
 Complex temp; 
 double x;
 if (Abs(s/4./m/m) < 0.01)  // use large mass expansion here:
     { return s/6./m/m * ( 1. + 1./10. * s/m/m * (
			   1. + 1./7. * s/m/m * ( 
			   1. + 1./6. * s/m/m ))) ; }
 if (Abs(s/4./m/m) > 300.0) // use small mass expansion here:
   { if (s < 0) { return 2. - log(-s/m/m)
		    + 2.*m*m/s * ( 1. + log(-s/m/m) )
                    + m*m*m*m/s/s * ( -1. + 2. * log(-s/m/m) ) ; }
     else { return 2. - log(s/m/m)
		    + 2.*m*m/s * ( 1. + log(s/m/m) )
	            + m*m*m*m/s/s * ( -1. + 2. * log(s/m/m) ) 
	    + I * PI_DEF * ( 1. - 2.*m*m/s - 2.*m*m*m*m/s/s ) ; }
   }
 if (s < 0) 
    { x = sqrt(1.-4.*m*m/s); 
      temp = 2. + x * log((x-1.)/(x+1.)) ; }
 if (s >= 0 && s < 4.*m*m) 
      { x = sqrt(-1.+4.*m*m/s) ;
        temp = real( 2. + I * x * log((I*x-1.)/(I*x+1.)) ) ; }
 if (s >= 4.*m*m)
    { x = sqrt(1.-4.*m*m/s); 
      temp = 2. + x * ( log((1.-x)/(x+1.)) + I * PI_DEF ) ; }
 return temp; 
}

Complex C0(double m, double s){
 Complex temp, temp1; 
 double x;
 if (Abs(s/4./m/m) < 0.01)    // use large mass expansion here:
     { return -0.5/m/m * ( 1. + 1./12. * s/m/m * (
			   1. + 2./15. * s/m/m * ( 
			   1. + 9./56. * s/m/m ))) ; }
 if (Abs(s/4./m/m) > 300.0) // use small mass expansion here:
   { if (s < 0) { temp1 = log(-s/m/m) - 2.*m*m/s - 3.*m*m*m*m/s/s 
			       - 20./3.*m*m*m*m*m*m/s/s/s ; 
                  return 0.5/s * temp1 * temp1 ; }
     else { temp1 = log(s/m/m) - 2.*m*m/s - 3.*m*m*m*m/s/s 
			       - 20./3.*m*m*m*m*m*m/s/s/s - I * PI_DEF ;
            return 0.5/s * temp1 * temp1 ; }
   }
 if (s < 0) 
    { x = sqrt(1.-4.*m*m/s); 
      temp1 = log((x+1.)/(x-1.));
      temp = 0.5/s * temp1 * temp1 ; }
 if (s > 0 && s < 4.*m*m) 
      { x = sqrt(-1.+4.*m*m/s) ;
        temp1 = log((I*x-1.)/(I*x+1.)) ;
        temp = 0.5/s * temp1 * temp1 ; }
 if (s >= 4.*m*m)
    { x = sqrt(1-4.*m*m/s); 
      temp1 = log((1.-x)/(x+1.)) + I * PI_DEF ;
      temp = 0.5/s * temp1 * temp1 ; }
 return temp; 
}

Complex D0(double m, double s, double t){
 Complex bu, temp; 
 double u, v, a, bv, temp_re, temp_im;
 // Complex u, v, a, bv, temp_re;
 if (t > s) { return D0(m,t,s); }
 u = s/4./m/m;   
 v = t/4./m/m;
 if ((Abs(u) < 0.01) && (Abs(v) < 0.01))   // use large mass expansion here:
   { return 1./m/m/m/m * ( 1./6. + 1./15.*(u+v) 
			   + 2./105. * (2.*(u*u+v*v)+u*v)
                           + 8./945. * (u+v)*(3*(u*u+v*v)-2*u*v) ); }
 if ((Abs(u) > 10000.0) && (Abs(v) > 10000.0)) // small mass expansion here:
 { if (u < 0.) { return 
         ( - (log(s/t)*log(s/t) + PI_DEF*PI_DEF) 
           + 2.*s * C0(m,s) + 2.*t * C0(m,t) )/ (s*t-4.*m*m*(s+t)) ; }
   else return ( - ( log(-s/t)*log(-s/t) - 2. * I * PI_DEF * log(-s/t) ) 
		 + 2.*s * C0(m,s) + 2.*t * C0(m,t) )/ (s*t-4.*m*m*(s+t)) ; 
 }
 a = sqrt(1.-(u+v)/u/v);
 if ((u<0.) || (u>1.)) { bu = sqrt((u-1.)/u); }
 else { bu = I * sqrt((1.-u)/u); }
 bv = sqrt((v-1.)/v);
 temp_re = 2./s/t/a * ( 
              - ReLi2((a+1.)/(a+bu)) - ReLi2((a+1.)/(a-bu)) 
              + ReLi2((a-1.)/(a+bu)) + ReLi2((a-1.)/(a-bu)) 
              - ReLi2((a+1.)/(a+bv)) - ReLi2((a+1.)/(a-bv)) 
              + ReLi2((a-1.)/(a+bv)) + ReLi2((a-1.)/(a-bv)) ); 
 if (s > 4.*m*m) { 
   temp_im = real( -2./s/t * PI_DEF/a * log(v*(a+bu)*(a+bu)) );
   temp = temp_re + I * temp_im ; 
   return temp; }
 else return temp_re;
}

/*=====================================================================
   One-loop light-by-light helicity amplitudes, 
   from Gounaris, Porfyriadis, Renard, Eur. Phys. J C9, 673 (1999).   
   We use our all-outgoing helicity convention, though, and note that
   t and u should be exchanged vs. GPR conventions.
   Here x = cos(theta).  t = -s/2*(1-x), u = -s/2*(1+x). 
======================================================================*/

// W loop: (old hel. label used)

Complex lbyl_W_loop(int hel, double x, double m_W, double s){
 double t, u, m_W_2, m_W_4;
 Complex temp;
 t = -s/2.*(1.-x);    
 u = -s/2.*(1.+x);
 m_W_2 = m_W*m_W;
 m_W_4 = m_W_2*m_W_2;
 if (hel==1) { temp = // pppp
       - 12. + 24. * m_W_4 * ( D0(m_W,s,t) + D0(m_W,s,u) + D0(m_W,t,u) ) ; }
 else if (hel==2) { temp = 
       - 12. + 24. * m_W_4 * ( D0(m_W,s,t) + D0(m_W,s,u) + D0(m_W,t,u) )
           + 12. * m_W_2 * s * t * u 
                * ( D0(m_W,s,t)/u/u + D0(m_W,s,u)/t/t + D0(m_W,t,u)/s/s )
           - 24. * m_W_2 * ( 1./s + 1./t + 1./u )
                      * ( t * C0(m_W,t) + u * C0(m_W,u) + s * C0(m_W,s) ) ; } 
 else if (hel==3) { temp = // mmpp
         12. - 12. * (1. + 2.*u/s) * B0(m_W,u) 
             - 12. * (1. + 2.*t/s) * B0(m_W,t)
         + 24. * m_W_2 * t * u/s * D0(m_W,u,t)
         + 16. * (1. - 1.5 * m_W_2/s - 0.75 * t * u/s/s )
              * ( 2. * t * C0(m_W,t) + 2. * u * C0(m_W,u) 
                - t * u * D0(m_W,t,u) )
         + 8. * (s-m_W_2) * (s-3.*m_W_2)
                      * (D0(m_W,s,t) + D0(m_W,s,u) + D0(m_W,t,u)) ; }
 else if (hel==4)      // cross s <-> u from mmpp
   { temp = 
          12. - 12. * (1. + 2.*s/u) * B0(m_W,s) 
              - 12. * (1. + 2.*t/u) * B0(m_W,t)
         + 24. * m_W_2 * t * s/u * D0(m_W,s,t)
         + 16. * (1. - 1.5 * m_W_2/u - 0.75 * t * s/u/u )
              * ( 2. * t * C0(m_W,t) + 2. * s * C0(m_W,s) 
                - t * s * D0(m_W,t,s) )
         + 8. * (u-m_W_2) * (u-3.*m_W_2)
                      * (D0(m_W,u,t) + D0(m_W,u,s) + D0(m_W,t,s)) ; }
 else if (hel==5)      // exchange t <-> u, i.e. x <-> -x, in mpmp=4.
   { return lbyl_W_loop(4,-x,m_W,s) ; }
 return temp;
}

// fermion loop: (old hel. label used)

Complex lbyl_f_loop(int hel, double x, double m_f, double s){
 double t, u, m_f_2, m_f_4;
 Complex temp;
 t = -s/2.*(1.-x);    
 u = -s/2.*(1.+x);
 m_f_2 = m_f*m_f;
 m_f_4 = m_f_2*m_f_2;
 if (hel==1) { return - 2./3. * lbyl_W_loop(1,x,m_f,s) ; } // pppp
 else if (hel==2) { return - 2./3. * lbyl_W_loop(2,x,m_f,s) ; }
 else if (hel==3) { temp = // mmpp
           - 8. + 8. * (1. + 2.*u/s) * B0(m_f,u)  
                + 8. * (1. + 2.*t/s) * B0(m_f,t)
           - 8. * ( (t*t+u*u)/s/s - 4.*m_f_2/s ) 
               * ( t * C0(m_f,t) + u * C0(m_f,u) )
           + 8. * m_f_2 * (s-2.*m_f_2) * (D0(m_f,s,t) + D0(m_f,s,u))
           - 4. * ( 4. * m_f_4 - (2.*s*m_f_2 + t*u) * (t*t+u*u)/s/s 
                  + 4. * m_f_2 * t*u/s ) * D0(m_f,t,u) ; }
 else if (hel==4)   // cross s <-> u from mmpp
   { temp = - 8. + 8. * (1. + 2.*s/u) * B0(m_f,s) 
                + 8. * (1. + 2.*t/u) * B0(m_f,t)
           - 8. * ( (t*t+s*s)/u/u - 4.*m_f_2/u ) 
               * ( t * C0(m_f,t) + s * C0(m_f,s) )
           + 8. * m_f_2 * (u-2.*m_f_2) * (D0(m_f,u,t) + D0(m_f,u,s))
           - 4. * ( 4. * m_f_4 - (2.*u*m_f_2 + t*s) * (t*t+s*s)/u/u 
                    + 4. * m_f_2 * t*s/u ) * D0(m_f,t,s) ; }
 else if (hel==5)      // exchange t <-> u, i.e. x <-> -x, in mpmp=4.
   { return lbyl_f_loop(4,-x,m_f,s) ; }
 return temp;
}

// scalar loop: (old hel. label used)

Complex lbyl_s_loop(int hel, double x, double m_sc, double s){
 double t, u, m_s_2, m_s_4;
 Complex temp;
 t = -s/2.*(1.-x);    
 u = -s/2.*(1.+x);
 m_s_2 = m_sc*m_sc;
 m_s_4 = m_s_2*m_s_2;
 if (hel==1) { return 1./3. * lbyl_W_loop(1,x,m_sc,s) ; } // pppp
 else if (hel==2) { return 1./3. * lbyl_W_loop(2,x,m_sc,s) ; }
 else if (hel==3) { temp = // mmpp
             4. - 4. * (1. + 2.*u/s) * B0(m_sc,u) 
                - 4. * (1. + 2.*t/s) * B0(m_sc,t)
           + 8. * m_s_2 * t*u/s * D0(m_sc,t,u)
           - 8. * m_s_2/s * ( 1. + u*t/2./m_s_2/s )
               * ( 2. * t * C0(m_sc,t) + 2. * u * C0(m_sc,u) 
                 - t * u * D0(m_sc,t,u) )
          + 8. * m_s_4 * (D0(m_sc,s,t) + D0(m_sc,s,u) + D0(m_sc,t,u)) ; }
 else if (hel==4)    // cross s <-> u from mmpp
   { temp = 4. - 4. * (1. + 2.*s/u) * B0(m_sc,s) 
               - 4. * (1. + 2.*t/u) * B0(m_sc,t)
           + 8. * m_s_2 * t*s/u * D0(m_sc,t,s)
           - 8. * m_s_2/u * ( 1. + s*t/2./m_s_2/u )
               * ( 2. * t * C0(m_sc,t) + 2. * s * C0(m_sc,s) 
                 - t * s * D0(m_sc,t,s) )
         + 8. * m_s_4 * (D0(m_sc,u,t) + D0(m_sc,u,s) + D0(m_sc,t,s)) ; }
 else if (hel==5)      // exchange t <-> u, i.e. x <-> -x, in mpmp=4.
   { return lbyl_s_loop(4,-x,m_sc,s) ; }
 return temp;
}

/*==================================================
 Add up all Standard Model one-loop contributions to 
 gg -> gamma gamma amplitude.  s is in GeV^2
 (Formula is not reliable at hadronic thresholds, of course.) */

// old hel. label used

Complex A_cont_1l(int hel, double x, double s, double muR){
 // for the heaviest quark we'll use running masses (because we have them);
 // for u,d,s it doesn't matter, so we'll just leave them as m_q(m_q).
 double mt = M_t(muR);  double mb = M_b(muR); double mc = M_c(muR);
 return 1./2. * alpha0 * alpha_s(muR) 
     * ( 4./9. * ( lbyl_f_loop(hel,x,m_u,s) + lbyl_f_loop(hel,x,mc,s)
                 + lbyl_f_loop(hel,x,mt,s) )
       + 1./9. * ( lbyl_f_loop(hel,x,m_d,s) + lbyl_f_loop(hel,x,m_s,s)
	        + lbyl_f_loop(hel,x,mb,s) ) ) ;
}

/* ================================================================
   The continuum g g -> gamma gamma amplitude with only
   two interference relevent helicity configurations */

// new hel. label used

Complex A_c_1l(double M_H, double cth, int hel, double muR){
  int oldhel;
  if (hel==2) oldhel=1; else if (hel==1) oldhel=3;
  return A_cont_1l(oldhel, cth, M_H*M_H, muR);
}
