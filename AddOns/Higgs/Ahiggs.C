
/*=================================================================
  LO production cross section, sigma_0(g g -> H) (see CdFG, hep-ph/0102227).
  We also use some of these functions in Gamma(H -> g g) below.  */

Complex f_AQ(double x) {
  if (x >= 1.) { return asin(1./sqrt(x)) * asin(1./sqrt(x)); }
  else { return - 0.25 * ( log((1.+sqrt(1.-x))/(1.-sqrt(1.-x))) - I*PI )
                       * ( log((1.+sqrt(1.-x))/(1.-sqrt(1.-x))) - I*PI ) ; }
}

Complex A_Q(double M_Q, double M_H){
  double x = 4.*M_Q*M_Q/M_H/M_H;
  return 1.5 * x * ( 1. + (1.-x) * f_AQ(x) );
}

// quark loop contribution to H -> gamma gamma:

inline Complex A_QH(double M_Q, double M_H){ return 4./3. * A_Q(M_Q,M_H); }

// W loop contribution to H -> gamma gamma:

Complex A_WH(double M_H){
 double x = 4.*m_W*m_W/M_H/M_H;
 return - x * ( 3. + 2./x + 3. * (2.-x) * f_AQ(x) );
}

/* ================================================================
   The one-loop gg -> H production amplitude  */

Complex A_prod_1l(double M_H, double mu_R){
 // A_Q(M_q=infinity,M_H) = 1
 // A_Q(M_q=0,M_H) = 0 
  // return alpha_s(mu_R)/3. * M_H*M_H/PI * sqrt(G_F/2./SQRT2);// EFT
 return alpha_s(mu_R)/3. * M_H*M_H/PI * sqrt(G_F/2./SQRT2)
   * ( A_Q(M_t(mu_R),M_H) + A_Q(M_b(mu_R),M_H) + A_Q(M_c(mu_R),M_H) ) ;
}

/* The two-loop QCD correction to A_prod for a heavy top loop
  (ignoring the singular part which cancels against A_c_2l in the
  perturbative expansion of Rinterf). */

Complex A_prod_2l(double M_H, double mu_R){
 return alpha_s(mu_R)/3. * M_H*M_H/PI * sqrt(G_F/2./SQRT2) 
   * ( A_Q(M_t(mu_R),M_H) + A_Q(M_b(mu_R),M_H) + A_Q(M_c(mu_R),M_H) )
   * 11./6. * C_A * alpha_s(mu_R)/2./PI ;
}

/* ================================================================
   The one-loop H -> gamma gamma decay amplitude */

Complex A_dec_1l(double M_H, double mu_R){
 // A_QH(...) == 4/3 A_Q(...)
 // A_WH(m_W=infinity,M_H) = -7
 // A_WH(m_W=0,M_H) = 0
  // return alpha0 * M_H*M_H/PI * sqrt(G_F/2./SQRT2) * 47./18.;// EFT
 return alpha0/2. * M_H*M_H/PI * sqrt(G_F/2./SQRT2) 
             * ( 3. * e_t*e_t * A_QH(M_t(mu_R),M_H) 
               + 3. * e_b*e_b * A_QH(M_b(mu_R),M_H) 
               + 3. * e_c*e_c * A_QH(M_c(mu_R),M_H)
		                 + A_QH(m_tau,M_H)
		                       + A_WH(M_H) ) ;
}

/* The two-loop QCD correction to A_dec for a heavy top loop,
   and for a light bottom loop. Light bottom loop is from Akhoury, Wang
   and Yakovlev, hep-ph/0102105, eq. (41), evaluated for mu = M_H.
   No resummation of ln(M_H/m_b) yet, though... */ 

Complex A_dec_2l(double M_H, double mu_R){
  double lnb = 2. * log(M_H/M_b(mu_R));
  return alpha0/2. * M_H*M_H/PI * sqrt(G_F/2./SQRT2) 
             * ( 3. * e_t*e_t * A_QH(M_t(mu_R),M_H) 
                    * (-1.) * alpha_s(mu_R)/PI
               + 3. * e_b*e_b * A_QH(M_b(mu_R),M_H) 
		    * C_F * alpha_s(mu_R)/2./PI 
                          * ( -1./12. * lnb*lnb + 2. * lnb ) );
}
