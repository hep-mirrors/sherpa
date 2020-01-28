namespace CFPSHOWER {
  class Functions {
  public:
    /////////////////////////////////////////////////////////////////////////////////////
    // Simple splitting kernels
    /////////////////////////////////////////////////////////////////////////////////////
    inline static double Pqq(const double & z) {
      return (1.+ATOOLS::sqr(z))/(1.-z);
    }
    inline static double Pgq(const double & z) {
      return ATOOLS::sqr(z)+ATOOLS::sqr(1.-z);
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // Real and (collinear) subtraction terms
    /////////////////////////////////////////////////////////////////////////////////////
    // The funtions below are based on
    // [CG] Catani & Grazzini, "Infrared factorization ...", 9908523,
    //      in the limit of epsilon->0
    // [HP] Hoeche & Prestel, "Triple Collinear Emissions in Parton Showers", 1705.00742
    /////////////////////////////////////////////////////////////////////////////////////
    //
    // Eqs (CG 57) & (CG 59) / Eq (HP, 34), 
    // excluding colour factors CF*TR for Rqqprime and for the interference term
    // we assume q(aij) -> qbar'(a) + q'(i) + q(j) and q(aij) -> qbar(a) + q(i) + q(j)
    inline static double
    RqqprimeF(const double & za,const double & zi,const double & zj,
	      const double & sai, const double & saj, const double & sij,
	      const double & saij) {
      double tai_j = (2.*(za*sij-zi*saj)+(za-zi)*sai)/(za+zi);
      return 1./2. * saij/sai * (-ATOOLS::sqr(tai_j)/(sai*saij) +
				 (4.*zj+ATOOLS::sqr(za-zi))/(za+zi) +
				 (za+zi-sai/saij) );
    }
    inline static double
    Softqqprime(const double & za,const double & zi,const double & zj,
		const double & sai, const double & saj, const double & sij,
		const double & saij) {
      double tai_j = (2.*(za*sij-zi*saj)+(za-zi)*sai)/(za+zi);
      return 1./2. * (-ATOOLS::sqr(tai_j/sai - (za-zi)/(za+zi)) -
		      4.*zj*(sai-saij)/((za+zi)*sai) );
    }
    inline static double
    RqqbarF(const double & za,const double & zi,const double & zj,
	    const double & sai, const double & saj, const double & sij,
	    const double & saij) {
      return ( RqqprimeF(za,zi,zj,sai,saj,sij,saij) + helpRqqbarF(za,zi,zj,sai,saj,sij,saij) +
	       RqqprimeF(za,zj,zi,saj,sai,sij,saij) + helpRqqbarF(za,zj,zi,saj,sai,sij,saij) );
      
    }
    inline static double
    helpRqqbarF(const double & za,const double & zi,const double & zj,
		const double & sai, const double & saj, const double & sij,
		const double & saij) {
      double Nc = 3.;
      return -1./Nc * saij/sai * (2.*sij/saij + (1.+ATOOLS::sqr(za))/(1.-zi) -
				  2.*zi/(1.-zj) -
				  saij/saj*za/2.*(1.+ATOOLS::sqr(za))/((1.-zi)*(1.-zj)));
    }
    // Eqs (61) & (62), Catani & Grazzini, "Infrared factorization ...", 9908523,
    //                  in the limit of epsilon->0 
    // excluding colour factors CF*CF for Rqqgg_Abelian and CF*CA for Rqqgg_nonAbelian
    // we assume q -> g(a) + g(i) + q(j)
    inline static double
    RqqggF_Abelian(const double & za,const double & zi,const double & zj,
		   const double & sai, const double & saj, const double & sij,
		   const double & saij) {
      return ( ATOOLS::sqr(saij)/(2.*saj*sij) * zj * (1.+ATOOLS::sqr(zj))/(za*zi) +
	       saij/saj * (zj*(1.-za)+pow(1.-zi,3))/(zi*zj) - sij/saj );
    }
    inline static double
    RqqggF_nonAbelian(const double & za,const double & zi,const double & zj,
		      const double & sai, const double & saj, const double & sij,
		      const double & saij) {
      double tai_j = (2.*(za*sij-zi*saj)+(za-zi)*sai)/(za+zi);
      double saij2 = ATOOLS::sqr(saij), za2 = ATOOLS::sqr(za), zj2 = ATOOLS::sqr(zi);
      double one_min_zi3 = pow(1.-zi,3);
      return ( 1./4. * (ATOOLS::sqr(tai_j/sij) + 1.) +
	       saij2/(saj*sij) * ( 1./2.*((1.+zj2)/zi + (1.+ATOOLS::sqr(1.-zi))/(1.-zj)) -
				   zj/4.* (1.+zj2)/(za*zi) ) +
	       saij/(2.*sai) * ( ( za*(1.+ATOOLS::sqr(1.-za))-
				   zi*(1.+ATOOLS::sqr(1.-zi)+4.*(1.-zi)) )/(zi*(1.-zj))) +
	       saij/(2.*saj) * ( (one_min_zi3+zj2-zi)/(zi*(1.-zj)) -
				 (one_min_zi3+zj*(1.-za))/(za*zi) ) );
    }
    // Eqs (68) & (69), Catani & Grazzini, "Infrared factorization ...", 9908523,
    //     in the limit of epsilon->0 
    // excluding colour factors CF*TR for Rggqq_Abelian and CA*TR for Rggqq_nonAbelian
    // we assume g -> qbar(a) + q(i) + g(j)
    inline static double
    RggqqF_Abelian(const double & za,const double & zi,const double & zj,
		   const double & sai, const double & saj, const double & sij,
		   const double & saij) {
      return ( -2 - sai*(saj+sij)/(saj*sij) +
	       2.*ATOOLS::sqr(saij)/(saj*sij)*(1.-zj+ATOOLS::sqr(zj)+2.*za*zi) -
	       saij/saj*(1.+2.*zi) - saij/sij*(1.-2.*za) );
    }
    inline static double
    RggqqF_nonAbelian(const double & za,const double & zi,const double & zj,
		      const double & sai, const double & saj, const double & sij,
		      const double & saij) {
      double tai_j = (2.*(za*sij-zi*saj)+(za-zi)*sai)/(za+zi);
      return (helpggqqF_nonAbelian(za,zi,zj,sai,saj,sij,saij,tai_j) +
	      helpggqqF_nonAbelian(zi,za,zj,sai,sij,saj,saij,tai_j));
    }
    inline static double
    helpggqqF_nonAbelian(const double & za,const double & zi,const double & zj,
			 const double & sai, const double & saj, const double & sij,
			 const double & saij,const double tai_j) {
      double saij2 = ATOOLS::sqr(saij);
      double zj3 = pow(zj,3), one_min_zj3 = pow(1.-zj,3), denom = 1./(zj*(1.-zj));
      return ( -1./4.*ATOOLS::sqr(tai_j/sai) +
	       saij2/(2.*saj*sij)*zi*(one_min_zj3-zj3-2.*zi*(1.-zi-2.*zj*za))*denom +
	       saij/(2.*sij)*(1.-za)*(1.-2*za*(1.-za))*denom +
	       saij/(2.*sai)*(1.+zj3+zj*ATOOLS::sqr(zi-za)-2.*za*zi*(1.+zj))*denom -
	       1./4. - saij2/(2.*saj*sij)*(1.-zj+ATOOLS::sqr(zj)+2.*za*zi) );
    }
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Subtractions from iteraed splitting functions
    /////////////////////////////////////////////////////////////////////////////////////
    // Eq (HP 35), excluding colour factors CF*TR for Sqqprime
    inline static double
    SqqprimeF(const double & za,const double & zi,const double & zj,
	      const double & sai, const double & saj, const double & sij,
	      const double & saij,const double & cosphi_aj) {
      return saij/sai * (Pqq(zj)*(ATOOLS::sqr(za)+ATOOLS::sqr(zi))/ATOOLS::sqr(za+zi) +
			 (4.*za*zi*zj)/pow(1.-zj,3)*(1.-2.*ATOOLS::sqr(cosphi_aj)) );
    }
    inline static double
    SqqbarF(const double & za,const double & zi,const double & zj,
	    const double & sai, const double & saj, const double & sij,
	    const double & saij,
	    const double & cosphi_aj,const double & cosphi_ai) {
      return ( SqqprimeF(za,zi,zj,sai,saj,sij,saij,cosphi_aj) +
	       SqqprimeF(za,zj,zi,saj,sai,sij,saij,cosphi_ai) ); 
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // End-point contributions
    /////////////////////////////////////////////////////////////////////////////////////
    // Eq (HP 49), (HP 51)
    inline static double
    DeltaIqqprimeF(const double & za,const double & zi,const double & zj) {
      return Itildeqqprime(za,zi,zj,za) - itildeqqprime(za,zi,zj,za+zi);
    }
    inline static double
    DeltaIqqbarF(const double & za,const double & zi,const double & zj) {
      return DeltaIqqprimeF(za,zi,zj)+DeltaIqqprimeF(za,zj,zi);
    }
    // Eq (HP 50), omitting colour factors CF*TR for Itilde and 2*CF for itilde
    inline static double
    Itildeqqprime(const double & za,const double & zi,
		  const double & zj,const double & x) {
      return ( Pqq(zj) +
	       (ATOOLS::sqr(za)+ATOOLS::sqr(zi))/ATOOLS::sqr(za+zi)*
	       (1-zj+Pqq(zj)*(log(x*zi*zj)-1.)) );
    }
    inline static double
    itildeqqprime(const double & za,const double & zi,
		  const double & zj,const double & x) {
      return 2. * ( Pqq(zj)*log(x*zj) + (1.-zj) ) * Pgq(za/(za+zi));
    }
  };
}
