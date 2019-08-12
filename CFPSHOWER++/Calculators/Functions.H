#ifndef CFPSHOWER_Calculators_Functions_H
#define CFPSHOWER_Calculators_Functions_H

#include "ATOOLS/Math/Vector.H"
#include <vector>
#include <algorithm>


namespace CFPSHOWER {
  class Functions {
  private:
    static double s_CF, s_CA, s_TR, s_Nc;
    inline static double P_qq(const double & z) {
      return (1.+ATOOLS::sqr(z))/(1.-z);
    }
    inline static double P_gq(const double & z) {
      return ATOOLS::sqr(z)+ATOOLS::sqr(1.-z);
    }
    inline static double Iqqp(const double & zi, const double & zj,
			      const double & zk, const double & x) {
      double P_qq_zk = P_qq(zk);
      return s_TR * ( P_qq_zk +
		      (1.-(2.*zi*zj)/ATOOLS::sqr(zi+zj)) *
		      (1.-zk+P_qq_zk) * (log(x*zj*zk)-1.) );
    }
    inline static double iqqp(const double & zi, const double & zj,
			      const double & zk, const double & x) {
      return 2 * s_TR * ( P_qq(zk)*log(x*zk) + (1.-zk) ) * P_gq(zi/(zi+zj));
      
    }
  public:
    inline static double DeltaI_qqp_F(const double & zi, const double & zj,
				      const double & zk) {
      return Iqqp(zi,zj,zk,zi) - iqqp(zi,zj,zk,zi+zj);
    }
    inline static double DeltaI_qqbar_F(const double & zi, const double & zj,
					const double & zk) {
      return DeltaI_qqp_F(zi,zj,zk) + DeltaI_qqp_F(zi,zk,zj);
    }
  };
}

#endif
