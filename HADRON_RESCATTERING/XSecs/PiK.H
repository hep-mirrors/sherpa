#ifndef HADRON_RESCATTERING_XSecs_PiK_H
#define HADRON_RESCATTERING_XSecs_PiK_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"

//////////////////////////////////////////////////////////////////////////////
//
// Parametrizations taken from hep-ph/0411334, Appendix A, and
// ignoring the relatively small F and G-waves (with l>3)
// (this was the simplest case that looks more or less reasonable).
//
// Regge parametrizations taken from hep-ph/0312187, Eqs (1-4, 11, 18).
//
// Todo: We should update this - at least in parts - with more
//       recent parametrizations from 1102.2183 or 1907.13162.
//
//////////////////////////////////////////////////////////////////////////////

namespace HADRON_RESCATTERING {  
  class PiK {
  public:
  private:
    double  m_mpi, m_mpi2, m_mK, m_mK2, m_meta, m_meta2;
    double  m_sum2, m_Delta2, m_sKeta;
    double  m_beta;
    double  m_gKpi;
    double  m_fpi_P, m_fpi_Pp, m_fpi_rho, m_fpi_2, m_s0, m_b;
    Complex m_xi_Pom, m_xi_PomP, m_xi_rho, m_xi_2;
    
    inline const double sigma(const double & s) const {
      return sqrt(ATOOLS::sqr(s-m_mK2-m_mpi2)-4.*m_mK2*m_mpi2)/s;
    }
    inline const double k2(const double & s) const {
      return (ATOOLS::sqr(s-m_mK2-m_mpi2)-4.*m_mK2*m_mpi2)/(4.*s);
    }
    inline const double keta2(const double & s) const {
      return (ATOOLS::sqr(s-m_mK2-m_meta2)-4.*m_mK2*m_meta2)/(4.*s);
    }



    inline const double omega(const double & s,const double & y0,
			      const double & alpha=1.) const {
      double sqrt_y = (s-m_Delta2)/(s+m_Delta2); 
      double y      = ATOOLS::sqr(sqrt_y), sqrt_y0y = sqrt(y0-y);
      return (sqrt_y-alpha*sqrt_y0y)/(sqrt_y+alpha*sqrt_y0y);
    }
    inline const Complex Sn(const double & s,const double & srn,
			    const double & Pn, const double & Qn) const {
      // Eq (10) from 1602.08404
      double sterm = srn-s, PQsum = Pn+Qn, PQdiff = Pn-Qn;
      return ( Complex(ATOOLS::sqr(sterm)-(ATOOLS::sqr(Pn)-ATOOLS::sqr(Qn)),
		       sterm*(PQsum+PQdiff))/
	       (ATOOLS::sqr(sterm)+ATOOLS::sqr(PQsum)) );
    }
    inline const double p1(const double & x) const {
      return 1. - 5.22*ATOOLS::sqr(x) + 7.57*pow(x,4);
    }
    inline const double p2(const double & x) const {
      return 1. - 1.72*ATOOLS::sqr(x);
    }
    
    const Complex ampl(const double & s,const int & I, const int & l) const;
    const Complex ampl_10(const double & s) const;
    inline const double xsel_partial(const int & I,const int & l,
				     const double & s) const {
      return 0.;
    }
    inline const double xsinel_partial(const int & I,const int & l,
				       const double & s) const {
      return 0.;
    }
    
    const double eta(const size_t & I,const size_t & l,
		     const double & s) const;
    const double delta(const size_t & I,const size_t & l,
		       const double & s) const;

    
    
    const Complex sigma_I0(const double & s,const double & t=0.) const; 
    const Complex sigma_I1(const double & s,const double & t=0.) const; 
    

    void Test();
  public:
    PiK();
    ~PiK() {}
    const double XStot(const double & s,const int & I1,const int & I2) const;
    const double hatXSel(const double & s,const int & I1,const int & I2) const;
    const double hatXSinel(const double & s,const int & I1,const int & I2) const;
    const double hatXStot(const double & s,const int & I1,const int & I2) const;
  };
}

#endif
