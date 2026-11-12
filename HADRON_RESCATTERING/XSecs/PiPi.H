#ifndef HADRON_RESCATTERING_XSecs_PiPi_H
#define HADRON_RESCATTERING_XSecs_PiPi_H

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
  class PiPi {
  public:
  private:
    double  m_mpi, m_mpi2, m_mK, m_mK2;
    
    inline const double sigma(const double & s) const {
      return sqrt(1.-4.*m_mpi2/s);
    }
    inline const double k(const double & s) const {
      return sqrt(s/4.-m_mpi2);
    }
    inline const Complex t(const size_t & I,const size_t & l,
			   const double & s) const {
      return ( sqrt(s)/(M_PI*k(s)) *
	       ( eta(I,l,s) * exp(Complex(2.,0.)*delta(I,l,s)) -
		 Complex(1.,0.))/Complex(0.,1.)
	       );
    }
    inline const double xsel_partial(const int & I,const int & l,
				     const double & s) const {
      double etaIl = eta(I,l,s), deltaIl = delta(I,l,s);
      return (2.*l+1)*((1.+etaIl*etaIl)/2.-etaIl*cos(2.*deltaIl) )/2.;
    }
    inline const double xsinel_partial(const int & I,const int & l,
				       const double & s) const {
      double etaIl = eta(I,l,s);
      return (2.*l+1)*(1.-etaIl*etaIl)/4.;
    }
    
    const double eta(const size_t & I,const size_t & l,
		     const double & s) const;
    const double delta(const size_t & I,const size_t & l,
		       const double & s) const;

    
    const Complex sigma_I0(const double & s,const double & t=0.) const; 
    const Complex sigma_I1(const double & s,const double & t=0.) const; 
    const Complex sigma_I2(const double & s,const double & t=0.) const; 
    

    void Test();
  public:
    PiPi();
    ~PiPi() {}
    const double XStot(const double & s,const int & I1,const int & I2) const;
    const double hatXSel(const double & s,const int & I1,const int & I2) const;
    const double hatXSinel(const double & s,const int & I1,const int & I2) const;
    const double hatXStot(const double & s,const int & I1,const int & I2) const;
  };
}

#endif
