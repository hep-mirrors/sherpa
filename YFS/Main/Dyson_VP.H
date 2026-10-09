#ifndef YFS_Main_Dyson_VP_H
#define YFS_Main_Dyson_VP_H

/*!
  \file Dyson_VP.H

  YFS: DYSON_VP - the Dyson resummation of the photon vacuum polarisation
  beyond the one loop of the provider

*/

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

#include <string>
#include <vector>

namespace YFS {

  //! One massless helicity configuration of the 2 -> 2 Born.
  struct VP_Helicity {
    int l[4];        //!< physical helicities of e-, e+, f, fbar (+-1)
    double w;        //!< |A(h)|^2 / sum_h |A(h)|^2
    Complex r;       //!< sum_x D(x) A_gamma,x(h) / A(h)
  };

  //! The event's dressing.
  struct VP_Factors {
    bool ok = false;
    double s = 0., t = 0.;
    Complex Ds = Complex(0., 0.), Dt = Complex(0., 0.);
    std::vector<VP_Helicity> hel;    //!< the configurations with A(h) != 0
    Complex rbar = Complex(0., 0.);  //!< sum_h w r, for mass-suppressed labels
    //! squared level, relative to the Born: sum_h w (|1 + r|^2 - 1)
    double v2 = 0.;
    //! r for the physical helicities l1..l4, rbar if A(h) = 0 there
    Complex Share(int l1, int l2, int l3, int l4) const;
  };

  class Dyson_VP {
  public:
    //! complete: the provider has no photon self-energy, so the whole
    //! (1 - Delta r)/(1 - Pi) - 1 is added rather than its part beyond O(alpha)
    explicit Dyson_VP(bool complete = false);
    //! complex Delta alpha(x) = Pi(x), on shell, alpha(0); x the signed invariant
    Complex DeltaAlpha(double x) const;
    //! D(x) = (1 - Delta r)/(1 - Pi(x)) - 1 - (Pi(x) - Delta r)
    Complex DysonFactor(double x) const;
    //! the helicity factors at the 2 -> 2 Born point p (e-, e+, f, fbar)
    VP_Factors Factors(const ATOOLS::Flavour_Vector &fl,
                       const ATOOLS::Vec4D_Vector &p) const;
    inline double DeltaR() const { return m_deltar; }
    std::string Describe() const;
  private:
    double m_alpha0, m_deltar, m_mz, m_gz, m_alphas;
    bool m_complete;
    Complex m_sw2;
    //! Delta alpha of one fermion loop, exact one-loop in its mass
    Complex FermionLoop(double x, double m, double nc, double q) const;
    //! Delta alpha of the five light quarks
    Complex Hadronic(double x) const;
    //! the Z coupling of a fermion of physical helicity l (e^2 units)
    Complex ZCoupling(const ATOOLS::Flavour &f, int l) const;
  };

}

#endif
