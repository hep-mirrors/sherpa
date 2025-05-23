#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include <array>

namespace MODEL{ class Model_Base; }

namespace AMISIC {
  enum class xsec_type {
    nucleon_nucleon = 0,
    photon_nucleon  = 1,
    nucleon_photon  = 2,
    photon_photon   = 3,
    hadron_hadron   = 1,
    none            = 99
  };

  class Hadronic_XSec_Calculator {
  private:
    static double s_slopes[4], s_X[4][4], s_Y[4][4], s_c[4][4][8], s_d[4][4][9];
    // all cross sections in mb.
    static std::map<kf_code, size_t> s_indexmap;
    static std::map<kf_code, double> s_fVs;
    double m_beta0[4];

    std::array<ATOOLS::Flavour,2> m_flavs;
    xsec_type                      m_type;
    double m_mass_proton, m_mass_proton2, m_mass_pi;
    std::array<double, 2> m_masses, m_masses2;
    double m_mres, m_cres, m_s1, m_Ypp, m_c0, m_c1, m_prefElastic, m_prefSD, m_prefDD;
    double m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_eta_reggeon;
    double m_alphaQED;
    double m_s, m_xstot, m_xsel, m_xssdA, m_xssdB, m_xsdd, m_xsnd, m_xsnd_norm;
    size_t m_testmode;

    void   FixType();
    void   TestXSecs();
    void   CalculateHHXSecs();
    void   CalculateHGammaXSecs(const size_t photon);
    void   CalculatePhotonPhotonXSecs();
    double TotalXSec(const size_t hadtags[2]) const;
    double IntElXSec(const size_t hadtags[2], const double & xstot) const;
    double IntSDXSec(const size_t hadtags[2], const size_t & diff,const double masses[2]) const;
    double IntDDXSec(const size_t hadtags[2], const double masses[2]) const;
    double CalculateDoubleDXSec() const;
  public:
    explicit Hadronic_XSec_Calculator(MODEL::Model_Base *model,
				      const ATOOLS::Flavour & fl1,
				      const ATOOLS::Flavour & fl2);
    ~Hadronic_XSec_Calculator() = default;

    void operator()(double s);
    void Output() const;
    inline double XSnd(const double & s=-1.) {
      if (s>0.) (*this)(s);
      return m_xsnd;
    }
    inline const double & XSndNorm() const { return m_xsnd_norm; }
  };
}

#endif
