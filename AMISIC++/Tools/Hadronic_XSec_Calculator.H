#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include <array>
#include <list>
#include <algorithm>

namespace MODEL{ class Model_Base; }

namespace AMISIC {
  enum class xsec_type {
    nucleon_nucleon = 0,
    photon_nucleon  = 1,
    nucleon_photon  = 2,
    photon_photon   = 3,
    hadron_hadron   = 10,
    none            = 99
  };

  class Hadronic_XSec_Calculator {
  public:
    static double s_slopes[4], s_X[4][4], s_Y[4][4], s_c[4][4][8], s_d[4][4][9];
  private:
    // all cross sections in mb.
    std::map<ATOOLS::Flavour, std::list<ATOOLS::Flavour> > m_hadroncomponents;
    std::map<ATOOLS::Flavour, size_t> m_indexmap;
    std::map<ATOOLS::Flavour, double> m_fVs;
    double m_beta0[4];
    void   FixTables();

    std::array<ATOOLS::Flavour,2> m_flavs;
    xsec_type                      m_type;
    double m_mmin, m_smin, m_mass_proton, m_mass_proton2, m_mass_pi;
    std::array<double, 2> m_masses, m_masses2;
    double m_mres, m_cres, m_s0, m_s1, m_Ypp, m_c0, m_c1, m_prefElastic, m_prefSD, m_prefDD;
    double m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_eta_reggeon;
    double m_alphaQED;
    double m_s, m_xstot, m_xsel, m_xssdA, m_xssdB, m_xsdd, m_xsnd, m_xsnd_norm;
    std::map<std::pair<ATOOLS::Flavour, ATOOLS::Flavour>, double> m_sigmaEl, m_sigmaSDA, m_sigmaSDB, m_sigmaDD;
    size_t m_testmode;

    void   FixType();
    void   TestXSecs();
    void   CalculateXSecs();
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
    inline const double & Smin()   const { return m_smin; }
    inline double XSnd(const double & s=-1.) {
      if (s>0.) (*this)(s);
      return m_xsnd;
    }
    inline const double & XSndNorm()      const { return m_xsnd_norm; }
    inline const double & EpsPomeron()    const { return m_eps_pomeron; }
    inline const double & AlphaPPomeron() const { return m_alphaP_pomeron; }
    inline const double & TriplePomeron() const { return m_triple_pomeron; }
    inline const double & AlphaQED()      const { return m_alphaQED; }

    inline const double & SigmaEl()  const { return m_xsel; }
    inline const double & SigmaSDA() const { return m_xssdA; }
    inline const double & SigmaSDB() const { return m_xssdB; }
    inline const double & SigmaDD()  const { return m_xsdd; }

    bool SelectEl(std::array<ATOOLS::Flavour, 2> & flavs);
    bool SelectSDA(std::array<ATOOLS::Flavour, 2> & flavs);
    bool SelectSDB(std::array<ATOOLS::Flavour, 2> & flavs);
    bool SelectDD(std::array<ATOOLS::Flavour, 2> & flavs);

    inline const std::array<ATOOLS::Flavour,2> GetFlavs() const { return m_flavs; }

    inline size_t Index(ATOOLS::Flavour & flav,size_t & beam) {
      if (std::find(m_hadroncomponents[m_flavs[beam]].begin(),
		    m_hadroncomponents[m_flavs[beam]].end(), flav) !=
	  m_hadroncomponents[m_flavs[beam]].end() &&
	  m_indexmap.find(flav)!=m_indexmap.end()) return m_indexmap[flav];
      return 0;
    }
  };
}

#endif
