#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

#include "AMISIC++/Tools/MI_Parameters.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
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

  class MI_Processes;

  class Hadronic_XSec_Calculator {
  public:
    static double s_slopes[4], s_X[4][4], s_Y[4][4];
    static double s_c[4][4][8], s_d[4][4][9];
  private:
    // all cross sections in mb.
    std::map<ATOOLS::Flavour, std::list<ATOOLS::Flavour> > m_hadroncomponents;
    static std::map<kf_code, size_t> s_indexmap;
    static std::map<kf_code, double> s_fVs;
    double m_beta0[4];

    std::array<ATOOLS::Flavour,2> m_flavs;
    xsec_type      m_type;
    evt_type::code m_evttype;
    double m_mmin, m_smin, m_mass_proton, m_mass_proton2, m_mass_pi;
    std::array<double, 2> m_masses, m_masses2;
    double m_mres, m_cres, m_s0, m_s1, m_Ypp, m_c0, m_c1;
    double m_prefElastic, m_prefSD, m_prefDD;
    double m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_eta_reggeon;
    double m_alphaQED;
    double m_s, m_xstot, m_xsel, m_xssdA, m_xssdB, m_xsdd, m_xsnd, m_xsnd_norm;
    std::map<std::pair<ATOOLS::Flavour, ATOOLS::Flavour>, double> m_sigmaEl;
    std::map<std::pair<ATOOLS::Flavour, ATOOLS::Flavour>, double> m_sigmaSDA;
    std::map<std::pair<ATOOLS::Flavour, ATOOLS::Flavour>, double> m_sigmaSDB;
    std::map<std::pair<ATOOLS::Flavour, ATOOLS::Flavour>, double> m_sigmaDD;
    double m_GeV2mb;
    OneDim_Table * p_xsratio, * p_xshard;
    size_t m_testmode;

    void   FixType();
    void   FixTables();
    void   TestXSecs();
    void   CalculateXSecs();
    double TotalXSec(const size_t hadtags[2]) const;
    double IntElXSec(const size_t hadtags[2], const double & xstot) const;
    double IntSDXSec(const size_t hadtags[2], const size_t & diff,
		     const double masses[2]) const;
    double IntDDXSec(const size_t hadtags[2], const double masses[2]) const;
    double CalculateDoubleDXSec() const;
    void OutputXSratios(axis * sbins);
  public:
    explicit Hadronic_XSec_Calculator();
    ~Hadronic_XSec_Calculator();

    void Initialize(const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2,
		    MODEL::Model_Base * model,const evt_type::code & type);
    void operator()(double s);
    void CalculateXSratios(MI_Processes * processes,axis * sbins);
    void Output() const;

    inline const double & Smin()   const { return m_smin; }
    //////////////////////////////////////////////////////////////////////////
    // The ratio of parton-level cross section obtained by integration over
    // the partonic phase space and the (renormalised) non-diffractive cross
    // section
    //////////////////////////////////////////////////////////////////////////
    inline OneDim_Table * GetXSratio() { return p_xsratio; }
    inline double XSratio(const double & s=-1.) { return (*p_xsratio)(s); }
    //////////////////////////////////////////////////////////////////////////
    // Hard cross section
    //////////////////////////////////////////////////////////////////////////
    inline double XShard(const double & s=-1.)  { return (*p_xshard)(s); }
    //////////////////////////////////////////////////////////////////////////
    // Non-diffractive cross section, given in GeV^-2
    //////////////////////////////////////////////////////////////////////////
    inline double XSnd(const double & s=-1.) {
      if (s>0.) (*this)(s); return m_xsnd;
    }
    //////////////////////////////////////////////////////////////////////////
    // Total cross section, given in mb
    //////////////////////////////////////////////////////////////////////////
    inline double XStot(const double & s=-1.) {
      if (s>0.) (*this)(s); return m_xstot;
    }
    //////////////////////////////////////////////////////////////////////////
    // Total cross section, given in mb
    //////////////////////////////////////////////////////////////////////////
    inline const double SimulatedXSec(const double & s=-1.) {
      if (s>0. && m_evttype!=evt_type::Perturbative) (*this)(s);
      double xsec = 0;
      if (m_evttype==evt_type::Perturbative) xsec += (*p_xshard)(s);
      if (m_evttype==evt_type::Elastic       ||
	  m_evttype==evt_type::QuasiElastic) xsec += SigmaEl();
      if (m_evttype==evt_type::DiffractiveA  ||
	  m_evttype==evt_type::QuasiElastic) xsec += SigmaSDA();
      if (m_evttype==evt_type::DiffractiveB  ||
	  m_evttype==evt_type::QuasiElastic) xsec += SigmaSDB();
      if (m_evttype==evt_type::DiffractiveAB ||
	  m_evttype==evt_type::QuasiElastic) xsec += SigmaDD();
      return xsec;
    }
    //////////////////////////////////////////////////////////////////////////
    // Tuning value to (re-)normalise the non-diffractive cross section,
    // default is 1.
    //////////////////////////////////////////////////////////////////////////
    inline const double & XSndNorm()         const { return m_xsnd_norm; }
    inline const double & EpsPomeron()       const { return m_eps_pomeron; }
    inline const double & AlphaPPomeron()    const { return m_alphaP_pomeron; }
    inline const double & TriplePomeron()    const { return m_triple_pomeron; }
    inline const double & AlphaQED()         const { return m_alphaQED; }
    inline const double & Diffractive_cres() const { return m_cres; }
    inline const double & Diffractive_Mres() const { return m_mres; }
    inline const double & Diffractive_s1()   const { return m_s1; }

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
	  s_indexmap.find(flav)!=s_indexmap.end()) return s_indexmap[flav];
      return 0;
    }
  };
}

#endif
