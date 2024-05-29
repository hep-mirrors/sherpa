#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

#include "AMISIC++/Tools/MI_Parameters.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
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

  class MI_Processes;
  
  class Hadronic_XSec_Calculator {
  public:
    static double s_slopes[4], s_X[4][4], s_Y[4][4];
    static double s_c[4][4][8], s_d[4][4][9];
  private:
    // all cross sections in mb.
    std::map<ATOOLS::Flavour, size_t> m_indexmap;
    std::map<ATOOLS::Flavour, double> m_fVs;
    double m_beta0[4];
    void   FixTables();

    std::array<ATOOLS::Flavour,2> m_flavs;
    xsec_type                      m_type;
    double m_mass_proton, m_mass_proton2, m_mass_pi;
    std::array<double, 2> m_masses, m_masses2;
    double m_mres, m_cres, m_s1, m_Ypp, m_c0, m_c1;
    double m_prefElastic, m_prefSD, m_prefDD;
    double m_eps_pomeron, m_alphaP_pomeron, m_triple_pomeron, m_eta_reggeon;
    double m_alphaQED;
    double m_s, m_xstot, m_xsel, m_xssdA, m_xssdB, m_xsdd, m_xsnd, m_xsnd_norm;
    double m_GeV2mb;
    OneDim_Table * p_xsratio, * p_xshard; 
    size_t m_testmode;

    void   FixType();
    void   TestXSecs();
    void   CalculateHHXSecs();
    void   CalculateHGammaXSecs(const size_t photon);
    void   CalculatePhotonPhotonXSecs();
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
		    MODEL::Model_Base * model);
    void operator()(double s);
    void CalculateXSratios(MI_Processes * processes,axis * sbins);
    void Output() const;
    //////////////////////////////////////////////////////////////////////////
    // The ratio of parton-level cross section obtained by integration over 
    // the partonic phase space and the (renormalised) non-diffractive cross
    // section
    //////////////////////////////////////////////////////////////////////////
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
    // Tuning value to (re-)normalise the non-diffractive cross section,
    // default is 1.
    //////////////////////////////////////////////////////////////////////////
    inline const double & XSndNorm() const { return m_xsnd_norm; }
  };
}

#endif
