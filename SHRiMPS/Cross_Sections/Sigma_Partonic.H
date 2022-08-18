#ifndef SHRIMPS_Cross_Sections_Sigma_Partonic_H
#define SHRIMPS_Cross_Sections_Sigma_Partonic_H
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  class Remnant_Handler;
  
  struct xs_mode {
    enum code {
      perturbative = 2,
      integrated   = 1,
      Regge        = 0
    };
  };
  
  class Sigma_Partonic {
    Continued_PDF          * p_pdf[2];
    MODEL::Strong_Coupling * p_alphaS;
    xs_mode::code            m_mode;
    bool                     m_fixflavour;
    double          m_Ymax, m_S, m_eta, m_smin, m_tmin;
    double          m_accu, m_sigma, m_maxdsigma;
    int             m_Nmaxtrials;
    long int        m_kinX_fails;
    double          m_y[2], m_yhat, m_dy, m_coshdy, m_pt2, m_phi, m_shat, m_that;
    double          m_dsigma, m_x[2], m_xpdf[2];
    ATOOLS::Flavour m_flavs[2];

    bool                                       m_ana;
    std::map<std::string,ATOOLS::Histogram * > m_histos;


    inline const double ColourFactor(const ATOOLS::Flavour & flav) const {
      if (flav.IsQuark()) return 4./3.;
      if (flav.IsGluon()) return 3.;
      return 0.;
    }
    const double ME2(const double & shat,const double & that,const double & scale);
    const double MakePoint();
    const bool   Calculate();
    void  SelectFlavours(const bool & fixflavour);
  public:
    Sigma_Partonic(const xs_mode::code & mode=xs_mode::Regge);
    ~Sigma_Partonic();

    void Initialise(Remnant_Handler * remnants); 
    const double MakeEvent();
    const double dSigma();
    const double dSigma(const double & shat,const double & yhat);
    const double PDF(const size_t beam,const double & x,const double & Q2,
		     const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_gluon));
    void SetAlphaS(MODEL::Strong_Coupling * as)          { p_alphaS = as; }
    void SetPDF(const size_t & beam,Continued_PDF * pdf) { p_pdf[beam] = pdf; }
    void SetEta(const double & eta=0.)   { m_eta  = eta; }
    void SetTmin(const double & tmin=1.) { m_tmin = tmin; }
    void SetSmin(const double & smin=1.) { m_smin = smin; }
    const double          & SigmaTot() const { return m_sigma; }
    const double          & SigmaHat() const { return m_dsigma; }
    const double          & SHat()     const { return m_shat; }
    const double          & YHat()     const { return m_yhat; }
    const double          & Smin()     const { return m_smin; }
    const double          & Tmin()     const { return m_tmin; }
    const double          & PT2()      const { return m_pt2; }
    const double          & Phi()      const { return m_phi; }
    const double          & Y(const size_t & beam)    const { return m_y[beam]; }
    const double          & X(const size_t & beam)    const { return m_x[beam]; }
    const ATOOLS::Flavour & Flav(const size_t & beam) const { return m_flavs[beam]; }
  };
}
#endif
