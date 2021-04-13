#ifndef SHRIMPS_Cross_Sections_Sigma_Partonic_H
#define SHRIMPS_Cross_Sections_Sigma_Partonic_H
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>

namespace SHRIMPS {
  struct xs_mode {
    enum code {
      perturbative = 1,
      Regge        = 0
    };
  };
  
  class Sigma_Partonic {
    Continued_PDF          * p_pdf[2];
    MODEL::Strong_Coupling * p_alphaS;
    xs_mode::code   m_mode;
    double          m_Ymax, m_S, m_eta, m_smin, m_tmin;
    double          m_accu, m_sigma, m_maxdsigma;
    long int        m_Nmaxtrials;
    double          m_x[2], m_xpdf[2];
    ATOOLS::Flavour m_flavs[2];
    

    inline const double ColourFactor(const ATOOLS::Flavour & flav) const {
      if (flav.IsQuark()) return 4./3.;
      if (flav.IsGluon()) return 3.;
      return 0.;
    }
    const double ME2(const double & s);
    const double MakePoint(double & s,double & y);
    const double dSigma(const double & s,const double & y);
    const bool   Calculate();
    void  SelectFlavours(const bool & fixflavour);
  public:
    Sigma_Partonic(const xs_mode::code & mode=xs_mode::Regge);
    ~Sigma_Partonic();

    void       Initialise(); 
    const double MakeEvent(const bool & fixflavour);
    void SetPDF(const size_t & beam,Continued_PDF * pdf) { p_pdf[beam] = pdf; }
    void SetEta(const double & eta=0.)   { m_eta  = eta; }
    void SetTmin(const double & tmin=1.) { m_tmin = tmin; }
    void SetSmin(const double & smin=1.) { m_smin = smin; }
    const double & SigmaTot() const { return m_sigma; }
    const double          & X(const size_t beam)    const { return m_x[beam]; }
    const ATOOLS::Flavour & Flav(const size_t beam) const { return m_flavs[beam]; }
  };
}
#endif
