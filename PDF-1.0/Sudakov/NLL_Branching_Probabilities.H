#ifndef NLL_Branching_Probabilities_H
#define NLL_Branching_Probabilities_H

#include "NLL_Branching_Probability_Base.H"
//#include "PDF_Base.H"
#include "Flavour.H"

namespace MODEL {
  class Running_AlphaS;
}

namespace SHERPA {
  class BPType {
  public:
    enum code {
      gamma_q2qg  = 0,
      gamma_q2gq  = 1,
      gamma_g2gg  = 2,
      gamma_g2qq  = 3
    };
  };

  class BPMode {
  public:
    enum code {
      none        = 0,
      fs          = 1,
      is          = 2,
      linear_term = 4,
      power_corrs = 8,
      soft_kfac   = 16
    };
  };

  
  class Gamma_Lambda_Base : public NLL_Branching_Probability_Base {
  protected:
    BPType::code m_type;
    BPMode::code m_mode;
    double m_colfac, m_dlog, m_slog, m_power;
    double m_lambda;
    MODEL::Running_AlphaS * p_runas;
    //PDF::PDF_Base         * p_pdf;
    //ATOOLS::Flavour         m_flav;
    double m_kfac;
  public:
    Gamma_Lambda_Base(BPType::code type, BPMode::code mode, double lambda, 
		      MODEL::Running_AlphaS * =0, int nf=-1);
    double AlphaS(double q2);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
    //void   SetPDF(PDF::PDF_Base * pdf,ATOOLS::Flavour & flav) {
    //  p_pdf = pdf; m_flav = flav; 
    //}
  };
  
  class GammaQ_QG_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaQ_QG_Lambda(BPMode::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };

  class GammaQ_GQ_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaQ_GQ_Lambda(BPMode::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };

  class GammaG_GG_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaG_GG_Lambda(BPMode::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };

  class GammaG_QQ_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaG_QQ_Lambda(BPMode::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };


  inline BPMode::code operator|(BPMode::code a, BPMode::code b) { return (BPMode::code)(int(a)|int(b)); }
  inline bool         operator&(BPMode::code a, BPMode::code b) { return (int(a)&int(b)); }
}

#endif
