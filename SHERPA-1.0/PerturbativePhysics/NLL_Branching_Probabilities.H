#ifndef NLL_Branching_Probabilities_H
#define NLL_Branching_Probabilities_H

#include "NLL_Branching_Probability_Base.H"
#include "Flavour.H"

namespace MODEL {
  class Running_AlphaS;
}

namespace SHERPA {
  class BP {
  public:
    enum code {
      gamma  = 0,
      gammaq = 1,
      gammag = 2,
      gammaf = 3,
      gamma_corr  = 4,
      gammaq_corr = 5,
      gammag_corr = 6,
      gammaf_corr = 7,
      gamma_powercorr  = 8,
      gammaq_powercorr = 9,
      gammag_powercorr = 10,
      gammaf_powercorr = 11,
      gamma_kfac  = 16,
      gammaq_kfac = 17,
      gammag_kfac = 18,
      gammaf_kfac = 19,
      gamma_cut  = 32,
      gammaq_cut = 33,
      gammag_cut = 34,
      gammaf_cut = 35,
      gamma_kinlim  = 64,
      gammaq_kinlim = 65,
      gammag_kinlim = 66,
      gammaf_kinlim = 67
    };
  };

  BP::code operator|(BP::code a, BP::code b);

  class Gamma_Lambda_Base : public NLL_Branching_Probability_Base {
  protected:
    BP::code  m_mode;
    double m_cc;
    double m_f1, m_f2, m_f3, m_f4;
    double m_lambda;
    double m_qlimit;
    int    m_powercorr;
    int    m_nlo;
    double m_kfac;
    MODEL::Running_AlphaS * p_runas;
  public:
    Gamma_Lambda_Base(BP::code mode, double lambda, MODEL::Running_AlphaS * =0, int nf=-1);
    double AlphaS(double q2);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
    double GammaF(double q, double Q);
    double IntGammaF(double q, double Q);
  };

  class GammaQ_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaQ_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };

  class GammaG_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaG_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * =0);
  };

  class GammaF_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaF_Lambda(BP::code mode, double lambda, MODEL::Running_AlphaS * =0, int nf=-1);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
  };

  class Gamma_AlphaS_Base : public NLL_Branching_Probability_Base {
    BP::code  m_mode;
    double m_cc;
    double m_f1, m_f2, m_f3;

    int    m_nlo;
    double m_kfac;
    double m_asmu, m_mu2;
  public:
    Gamma_AlphaS_Base(BP::code mode, double asmu, double mu2, int nf=-1);
    double AlphaS(double q2);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
    double GammaF(double q, double Q);
    double IntGammaF(double q, double Q);
  };

  class GammaQ_AlphaS :  public Gamma_AlphaS_Base {
  public:
    GammaQ_AlphaS(BP::code mode, double asmu, double mu2);
  };

  class GammaG_AlphaS :  public Gamma_AlphaS_Base {
  public:
    GammaG_AlphaS(BP::code mode, double asmu, double mu2);
  };

  class GammaF_AlphaS :  public Gamma_AlphaS_Base {
  public:
    GammaF_AlphaS(BP::code mode, double asmu, double mu2, int nf);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
  };



  class Gamma_Lambda_Massive : public Gamma_Lambda_Base {
    double m_mass;
  public:
    Gamma_Lambda_Massive(BP::code mode, double lambda, MODEL::Running_AlphaS *, ATOOLS::Flavour fl);
    double GammaQ(double q, double Q);
    double IntGammaQ(double q, double Q);
    double GammaF(double q, double Q);
    double IntGammaF(double q, double Q);
  };

  class GammaQ_Lambda_Massive :  public Gamma_Lambda_Massive {
  public:
    GammaQ_Lambda_Massive(BP::code mode, double lambda, MODEL::Running_AlphaS *, ATOOLS::Flavour fl);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
  };

  class GammaG_Lambda_Massive :  public Gamma_Lambda_Massive {
  public:
    GammaG_Lambda_Massive(BP::code mode, double lambda, MODEL::Running_AlphaS *);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
  };

  class GammaF_Lambda_Massive :  public Gamma_Lambda_Massive {
  public:
    GammaF_Lambda_Massive(BP::code mode, double lambda, MODEL::Running_AlphaS * , ATOOLS::Flavour fl);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);
  };


  inline BP::code operator|(BP::code a, BP::code b) 
  {
    return (BP::code)(int(a)|int(b));
  }

}

#endif
