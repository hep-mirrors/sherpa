#ifndef NLL_Branching_Probabilities_H
#define NLL_Branching_Probabilities_H

#include "NLL_Branching_Probability_Base.H"

namespace MODEL { class Running_AlphaS; }

namespace SHERPA {

  class bpt {
  public:
    enum code {
      gamma_q2qg  = 0,
      gamma_q2gq  = 1,
      gamma_g2gg  = 2,
      gamma_g2qq  = 3
    };
  };

  class bpm {
  public:
    enum code {
      none        = 0,
      fs          = 1,
      is          = 2,
      z_shift     = 4,
      power_corrs = 8,
      soft_kfac   = 16,
      massive     = 32,
      dead_cone   = 64
    };
  };

  
  class Gamma_Lambda_Base : public NLL_Branching_Probability_Base {
  protected:

    MODEL::Running_AlphaS *p_runas;

    bpt::code m_type;
    bpm::code m_mode;

    double m_colfac, m_dlog, m_slog, m_power[5], m_qmass;
    double m_lambda, m_as_factor, m_kfac, m_lastas, m_z[2];

    int m_pca;

    void SetZBounds(const double &q, const double &Q,
		    const double &m,const double &m1,const double &m2);

  public:

    Gamma_Lambda_Base(bpt::code type, bpm::code mode, double lambda, 
		      MODEL::Running_AlphaS * =0, double qmass=0.0, 
		      double asfac=1.0);
    double AlphaS(double q2);
    double Gamma(double q, double Q);
    double IntGamma(double q, double Q);

  };
  
  class GammaQ_QG_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaQ_QG_Lambda(bpm::code mode, double lambda, 
		     MODEL::Running_AlphaS * =0,double qmass=0.0,
		     double asfac=1.0);
    double Gamma(double q, double Q);
  };

  class GammaQ_GQ_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaQ_GQ_Lambda(bpm::code mode, double lambda, 
		     MODEL::Running_AlphaS * =0,double qmass=0.0,
		     double asfac=1.0);
    double Gamma(double q, double Q);
  };

  class GammaG_GG_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaG_GG_Lambda(bpm::code mode, double lambda, 
		     MODEL::Running_AlphaS * =0,double asfac=1.0);
    double Gamma(double q, double Q);
  };

  class GammaG_QQ_Lambda :  public Gamma_Lambda_Base {
  public:
    GammaG_QQ_Lambda(bpm::code mode, double lambda, 
		     MODEL::Running_AlphaS * =0,double qmass=0.0,
		     double asfac=1.0);
    double Gamma(double q, double Q);
  };

  inline bpm::code operator|(bpm::code a, bpm::code b) 
  { return (bpm::code)(int(a)|int(b)); }
  inline bool operator&(bpm::code a, bpm::code b) 
  { return (int(a)&int(b)); }

}// end of namespace SHERPA

#endif
