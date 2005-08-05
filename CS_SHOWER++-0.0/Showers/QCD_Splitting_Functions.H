#ifndef QCD_Splitting_Functions_h
#define QCD_Splitting_Functions_h

#include "Splitting_Function_Base.H"
#include "Flavour.H"

namespace CS_SHOWER {
  struct splittype {
    enum code {
      FF   = 11,
      FI   = 12,
      IF   = 21,
      II   = 22,
      none = 0
    };
  };

  class q_qg_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double);
  public:
    q_qg_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double);
    double OverIntegrated(const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double);
    double Z();
  };

  class q_gq_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double);
  public:
    q_gq_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double);
    double OverIntegrated(const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double);
    double Z();
  };

  class g_gg_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double);
  public:
    g_gg_FF();
    double operator() (const double,const double);
    double OverIntegrated(const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double);
    double Z();
  };

  class g_qq_FF : public Splitting_Function_Base {
    double m_Jmax;
    double J(const double);
  public:
    g_qq_FF(const ATOOLS::Flavour &);
    double operator() (const double,const double);
    double OverIntegrated(const double,const double,const double);
    double Overestimated(const double,const double);
    double RejectionWeight(const double,const double);
    double Z();
  };
}

#endif
