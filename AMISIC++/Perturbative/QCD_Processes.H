#ifndef AMISIC_Perturbative_QED_Processes_H
#define AMISIC_Perturbative_QED_Processes_H

#include "AMISIC++/Perturbative/MI_Process.H"
#include "ATOOLS/Phys/Flavour.H"

namespace AMISIC {
  class gg_gg: public XS_Base {
  public:
    gg_gg();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };
  
  class gg_qqbar: public XS_Base {
  public:
    gg_qqbar();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class qqbar_gg: public XS_Base {    
  public:
    qqbar_gg();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class qg_qg: public XS_Base {
  public:
    qg_qg();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class qq_qq: public XS_Base {
  public:
    qq_qq();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class qqbar_qqbar: public XS_Base {
  public:
    qqbar_qqbar();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class q1q2_q1q2: public XS_Base {
  public:
    q1q2_q1q2();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class q1q1bar_q2q2bar: public XS_Base {
  public:
    q1q1bar_q2q2bar();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector & flavs);
  };

  class XS_gg_g3S1_oct : public XS_Base {
  private:
    size_t m_S;
    double m_alphaS, LDME, m_mass, m_mass2, m_pref;
    kf_code kfc;
  public:
    XS_gg_g3S1_oct();
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);

  };

  class XS_gg_g3S1 : public XS_Base {
  private:
    size_t m_S;
    double m_alphaS, m_R02, m_mass, m_mass2, m_pref;
    kf_code kfc;
  public:
    XS_gg_g3S1(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);

  };
}

#endif
