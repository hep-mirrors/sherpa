#ifndef AMISIC_Perturbative_QED_Processes_H
#define AMISIC_Perturbative_QED_Processes_H

#include "AMISIC++/Perturbative/MI_Process.H"
#include "ATOOLS/Phys/Flavour.H"
#include "EXTRA_XS/Two2Two/XS_Quarkonia.H"

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


  class XS_gg_g3S1 : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g3S1(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  //

  class XS_qg_q1S0_oct : public XS_Base {
  private:
    size_t m_a, m_q, m_g;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qg_q1S0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  class XS_qqbar_g1S0_oct : public XS_Base {
  private:
    size_t m_a;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qqbar_g1S0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_gg_g1S0_oct : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g1S0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_qg_q3S1_oct : public XS_Base {
  private:
    size_t m_a, m_q, m_g;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qg_q3S1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  class XS_qqbar_g3S1_oct : public XS_Base {
  private:
    size_t m_a;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qqbar_g3S1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_gg_g3S1_oct : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g3S1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_qg_q3P0_oct : public XS_Base {
  private:
    size_t m_a, m_q, m_g;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qg_q3P0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  class XS_qqbar_g3P0_oct : public XS_Base {
  private:
    size_t m_a;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qqbar_g3P0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_gg_g3P0_oct : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g3P0_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_qg_q3P1_oct : public XS_Base {
  private:
    size_t m_a, m_q, m_g;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qg_q3P1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  class XS_qqbar_g3P1_oct : public XS_Base {
  private:
    size_t m_a;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qqbar_g3P1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_gg_g3P1_oct : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g3P1_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_qg_q3P2_oct : public XS_Base {
  private:
    size_t m_a, m_q, m_g;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qg_q3P2_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  class XS_qqbar_g3P2_oct : public XS_Base {
  private:
    size_t m_a;
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_qqbar_g3P2_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };

  class XS_gg_g3P2_oct : public XS_Base {
  private:
    double m_alphaS, m_pref;
    QuarkoniumState Q;
  public:
    XS_gg_g3P2_oct(const ATOOLS::Flavour_Vector & flavs);
    void Calc(const double & s,const double & t,const double & u);
    bool SetColours(const ATOOLS::Flavour_Vector& flavs);
  };
  

}

#endif
