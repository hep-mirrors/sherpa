#ifndef CFPSHOWER_Calculators_Gauge_Base_H
#define CFPSHOWER_Calculators_Gauge_Base_H

#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "CFPSHOWER++/Tools/Color.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace MODEL {
  class Running_AlphaS;
}

namespace CFPSHOWER {
  class Gauge_Base {
  protected:
    kernel_type::code         m_type;
    MODEL::Running_AlphaS   * p_alphaS;
    MODEL::Running_AlphaQED * p_alpha;
    double m_CF, m_CA, m_TR, m_zeta3;
    double m_cplmax, m_K1max, m_K2max, m_charge;
    double m_muR2factor, m_asfactor;
    int    m_orderA, m_orderB;
    bool   m_swap;
    
    std::string  m_name;
    Color_Vector m_colors;
    
    virtual const double Beta0(const double & NF) const;
    virtual const double Beta1(const double & NF) const;
    virtual const double NF(const double & q2) const;
    virtual const double K(const double & NF) const;
    virtual const double K1(const double & NF) const;
    virtual const double K2(const double & NF) const;
  public:
    Gauge_Base(const Kernel_Info & info);
    virtual double operator()(const Splitting & split);
    virtual const double K(const Splitting & split) const;
    virtual const double KMax(const Splitting & split) const;
    virtual inline const double CplMax() { return m_cplmax; }
    virtual const double Scale(const Splitting & split) const;
    virtual const double OverEstimate(const Splitting & split) const;
    virtual inline const double Charge(const double & scale=0.) const {
      return m_charge;
    }
    
    virtual const std::string & Name()       const { return m_name; }
    virtual void SetName(const std::string & name) { m_name = name; }
    virtual const Color & GetColor(size_t i) const { return m_colors[i]; }
    
    virtual bool SetColours(Splitting & split)=0;
  };

  typedef ATOOLS::Getter_Function<Gauge_Base,Kernel_Info,
				  std::less<std::string> > GP_Getter;
}

#endif
