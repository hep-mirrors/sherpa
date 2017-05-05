#ifndef DIRE__Shower__Alpha_QCD_H
#define DIRE__Shower__Alpha_QCD_H

#include "DIRE/Shower/Gauge.H"

namespace MODEL {
  class Running_AlphaS;
}

namespace DIRE {
  
  class Alpha_QCD: public Gauge {
  protected:
    
    MODEL::Running_AlphaS *p_cpl;

    double m_fac, m_max;
    double m_Nc, m_CF, m_CA, m_TR;

    int m_lc;

    bool  m_mab, m_mba, m_maa, m_mbb;
    Color m_ca, m_cb;

    double B0(const double &nf) const;
    double B1(const double &nf) const;

    double G2(const double &nf) const;
    double G3(const double &nf) const;

    double CplFac(const double &scale) const;

  public:

    Alpha_QCD(const Kernel_Key &key);

    double Coupling(const Splitting &s) const;
    double CplMax(const Splitting &s) const;

    double K(const Splitting &s) const;
    double KMax(const Splitting &s) const;

    double Nf(const Splitting &s) const;

    void SetLimits();

    double Solve(const double &cpl) const;

  };// end of class Alpha_QCD

}// end of namespace DIRE

#endif
