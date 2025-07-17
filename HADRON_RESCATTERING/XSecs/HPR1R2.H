#ifndef HADRON_RESCATTERING_XSecs_HPR1R2_H
#define HADRON_RESCATTERING_XSecs_HPR1R2_H

#include <vector>
#include <map>

namespace HADRON_RESCATTERING {
  struct hpr1r2 {
    enum code {
      unknown = 0,
      pp      = 1,
      pn      = 2,
      ppbar   = 3
    };
  };
    
  class HPR1R2 {
    double m_eta1, m_eta2;
    std::map<hpr1r2::code, std::vector<double> > m_params;

    void InitParams();
  public:
    HPR1R2();
    ~HPR1R2();

    double xs_tot(hpr1r2::code tag, const double & s);
    double xs_el(hpr1r2::code tag, const double & s);
  };
}

#endif

