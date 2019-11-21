#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

namespace AMISIC {
  class Hadronic_XSec_Calculator {
  private:
    // all cross sections in mb.
    double m_Ecms, m_s, m_massp, m_masspi;
    double m_pomeron, m_reggeon, m_slope, m_xsecpom, m_xsecregge;
    double m_xstot, m_xsel, m_xssd, m_xsdd, m_xsnd;
    
    double CalculateTotalXSec();
    double CalculateElasticXSec(const double & xstot);
    double CalculateSingleDXSec();
    double CalculateDoubleDXSec();
  public:
    Hadronic_XSec_Calculator();
    ~Hadronic_XSec_Calculator() {}

    void operator()();

    const double XStot() const { return m_xstot; }
    const double XSel()  const { return m_xsel; }
    const double XSnd()  const { return m_xsnd; }
  };
}

#endif
