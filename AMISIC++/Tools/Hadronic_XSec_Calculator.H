#ifndef AMISIC_Tools_Hadronic_XSec_Calculator_H
#define AMISIC_Tools_Hadronic_XSec_Calculator_H

namespace AMISIC {
  class Hadronic_XSec_Calculator {
  private:
    // all cross sections in mb.
    double m_s, m_massp, m_masspi;
    double m_pomeron, m_reggeon, m_slope, m_xsecpom, m_xsecregge;
    double m_xstot, m_xsel, m_xssd, m_xsdd, m_xsnd;
    
    double CalculateTotalXSec() const;
    double CalculateElasticXSec(const double & xstot) const;
    double CalculateSingleDXSec() const;
    double CalculateDoubleDXSec() const;
  public:
    explicit Hadronic_XSec_Calculator(int type);
    ~Hadronic_XSec_Calculator() = default;

    void operator()(double s);
    void Output() const;

    double XSnd()  const { return m_xsnd; }
  };
}

#endif
