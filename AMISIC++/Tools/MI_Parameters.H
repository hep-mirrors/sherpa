#ifndef AMISIC_Tools_MI_Parameters_H
#define AMISIC_Tools_MI_Parameters_H

#include <map>
#include <string>

namespace AMISIC {
  class axis;
}

namespace AMISIC {
  enum class mitype {
    none          = 0,
    hadron_hadron = 1,
    gamma_hadron  = 2,
    gamma_gamma   = 3
  };

  enum class scale_scheme {
    PT           = 1,
    PT_with_Raps = 2
  };
  std::ostream &operator<<(std::ostream&, const scale_scheme&);
  std::istream &operator>>(std::istream&, scale_scheme&);

  class MI_Parameters {
  private:
    double m_pt02ref, m_pt02IR, m_ptmin2ref, m_ptmin2IR;
    double m_Eref, m_Sref, m_Ecms, m_Scms, m_eta;
    scale_scheme m_scaleRscheme, m_scaleFscheme;
    std::map<std::string,double> m_parameters;

    axis * p_sbins;
  public:
    MI_Parameters();
    ~MI_Parameters() = default;

    double operator()(const std::string& keyword) const;
    double CalculatePT02(const double & s=-1) const;

    const scale_scheme & GetScaleRScheme() const { return m_scaleRscheme; }
    const scale_scheme & GetScaleFScheme() const { return m_scaleFscheme; }
  };

  extern const MI_Parameters * mipars;
}
#endif
