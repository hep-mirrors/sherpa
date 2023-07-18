#ifndef AMISIC_Tools_MI_Parameters_H
#define AMISIC_Tools_MI_Parameters_H

#include <map>
#include <string>

namespace AMISIC {
  struct mitype {
    enum code {
      none           = 0,
      hadron_hadron  = 1,
      gamma_hadron   = 2,
      gamma_gamma    = 3
    };
  };
  struct overlap_form {
    enum code {
      Single_Gaussian = 1,
      Double_Gaussian = 2,
      unknown         = -1
    };
  };
  std::ostream &operator<<(std::ostream&, const overlap_form::code&);
  std::istream &operator>>(std::istream&, overlap_form::code&);
  struct scale_scheme {
    enum class code {
      PT,
      PT_with_Raps
    };
  };
  std::ostream &operator<<(std::ostream&, const scale_scheme::code&);
  std::istream &operator>>(std::istream&, scale_scheme::code&);

  class MI_Parameters {
  private:
    double m_pt02ref, m_pt02IR, m_ptmin2ref, m_ptmin2IR;
    double m_Eref, m_Sref, m_Ecms, m_Scms, m_eta;
    overlap_form::code m_overlapform;
    scale_scheme::code m_scalescheme;
    std::map<std::string,double> m_parameters;

  public:
    MI_Parameters();
    ~MI_Parameters() = default;

    double operator()(const std::string& keyword) const;
    double CalculatePT02(const double & s=-1) const;

    const overlap_form::code & GetOverlapForm() const { return m_overlapform; }
    const scale_scheme::code & GetScaleScheme() const { return m_scalescheme; }
  };
  
  extern const MI_Parameters * mipars;
}
#endif
