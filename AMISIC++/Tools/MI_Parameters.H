#ifndef AMISIC_Tools_MI_Parameters_H
#define AMISIC_Tools_MI_Parameters_H

#include <map>
#include <string>

namespace AMISIC {
  struct evt_type {
    enum code {
      Perturbative    = 0,
      Elastic         = 1,
      DiffractiveA    = 2,
      DiffractiveB    = 3,
      DiffractiveAB   = 4,
      QuasiElastic    = 5,
      AllMinimumBias  = 10
    };
  };
  std::ostream &operator<<(std::ostream&, const evt_type::code&);
  std::istream &operator>>(std::istream&, evt_type::code&);

  enum class mitype {
    none          = 0,
    hadron_hadron = 1,
    gamma_hadron  = 2,
    gamma_gamma   = 3
  };
  enum class overlap_form {
    Single_Gaussian = 1,
    Double_Gaussian = 2,
    unknown         = -1
  };
  std::ostream& operator<<(std::ostream&, const overlap_form&);
  std::istream& operator>>(std::istream&, overlap_form&);
  enum class scale_scheme { PT, PT_with_Raps };
  std::ostream& operator<<(std::ostream&, const scale_scheme&);
  std::istream& operator>>(std::istream&, scale_scheme&);

  class MI_Parameters {
  private:
    double m_pt02ref, m_pt02IR, m_ptmin2ref, m_ptmin2IR;
    double m_Eref, m_Sref, m_Ecms, m_Scms, m_eta;
    evt_type::code     m_evttype;
    overlap_form m_overlapform;
    scale_scheme m_scalescheme;
    std::map<std::string,double> m_parameters;

  public:
    MI_Parameters();
    ~MI_Parameters() = default;

    double operator()(const std::string& keyword) const;
    double CalculatePT02(const double & s=-1) const;

    const evt_type::code     & GetEvtType()     const { return m_evttype; }
    const overlap_form& GetOverlapForm() const { return m_overlapform; }
    const scale_scheme& GetScaleScheme() const { return m_scalescheme; }
  };

  extern const MI_Parameters * mipars;
}
#endif
