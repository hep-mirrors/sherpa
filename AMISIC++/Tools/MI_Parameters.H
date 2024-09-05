#ifndef AMISIC_Tools_MI_Parameters_H
#define AMISIC_Tools_MI_Parameters_H

#include <map>
#include <string>

namespace AMISIC {
  class axis;
}

namespace AMISIC {
  struct evt_type {
    enum code {
      Perturbative    = 0,
      NonPerturbative = 1,
      Elastic         = 11,
      DiffractiveA    = 12,
      DiffractiveB    = 13,
      DiffractiveAB   = 14,
      QuasiElastic    = 15,
      AllMinimumBias  = 20
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

  enum class scale_scheme {
    PT           = 1,
    PT_with_Raps = 2
  };
  std::ostream &operator<<(std::ostream&, const scale_scheme&);
  std::istream &operator>>(std::istream&, scale_scheme&);
  struct two_pions {
    enum code {
      none            = 0,
      rho_only        = 1,
      rho_omega       = 2,
      rho_omega_cont  = 3,
      cont_only       = 4
    };
  };

  class MI_Parameters {
  private:
    double m_pt02ref, m_pt02IR, m_ptmin2ref, m_ptmin2IR;
    double m_Eref, m_Sref, m_Ecms, m_Scms, m_eta;
    evt_type::code  m_evttype;
    scale_scheme    m_scaleRscheme, m_scaleFscheme;
    two_pions::code m_twopions;
    std::map<std::string,double> m_parameters;

    axis * p_sbins;
  public:
    MI_Parameters();
    ~MI_Parameters() = default;

    double operator()(const std::string& keyword) const;
    double CalculatePT02(const double & s=-1) const;

    const evt_type::code  & GetEvtType()          const { return m_evttype; }
    const scale_scheme    & GetScaleRScheme()     const { return m_scaleRscheme; }
    const scale_scheme    & GetScaleFScheme()     const { return m_scaleFscheme; }
    const two_pions::code & GetTwoPionTreatment() const { return m_twopions; }
  };

  extern const MI_Parameters * mipars;
}
#endif
