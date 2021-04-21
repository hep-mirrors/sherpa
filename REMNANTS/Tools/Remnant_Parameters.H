#ifndef REMNANTS_Tools_Remnant_Parameters_H
#define REMNANTS_Tools_Remnant_Parameters_H

#include <map>
#include <string>

namespace REMNANTS {
  struct prim_kperp_form {
    enum code {
      undefined      = 0,
      gauss          = 1,
      gauss_limited  = 2,
      dipole         = 3,
      dipole_limited = 4,
    };
  };

  class Remnant_Parameters {
  private:
    std::string           m_defform;
    prim_kperp_form::code m_form[2];
    double m_defmean, m_defsigma, m_refE, m_scaleexpo, m_defQ2, m_defktmax, m_defeta;
    overlap_form::code    m_overlapform;
    std::map<std::string,double> m_parameters;

    void                  ReadParameters(ATOOLS::Default_Reader *const defaultreader);
    prim_kperp_form::code SelectForm(const std::string & form);
  public:
    Remnant_Parameters();
    ~Remnant_Parameters() {}

    bool Init(A);
    double operator()(std::string keyword);
    const prim_kperp_form::code & GetKTForm(const size_t & beam) { return m_form[beam]; }
    const overlap_form::code    & GetBTForm()                    { return m_overlapform; }
  };
  
  extern Remnant_Parameters * repars;
}
#endif
