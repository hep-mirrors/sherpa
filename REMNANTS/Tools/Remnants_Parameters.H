#ifndef REMNANTS_Tools_Remnants_Parameters_H
#define REMNANTS_Tools_Remnants_Parameters_H

#include "ATOOLS/Phys/Flavour.H"
#include <map>
#include <string>

namespace REMNANTS {
  enum class primkT_form {
    none           = 0,
    gauss          = 1,
    gauss_limited  = 2,
    dipole         = 3,
    dipole_limited = 4,
    undefined      = 99
  };
  enum class primkT_recoil { democratic = 1, beam_vs_shower = 2, undefined = 99 };

  enum class matter_form {
    none            = 0,
    single_gaussian = 1,
    double_gaussian = 2,
    unknown         = 99
  };

  std::ostream& operator<<(std::ostream& os, const primkT_form & form);
  std::ostream& operator<<(std::ostream& os, const primkT_recoil & recoil);
  std::ostream &operator<<(std::ostream& os, const matter_form & matter);
  std::istream& operator>>(std::istream& is, primkT_form& form);
  std::istream& operator>>(std::istream& is, primkT_recoil& recoil);
  std::istream &operator>>(std::istream& is, matter_form & matter);

  struct remnant_parameters {
  public:
    primkT_form   kT_form;
    primkT_recoil kT_recoil;
    matter_form   m_form;
    std::map<std::string,double> params;

    remnant_parameters() :
      kT_form(primkT_form::undefined), kT_recoil(primkT_recoil::undefined),
      m_form(matter_form::unknown) {}
    remnant_parameters(const remnant_parameters & parms);
  };

  std::ostream& operator<<(std::ostream& os, const remnant_parameters & parms);

  class Remnants_Parameters {
    std::map<ATOOLS::Flavour, remnant_parameters *> m_defaults, m_actuals;

    void SetNucleonDefaults();
    void SetMesonDefaults();
    void SetPhotonDefaults();
    void SetLeptonDefaults();
  public:
    Remnants_Parameters();
    ~Remnants_Parameters();
    void Init();
    double Get(const ATOOLS::Flavour & flav, std::string keyword);
    primkT_form   KT_Form(const ATOOLS::Flavour & flav);
    primkT_recoil KT_Recoil(const ATOOLS::Flavour & flav);
    matter_form   Matter_Form(const ATOOLS::Flavour & flav);
    void Output();
  };

  extern Remnants_Parameters * rempars;
}

#endif
