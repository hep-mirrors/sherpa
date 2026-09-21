#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"

using namespace REMNANTS;
using namespace ATOOLS;

Remnants_Parameters* REMNANTS::rempars = nullptr;

remnant_parameters::remnant_parameters(const remnant_parameters& parms)
{
  kT_form   = parms.kT_form;
  kT_recoil = parms.kT_recoil;
  m_form    = parms.m_form;
  for (std::map<std::string, double>::const_iterator pit = parms.params.begin();
       pit != parms.params.end(); pit++) {
    params[pit->first] = pit->second;
  }
  for (std::map<std::string, std::vector<double> >::const_iterator vit =
         parms.param_variations.begin();
       vit != parms.param_variations.end(); vit++) {
    param_variations[vit->first] = vit->second;
  }
}

Remnants_Parameters::Remnants_Parameters()
{
  SetNucleonDefaults();
  SetMesonDefaults();
  SetPhotonDefaults();
  SetLeptonDefaults();
}

Remnants_Parameters::~Remnants_Parameters()
{
  for (std::map<Flavour, remnant_parameters*>::iterator flit = m_defaults.begin();
       flit != m_defaults.end(); flit++)
    delete flit->second;
  m_defaults.clear();
}

void Remnants_Parameters::SetNucleonDefaults()
{
  remnant_parameters* parmsP                = new remnant_parameters;
  parmsP->kT_form                           = primkT_form::gauss_limited;
  parmsP->kT_recoil                         = primkT_recoil::beam_vs_shower;
  parmsP->params["SHOWER_INITIATOR_MEAN"]   = 1.00;
  parmsP->params["SHOWER_INITIATOR_SIGMA"]  = 1.10;
  parmsP->params["SHOWER_INITIATOR_Q2"]     = 0.77;
  parmsP->params["SHOWER_INITIATOR_KTMAX"]  = 2.70;
  parmsP->params["SHOWER_INITIATOR_KTEXPO"] = 5.12;
  parmsP->params["BEAM_SPECTATOR_MEAN"]     = 0.00;
  parmsP->params["BEAM_SPECTATOR_SIGMA"]    = 0.25;
  parmsP->params["BEAM_SPECTATOR_Q2"]       = 0.77;
  parmsP->params["BEAM_SPECTATOR_KTMAX"]    = 1.00;
  parmsP->params["BEAM_SPECTATOR_KTEXPO"]   = 5.00;
  parmsP->params["REFERENCE_ENERGY"]        = 7000.;
  parmsP->params["ENERGY_SCALING_EXPO"]     = 0.08;
  parmsP->m_form                            = matter_form::double_gaussian;
  parmsP->params["MATTER_RADIUS_1"]         = 0.85;
  parmsP->params["MATTER_RADIUS_2"]         = 1.00;
  parmsP->params["MATTER_FRACTION_1"]       = 0.65;
  parmsP->param_variations["MATTER_RADIUS_1"]   = {parmsP->params["MATTER_RADIUS_1"]};
  parmsP->param_variations["MATTER_RADIUS_2"]   = {parmsP->params["MATTER_RADIUS_2"]};
  parmsP->param_variations["MATTER_FRACTION_1"] = {parmsP->params["MATTER_FRACTION_1"]};
  parmsP->params["SOFT_EXPONENT"]           = 0.08;
  parmsP->param_variations["SOFT_EXPONENT"]     = {parmsP->params["SOFT_EXPONENT"]};
  m_defaults[Flavour(kf_p_plus)]            = parmsP;
  m_defaults[Flavour(kf_p_plus).Bar()]      = new remnant_parameters(*parmsP);
  m_defaults[Flavour(kf_n)]                 = new remnant_parameters(*parmsP);
  m_defaults[Flavour(kf_n).Bar()]           = new remnant_parameters(*parmsP);
}

void Remnants_Parameters::SetMesonDefaults()
{
  remnant_parameters* parmsP                = new remnant_parameters;
  parmsP->kT_form                           = primkT_form::gauss_limited;
  parmsP->kT_recoil                         = primkT_recoil::beam_vs_shower;
  parmsP->params["SHOWER_INITIATOR_MEAN"]   = 1.00;
  parmsP->params["SHOWER_INITIATOR_SIGMA"]  = 1.10;
  parmsP->params["SHOWER_INITIATOR_Q2"]     = 0.77;
  parmsP->params["SHOWER_INITIATOR_KTMAX"]  = 2.70;
  parmsP->params["SHOWER_INITIATOR_KTEXPO"] = 5.12;
  parmsP->params["BEAM_SPECTATOR_MEAN"]     = 0.00;
  parmsP->params["BEAM_SPECTATOR_SIGMA"]    = 0.25;
  parmsP->params["BEAM_SPECTATOR_Q2"]       = 0.77;
  parmsP->params["BEAM_SPECTATOR_KTMAX"]    = 1.00;
  parmsP->params["BEAM_SPECTATOR_KTEXPO"]   = 5.00;
  parmsP->params["REFERENCE_ENERGY"]        = 7000.;
  parmsP->params["ENERGY_SCALING_EXPO"]     = 0.08;
  parmsP->m_form                            = matter_form::single_gaussian;
  parmsP->params["MATTER_RADIUS_1"]         = 0.75;
  parmsP->params["MATTER_RADIUS_2"]         = 0.00;
  parmsP->params["MATTER_FRACTION_1"]       = 1.00;
  parmsP->param_variations["MATTER_RADIUS_1"]   = {parmsP->params["MATTER_RADIUS_1"]};
  parmsP->param_variations["MATTER_RADIUS_2"]   = {parmsP->params["MATTER_RADIUS_2"]};
  parmsP->param_variations["MATTER_FRACTION_1"] = {parmsP->params["MATTER_FRACTION_1"]};
  parmsP->params["SOFT_EXPONENT"]           = 0.;
  parmsP->param_variations["SOFT_EXPONENT"]     = {parmsP->params["SOFT_EXPONENT"]};
  m_defaults[Flavour(kf_pi)]                = parmsP;
  m_defaults[Flavour(kf_pi_plus)]           = new remnant_parameters(*parmsP);
  m_defaults[Flavour(kf_pi_plus).Bar()]     = new remnant_parameters(*parmsP);
}

void Remnants_Parameters::SetPhotonDefaults()
{
  remnant_parameters* parmsP                = new remnant_parameters;
  parmsP->kT_form                           = primkT_form::gauss_limited;
  parmsP->kT_recoil                         = primkT_recoil::beam_vs_shower;
  parmsP->params["SHOWER_INITIATOR_MEAN"]   = 1.00;
  parmsP->params["SHOWER_INITIATOR_SIGMA"]  = 1.10;
  parmsP->params["SHOWER_INITIATOR_Q2"]     = 0.77;
  parmsP->params["SHOWER_INITIATOR_KTMAX"]  = 2.70;
  parmsP->params["SHOWER_INITIATOR_KTEXPO"] = 5.12;
  parmsP->params["BEAM_SPECTATOR_MEAN"]     = 0.00;
  parmsP->params["BEAM_SPECTATOR_SIGMA"]    = 0.25;
  parmsP->params["BEAM_SPECTATOR_Q2"]       = 0.77;
  parmsP->params["BEAM_SPECTATOR_KTMAX"]    = 1.00;
  parmsP->params["BEAM_SPECTATOR_KTEXPO"]   = 5.00;
  parmsP->params["REFERENCE_ENERGY"]        = 7000.;
  parmsP->params["ENERGY_SCALING_EXPO"]     = 0.08;
  parmsP->m_form                            = matter_form::single_gaussian;
  parmsP->params["MATTER_RADIUS_1"]         = 0.75;
  parmsP->params["MATTER_RADIUS_2"]         = 0.00;
  parmsP->params["MATTER_FRACTION_1"]       = 1.00;
  parmsP->param_variations["MATTER_RADIUS_1"]   = {parmsP->params["MATTER_RADIUS_1"]};
  parmsP->param_variations["MATTER_RADIUS_2"]   = {parmsP->params["MATTER_RADIUS_2"]};
  parmsP->param_variations["MATTER_FRACTION_1"] = {parmsP->params["MATTER_FRACTION_1"]};
  parmsP->params["SOFT_EXPONENT"]           = 0.;
  parmsP->param_variations["SOFT_EXPONENT"]     = {parmsP->params["SOFT_EXPONENT"]};
  m_defaults[Flavour(kf_photon)]            = parmsP;
}

void Remnants_Parameters::SetLeptonDefaults()
{
  remnant_parameters* parmsE                = new remnant_parameters;
  parmsE->kT_form                           = primkT_form::none;
  parmsE->kT_recoil                         = primkT_recoil::beam_vs_shower;
  parmsE->params["SHOWER_INITIATOR_MEAN"]   = 0.00;
  parmsE->params["SHOWER_INITIATOR_SIGMA"]  = 0.00;
  parmsE->params["SHOWER_INITIATOR_Q2"]     = 0.00;
  parmsE->params["SHOWER_INITIATOR_KTMAX"]  = 0.00;
  parmsE->params["SHOWER_INITIATOR_KTEXPO"] = 0.00;
  parmsE->params["BEAM_SPECTATOR_MEAN"]     = 0.00;
  parmsE->params["BEAM_SPECTATOR_SIGMA"]    = 0.00;
  parmsE->params["BEAM_SPECTATOR_Q2"]       = 0.00;
  parmsE->params["BEAM_SPECTATOR_KTMAX"]    = 0.00;
  parmsE->params["BEAM_SPECTATOR_KTEXPO"]   = 0.00;
  parmsE->params["REFERENCE_ENERGY"]        = 0.00;
  parmsE->params["ENERGY_SCALING_EXPO"]     = 0.00;
  parmsE->m_form                            = matter_form::none;
  parmsE->params["MATTER_RADIUS_1"]         = 1.e-12;
  parmsE->params["MATTER_RADIUS_2"]         = 0.;
  parmsE->params["MATTER_FRACTION_1"]       = 1.00;
  parmsE->param_variations["MATTER_RADIUS_1"]   = {parmsE->params["MATTER_RADIUS_1"]};
  parmsE->param_variations["MATTER_RADIUS_2"]   = {parmsE->params["MATTER_RADIUS_2"]};
  parmsE->param_variations["MATTER_FRACTION_1"] = {parmsE->params["MATTER_FRACTION_1"]};
  parmsE->params["SOFT_EXPONENT"]           = 0.;
  parmsE->param_variations["SOFT_EXPONENT"]     = {parmsE->params["SOFT_EXPONENT"]};
  m_defaults[Flavour(kf_e)]                 = parmsE;
  m_defaults[Flavour(kf_e).Bar()]           = new remnant_parameters(*parmsE);
  m_defaults[Flavour(kf_mu)]                = new remnant_parameters(*parmsE);
  m_defaults[Flavour(kf_mu).Bar()]          = new remnant_parameters(*parmsE);
}

remnant_parameters* Remnants_Parameters::Defaults(const ATOOLS::Flavour& flav)
{
  std::map<Flavour, remnant_parameters*>::iterator it = m_defaults.find(flav);
  if (it != m_defaults.end()) return it->second;
  // Fall back by hadron class. Use find() (not operator[]) so a missing
  // fallback fails fast instead of default-inserting a null and dereferencing
  // it in the callers (Get/KT_Form/...).
  Flavour key = flav.IsBaryon() ? Flavour(kf_p_plus)
              : flav.IsMeson()  ? Flavour(kf_pi_plus)
                                : Flavour(kf_e);
  it = m_defaults.find(key);
  if (it == m_defaults.end())
    THROW(fatal_error, "no remnant defaults for " + ToString(flav));
  return it->second;
}

std::vector<double>
Remnants_Parameters::GetVariationVector(const ATOOLS::Flavour& flav,
                                        std::string            keyword)
{
  remnant_parameters* defaults = Defaults(flav);
  std::map<std::string, std::vector<double> >::const_iterator vit =
    defaults->param_variations.find(keyword);
  if (vit == defaults->param_variations.end()) return std::vector<double>();
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  return data[ToString(flav.Kfcode())][keyword]
          .SetDefault(vit->second)
          .GetVector<double>();
}

double Remnants_Parameters::Get(const ATOOLS::Flavour& flav, std::string keyword)
{
  remnant_parameters* defaults = Defaults(flav);
  // Parameters that can be varied on the fly may be given as a list, with the
  // nominal value first - read them as a vector and return the nominal value.
  if (defaults->param_variations.find(keyword) !=
      defaults->param_variations.end()) {
    const std::vector<double> values = GetVariationVector(flav, keyword);
    if (values.empty())
      THROW(fatal_error, "Empty value list for REMNANTS:" +
                             ToString(flav.Kfcode()) + ":" + keyword);
    return values.front();
  }
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  return data[ToString(flav.Kfcode())][keyword]
          .SetDefault(defaults->params.at(keyword))
          .Get<double>();
}

primkT_form Remnants_Parameters::KT_Form(const ATOOLS::Flavour& flav)
{
  if (flav == Flavour(kf_none)) return primkT_form::none;
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  return data[ToString(flav.Kfcode())]["KT_FORM"]
          .SetDefault(Defaults(flav)->kT_form)
          .Get<primkT_form>();
}

primkT_recoil Remnants_Parameters::KT_Recoil(const ATOOLS::Flavour& flav)
{
  if (flav == Flavour(kf_none)) return primkT_recoil::none;
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  return data[ToString(flav.Kfcode())]["KT_RECOIL"]
          .SetDefault(Defaults(flav)->kT_recoil)
          .Get<primkT_recoil>();
}

matter_form Remnants_Parameters::Matter_Form(const ATOOLS::Flavour& flav)
{
  if (flav == Flavour(kf_none)) return matter_form::none;
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  return data[ToString(flav.Kfcode())]["MATTER_FORM"]
          .SetDefault(Defaults(flav)->m_form)
          .Get<matter_form>();
}

std::ostream& REMNANTS::operator<<(std::ostream&                os,
                                   const REMNANTS::primkT_form& form)
{
  switch (form) {
    case primkT_form::none: return os << "None";
    case primkT_form::gauss: return os << "Gauss";
    case primkT_form::gauss_limited: return os << "Gauss_Limited";
    case primkT_form::dipole: return os << "Dipole";
    case primkT_form::dipole_limited: return os << "Dipole_Limited";
    default: break;
  }
  return os << "Undefined";
}

std::ostream& REMNANTS::operator<<(std::ostream&                  os,
                                   const REMNANTS::primkT_recoil& recoil)
{
  switch (recoil) {
    case primkT_recoil::democratic: return os << "Democratic";
    case primkT_recoil::beam_vs_shower: return os << "Beam_vs_Shower";
    default: break;
  }
  return os << "Undefined";
}

std::ostream& REMNANTS::operator<<(std::ostream&                os,
                                   const REMNANTS::matter_form& f)
{
  switch (f) {
    case matter_form::none:                 return os << "None";
    case matter_form::single_gaussian:      return os << "Single_Gaussian";
    case matter_form::double_gaussian:      return os << "Double_Gaussian";
    case matter_form::x_dependent_gaussian: return os << "X-Dependent_Gaussian";
    case matter_form::unknown:              return os << "Unknown";
    default: break;
  }
  return os << "Undefined";
}

std::istream& REMNANTS::operator>>(std::istream&          is,
                                   REMNANTS::primkT_form& form)
{
  std::string tag;
  is >> tag;
  if (tag == "None") form = primkT_form::none;
  else if (tag == "Gauss")
    form = primkT_form::gauss;
  else if (tag == "Gauss_Limited")
    form = primkT_form::gauss_limited;
  else if (tag == "Dipole")
    form = primkT_form::dipole;
  else if (tag == "Dipole_Limited")
    form = primkT_form::dipole_limited;
  else
    form = primkT_form::undefined;
  return is;
}

std::istream& REMNANTS::operator>>(std::istream&            is,
                                   REMNANTS::primkT_recoil& recoil)
{
  std::string tag;
  is >> tag;
  if (tag == "Democratic") recoil = primkT_recoil::democratic;
  else if (tag == "Beam_vs_Shower")
    recoil = primkT_recoil::beam_vs_shower;
  else
    recoil = primkT_recoil::undefined;
  return is;
}

std::istream& REMNANTS::operator>>(std::istream& is, REMNANTS::matter_form& f)
{
  std::string tag;
  is >> tag;
  if (tag == "None") f = matter_form::none;
  else if (tag == "Single_Gaussian") f = matter_form::single_gaussian;
  else if (tag == "Double_Gaussian")
    f = matter_form::double_gaussian;
  else if (tag == "x_Dependent_Gaussian" || tag == "X-Dependent_Gaussian")
    f = matter_form::x_dependent_gaussian;
  else if (tag == "Unknown") f = matter_form::unknown;
  else
    THROW(fatal_error, "Unknown matter form \"" + tag + "\"");
  return is;
}
