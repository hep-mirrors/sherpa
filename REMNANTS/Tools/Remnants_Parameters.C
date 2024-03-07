#include "ATOOLS/Org/Message.H"
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
  if (!m_defaults.empty()) {
    for (std::map<Flavour, remnant_parameters*>::iterator flit =
                 m_defaults.begin();
         flit != m_defaults.end(); flit++)
      delete flit->second;
    m_defaults.clear();
  }
  if (!m_actuals.empty()) {
    for (std::map<Flavour, remnant_parameters*>::iterator flit =
                 m_actuals.begin();
         flit != m_actuals.end(); flit++)
      delete flit->second;
    m_actuals.clear();
  }
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
  parmsP->m_form                            = matter_form::single_gaussian;
  parmsP->params["RADIUS_1"]                = 0.86;
  parmsP->params["RADIUS_2"]                = 0.00;
  parmsP->params["FRACTION_1"]              = 1.00;
  m_defaults[Flavour(kf_p_plus)]            = parmsP;
  m_defaults[Flavour(kf_n)]                 = new remnant_parameters(*parmsP);
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
  parmsP->params["RADIUS_1"]                = 0.75;
  parmsP->params["RADIUS_2"]                = 0.00;
  parmsP->params["FRACTION_1"]              = 1.00;
  m_defaults[Flavour(kf_pi)]                = parmsP;
  m_defaults[Flavour(kf_pi_plus)]           = new remnant_parameters(*parmsP);
}

void Remnants_Parameters::SetPhotonDefaults()
{
  remnant_parameters* parmsP                = new remnant_parameters;
  m_defaults[Flavour(kf_photon)]            = parmsP;
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
  parmsP->params["RADIUS_1"]                = 0.75;
  parmsP->params["RADIUS_2"]                = 0.00;
  parmsP->params["FRACTION_1"]              = 1.00;
  m_defaults[Flavour(kf_photon)]            = parmsP;
}

void Remnants_Parameters::SetLeptonDefaults()
{
  remnant_parameters* parmsE                = new remnant_parameters;
  m_defaults[Flavour(kf_e)]                 = parmsE;
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
  parmsE->params["RADIUS_1"]                = 1.e-12;
  parmsE->params["RADIUS_2"]                = 0.;
  parmsE->params["FRACTION_1"]              = 1.00;
  m_defaults[Flavour(kf_e)]                 = parmsE;
  m_defaults[Flavour(kf_mu)]                = new remnant_parameters(*parmsE);
}

void Remnants_Parameters::Init()
{
  msg_Out() << METHOD << "\n";
  Scoped_Settings data = Settings::GetMainSettings()["REMNANTS"];
  for (const auto& pid : data.GetKeys()) {
    kf_code             kf = ToType<kf_code>(pid);
    Flavour             flav(kf);
    remnant_parameters* defaults;
    if (m_defaults.find(flav) != m_defaults.end()) defaults = m_defaults[flav];
    else {
      msg_Error() << "Warning in " << METHOD
                  << ": did not find default settings for " << flav << "\n"
                  << "   Will continue with defaults for a proton and hope for "
                     "the best.\n";
      defaults = m_defaults[Flavour(kf_p_plus)];
    }
    remnant_parameters* actuals = new remnant_parameters;
    ////////////////////////////////////////////////////////////////////////////////////
    // Fix the intrinsic kT parametrization: form and parameters
    ///////////////////////////////////////////////////////////////////////////////////
    actuals->kT_form   =
      (data[pid]["KT_FORM"]
       .SetDefault(defaults->kT_form)
       .Get<primkT_form>());
    actuals->kT_recoil =
      (data[pid]["KT_RECOIL"]
       .SetDefault(defaults->kT_recoil)
       .Get<primkT_recoil>());
    actuals->params["SHOWER_INITIATOR_MEAN"] =
      (data[pid]["SHOWER_INITIATOR_MEAN"]
       .SetDefault(defaults->params["SHOWER_INITIATOR_MEAN"])
       .Get<double>());
    actuals->params["SHOWER_INITIATOR_SIGMA"] =
      (data[pid]["SHOWER_INITIATOR_SIGMA"]
       .SetDefault(defaults->params["SHOWER_INITIATOR_SIGMA"])
       .Get<double>());
    actuals->params["SHOWER_INITIATOR_Q2"] =
      (data[pid]["SHOWER_INITIATOR_Q2"]
       .SetDefault(defaults->params["SHOWER_INITIATOR_Q2"])
       .Get<double>());
    actuals->params["SHOWER_INITIATOR_KTMAX"] =
      (data[pid]["SHOWER_INITIATOR_KTMAX"]
       .SetDefault(defaults->params["SHOWER_INITIATOR_KTMAX"])
       .Get<double>());
    actuals->params["SHOWER_INITIATOR_KTEXPO"] =
      (data[pid]["SHOWER_INITIATOR_KTEXPO"]
       .SetDefault(defaults->params["SHOWER_INITIATOR_KTEXPO"])
       .Get<double>());
    actuals->params["REFERENCE_ENERGY"] =
      (data[pid]["REFERENCE_ENERGY"]
       .SetDefault(defaults->params["REFERENCE_ENERGY"])
       .Get<double>());
    actuals->params["ENERGY_SCALING_EXPO"] =
      (data[pid]["ENERGY_SCALING_EXPO"]
       .SetDefault(defaults->params["ENERGY_SCALING_EXPO"])
       .Get<double>());
    actuals->params["BEAM_SPECTATOR_MEAN"] =
      (data[pid]["BEAM_SPECTATOR_MEAN"]
       .SetDefault(defaults->params["BEAM_SPECTATOR_MEAN"])
       .Get<double>());
    actuals->params["BEAM_SPECTATOR_SIGMA"] =
      (data[pid]["BEAM_SPECTATOR_SIGMA"]
       .SetDefault(defaults->params["BEAM_SPECTATOR_SIGMA"])
       .Get<double>());
    actuals->params["BEAM_SPECTATOR_Q2"] =
      (data[pid]["BEAM_SPECTATOR_Q2"]
       .SetDefault(defaults->params["BEAM_SPECTATOR_Q2"])
       .Get<double>());
    actuals->params["BEAM_SPECTATOR_KTMAX"] =
      (data[pid]["BEAM_SPECTATOR_KTMAX"]
       .SetDefault(defaults->params["BEAM_SPECTATOR_KTMAX"])
       .Get<double>());
    actuals->params["BEAM_SPECTATOR_KTEXPO"] =
      (data[pid]["BEAM_SPECTATOR_KTEXPO"]
       .SetDefault(defaults->params["BEAM_SPECTATOR_KTEXPO"])
       .Get<double>());
    ////////////////////////////////////////////////////////////////////////////////////
    // Fix the matter distribution: form and parameters
    ///////////////////////////////////////////////////////////////////////////////////
    actuals->m_form =
      (data[pid]["MATTER_FORM"]
       .SetDefault(defaults->m_form)
       .Get<matter_form>());
    actuals->params["RADIUS_1"] =
      (data[pid]["RADIUS_1"]
       .SetDefault(defaults->params["RADIUS_1"])
       .Get<double>());
    actuals->params["RADIUS_2"] =
      (data[pid]["RADIUS_2"]
       .SetDefault(defaults->params["RADIUS_2"])
       .Get<double>());
    actuals->params["FRACTION_1"] =
      (data[pid]["FRACTION_1"]
       .SetDefault(defaults->params["FRACTION_1"])
       .Get<double>());
    m_actuals[flav] = actuals;
    msg_Out() << "Reading in parameters for " << flav << " yields:\n"
              << (*m_actuals[flav]) << "\n";
  }
  rempars->Output();
}

double Remnants_Parameters::Get(const ATOOLS::Flavour& flav,
                                std::string            keyword)
{
  if (m_actuals.find(flav) != m_actuals.end() &&
      m_actuals[flav]->params.find(keyword) != m_actuals[flav]->params.end())
    return m_actuals[flav]->params[keyword];
  else if (m_defaults.find(flav) != m_defaults.end() &&
           m_defaults[flav]->params.find(keyword) !=
                   m_defaults[flav]->params.end())
    return m_defaults[flav]->params[keyword];
  THROW(fatal_error, "Keyword " + keyword + " not found for " + flav.IDName());
  return 0.;
}

primkT_form Remnants_Parameters::KT_Form(const ATOOLS::Flavour& flav)
{
  if (m_actuals.find(flav) != m_actuals.end()) return m_actuals[flav]->kT_form;
  else if (m_defaults.find(flav) != m_defaults.end())
    return m_defaults[flav]->kT_form;
  THROW(fatal_error, "Keyword kT_form not found for " + flav.IDName());
  return primkT_form::undefined;
}

primkT_recoil Remnants_Parameters::KT_Recoil(const ATOOLS::Flavour& flav)
{
  if (m_actuals.find(flav) != m_actuals.end())
    return m_actuals[flav]->kT_recoil;
  else if (m_defaults.find(flav) != m_defaults.end())
    return m_defaults[flav]->kT_recoil;
  THROW(fatal_error, "Keyword kT_recoil not found for " + flav.IDName());
  return primkT_recoil::undefined;
}

matter_form Remnants_Parameters::Matter_Form(const ATOOLS::Flavour& flav)
{
  if (m_actuals.find(flav) != m_actuals.end()) return m_actuals[flav]->m_form;
  else if (m_defaults.find(flav) != m_defaults.end())
    return m_defaults[flav]->m_form;
  THROW(fatal_error, "Keyword matter_form not found for " + flav.IDName());
  return matter_form::unknown;
}

void Remnants_Parameters::Output()
{
  if (!msg->LevelIsInfo()) return;
  msg_Out() << "==============================================================="
               "========\n";
  for (std::map<Flavour, remnant_parameters*>::iterator flrpit =
               m_defaults.begin();
       flrpit != m_defaults.end(); flrpit++) {
    bool act = (m_actuals.find(flrpit->first) != m_actuals.end());
    msg_Out() << "-------------------------------------------------------------"
                 "----------\n"
              << "Remnant default (actuals) for " << flrpit->first << ":\n";
    msg_Out() << "   Primordial KT Form:   " << flrpit->second->kT_form;
    if (act) msg_Out() << " (" << m_actuals[flrpit->first]->kT_form << ")";
    msg_Out() << "\n";
    msg_Out() << "   Primordial KT Recoil: " << flrpit->second->kT_recoil;
    if (act) msg_Out() << " (" << m_actuals[flrpit->first]->kT_recoil << ")";
    msg_Out() << "\n";
    msg_Out() << "   Matter Form:          " << flrpit->second->m_form;
    if (act) msg_Out() << " (" << m_actuals[flrpit->first]->m_form << ")";
    msg_Out() << "\n";
    for (std::map<std::string, double>::iterator pit =
                 m_defaults[flrpit->first]->params.begin();
         pit != m_defaults[flrpit->first]->params.end(); pit++) {
      msg_Out() << "   " << pit->first << ": " << pit->second;
      if (act)
        msg_Out() << " (" << m_actuals[flrpit->first]->params[pit->first]
                  << ")";
      msg_Out() << "\n";
    }
    msg_Out() << "-------------------------------------------------------------"
                 "----------\n";
  }
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
    case matter_form::single_gaussian: return os << "Single_Gaussian";
    case matter_form::double_gaussian: return os << "Double_Gaussian";
    case matter_form::unknown: return os << "Unknown";
  }
  return os;
}

std::ostream& REMNANTS::operator<<(std::ostream&             os,
                                   const remnant_parameters& parms)
{
  os << "   Primordial k_T Form   = " << parms.kT_form << "\n"
     << "   Primordial k_T Recoil = " << parms.kT_recoil << "\n"
     << "   Matter Form           = " << parms.m_form << "\n";
  for (std::map<std::string, double>::const_iterator pit = parms.params.begin();
       pit != parms.params.end(); pit++)
    os << "   " << pit->first << " = " << pit->second << "\n";
  return os;
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
  if (tag == "Single_Gaussian") f = matter_form::single_gaussian;
  else if (tag == "Double_Gaussian")
    f = matter_form::double_gaussian;
  else
    THROW(fatal_error, "Unknown matter form \"" + tag + "\"");
  return is;
}
