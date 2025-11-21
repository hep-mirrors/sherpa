#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon_Kelly.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon_Dipole.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include <cmath>

using namespace ATOOLS;

FormFactor_EMnucleon::FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon) 
      : m_boson_type(boson), m_nucleon_type(nucleon), p_model_impl(nullptr)
{
  msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Constructor called"<<std::endl;
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_formfactor_model = s["Model"].SetDefault(ffmodel::off).Get<ffmodel::code>();
  
  // Global params
  m_gA = s["gA"].SetDefault(1.267).Get<double>();
  m_sin2thetaW = s["sin2thetaW"].SetDefault(0.231).Get<double>();
  
  msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Using model: "<<m_formfactor_model<<std::endl;
  
  // Create model-specific implementation (ONE check, then pointer handles everything!)
  switch (m_formfactor_model)
  {
    case ffmodel::off:
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Off model - point-like form factors"<<std::endl;
      RegisterDefaultsOff();
      p_model_impl = nullptr; // No model implementation, use fallback values given in this file
      break;
    case ffmodel::kelly:
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Kelly model - creating Kelly implementation"<<std::endl;
      p_model_impl = std::make_unique<FormFactor_EMnucleon_Kelly>(boson, nucleon, m_gA, m_sin2thetaW);
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Kelly model implementation created"<<std::endl;
      break;
    case ffmodel::dipole:
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Dipole model - creating Dipole implementation"<<std::endl;
      p_model_impl = std::make_unique<FormFactor_EMnucleon_Dipole>(boson, nucleon, m_gA, m_sin2thetaW);
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Dipole model implementation created"<<std::endl;
      break;
    default:
      THROW(fatal_error, "Unknown form factor model: " + ToString(m_formfactor_model));
  }
  
  msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Constructor complete"<<std::endl;
}

FormFactor_EMnucleon::~FormFactor_EMnucleon() {}

void FormFactor_EMnucleon::RegisterDefaultsOff() {
  // Set form factors to point-like values
  m_F1p = 1.0;
  m_F1n = 0.0;
  m_F1W = 1.0;
  m_F1Zp = 0.5 - m_sin2thetaW;
  m_F1Zn = -0.5;
  m_FAW = 1.0;
  m_FAZp = 0.5 * m_gA;
  m_FAZn = -0.5 * m_gA;

  m_F2p = 0.0;
  m_F2n = 0.0;
  m_F2W = 0.0;
  m_F2Zp = 0.0;
  m_F2Zn = 0.0;
  m_FPW = 0.0;
  m_FPZp = 0.0;
  m_FPZn = 0.0;
}

double FormFactor_EMnucleon::F1p(const double &q2) {
  if (p_model_impl) return p_model_impl->F1p(q2);  // call model-specific implementation
  return m_F1p; // off model fallback
}

double FormFactor_EMnucleon::F1n(const double &q2) {
  if (p_model_impl) return p_model_impl->F1n(q2);  
  return m_F1n; 
}

double FormFactor_EMnucleon::F2p(const double &q2) {
  if (p_model_impl) return p_model_impl->F2p(q2);  
  return m_F2p; 
}

double FormFactor_EMnucleon::F2n(const double &q2) {
  if (p_model_impl) return p_model_impl->F2n(q2);  
  return m_F2n; 
}

double FormFactor_EMnucleon::F1W(const double &q2) {
  if (p_model_impl) return p_model_impl->F1W(q2);
  return m_F1W;
}

double FormFactor_EMnucleon::F2W(const double &q2) {
  if (p_model_impl) return p_model_impl->F2W(q2);
  return m_F2W;
}

double FormFactor_EMnucleon::FAW(const double &q2) {
  if (p_model_impl) return p_model_impl->FAW(q2);  
  return m_FAW; 
}

double FormFactor_EMnucleon::FPW(const double &q2) {
  if (p_model_impl) return p_model_impl->FPW(q2);  
  return m_FPW; 
}

double FormFactor_EMnucleon::F1Zp(const double &q2) {
  if (p_model_impl) return p_model_impl->F1Zp(q2);  
  return m_F1Zp; 
}

double FormFactor_EMnucleon::F2Zp(const double &q2) {
  if (p_model_impl) return p_model_impl->F2Zp(q2);  
  return m_F2Zp; 
}

double FormFactor_EMnucleon::FAZp(const double &q2) {
  if (p_model_impl) return p_model_impl->FAZp(q2);  
  return m_FAZp; 
}

double FormFactor_EMnucleon::FPZp(const double &q2) {
  if (p_model_impl) return p_model_impl->FPZp(q2);  
  return m_FPZp; 
}

double FormFactor_EMnucleon::F1Zn(const double &q2) {
  if (p_model_impl) return p_model_impl->F1Zn(q2);  
  return m_F1Zn; 
} 

double FormFactor_EMnucleon::F2Zn(const double &q2) {
  if (p_model_impl) return p_model_impl->F2Zn(q2);  
  return m_F2Zn; 
}

double FormFactor_EMnucleon::FAZn(const double &q2) {
  if (p_model_impl) return p_model_impl->FAZn(q2);  
  return m_FAZn; 
}

double FormFactor_EMnucleon::FPZn(const double &q2) {
  if (p_model_impl) return p_model_impl->FPZn(q2);  
  return m_FPZn; 
}

NucleonFormFactors FormFactor_EMnucleon::Photon_Proton(const double &q2) {
  return NucleonFormFactors(F1p(q2), F2p(q2), 0.0, 0.0);
}

NucleonFormFactors FormFactor_EMnucleon::Photon_Neutron(const double &q2) {
  return NucleonFormFactors(F1n(q2), F2n(q2), 0.0, 0.0);
}

NucleonFormFactors FormFactor_EMnucleon::W_Boson(const double &q2) {
  return NucleonFormFactors(F1W(q2), F2W(q2), FAW(q2), FPW(q2));
}

NucleonFormFactors FormFactor_EMnucleon::Z_Proton(const double &q2) {
  return NucleonFormFactors(F1Zp(q2), F2Zp(q2), FAZp(q2), FPZp(q2));
}

NucleonFormFactors FormFactor_EMnucleon::Z_Neutron(const double &q2) {
  return NucleonFormFactors(F1Zn(q2), F2Zn(q2), FAZn(q2), FPZn(q2));
}

NucleonFormFactors FormFactor_EMnucleon::GetFormFactors(const double &q2) {
  switch (m_boson_type) {
    case incomingboson::photon:
      if (m_nucleon_type == incomingnucleon::proton)
        return Photon_Proton(q2);
      else if (m_nucleon_type == incomingnucleon::neutron)
        return Photon_Neutron(q2);
      break;
    case incomingboson::W:
      return W_Boson(q2);
    case incomingboson::Z:
      if (m_nucleon_type == incomingnucleon::proton)
        return Z_Proton(q2);
      else if (m_nucleon_type == incomingnucleon::neutron)
        return Z_Neutron(q2);
      break;
    default:
      THROW(fatal_error, "Unknown boson type in GetFormFactors");
  }
  // Fallback (shouldn't reach here)
  return NucleonFormFactors(0.0, 0.0, 0.0, 0.0);
}

// Stream operators for BosonNucleonType 
std::ostream &ATOOLS::operator<<(std::ostream &str, const BosonNucleonType &type)
{
  if (type.boson == incomingboson::photon) {
    if (type.nucleon == incomingnucleon::proton)
      return str << "Photon-Proton";
    else if (type.nucleon == incomingnucleon::neutron)
      return str << "Photon-Neutron";
  } else if (type.boson == incomingboson::W) {
    return str << "W-Boson";
  } else if (type.boson == incomingboson::Z) {
    if (type.nucleon == incomingnucleon::proton)
      return str << "Z-Proton";
    else if (type.nucleon == incomingnucleon::neutron)
      return str << "Z-Neutron";
  }
  return str << "Unknown FormFactor process";
}

std::istream &ATOOLS::operator>>(std::istream &str, BosonNucleonType &type)
{
  std::string tag;
  str >> tag;
  if (tag.find("Photon-Proton") != std::string::npos) {
    type.boson = incomingboson::photon;
    type.nucleon = incomingnucleon::proton;
  }
  else if (tag.find("Photon-Neutron") != std::string::npos) {
    type.boson = incomingboson::photon;
    type.nucleon = incomingnucleon::neutron;
  }
  else if (tag.find("W-Boson") != std::string::npos) {
    type.boson = incomingboson::W;
    type.nucleon = incomingnucleon::off; 
  }
  else if (tag.find("Z-Proton") != std::string::npos) {
    type.boson = incomingboson::Z;
    type.nucleon = incomingnucleon::proton;
  }
  else if (tag.find("Z-Neutron") != std::string::npos) {
    type.boson = incomingboson::Z;
    type.nucleon = incomingnucleon::neutron;
  }
  else
    THROW(fatal_error, "Unknown Form_Factor: Mode ");
  return str;
}

std::ostream &ATOOLS::operator<<(std::ostream &str,
                                 const ffmodel::code &ff)
{
  if (ff == ffmodel::kelly)
    return str << "Kelly";
  else if (ff == ffmodel::dipole)
    return str << "Dipole";
  else if (ff == ffmodel::off)
    return str << "Off";
  return str << "unknown";
}

std::istream &ATOOLS::operator>>(std::istream &str, ffmodel::code &mode)
{
  std::string tag;
  str >> tag;
  if (tag.find("Kelly") != std::string::npos)
    mode = ffmodel::kelly;
  else if (tag.find("Dipole") != std::string::npos)
    mode = ffmodel::dipole;
  else if (tag.find("Off") != std::string::npos)
    mode = ffmodel::off;
  else if (tag.find("None") != std::string::npos)
    mode = ffmodel::off;
  else
    THROW(fatal_error, "Unknown Form_Factor: Mode ");
  return str;
}
