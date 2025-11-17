#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include <cmath>

using namespace ATOOLS;

FormFactor_EMnucleon::FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon) 
      : m_boson_type(boson), m_nucleon_type(nucleon)
{
  msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Constructor called"<<std::endl;
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_formfactor_model = s["Model"].SetDefault(ffmodel::kelly).Get<ffmodel::code>();
  
  // Global params
  m_gA = s["gA"].SetDefault(1.267).Get<double>();
  m_sin2thetaW = s["sin2thetaW"].SetDefault(0.231).Get<double>();
  
  msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Using model: "<<m_formfactor_model<<std::endl;
  
  // Load model-specific parameters based on m_formfactor_model
  switch (m_formfactor_model)
  {
    case ffmodel::off:
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Off model - point-like form factors"<<std::endl;
      RegisterDefaultsOff(); 
      break;
    case ffmodel::kelly:
      msg_Out()<<"FormFactor_EMnucleon::FormFactor_EMnucleon(): Kelly model - loading parameters"<<std::endl;
      RegisterDefaultsKelly();
      if (boson == incomingboson::W || boson == incomingboson::Z)
      {
        msg_Out() << "FormFactor_EMnucleon::FormFactor_EMnucleon(): Loading axial parameters" << std::endl;
        RegisterDefaultsAxial();
      }
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

  m_F2p, m_F2n, m_F2W, m_F2Zp, m_F2Zn = 0.0;
  m_FPW, m_FPZp, m_FPZn = 0.0;
}

void FormFactor_EMnucleon::RegisterDefaultsKelly() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_massp = s["mass_p"].SetDefault(0.938272081).Get<double>();
  m_mup = s["mu_p"].SetDefault(2.79284734463).Get<double>();
  // Kelly parametrisation for proton EM form factors
  m_a0pE = s["a0"].SetDefault(1.).Get<double>();
  m_a1pE = s["a1"].SetDefault(-0.24).Get<double>();
  m_b1pE = s["b1"].SetDefault(10.98).Get<double>();
  m_b2pE = s["b2"].SetDefault(12.82).Get<double>();
  m_b3pE = s["b3"].SetDefault(21.97).Get<double>();
  m_a0pM = s["a0_m"].SetDefault(1.).Get<double>();
  m_a1pM = s["a1_m"].SetDefault(0.12).Get<double>();
  m_b1pM = s["b1_m"].SetDefault(10.97).Get<double>();
  m_b2pM = s["b2_m"].SetDefault(18.86).Get<double>();
  m_b3pM = s["b3_m"].SetDefault(6.55).Get<double>();
  m_massn = s["mass_n"].SetDefault(0.9395654133).Get<double>();
  m_mun = s["mu_n"].SetDefault(-1.9130427).Get<double>();
  // Galster parametrisation for neutron EM form factors
  m_An = s["A_n"].SetDefault(1.7).Get<double>();
  m_Bn = s["B_n"].SetDefault(3.3).Get<double>();
  m_Deltasq = s["Delta_sq"].SetDefault(0.71).Get<double>();
  // Kelly parametrisation for neutron EM form factors
  m_a0nE = s["a0_n"].SetDefault(0.).Get<double>();
  m_a1nE = s["a1_n"].SetDefault(m_An).Get<double>();
  m_b1nE = s["b1_n"].SetDefault(m_Bn + (8 * m_massn * m_massn) / m_Deltasq).Get<double>();
  m_b2nE = s["b2_n"].SetDefault((16 * m_massn * m_massn * m_massn * m_massn) / (m_Deltasq * m_Deltasq) + m_Bn * (8 * m_massn * m_massn) / m_Deltasq).Get<double>();
  m_b3nE = s["b3_n"].SetDefault(m_Bn * (16 * m_massn * m_massn * m_massn * m_massn) / (m_Deltasq * m_Deltasq)).Get<double>();
  m_a0nM = s["a0_nm"].SetDefault(1.).Get<double>();
  m_a1nM = s["a1_nm"].SetDefault(2.33).Get<double>();
  m_b1nM = s["b1_nm"].SetDefault(14.72).Get<double>();
  m_b2nM = s["b2_nm"].SetDefault(24.2).Get<double>();
  m_b3nM = s["b3_nm"].SetDefault(84.1).Get<double>();
}

void FormFactor_EMnucleon::RegisterDefaultsAxial() {
  // adjustable parameters for axial form factors
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_massA = s["mass_A"].SetDefault(1.0).Get<double>();
  m_fA = s["f_A"].SetDefault(0.0).Get<double>();
  m_masspi = s["mass_pi"].SetDefault(0.13957018).Get<double>();
}

double FormFactor_EMnucleon::Q2_check(const double &Q2) {
  // ensure Q2 is positive
  if (Q2 < 0) {
    msg_Out() << "FormFactor_EMnucleon::Q2_check(): Q2 is negative (" << Q2 << "). Using its absolute value." << std::endl;
    return -Q2;
  }
  else {
    return Q2;
  }
}

double FormFactor_EMnucleon::tau_eval(const double &Q2, const double &mass) {
  return Q2 / (4. * mass * mass);
}

double FormFactor_EMnucleon::Kelly_func(const double &Q2, const double &a0, const double &a1, const double &b1, const double &b2, const double &b3) {
  double tau = tau_eval(Q2, m_massp); // assuming proton mass for Kelly function
  double num = a0 + a1 * tau;
  double den = 1. + b1 * tau + b2 * tau * tau + b3 * tau * tau * tau;
  return num / den;
}

double FormFactor_EMnucleon::F1p(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F1p;
    case ffmodel::kelly: {
      double tau = tau_eval(Q2, m_massp);
      double GE = Kelly_func(Q2, m_a0pE, m_a1pE, m_b1pE, m_b2pE, m_b3pE);
      double GM = m_mup * Kelly_func(Q2, m_a0pM, m_a1pM, m_b1pM, m_b2pM, m_b3pM);
      return (GE + tau * GM) / (1. + tau);
    }
    default:
      THROW(fatal_error, "Unknown form factor model in F1p");
  }
}

double FormFactor_EMnucleon::F1n(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F1n;
    case ffmodel::kelly: {
      double tau = tau_eval(Q2, m_massn);
      double GE = Kelly_func(Q2, m_a0nE, m_a1nE, m_b1nE, m_b2nE, m_b3nE);
      double GM = m_mun * Kelly_func(Q2, m_a0nM, m_a1nM, m_b1nM, m_b2nM, m_b3nM);
      return (GE + tau * GM) / (1. + tau);
    }
    default:
      THROW(fatal_error, "Unknown form factor model in F1n");
  }
}

double FormFactor_EMnucleon::F2p(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F2p;
    case ffmodel::kelly: {
      double tau = tau_eval(Q2, m_massp);
      double GE = Kelly_func(Q2, m_a0pE, m_a1pE, m_b1pE, m_b2pE, m_b3pE);
      double GM = m_mup * Kelly_func(Q2, m_a0pM, m_a1pM, m_b1pM, m_b2pM, m_b3pM);
      return (GM - GE) / (1. + tau);
    }
    default:
      THROW(fatal_error, "Unknown form factor model in F2p");
  }
}

double FormFactor_EMnucleon::F2n(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F2n;
    case ffmodel::kelly: {
      double tau = tau_eval(Q2, m_massn);
      double GE = Kelly_func(Q2, m_a0nE, m_a1nE, m_b1nE, m_b2nE, m_b3nE);
      double GM = m_mun * Kelly_func(Q2, m_a0nM, m_a1nM, m_b1nM, m_b2nM, m_b3nM);
      return (GM - GE) / (1. + tau);
    }
    default:
      THROW(fatal_error, "Unknown form factor model in F2n");
  }
}

double FormFactor_EMnucleon::F1W(const double &Q2) {
  return F1p(Q2) - F1n(Q2);
}

double FormFactor_EMnucleon::F2W(const double &Q2) {
  return F2p(Q2) - F2n(Q2);
}

double FormFactor_EMnucleon::FAW(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FAW;
    case ffmodel::kelly:
      return m_gA / pow(1. + Q2 / (m_massA * m_massA), 2);
    default:
      THROW(fatal_error, "Unknown form factor model in FAW");
  }
}

double FormFactor_EMnucleon::FPW(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FPW;
    case ffmodel::kelly: {
      // PCAC relation and pion-pole dominance
      double FA = FAW(Q2);
      return (4 * m_massp * m_massp * FA) / (Q2 + m_masspi * m_masspi);
    }
    default:
      THROW(fatal_error, "Unknown form factor model in FPW");
  }
}

double FormFactor_EMnucleon::F1Zp(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F1Zp;
    case ffmodel::kelly:
      return (0.5 - 2 * m_sin2thetaW) * F1p(Q2) - 0.5 * F1n(Q2);
  }
}

double FormFactor_EMnucleon::F2Zp(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F2Zp;
    case ffmodel::kelly:
      return (0.5 - 2 * m_sin2thetaW) * F2p(Q2) - 0.5 * F2n(Q2);
  }
}

double FormFactor_EMnucleon::FAZp(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FAZp;
    case ffmodel::kelly:
      return 0.5 * FAW(Q2);
  }
}

double FormFactor_EMnucleon::FPZp(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FPZp;
    case ffmodel::kelly:
      return 0.5 * FPW(Q2); //check
  }
}

double FormFactor_EMnucleon::F1Zn(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F1Zn;
    case ffmodel::kelly:
      return (0.5 - 2 * m_sin2thetaW) * F1n(Q2) - 0.5 * F1p(Q2);
  }
} 

double FormFactor_EMnucleon::F2Zn(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_F2Zn;
    case ffmodel::kelly:
      return (0.5 - 2 * m_sin2thetaW) * F2n(Q2) - 0.5 * F2p(Q2);
  }
}

double FormFactor_EMnucleon::FAZn(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FAZn;
    case ffmodel::kelly:
      return -0.5 * FAW(Q2);
  }
}

double FormFactor_EMnucleon::FPZn(const double &Q2) {
  switch (m_formfactor_model) {
    case ffmodel::off:
      return m_FPZn;
    case ffmodel::kelly:
      return -0.5 * FPW(Q2); //check
  }
}

NucleonFormFactors FormFactor_EMnucleon::Photon_Proton(const double &Q2) {
  return NucleonFormFactors(F1p(Q2), F2p(Q2), 0.0, 0.0);
}

NucleonFormFactors FormFactor_EMnucleon::Photon_Neutron(const double &Q2) {
  return NucleonFormFactors(F1n(Q2), F2n(Q2), 0.0, 0.0);
}

NucleonFormFactors FormFactor_EMnucleon::W_Boson(const double &Q2) {
  return NucleonFormFactors(F1W(Q2), F2W(Q2), FAW(Q2), FPW(Q2));
}

NucleonFormFactors FormFactor_EMnucleon::Z_Proton(const double &Q2) {
  return NucleonFormFactors(F1Zp(Q2), F2Zp(Q2), FAZp(Q2), FPZp(Q2));
}

NucleonFormFactors FormFactor_EMnucleon::Z_Neutron(const double &Q2) {
  return NucleonFormFactors(F1Zn(Q2), F2Zn(Q2), FAZn(Q2), FPZn(Q2));
}

NucleonFormFactors FormFactor_EMnucleon::GetFormFactors(const double &Q2) {
  switch (m_boson_type) {
    case incomingboson::photon:
      if (m_nucleon_type == incomingnucleon::proton)
        return Photon_Proton(Q2);
      else if (m_nucleon_type == incomingnucleon::neutron)
        return Photon_Neutron(Q2);
      break;
    case incomingboson::W:
      return W_Boson(Q2);
    case incomingboson::Z:
      if (m_nucleon_type == incomingnucleon::proton)
        return Z_Proton(Q2);
      else if (m_nucleon_type == incomingnucleon::neutron)
        return Z_Neutron(Q2);
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
  else if (ff == ffmodel::off)
    return str << "Off";
  return str << "unknown";
}

std::istream &ATOOLS::operator>>(std::istream &str, ffmodel::code &mode)
{
  std::string tag;
  str >> tag;
  // mode=wgt::off;
  if (tag.find("Kelly") != std::string::npos)
    mode = ffmodel::kelly;
  else if (tag.find("Off") != std::string::npos)
    mode = ffmodel::off;
  else if (tag.find("None") != std::string::npos)
    mode = ffmodel::off;
  else
    THROW(fatal_error, "Unknown Form_Factor: Mode ");
  return str;
}
