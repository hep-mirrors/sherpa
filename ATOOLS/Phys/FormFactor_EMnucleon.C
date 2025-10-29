#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace ATOOLS;

FormFactor_EMnucleon::FormFactor_EMnucleon(incomingboson::code boson, incomingnucleon::code nucleon) 
      : m_boson_type(boson), m_nucleon_type(nucleon)
{
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  if (boson == incomingboson::photon) 
  {
    if (nucleon == incomingnucleon::proton) 
    {
      RegisterDefaultsProton();
    } 
    else if (nucleon == incomingnucleon::neutron) 
    {
      RegisterDefaultsNeutron();
    }
  }
  else if (boson == incomingboson::W) 
  {
    RegisterDefaultsProton();
    RegisterDefaultsNeutron();
    RegisterDefaultsAxial();
  }
  else if (boson == incomingboson::Z)
  {
    RegisterDefaultsProton();
    RegisterDefaultsNeutron();
    RegisterDefaultsAxial();
    RegisterDefaultsZ();
  }
}

FormFactor_EMnucleon::~FormFactor_EMnucleon() {}

void FormFactor_EMnucleon::RegisterDefaultsProton() {
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
}

void FormFactor_EMnucleon::RegisterDefaultsNeutron() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
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
  m_gA = s["g_A"].SetDefault(1.2695).Get<double>();
  m_fA = s["f_A"].SetDefault(0.0).Get<double>();
  m_masspi = s["mass_pi"].SetDefault(0.13957018).Get<double>();
}

void FormFactor_EMnucleon::RegisterDefaultsZ() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_sin2thetaW = s["sin2thetaW"].SetDefault(0.231).Get<double>();
}

double FormFactor_EMnucleon::Q2_eval(const double &q2) {
  return sqr(q2*q2); // ensure Q2 is positive
}

double FormFactor_EMnucleon::tau_eval(const double &q2, const double &mass) {
  double Q2 = Q2_eval(q2);
  return Q2 / (4. * mass * mass);
}

double FormFactor_EMnucleon::GEp(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massp);
  double num = m_a0pE + m_a1pE * tau;
  double den = 1. + m_b1pE * tau + m_b2pE * tau * tau + m_b3pE * tau * tau * tau;
  return num / den;
}

double FormFactor_EMnucleon::GMp(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massp);
  double num = m_a0pM + m_a1pM * tau;
  double den = 1. + m_b1pM * tau + m_b2pM * tau * tau + m_b3pM * tau * tau * tau;
  return m_mup * num / den;
}

double FormFactor_EMnucleon::GEn(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massn);
  double num = m_a0nE + m_a1nE * tau;
  double den = 1. + m_b1nE * tau + m_b2nE * tau * tau + m_b3nE * tau * tau * tau;
  return num / den;
}

double FormFactor_EMnucleon::GMn(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massn);
  double num = m_a0nM + m_a1nM * tau;
  double den = 1. + m_b1nM * tau + m_b2nM * tau * tau + m_b3nM * tau * tau * tau;
  return m_mun * num / den;
}

double FormFactor_EMnucleon::F1p(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massp);
  double GE = GEp(Q2);
  double GM = GMp(Q2);
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon::F1n(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massn);
  double GE = GEn(Q2);
  double GM = GMn(Q2);
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon::F2p(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massp);
  double GE = GEp(Q2);
  double GM = GMp(Q2);
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon::F2n(const double &q2) {
  double Q2 = Q2_eval(q2);
  double tau = tau_eval(Q2, m_massn);
  double GE = GEn(Q2);
  double GM = GMn(Q2);
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon::F1W(const double &q2) {
  return F1p(q2) - F1n(q2);
}

double FormFactor_EMnucleon::F2W(const double &q2) {
  return F2p(q2) - F2n(q2);
}

double FormFactor_EMnucleon::FAW(const double &q2) {
  double Q2 = Q2_eval(q2);
  return m_gA / pow(1. + Q2 / (m_massA * m_massA), 2);
}

double FormFactor_EMnucleon::FPW(const double &q2) {
  // PCAC relation
  double Q2 = Q2_eval(q2);
  return (2. * m_massp * m_gA) / (Q2 + m_masspi * m_masspi) * pow(1. + Q2 / (m_massA * m_massA), -2);
}

double FormFactor_EMnucleon::F1Zp(const double &q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1p(q2) - 0.5 * F1n(q2);
}

double FormFactor_EMnucleon::F2Zp(const double &q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2p(q2) - 0.5 * F2n(q2);
}

double FormFactor_EMnucleon::FAZp(const double &q2) {
  return 0.5 * FAW(q2);
}

double FormFactor_EMnucleon::FPZp(const double &q2) {
  // TODO: check 
  return 0.5 * FPW(q2);
}

double FormFactor_EMnucleon::F1Zn(const double &q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1n(q2) - 0.5 * F1p(q2);
} 

double FormFactor_EMnucleon::F2Zn(const double &q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2n(q2) - 0.5 * F2p(q2);
}

double FormFactor_EMnucleon::FAZn(const double &q2) {
  return -0.5 * FAW(q2);
}

double FormFactor_EMnucleon::FPZn(const double &q2) {
  // TODO: check
  return -0.5 * FPW(q2);
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