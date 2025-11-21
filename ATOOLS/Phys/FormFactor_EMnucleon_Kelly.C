#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon_Kelly.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include <cmath>

using namespace ATOOLS;

FormFactor_EMnucleon_Kelly::FormFactor_EMnucleon_Kelly(incomingboson::code boson, incomingnucleon::code nucleon, double gA, double sin2thetaW) 
      : m_boson_type(boson), m_nucleon_type(nucleon), m_gA(gA), m_sin2thetaW(sin2thetaW), m_formfactor_model(ffmodel::kelly)
{
  msg_Out()<<"FormFactor_EMnucleon_Kelly::FormFactor_EMnucleon_Kelly(): Constructor called"<<std::endl;    

  // Load model-specific parameters based on m_formfactor_model
  RegisterDefaultsKelly();
  if (boson == incomingboson::W || boson == incomingboson::Z)
  {
    RegisterDefaultsAxial();
  }

  msg_Out()<<"FormFactor_EMnucleon_Kelly::FormFactor_EMnucleon_Kelly(): Constructor complete"<<std::endl;
}

FormFactor_EMnucleon_Kelly::~FormFactor_EMnucleon_Kelly() {}

void FormFactor_EMnucleon_Kelly::RegisterDefaultsKelly() {
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

void FormFactor_EMnucleon_Kelly::RegisterDefaultsAxial() {
  // adjustable parameters for axial form factors
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_massA = s["mass_A"].SetDefault(1.0).Get<double>();
  m_fA = s["f_A"].SetDefault(0.0).Get<double>();
  m_masspi = s["mass_pi"].SetDefault(0.13957018).Get<double>();
}

double FormFactor_EMnucleon_Kelly::Q2_check(const double &Q2) {
  // ensure Q2 is positive
  if (Q2 < 0) {
    msg_Out() << "FormFactor_EMnucleon_Kelly::Q2_check(): Q2 is negative (" << Q2 << "). Using its absolute value." << std::endl;
    return -Q2;
  }
  else {
    return Q2;
  }
}

double FormFactor_EMnucleon_Kelly::tau_eval(const double &Q2, const double &mass) {
  return Q2 / (4. * mass * mass);
}

double FormFactor_EMnucleon_Kelly::Kelly_func(const double &Q2, const double &a0, const double &a1, const double &b1, const double &b2, const double &b3) {
  double Q2_checked = Q2_check(Q2); // ensure Q2 is positive
  double tau = tau_eval(Q2_checked, m_massp); // assuming proton mass for Kelly function
  double num = a0 + a1 * tau;
  double den = 1. + b1 * tau + b2 * tau * tau + b3 * tau * tau * tau;
  return num / den;
}

double FormFactor_EMnucleon_Kelly::F1p(const double &Q2) {
  double tau = tau_eval(Q2, m_massp);
  double GE = Kelly_func(Q2, m_a0pE, m_a1pE, m_b1pE, m_b2pE, m_b3pE);
  double GM = m_mup * Kelly_func(Q2, m_a0pM, m_a1pM, m_b1pM, m_b2pM, m_b3pM);
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon_Kelly::F1n(const double &Q2) {
  double tau = tau_eval(Q2, m_massn);
  double GE = Kelly_func(Q2, m_a0nE, m_a1nE, m_b1nE, m_b2nE, m_b3nE);
  double GM = m_mun * Kelly_func(Q2, m_a0nM, m_a1nM, m_b1nM, m_b2nM, m_b3nM);
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon_Kelly::F2p(const double &Q2) {
  double tau = tau_eval(Q2, m_massp);
  double GE = Kelly_func(Q2, m_a0pE, m_a1pE, m_b1pE, m_b2pE, m_b3pE);
  double GM = m_mup * Kelly_func(Q2, m_a0pM, m_a1pM, m_b1pM, m_b2pM, m_b3pM);
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon_Kelly::F2n(const double &Q2) {
  double tau = tau_eval(Q2, m_massn);
  double GE = Kelly_func(Q2, m_a0nE, m_a1nE, m_b1nE, m_b2nE, m_b3nE);
  double GM = m_mun * Kelly_func(Q2, m_a0nM, m_a1nM, m_b1nM, m_b2nM, m_b3nM);
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon_Kelly::F1W(const double &Q2) {
  return F1p(Q2) - F1n(Q2);
}

double FormFactor_EMnucleon_Kelly::F2W(const double &Q2) {
  return F2p(Q2) - F2n(Q2);
}

double FormFactor_EMnucleon_Kelly::FAW(const double &Q2) {
  return m_gA / pow(1. + Q2 / (m_massA * m_massA), 2);
}

double FormFactor_EMnucleon_Kelly::FPW(const double &Q2) {
  // PCAC relation and pion-pole dominance
  double FA = FAW(Q2);
  return (4 * m_massp * m_massp * FA) / (Q2 + m_masspi * m_masspi);
}

double FormFactor_EMnucleon_Kelly::F1Zp(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1p(Q2) - 0.5 * F1n(Q2);
}

double FormFactor_EMnucleon_Kelly::F2Zp(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2p(Q2) - 0.5 * F2n(Q2);
}

double FormFactor_EMnucleon_Kelly::FAZp(const double &Q2) {
  return 0.5 * FAW(Q2);
}

double FormFactor_EMnucleon_Kelly::FPZp(const double &Q2) {
  return 0.5 * FPW(Q2);
}

double FormFactor_EMnucleon_Kelly::F1Zn(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1n(Q2) - 0.5 * F1p(Q2);
} 

double FormFactor_EMnucleon_Kelly::F2Zn(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2n(Q2) - 0.5 * F2p(Q2);
}

double FormFactor_EMnucleon_Kelly::FAZn(const double &Q2) {
  return -0.5 * FAW(Q2);
}

double FormFactor_EMnucleon_Kelly::FPZn(const double &Q2) {
  return -0.5 * FPW(Q2);
}


