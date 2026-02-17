#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon.H"
#include "ATOOLS/Phys/FormFactor_EMnucleon_Dipole.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include <cmath>

using namespace ATOOLS;

// Very simple experimental independent dipole form factor model based on https://arxiv.org/abs/nucl-ex/0307001v1 

FormFactor_EMnucleon_Dipole::FormFactor_EMnucleon_Dipole(incomingboson::code boson, incomingnucleon::code nucleon, double gA, double sin2thetaW) 
      : m_boson_type(boson), m_nucleon_type(nucleon), m_gA(gA), m_sin2thetaW(sin2thetaW), m_formfactor_model(ffmodel::dipole)
{
  msg_Out()<<"FormFactor_EMnucleon_Dipole::FormFactor_EMnucleon_Dipole(): Constructor called"<<std::endl;    

  // Load model-specific parameters based on m_formfactor_model
  RegisterDefaultsDipole();
  if (boson == incomingboson::W || boson == incomingboson::Z)
  {
    RegisterDefaultsAxial();
  }

  msg_Out()<<"FormFactor_EMnucleon_Dipole::FormFactor_EMnucleon_Dipole(): Constructor complete"<<std::endl;
}

FormFactor_EMnucleon_Dipole::~FormFactor_EMnucleon_Dipole() {}

void FormFactor_EMnucleon_Dipole::RegisterDefaultsDipole() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_massp = s["mass_p"].SetDefault(Flavour(kf_p_plus).Mass()).Get<double>(); //0.9382720813
  m_mup = s["mu_p"].SetDefault(2.79284734463).Get<double>();
  m_massn = s["mass_n"].SetDefault(Flavour(kf_n).Mass()).Get<double>(); //0.9395654133
  m_mun = s["mu_n"].SetDefault(-1.9130427).Get<double>();
  // Dipole parametrisation parameters
  m_Lambda2 = s["Lambda2"].SetDefault(0.71).Get<double>();
  m_norm = s["norm"].SetDefault(1.0).Get<double>();
}

void FormFactor_EMnucleon_Dipole::RegisterDefaultsAxial() {
  // adjustable parameters for axial form factors
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_massA = s["mass_A"].SetDefault(1.0).Get<double>();
  m_fA = s["f_A"].SetDefault(0.0).Get<double>();
  m_masspi = s["mass_pi"].SetDefault(Flavour(kf_pi_plus).Mass()).Get<double>(); //0.13957018
}

double FormFactor_EMnucleon_Dipole::Q2_check(const double &Q2) {
  // ensure Q2 is positive
  if (Q2 < 0) {
    msg_Out() << "FormFactor_EMnucleon_Dipole::Q2_check(): Q2 is negative (" << Q2 << "). Using its absolute value." << std::endl;
    return -Q2;
  }
  else {
    return Q2;
  }
}

double FormFactor_EMnucleon_Dipole::tau_eval(const double &Q2, const double &mass) {
  return Q2 / (4. * mass * mass);
}

double FormFactor_EMnucleon_Dipole::Dipole_func(const double &Q2, const double &Lambda2, const double &norm) {
  double Q2_checked = Q2_check(Q2); // ensure Q2 is positive
  double den = pow(1. + Q2_checked / Lambda2, 2);
  double num = norm;
  return num / den;
}

double FormFactor_EMnucleon_Dipole::F1p(const double &Q2) {
  double tau = tau_eval(Q2, m_massp);
  double GE = Dipole_func(Q2, m_Lambda2, m_norm);
  double GM = m_mup * Dipole_func(Q2, m_Lambda2, m_norm); //+ve mag mom
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon_Dipole::F1n(const double &Q2) {
  double tau = tau_eval(Q2, m_massn);
  double GE = Dipole_func(Q2, m_Lambda2, m_norm);
  double GM = m_mun * Dipole_func(Q2, m_Lambda2, m_norm); // -ve mag mom
  return (GE + tau * GM) / (1. + tau);
}

double FormFactor_EMnucleon_Dipole::F2p(const double &Q2) {
  double tau = tau_eval(Q2, m_massp);
  double GE = Dipole_func(Q2, m_Lambda2, m_norm);
  double GM = m_mup * Dipole_func(Q2, m_Lambda2, m_norm); 
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon_Dipole::F2n(const double &Q2) {
  double tau = tau_eval(Q2, m_massn);
  double GE = Dipole_func(Q2, m_Lambda2, m_norm);
  double GM = m_mun * Dipole_func(Q2, m_Lambda2, m_norm);
  return (GM - GE) / (1. + tau);
}

double FormFactor_EMnucleon_Dipole::F1W(const double &Q2) {
  return F1p(Q2) - F1n(Q2);
}

double FormFactor_EMnucleon_Dipole::F2W(const double &Q2) {
  return F2p(Q2) - F2n(Q2);
}

double FormFactor_EMnucleon_Dipole::FAW(const double &Q2) {
  return Dipole_func(Q2, m_massA * m_massA, m_gA);
}

double FormFactor_EMnucleon_Dipole::FPW(const double &Q2) {
  // PCAC relation and pion-pole dominance
  double FA = FAW(Q2);
  return (4 * m_massp * m_massp * FA) / (Q2 + m_masspi * m_masspi);
}

double FormFactor_EMnucleon_Dipole::F1Zp(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1p(Q2) - 0.5 * F1n(Q2);
}

double FormFactor_EMnucleon_Dipole::F2Zp(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2p(Q2) - 0.5 * F2n(Q2);
}

double FormFactor_EMnucleon_Dipole::FAZp(const double &Q2) {
  return 0.5 * FAW(Q2);
}

double FormFactor_EMnucleon_Dipole::FPZp(const double &Q2) {
  return 0.5 * FPW(Q2);
}

double FormFactor_EMnucleon_Dipole::F1Zn(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F1n(Q2) - 0.5 * F1p(Q2);
} 

double FormFactor_EMnucleon_Dipole::F2Zn(const double &Q2) {
  return (0.5 - 2 * m_sin2thetaW) * F2n(Q2) - 0.5 * F2p(Q2);
}

double FormFactor_EMnucleon_Dipole::FAZn(const double &Q2) {
  return -0.5 * FAW(Q2);
}

double FormFactor_EMnucleon_Dipole::FPZn(const double &Q2) {
  return -0.5 * FPW(Q2);
}
