#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/FormFactor.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace ATOOLS;

FormFactor::FormFactor() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_form_mode = s["Mode"].SetDefault(finalstate::off).Get<finalstate::code>();

  m_omega_mass = s["omega_mass"].SetDefault(0.78248).Get<double>();
  m_omega_g = s["omega_width"].SetDefault(0.00855).Get<double>();

  m_omegap_mass = s["omegap_mass"].SetDefault(1.410).Get<double>();
  m_omegap_g = s["omegap_width"].SetDefault(0.290).Get<double>();

  m_omegapp_mass = s["omegapp_mass"].SetDefault(1.67).Get<double>();
  m_omegapp_g = s["omegapp_width"].SetDefault(0.315).Get<double>();

  m_rho_mass = s["rho_mass"].SetDefault(0.77456).Get<double>();
  m_rho_g = s["rho_width"].SetDefault(0.14832).Get<double>();

  m_rhop_mass = s["rhop_mass"].SetDefault(1.4859).Get<double>();
  m_rhop_g = s["rhop_width"].SetDefault(0.37360).Get<double>();

  m_rhopp_mass = s["rhopp_mass"].SetDefault(1.8668).Get<double>();
  m_rhopp_g = s["rhopp_width"].SetDefault(0.30334).Get<double>();

  m_rhoppp_mass = s["rhoppp_mass"].SetDefault(2.2645).Get<double>();
  m_rhoppp_g = s["rhoppp_width"].SetDefault(0.11327).Get<double>();

  m_phi_mass = s["phi_mass"].SetDefault(1.01947).Get<double>();
  m_phi_g = s["phi_width"].SetDefault(0.00425).Get<double>();

  m_phip_mass = s["phip_mass"].SetDefault(1.680).Get<double>();
  m_phip_g = s["phip_width"].SetDefault(0.150).Get<double>();

  m_cphi = 0.00045;

  m_prho = 0;
  m_prhop = 3.7797;
  m_prhopp = 1.429;
  m_prhoppp = 0.921;

  m_pomega = 0.075;
  m_pphi = 2.888;

  m_I = Complex(0, 1);
  switch (m_form_mode) {
    case finalstate::pion:
      m_flv = Flavour(kf_pi_plus);
      RegisterDefaultsPion();
      break;
    case finalstate::kplus:
      m_flv = Flavour(kf_K_plus);
      RegisterDefaultsKaon();
      break;
    case finalstate::off:
      break;
    default:
      THROW(fatal_error, "Unknown final state");
  }
}

FormFactor::~FormFactor() {}

void FormFactor::RegisterDefaultsPion() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_comega = s["c_omega"].SetDefault(0.00158).Get<double>();
  m_crho = s["c_rho"].SetDefault(1.).Get<double>();
  m_crhop = s["c_rhop"].SetDefault(0.14104).Get<double>();
  m_crhopp = s["c_rhopp"].SetDefault(0.0614).Get<double>();
  m_crhoppp = s["c_rhoppp"].SetDefault(0.0047).Get<double>();
}

void FormFactor::RegisterDefaultsKaon() {
  Scoped_Settings s{Settings::GetMainSettings()["Form_Factor"]};
  m_comega = s["c_omega"].SetDefault(1.195).Get<double>();
  m_comegap = s["c_omegap"].SetDefault(-0.112).Get<double>();
  m_comegapp = s["c_omegapp"].SetDefault(-0.083).Get<double>();

  m_cphi = s["c_phi"].SetDefault(1.018).Get<double>();
  m_cphip = s["c_phip"].SetDefault(-0.018).Get<double>();

  m_crho = s["c_rho"].SetDefault(1.195).Get<double>();
  m_crhop = s["c_rhop"].SetDefault(-0.112).Get<double>();
  m_crhopp = s["c_rhopp"].SetDefault(-0.083).Get<double>();
}

void Initalize() {}

Complex FormFactor::ppi(const double &q2) {
  return 0.5 * csqrt(q2 - 4 * sqr(m_flv.Mass()));
}

Complex FormFactor::d(const double &m) {
  Complex L = log((m + 2. * ppi(m * m)) / 2. / m_flv.Mass());
  return 3. / M_PI * sqr(m_flv.Mass() / ppi(m * m)) * L +
         0.5 * m / ppi(m * m) / M_PI -
         m * sqr(m_flv.Mass()) / M_PI / pow(ppi(m * m), 3);
}

Complex FormFactor::h(const double &q2) {
  Complex L = log((sqrt(q2) + 2. * ppi(q2)) / 2. / m_flv.Mass());
  return 2. / M_PI * ppi(q2) / sqrt(q2) * L;
}

Complex FormFactor::dh(const double &q2) {
  return h(q2) / 8. * (1. / sqr(ppi(q2)) - 4. / q2) + 0.5 / M_PI / q2;
}

Complex FormFactor::f(const double &q2, const double &m, const double &g) {
  Complex t1 = sqr(ppi(q2)) * (h(q2) - h(m * m)) +
               sqr(ppi(m * m)) * (m * m - q2) * dh(m * m);
  return t1 * g * m * m / pow(ppi(m * m), 3);
}

Complex FormFactor::gamma(const double &q2, const double &m, const double &g) {
  return g * m / sqrt(q2) * pow(ppi(q2) / ppi(m * m), 3);
}

Complex FormFactor::b(const double &q2, const Flavour &fl) {
  double m = fl.Mass();
  double w = fl.Width();
  Complex I(0., 1.);

  return m * m / (m * m - q2 - m * w * I);
}

Complex FormFactor::BW_GS(const double &q2, const Flavour &fl) {
  //  Gounaris-Sakurai (GS) function https://inspirehep.net/literature/53152
  double mf = fl.Mass();
  double wf = fl.Width();
  Complex I(0, 1.);

  Complex num = mf * mf + d(mf) * mf * wf;
  Complex den = mf * mf - q2 + f(q2, mf, wf) - I * mf * gamma(q2, mf, wf);
  if (IsBad(num) || IsBad(den)) {
    msg_Error() << "NaN in " << METHOD << std::endl
                << "num = " << num << std::endl
                << "Flavour = " << fl << std::endl
                << "Final State = " << m_flv << std::endl
                << "Final State Mass = " << m_flv.Mass() << std::endl
                << "Final State Width = " << m_flv.Width() << std::endl
                << "den = " << num << std::endl
                << "mf = " << mf << std::endl
                << "wf = " << wf << std::endl
                << "d(mf) = " << d(mf) << std::endl
                << "f(q2,mf,wf) = " << f(q2, mf, wf) << std::endl
                << "gamma(q2,mf,wf) = " << gamma(q2, mf, wf) << std::endl;
  }
  return num / den;
}

Complex FormFactor::Pion(const double &q2) {
  Complex num1 =
      1. + (q2 / sqr(m_omega.Mass()) * m_comega * exp(m_I * m_pomega)) *
               b(q2, m_omega);
  num1 += q2 / (sqr(m_phi.Mass())) * m_cphi * exp(m_I * m_pphi) * b(q2, m_phi);
  num1 *= BW_GS(q2, m_rho);
  Complex den = 1. + m_crhop * exp(m_I * m_prhop) +
                m_crhopp * exp(m_I * m_prhopp) +
                m_crhoppp * exp(m_I * m_prhoppp);

  Complex num2 = BW_GS(q2, m_rhop) * m_crhop * exp(m_I * m_prhop);
  num2 += BW_GS(q2, m_rhopp) * m_crhopp * exp(m_I * m_prhopp);
  num2 += BW_GS(q2, m_rhoppp) * m_crhoppp * exp(m_I * m_prhoppp);
  if (IsBad(num1) || IsBad(num2) || IsBad(den)) {
    msg_Error() << "NaN in pion form-factor\n"
                << "Numertor 1 = " << num1 << std::endl
                << "Numertor 2 = " << num2 << std::endl
                << "Den 1 = " << den << std::endl;
  }
  Complex Form = (num1 + num2) / den;
  return Form;
}

Complex FormFactor::KPlus(const double &q2) {
  Complex Form =
      0.5 * (m_crho * BW_GS(q2, m_rho) + m_crhop * BW_GS(q2, m_rhop) +
             m_crhopp * BW_GS(q2, m_rhopp));
  Form += 1. / 6. *
          (m_comega * BW_GS(q2, m_omega) + m_comegap * BW_GS(q2, m_omegap) +
           m_comegapp * BW_GS(q2, m_omegapp));
  Form += 1. / 3. * (m_cphi * BW_GS(q2, m_phi) + m_cphip * BW_GS(q2, m_phip));
  return Form;
}

double FormFactor::Eval(const double &q2) {
  if (m_form_mode == finalstate::off) return 1;
  // https://gitlab.com/strong2020/monte-carlo-results/-/tree/root/pion-formfactor?ref_type=heads
  // for kplus see hep-ph/0409080 eq 64
  Particle_Info rhoi(kf_rho_770, m_rho_mass, 0.0, m_rho_g, 3, 0, 0, 0, 0, 1, 1,
                     "rho", "rho", "rho", "rho");

  Particle_Info rhopi(kf_rho_1450, m_rhop_mass, 0.0, m_rhop_g, 3, 0, 0, 0, 0, 1,
                      1, "rhopp", "rhopp", "rho", "rho");

  Particle_Info rhoppi(kf_rho_1700, m_rhopp_mass, 0.0, m_rhopp_g, 3, 0, 0, 0, 0,
                       1, 1, "rho ppp", "rho", "rho", "rho");

  Particle_Info rhopppi(kf_rho_2150, m_rhoppp_mass, 0.0, m_rhoppp_g, 3, 0, 0, 0,
                        0, 1, 1, "rho ppp", "rho", "rho", "rho");

  Particle_Info omegai(kf_omega_782, m_omega_mass, 0.0, m_omega_g, 3, 0, 0, 0,
                       0, 1, 1, "omega", "omega", "omega", "omega");
  Particle_Info omegapi(kf_omega_782, m_omegap_mass, 0.0, m_omegap_g, 3, 0, 0,
                        0, 0, 1, 1, "omegap", "omega", "omega", "omega");
  Particle_Info omegappi(kf_omega_782, m_omegapp_mass, 0.0, m_omegap_g, 3, 0, 0,
                         0, 0, 1, 1, "omegapp", "omega", "omega", "omega");

  Particle_Info phii(kf_phi_1020, m_phi_mass, 0.0, m_phi_g, 3, 0, 0, 0, 0, 1, 1,
                     "phi", "phi", "phi", "phi");
  Particle_Info phipi(kf_phi_1020, m_phip_mass, 0.0, m_phip_g, 3, 0, 0, 0, 0, 1,
                      1, "phi", "phi", "phi", "phi");

  m_rho = Flavour(rhoi);
  m_rhop = Flavour(rhopi);
  m_rhopp = Flavour(rhoppi);
  m_rhoppp = Flavour(rhopppi);

  m_omega = Flavour(omegai);
  m_omegap = Flavour(omegapi);
  m_omegapp = Flavour(omegappi);

  m_phi = Flavour(phii);
  m_phip = Flavour(phipi);
  Complex Form;
  switch (m_form_mode) {
    case finalstate::pion:
      Form = Pion(q2);
      break;
    case finalstate::kplus:
      Form = KPlus(q2);
      break;
    case finalstate::off:
      break;
    default:
      THROW(fatal_error, "Unknown final state");
  }
  // PRINT_VAR( (Form*conj(Form)).real());
  return (Form * conj(Form)).real();
}

std::ostream &ATOOLS::operator<<(std::ostream &str,
                                 const finalstate::code &fs) {
  if (fs == finalstate::pion)
    return str << "Pi+Pi-";
  else if (fs == finalstate::kplus)
    return str << "K+K-";
  else if (fs == finalstate::off)
    return str << "Off";
  return str << "unknown";
}

std::istream &ATOOLS::operator>>(std::istream &str, finalstate::code &mode) {
  std::string tag;
  str >> tag;
  // mode=wgt::off;
  if (tag.find("Pion") != std::string::npos)
    mode = finalstate::pion;
  else if (tag.find("Kplus") != std::string::npos)
    mode = finalstate::kplus;
  else if (tag.find("Off") != std::string::npos)
    mode = finalstate::off;
  else if (tag.find("None") != std::string::npos)
    mode = finalstate::off;
  else
    THROW(fatal_error, "Unknown Form_Factor: Mode ");
  return str;
}