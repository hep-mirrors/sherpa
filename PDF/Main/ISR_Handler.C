#include "PDF/Main/ISR_Handler.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Blob.H"
#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Base.H"
#include "REMNANTS/Main/Remnant_Base.H"

using namespace ATOOLS;
using namespace PDF;

static int s_nozeropdf = -1;

double Lambda2(double sp, double sp1, double sp2) {
  return (sp - sp1 - sp2) * (sp - sp1 - sp2) - 4.0 * sp1 * sp2;
}

struct isrmode {
  enum code {
    none = 0,
    hadron_hadron = 1,
    lepton_hadron = 2,
    hadron_lepton = 3,
    lepton_lepton = 4,
    unknown = 99
  };
};

std::ostream &PDF::operator<<(std::ostream &s, const PDF::isrmode::code mode) {
  switch (mode) {
  case isrmode::none:
    s << "none";
    break;
  case isrmode::hadron_hadron:
    s << "hadron_hadron";
    break;
  case isrmode::hadron_lepton:
    s << "hadron_lepton";
    break;
  case isrmode::lepton_hadron:
    s << "lepton_hadron";
    break;
  case isrmode::lepton_lepton:
    s << "lepton_lepton";
    break;
  case isrmode::unknown:
  default:
    s << "unknown";
    break;
  }
  return s;
}

ISR_Handler::ISR_Handler(ISR_Base **isrbase)
    : m_rmode(0), m_swap(0), m_info_lab(8), m_info_cms(8),
      m_freezePDFforLowQ(false) {
  p_isrbase[0] = isrbase[0];
  p_isrbase[1] = isrbase[1];
  Settings &s = Settings::GetMainSettings();
  m_freezePDFforLowQ = s["FREEZE_PDF_FOR_LOW_Q"].SetDefault(false).Get<bool>();
  if (s_nozeropdf < 0) {
    s_nozeropdf = s["NO_ZERO_PDF"].SetDefault(0).Get<int>();
  }
  m_xf1 = m_xf2 = 1.0;
  p_remnants[1] = p_remnants[0] = NULL;
  m_mode = 0;
  for (short int i = 0; i < 2; i++) {
    if (p_isrbase[i]->On())
      m_mode += i + 1;
    m_mass2[i] = sqr(p_isrbase[i]->Flavour().Mass());
    m_x[i] = 1.;
    m_mu2[i] = 0.;
  }
  FixType();
}

void ISR_Handler::FixType() {
  m_type = isrmode::unknown;
  isrtype::code type[2];
  for (size_t i = 0; i < 2; i++)
    type[i] = p_isrbase[i]->Type();
  if (type[0] == isrtype::hadron && type[1] == isrtype::hadron)
    m_type = isrmode::hadron_hadron;
  if ((type[0] == isrtype::lepton || type[0] == isrtype::intact) &&
      type[1] == isrtype::hadron)
    m_type = isrmode::lepton_hadron;
  if ((type[1] == isrtype::lepton || type[1] == isrtype::intact) &&
      type[0] == isrtype::hadron)
    m_type = isrmode::hadron_lepton;
  if ((type[0] == isrtype::lepton || type[0] == isrtype::intact) &&
      (type[1] == isrtype::lepton || type[1] == isrtype::intact))
    m_type = isrmode::lepton_lepton;
  if (type[0] == isrtype::intact && type[1] == isrtype::intact)
    m_type = isrmode::none;
}

void ISR_Handler::Output() {
  if (msg_LevelIsTracking() || msg_LevelIsInfo())
    msg_Out() << "ISR_Handler: type = " << m_type << "\n"
              << "    for " << p_isrbase[0]->Flavour()
              << " (internal structure = " << p_isrbase[0]->On() << ")\n"
              << "    and " << p_isrbase[1]->Flavour()
              << " (internal structure = " << p_isrbase[1]->On() << ")\n";
}

void ISR_Handler::Init(double *splimits) {
  double s = (p_beam[0]->OutMomentum() + p_beam[1]->OutMomentum()).Abs2();
  ATOOLS::rpa->gen.SetEcms(sqrt(s));

  m_splimits[0] = s * splimits[0];
  m_splimits[1] = ATOOLS::Min(s * splimits[1], s * Upper1() * Upper2());
  m_splimits[2] = s;
  m_fixed_smin = m_splimits[0];
  m_fixed_smax = m_splimits[1];
  m_ylimits[0] = -10.;
  m_ylimits[1] = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * p_isrbase[0]->Exponent() * p_isrbase[1]->Exponent();
  double E = ATOOLS::rpa->gen.Ecms();
  double x = 1. / 2. + (m_mass2[0] - m_mass2[1]) / (2. * E * E);
  double E1 = x * E;
  double E2 = E - E1;
  m_fixvecs[0] = Vec4D(E1, 0., 0., sqrt(sqr(E1) - m_mass2[0]));
  m_fixvecs[1] = Vec4D(E2, 0., 0., -m_fixvecs[0][3]);
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *bunches,
                                   ATOOLS::Flavour *partons) {
  bool fit = 1;
  for (int i = 0; i < 2; i++) {
    if (p_isrbase[i]->On()) {
      if (bunches[i] != PDF(i)->Bunch()) {
        fit = 0;
        break;
      }
      fit = PDF(i)->Contains(partons[i]);
      if (fit == 0)
        break;
    } else {
      bool found(false);
      if (p_isrbase[i]->Flavour().Includes(partons[i])) {
        found = true;
      }
      if (!found)
        return false;
    }
  }
  return fit;
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *partons) {
  bool fit = 1;
  for (int i = 0; i < 2; i++) {
    if (partons[i].Kfcode() == 0)
      continue;
    if (p_isrbase[i]->On()) {
      fit = PDF(i)->Contains(partons[i]);
      if (fit == 0) {
        for (size_t j(0); j < partons[i].Size(); ++j)
          fit |= PDF(i)->Contains(partons[i][j]);
      }
      if (fit == 0)
        break;
    } else {
      bool found(false);
      if (p_isrbase[i]->Flavour().Includes(partons[i])) {
        found = true;
      }
      if (!found)
        return false;
    }
  }
  return fit;
}

void ISR_Handler::SetMasses(const Flavour_Vector &fl) {
  m_mass2[0] = sqr(fl[0].Mass());
  m_mass2[1] = sqr(fl[1].Mass());
  double emin = 0.0;
  for (size_t i = 2; i < fl.size(); ++i)
    emin += fl[i].Mass();
  emin = ATOOLS::Max(emin, fl[0].Mass() + fl[1].Mass());
  m_splimits[0] = ATOOLS::Max(m_splimits[0], sqr(emin));
}

void ISR_Handler::SetPartonMasses(const Flavour_Vector &fl) {
  SetMasses(fl);
  double E = ATOOLS::rpa->gen.Ecms();
  double x = 1. / 2. + (m_mass2[0] - m_mass2[1]) / (2. * E * E);
  double E1 = x * E;
  double E2 = E - E1;
  m_fixvecs[0] = Vec4D(E1, 0., 0., sqrt(sqr(E1) - m_mass2[0]));
  m_fixvecs[1] = Vec4D(E2, 0., 0., -m_fixvecs[0][3]);
}

bool ISR_Handler::MakeISR(const double &sp, const double &y, Vec4D *p,
                          const Flavour_Vector &flavs) {
  if (m_mode == 0) {
    m_x[1] = m_x[0] = 1.;
    return true;
  }
  if (sp < m_splimits[0] || sp > m_splimits[1]) {
    msg_Error() << METHOD << "(..): " << om::red << "s' out of bounds.\n"
                << om::reset
                << "  s'_{min}, s'_{max 1,2} vs. s': " << m_splimits[0] << ", "
                << m_splimits[1] << ", " << m_splimits[2] << " vs. " << sp
                << std::endl;
    return false;
  }
  Vec4D p0 = p_beam[0]->OutMomentum();
  Vec4D p1 = p_beam[1]->OutMomentum();
  double gam = p0 * p1 + sqrt(sqr(p0 * p1) - p0.Abs2() * p1.Abs2());
  double bet = 1.0 / (1.0 - p0.Abs2() / gam * p1.Abs2() / gam);
  Vec4D pp = bet * (p0 - p0.Abs2() / gam * p1),
        pm = bet * (p1 - p1.Abs2() / gam * p0);
  double s = 2.0 * pp * pm;
  double tau = 0.5 / s * (sp - m_mass2[0] - m_mass2[1]);
  if (tau * tau < m_mass2[0] * m_mass2[1] / (s * s)) {
    msg_Error() << METHOD << "(): s' out of range." << std::endl;
    return false;
  }
  tau += sqrt(tau * tau - m_mass2[0] * m_mass2[1] / (s * s));
  if (m_mode == 1) {
    m_x[1] = m_xkey[5] = p1.PMinus() / pm.PMinus();
    m_x[0] = m_xkey[4] = tau / m_x[1];
  } else if (m_mode == 2) {
    m_x[0] = m_xkey[4] = p0.PPlus() / pp.PPlus();
    m_x[1] = m_xkey[5] = tau / m_x[0];
  } else if (m_mode == 3) {
    double yt =
        exp(y - 0.5 * log((tau + m_mass2[1] / s) / (tau + m_mass2[0] / s)) -
            (pp + pm).Y());
    tau = sqrt(tau);
    m_x[0] = m_xkey[4] = tau * yt;
    m_x[1] = m_xkey[5] = tau / yt;
  } else {
    THROW(fatal_error, "Invalid ISR mode");
  }
  if (PDF(0) && (m_x[0] < PDF(0)->XMin() || m_x[0] > PDF(0)->XMax()))
    return false;
  if (PDF(1) && (m_x[1] < PDF(1)->XMin() || m_x[1] > PDF(1)->XMax()))
    return false;
  p[0] = m_x[0] * pp + m_mass2[0] / s / m_x[0] * pm;
  p[1] = m_x[1] * pm + m_mass2[1] / s / m_x[1] * pp;
  if (m_swap) {
    std::swap<Vec4D>(p[0], p[1]);
  }
  m_cmsboost = Poincare(p[0] + p[1]);
  if (m_x[0] >= 1.0)
    m_x[0] = 1.0 - 1.0e-12;
  if (m_x[1] >= 1.0)
    m_x[1] = 1.0 - 1.0e-12;
  return true;
}

bool ISR_Handler::GenerateSwap(const ATOOLS::Flavour &f1,
                               const ATOOLS::Flavour &f2) {
  if (m_swap)
    m_swap = false;
  if (!AllowSwap(f1, f2))
    return false;
  m_swap = ran->Get() > 0.5;
  return true;
}

bool ISR_Handler::AllowSwap(const ATOOLS::Flavour &f1,
                            const ATOOLS::Flavour &f2) const {
  if (p_isrbase[0]->PDF() == NULL || p_isrbase[1]->PDF() == NULL)
    return false;
  int ok[2] = {0, 0};
  if (p_isrbase[0]->PDF()->Contains(f2))
    ok[0] = 1;
  else
    for (size_t j(0); j < f2.Size(); ++j)
      if (p_isrbase[0]->PDF()->Contains(f2[j])) {
        ok[0] = 1;
        break;
      }
  if (p_isrbase[1]->PDF()->Contains(f1))
    ok[1] = 1;
  else
    for (size_t j(0); j < f1.Size(); ++j)
      if (p_isrbase[1]->PDF()->Contains(f1[j])) {
        ok[1] = 1;
        break;
      }
  return ok[0] && ok[1];
}

void ISR_Handler::Reset() { m_splimits[1] = m_fixed_smax; }

void ISR_Handler::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign("ISR::s'", 5, 0, info);
  m_ykey.Assign("ISR::y", 3, 0, info);
  // Convention for m_xkey:
  // [x_{min,beam0}, x_{min,beam1}, x_{max,beam0}, x_{max,beam1}, x_{val,beam0},
  // x_{val,beam1}]. The limits, i.e. index 0,1,2,3 are saved as log(x), the
  // values are saved linearly.
  m_xkey.Assign("ISR::x", 6, 0, info);
  SetLimits();
}

void ISR_Handler::SetLimits() {
  for (short int i = 0; i < 3; ++i) {
    m_sprimekey[i] = m_splimits[i];
    if (i < 2)
      m_ykey[i] = m_ylimits[i];
  }
  m_xkey[0] = ((m_mass2[0] == 0.0)
                   ? -0.5 * std::numeric_limits<double>::max()
                   : log(m_mass2[0] / sqr(p_beam[0]->OutMomentum().PPlus())));
  m_xkey[1] = ((m_mass2[1] == 0.0)
                   ? -0.5 * std::numeric_limits<double>::max()
                   : log(m_mass2[1] / sqr(p_beam[1]->OutMomentum().PMinus())));
  double e1 = p_beam[0]->OutMomentum()[0];
  m_xkey[2] = ATOOLS::Min(e1 / p_beam[0]->OutMomentum().PPlus() *
                              (1.0 + sqrt(1.0 - m_mass2[0] / sqr(e1))),
                          Upper1());
  double e2 = p_beam[1]->OutMomentum()[0];
  m_xkey[3] = ATOOLS::Min(e2 / p_beam[1]->OutMomentum().PMinus() *
                              (1.0 + sqrt(1.0 - m_mass2[1] / sqr(e2))),
                          Upper2());
  m_sprimekey[1] = m_splimits[1] =
      Min(m_splimits[1], m_splimits[2] * m_xkey[2] * m_xkey[3]);
  m_xkey[2] = log(m_xkey[2]);
  m_xkey[3] = log(m_xkey[3]);
}

double ISR_Handler::PDFWeight(const int mode, Vec4D p1, Vec4D p2, double Q12,
                              double Q22, Flavour fl1, Flavour fl2, int warn) {
  // mode&2 -> override m_mode and only calc left beam
  // mode&4 -> override m_mode and only calc right beam
  // mode&8 -> return xf in mode 2 & 4
  if (m_mode == 0)
    return 1.;
  msg_IODebugging() << METHOD << "(mode = " << mode << ")\n";
  if (fl1.Size() > 1 || fl2.Size() > 1)
    THROW(fatal_error,
          "Do not try to calculate an ISR weight with containers.");
  double x1(0.), x2(0.);
  if (p1[3] < p2[3]) {
    std::swap<Flavour>(fl1, fl2);
    std::swap<Vec4D>(p1, p2);
    std::swap<double>(Q12, Q22);
  }
  x1 = CalcX(p1);
  x2 = CalcX(p2);
  if (IsBad(x1) || IsBad(x2) || IsBad(Q12) || IsBad(Q22)) {
    msg_Error() << "Bad PDF input: x1=" << x1 << ", x2=" << x2
                << ", Q12=" << Q12 << ", Q22=" << Q22 << std::endl;
    return 0.;
  }
  if (m_freezePDFforLowQ) {
    if (Q12 < PDF(0)->Q2Min())
      Q12 = 1.001 * PDF(0)->Q2Min();
    if (Q22 < PDF(1)->Q2Min())
      Q22 = 1.001 * PDF(1)->Q2Min();
  }
  msg_IODebugging() << "  " << p1 << " from " << p_beam[0]->OutMomentum()
                    << " -> " << p1.PPlus() << " / "
                    << p_beam[0]->OutMomentum().PPlus() << " = " << x1
                    << std::endl;
  msg_IODebugging() << "  " << p2 << " from " << p_beam[1]->OutMomentum()
                    << " -> " << p2.PMinus() << " / "
                    << p_beam[1]->OutMomentum().PMinus() << " = " << x2
                    << std::endl;
  if (warn && PDF(0) && (Q12 < PDF(0)->Q2Min() || Q12 > PDF(0)->Q2Max())) {
    msg_IODebugging() << "  Q_1^2 out of bounds" << std::endl;
    return 0.;
  }
  if (warn && PDF(1) && (Q22 < PDF(1)->Q2Min() || Q22 > PDF(1)->Q2Max())) {
    msg_IODebugging() << "  Q_2^2 out of bounds" << std::endl;
    return 0.;
  }
  m_mu2[mode & 1] = Q12;
  m_mu2[1 - (mode & 1)] = Q22;
  int cmode(((mode & 6) >> 1) ? ((mode & 6) >> 1) : m_mode);
  // cmode & 1 -> include first PDF; cmode & 2 -> include second PDF
  if ((cmode == 1 && PDF(0) == NULL) || (cmode == 2 && PDF(1) == NULL))
    return 1.0;
  const auto include_both_pdfs = (cmode == 3);
  // Checking remnant kinematics only makes sense when both PDFs are included;
  // NOTE: the check is done here because CheckRemnantKinematics internally
  // calls the PDF’s CalculateSpec, which overwrites its m_Q member; we don't
  // want this to happen after the following switch-statement, where this is
  // set in preparation of the p_isrbase[i]->Weight(fl) calls that get the PDF
  // values below
  const auto has_nonexistent_or_correct_remnant_kinematics =
      (!include_both_pdfs || (CheckRemnantKinematics(fl1, x1, 0, false) &&
                              CheckRemnantKinematics(fl2, x2, 1, false)));
  switch (cmode) {
  case 3:
    if (!p_isrbase[0]->PDF()->Contains(fl1) ||
        !p_isrbase[1]->PDF()->Contains(fl2)) {
      return 0.;
    }
    if (x1 > p_isrbase[0]->PDF()->RescaleFactor() ||
        x2 > p_isrbase[1]->PDF()->RescaleFactor()) {
      return 0.;
    }
    if (!(p_isrbase[0]->CalculateWeight(x1, 0.0, 0.0, Q12, warn) &&
          p_isrbase[1]->CalculateWeight(x2, 0.0, 0.0, Q22, warn))) {
      return 0.;
    }
    break;
  case 2:
    if (!p_isrbase[1]->PDF()->Contains(fl2)) {
      return 0.;
    }
    if (x2 > p_isrbase[1]->PDF()->RescaleFactor()) {
      return 0.;
    }
    if (!p_isrbase[1]->CalculateWeight(x2, 0.0, 0.0, Q22, warn)) {
      return 0.;
    }
    break;
  case 1:
    if (!p_isrbase[0]->PDF()->Contains(fl1)) {
      return 0.;
    }
    if (x1 > p_isrbase[0]->PDF()->RescaleFactor()) {
      return 0.;
    }
    if (!p_isrbase[0]->CalculateWeight(x1, 0.0, 0.0, Q12, warn)) {
      return 0.;
    }
    break;
  case 0:
    break;
  default:
    return 0.;
  }
  if (has_nonexistent_or_correct_remnant_kinematics) {
    double f1 = (cmode & 1) ? p_isrbase[0]->Weight(fl1) : 1.0;
    double f2 = (cmode & 2) ? p_isrbase[1]->Weight(fl2) : 1.0;
    m_xf1 = x1 * f1;
    m_xf2 = x2 * f2;
    msg_IODebugging() << "  PDF1: " << rpa->gen.Beam1() << " -> " << fl1
                      << " at (" << x1 << "," << sqrt(Q12) << ") -> "
                      << om::bold << f1 << om::reset << "\n";
    msg_IODebugging() << "  PDF2: " << rpa->gen.Beam2() << " -> " << fl2
                      << " at (" << x2 << "," << sqrt(Q22) << ") -> "
                      << om::bold << f2 << om::reset << "\n";
    msg_IODebugging() << "  Weight: " << f1 * f2 << std::endl;
    if (IsBad(f1 * f2))
      return 0.0;
    if (s_nozeropdf && f1 * f2 == 0.0)
      return pow(std::numeric_limits<double>::min(), 0.25);
    if (mode & 8) {
      if (mode & 2)
        return x1 * f1;
      if (mode & 4)
        return x2 * f2;
    }
    return f1 * f2;
  }
  return 0.;
}

double ISR_Handler::Flux(const Vec4D &p1, const Vec4D &p2) {
  return 0.25 / sqrt(sqr(p1 * p2) - p1.Abs2() * p2.Abs2());
}

double ISR_Handler::Flux(const Vec4D &p1) { return 0.5 / p1.Mass(); }

double ISR_Handler::CalcX(const ATOOLS::Vec4D &p) {
  if (p[3] > 0.)
    return Min(PDF(0) ? PDF(0)->XMax() : 1.,
               p.PPlus() / p_beam[0]->OutMomentum().PPlus());
  else
    return Min(PDF(1) ? PDF(1)->XMax() : 1.,
               p.PMinus() / p_beam[1]->OutMomentum().PMinus());
}

bool ISR_Handler::BoostInCMS(Vec4D *p, const size_t n) {
  for (size_t i = 0; i < n; ++i) {
    m_cmsboost.Boost(p[i]);
  }
  return true;
}

bool ISR_Handler::BoostInLab(Vec4D *p, const size_t n) {
  for (size_t i = 0; i < n; ++i) {
    m_cmsboost.BoostBack(p[i]);
  }
  return true;
}

bool ISR_Handler::CheckRemnantKinematics(const ATOOLS::Flavour &fl, double &x,
                                         int beam, bool swapped) {
  if (x > p_isrbase[beam]->PDF()->RescaleFactor())
    return false;
  if (m_rmode == 0)
    return true;
  p_remnants[beam]->Reset();
  double pp(beam == 0 ? x * p_beam[0]->OutMomentum().PPlus()
                      : x * p_beam[1]->OutMomentum().PMinus());
  double pm(sqr(fl.Mass()));
  pm /= pp;
  Vec4D mom((pp + pm) / 2.0, 0.0, 0.0,
            beam == 0 ? (pp - pm) / 2.0 : (pm - pp) / 2.0);
  return p_remnants[beam]->TestExtract(fl, mom);
}

void ISR_Handler::Reset(const size_t i) const {
  // if (p_isrbase[i]->PDF()!=NULL) p_isrbase[i]->Reset();
}

ATOOLS::Blob_Data_Base *ISR_Handler::Info(const int frame) const {
  if (frame == 0)
    return new ATOOLS::Blob_Data<std::vector<double>>(m_info_cms);
  return new ATOOLS::Blob_Data<std::vector<double>>(m_info_lab);
}
