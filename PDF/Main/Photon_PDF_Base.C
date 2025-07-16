#include "Photon_PDF_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Main/PDF_Base.H"

using namespace PDF;
using namespace ATOOLS;

Photon_PDF_Base::Photon_PDF_Base(const Flavour _bunch, const std::string _set, int nf) : PDF_Base() {
  m_d = m_u = m_s = m_c = m_b = m_g = m_t = m_ph = 0.;
  m_set = _set;
  m_bunch = _bunch;
  m_nf = nf;
  m_iset = 1;
  for (int i = 1; i < m_nf + 1; i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());

  // Insert the photon component
  Settings& s = Settings::GetMainSettings();
  m_include_photon_in_photon = s["INCLUDE_PHOTON_IN_PHOTON_PDF"].Get<bool>();
  if (m_include_photon_in_photon) {
    m_partons.insert(Flavour(kf_photon));
  }
}

double Photon_PDF_Base::GetPhotonCoefficient(double x,double Q2) {
  // The coefficient is proportional to a Dirac delta d(1-x)
  double dx = Max(1.-m_xmax, 1.e-6);
  if (x < 1. - dx)
    return m_ph = 0.;
  // Get the coefficient of the photon component in the photon, c.f. hep-ph/9605240
  // It is equal 1, subtracted by the VMD and anomalous (i.e. perturbative q-qbar) terms
  // The lepton component is neglected
  m_ph = 1.;
  const double alphaem = MODEL::s_model->ScalarFunction(std::string("alpha_QED"), 0);

  const double cutscale = 0.5; // cut-off scale for the perturbative contribution
  for (int i = 1; i < m_nf + 1; i++) {
    m_ph -= alphaem / 2 / M_PI * 2 * sqr(Flavour((kf_code)(i)).Charge()) * log(Q2 / sqr(cutscale));
  }

  // meson-photon couplings for rho, omega and phi: (might need updating)
  const std::array<double, 3> mesoncouplings = {2.2, 23.6, 18.4};
  for (auto mit = mesoncouplings.begin(); mit != mesoncouplings.end(); mit++) {
    m_ph -= alphaem / *mit;
  }
  msg_Debugging() << METHOD << ": Calculation photon->photon pdf, val = " << m_ph
            << "\n";
  if (m_ph < 0)
    msg_Error() << METHOD << ": Photon component is negative! Check the point, scale = "
                << Q2 << ", photon coefficient = " << m_ph << std::endl;
  return m_ph/dx;
}

double Photon_PDF_Base::GetXPDF(const ATOOLS::Flavour &infl) {
  return GetXPDF(infl.Kfcode(), false);
}

double Photon_PDF_Base::GetXPDF(const kf_code &kf, bool anti) {
  double value = 0.;

  if (kf == kf_gluon)
    value = m_g;
  else if (kf == kf_d)
    value = m_d;
  else if (kf == kf_u)
    value = m_u;
  else if (kf == kf_s)
    value = m_s;
  else if (kf == kf_c)
    value = m_c;
  else if (kf == kf_b)
    value = m_b;
  else if (kf == kf_t)
    value = m_t;
  else if (kf == kf_photon)
    value = m_ph;

  return m_rescale * value;
}