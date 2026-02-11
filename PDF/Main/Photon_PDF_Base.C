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
