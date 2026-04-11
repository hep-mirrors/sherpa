#ifndef SHERPA_PHOTON_PDF_BASE_H
#define SHERPA_PHOTON_PDF_BASE_H

#include "ATOOLS/Phys/Flavour.H"
#include "PDF_Base.H"

namespace PDF {
class Photon_PDF_Base : public PDF_Base {
protected:
  double m_d, m_u, m_s, m_c, m_b, m_g, m_t, m_ph;
  int m_iset;
public:
  Photon_PDF_Base(const ATOOLS::Flavour, const std::string _set, int nf);

  double GetXPDF(const ATOOLS::Flavour &flavour) override;
  double GetXPDF(const kf_code& kf, bool anti) override;

  double GetQ2Min() { return m_q2min; }
  double GetQ2Max() { return m_q2max; }
};
}
#endif
