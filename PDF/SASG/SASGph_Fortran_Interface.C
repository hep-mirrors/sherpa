#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/SASG/SASGph_Fortran_Interface.H"
#include <unistd.h>

#include <iostream>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
// TODO: interface the IP2 variable.
// CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
void sasgam_(int &, float &, float &, float &, int &, float &, float *);
}

SASGph_Fortran_Interface::SASGph_Fortran_Interface(const ATOOLS::Flavour _bunch,
                                                   const std::string _set) {
  m_xmin = 1.e-5;
  m_xmax = 1.;
  m_q2min = .25;
  m_q2max = 1.e6;
  m_nf = 6;

  m_set = _set;
  m_bunch = _bunch;
  m_d = m_u = m_s = m_c = m_b = m_g = m_t = 0.;

  if (m_set == std::string("SAS1D"))
    m_iset = 1;
  else if (m_set == std::string("SAS1M"))
    m_iset = 2;
  else if (m_set == std::string("SAS2D"))
    m_iset = 3;
  else if (m_set == std::string("SAS2M"))
    m_iset = 4;
  else {
    msg_Out() << METHOD
              << ": Cannot recognize the chosen PDF parametrization. Will "
                 "use the Leading Order parametrization. \n";
    m_iset = 1;
  }

  for (int i = 1; i < m_nf + 1; i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());

  rpa->gen.AddCitation(1,"The SaSg photon PDF is published under \\cite{Schuler:1995fk} and \\cite{Schuler:1996fc}.");
}

PDF_Base *SASGph_Fortran_Interface::GetCopy() {
  return new SASGph_Fortran_Interface(m_bunch, m_set);
}

void SASGph_Fortran_Interface::CalculateSpec(const double &_x,
                                             const double &_Q2) {
  float x = _x / m_rescale, Q2 = _Q2;

  float f2photon = 0;
  float pdf[13];

  // CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
  sasgam_(m_iset, x, Q2, P2, IP2, f2photon, pdf);
  m_g = pdf[6];
  m_d = pdf[7];
  m_u = pdf[8];
  m_s = pdf[9];
  m_c = pdf[10];
  m_b = pdf[11];
  m_t = pdf[12];
}

double SASGph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
  double value = 0.;

  if (infl.Kfcode() == kf_gluon)
    value = m_g;
  else if (infl.Kfcode() == kf_d)
    value = m_d;
  else if (infl.Kfcode() == kf_u)
    value = m_u;
  else if (infl.Kfcode() == kf_s)
    value = m_s;
  else if (infl.Kfcode() == kf_c)
    value = m_c;
  else if (infl.Kfcode() == kf_b)
    value = m_b;
  else if (infl.Kfcode() == kf_t)
    value = m_t;

  return m_rescale * value;
}

double SASGph_Fortran_Interface::GetXPDF(const kf_code &kf, bool anti) {
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

  return m_rescale * value;
}

DECLARE_PDF_GETTER(SASGph_Getter);

PDF_Base *SASGph_Getter::operator()(const Parameter_Type &args) const {
  if (!args.m_bunch.IsPhoton())
    return NULL;
  return new SASGph_Fortran_Interface(args.m_bunch, args.m_set);
}

void SASGph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
  str << "SASG photon PDF, see Z. Phys. C68 (1995) 607 and Phys. Lett. B376 "
         "(1996) 193.\n"
         "The following sets can be used: \n"
         " - SAS1D: SaS set 1D ('DIS',   Q0 = 0.6 GeV)\n"
         " - SAS1M: SaS set 1M ('MSbar', Q0 = 0.6 GeV)\n"
         " - SAS2D: SaS set 2D ('DIS',   Q0 =  2  GeV)\n"
         " - SAS2M: SaS set 2M ('MSbar', Q0 =  2  GeV)\n";
}

SASGph_Getter *p_get_sasg[4];

extern "C" void InitPDFLib() {
  p_get_sasg[0] = new SASGph_Getter("SAS1D");
  p_get_sasg[1] = new SASGph_Getter("SAS1M");
  p_get_sasg[2] = new SASGph_Getter("SAS2D");
  p_get_sasg[3] = new SASGph_Getter("SAS2M");
}

extern "C" void ExitPDFLib() {
  for (int i = 0; i < 4; i++)
    delete p_get_sasg[i];
}
