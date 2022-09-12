#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/Main/Photon_PDF_Base.H"
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
                                                   const std::string _set)
    : Photon_PDF_Base(_bunch, _set, 6) {
  m_xmin = 1.e-5;
  m_xmax = 1.;
  m_q2min = .25;
  m_q2max = 1.e6;

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

  rpa->gen.AddCitation(1, "The SaSg photon PDF is published under "
                          "\\cite{Schuler:1995fk} and \\cite{Schuler:1996fc}.");
}

PDF_Base *SASGph_Fortran_Interface::GetCopy() {
  return new SASGph_Fortran_Interface(m_bunch, m_set);
}

void SASGph_Fortran_Interface::CalculateSpec(const double &_x,
                                             const double &_Q2) {
  if (m_include_photon_in_photon)
    m_ph = GetPhotonCoefficient(_x, _Q2);

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
