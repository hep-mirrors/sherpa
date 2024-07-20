#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/H1/H1_Fortran_Interface.H"
#include <unistd.h>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
// H1: SUBROUTINE QCD_2006(Z,Q2,IFIT,XPQ,F2,FL,C2,CL)
void qcd_2006_(double&, double&, int&, double*, double*, double*, double*,
               double*);
}

H1_Fortran_Interface::H1_Fortran_Interface(const std::string& _set)
{
  m_x     = 0.;
  m_Q2    = 0.;
  m_xmin  = 0.001;
  m_xmax  = 0.99;
  m_q2min = 1;
  m_q2max = 30000.;
  m_nf    = 3;//< pdfs for c,b,t are always zero

  m_set   = _set;
  m_bunch = Flavour(kf_pomeron);
  // just call it once with dummy variables to initialize the lookup tables and
  // be able to call them faster, c.f. the comment in qcd_2006.f
  double xpdf[13];
  double f2[2], fl[2], c2[2], cl[2];
  if (_set == std::string("FitA")) {
    m_member = 1;
  } else if (_set == std::string("FitB")) {
    m_member = 2;
  } else {
    msg_Out() << "Warning: Unknown option for H1 pomeron PDF. Will fall back "
                 "to the FitB parametrizations.\n";
    m_member = 2;
  }
  double x = 0.5, q2 = 10.;
  qcd_2006_(x, q2, m_member, xpdf, f2, fl, c2, cl);

  for (int i = 1; i < m_nf + 1; i++) {
    m_partons.insert(Flavour((kf_code) (i)));
    m_partons.insert(Flavour((kf_code) (i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());

  rpa->gen.AddCitation(
          1, "H1 pomeron DPDF fit is published under \\cite{H1:2006zyl}.");
}

PDF_Base* H1_Fortran_Interface::GetCopy()
{
  return new H1_Fortran_Interface(m_set);
}

void H1_Fortran_Interface::CalculateSpec(const double& _x, const double& _Q2)
{
  m_x  = _x / m_rescale;
  m_Q2 = _Q2;
}

double H1_Fortran_Interface::GetXPDF(const ATOOLS::Flavour& flavour)
{
  return GetXPDF(flavour.Kfcode(), flavour.IsAnti());
}

double H1_Fortran_Interface::GetXPDF(const kf_code& kf, bool anti)
{
  double xpdf[13];
  double f2[2], fl[2], c2[2], cl[2];
  int    iset = 0;
  qcd_2006_(m_x, m_Q2, iset, xpdf, f2, fl, c2, cl);
  int index = static_cast<int>(kf) % 7 + 6;
  double result = xpdf[index] * m_rescale;
  return result;
}

DECLARE_PDF_GETTER(H1_Getter);

PDF_Base* H1_Getter::operator()(const Parameter_Type& args) const
{
  if (args.m_bunch.Kfcode() != kf_pomeron) return nullptr;
  return new H1_Fortran_Interface(args.m_set);
}

void H1_Getter::PrintInfo(std::ostream& str, const size_t width) const
{
  str << "H1 pomeron DPDF parametrizations, see Eur.Phys.J.C 48 (2006) "
         "715-748. Sets are FitA and FitB";
}

H1_Getter* p_get_H1[2];

extern "C" void InitPDFLib()
{
  p_get_H1[0] = new H1_Getter("FitA");
  p_get_H1[1] = new H1_Getter("FitB");
}

extern "C" void ExitPDFLib()
{
  for (int i(0); i < 2; ++i) delete p_get_H1[i];
}
