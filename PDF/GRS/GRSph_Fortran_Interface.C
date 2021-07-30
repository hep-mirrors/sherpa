#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/GRS/GRSph_Fortran_Interface.H"
#include <unistd.h>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
// GRSG99(ISET,X,Q2,UL,DL,SL,GL)
void grsg99_(int &, double &, double &, double &, double &, double &, double &);
}

extern "C" {
extern struct { char mfile[128]; } mrinput_;
#define input mrinput_
}

inline void MakeFortranString(char *output, std::string input,
                              unsigned int length) {
  for (unsigned int i = 0; i < length; ++i)
    output[i] = (char)32;
  for (size_t j = 0; j < input.length(); ++j)
    output[j] = (char)input[j];
}

GRSph_Fortran_Interface::GRSph_Fortran_Interface(const ATOOLS::Flavour _bunch) {
  /*
   * ISET = 1  LEADING ORDER SET (DATA FILE 'grsg99lo.grid')        *
   *           ISET = 2  NEXT-TO-LEADING ORDER MSbar SET            *
   *                     (DATA FILE 'grsg99m.grid')                 *
   *           ISET = 3  NEXT-TO-LEADING ORDER DISgamma SET         *
   *                    (DATA FILE 'grsg99d.grid')
   */
  iset = 1;
  m_path = rpa->gen.Variable("SHERPA_SHARE_PATH") + "/GRSGrid";

  m_xmin = 1.e-5;
  m_xmax = 1.;
  m_q2min = .4;
  m_q2max = 1.e6;
  m_nf = 3;

  m_set = "GRS";
  m_bunch = _bunch;
  m_d = m_u = m_s = m_g = 0.;

  for (int i = 1; i < 6; i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
  MakeFortranString(input.mfile, m_path + std::string("/grsg99lo.grid"), 128);
}

PDF_Base *GRSph_Fortran_Interface::GetCopy() {
  return new GRSph_Fortran_Interface(m_bunch);
}

void GRSph_Fortran_Interface::CalculateSpec(const double &_x,
                                            const double &_Q2) {
  double x = _x / m_rescale, Q2 = _Q2;

  grsg99_(iset, x, Q2, m_u, m_d, m_s, m_g);
}

double GRSph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
  double value = 0.;

  if (infl.Kfcode() == kf_gluon)
    value = m_g;
  else if (infl.Kfcode() == kf_d)
    value = m_d;
  else if (infl.Kfcode() == kf_u)
    value = m_u;
  else if (infl.Kfcode() == kf_s)
    value = m_s;

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                          sqr(rpa->gen.Ecms()));

  return m_rescale * value;
}

double GRSph_Fortran_Interface::GetXPDF(const kf_code &kf, bool anti) {
  double value = 0.;

  if (kf == kf_gluon)
    value = m_g;
  else if (kf == kf_d)
    value = m_d;
  else if (kf == kf_u)
    value = m_u;
  else if (kf == kf_s)
    value = m_s;

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                          sqr(rpa->gen.Ecms()));

  return m_rescale * value;
}

DECLARE_PDF_GETTER(GRSph_Getter);

PDF_Base *GRSph_Getter::operator()(const Parameter_Type &args) const {
  if (!args.m_bunch.IsPhoton())
    return NULL;
  return new GRSph_Fortran_Interface(args.m_bunch);
}

void GRSph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
  str << "GRS photon PDF, see Phys. Rev. D60(1999)054019";
}

GRSph_Getter *p_get_GRS;

extern "C" void InitPDFLib() { p_get_GRS = new GRSph_Getter("GRS"); }

extern "C" void ExitPDFLib() { delete p_get_GRS; }
