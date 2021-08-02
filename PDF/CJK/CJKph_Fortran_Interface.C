#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/CJK/CJKph_Fortran_Interface.H"
#include <unistd.h>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
// CJKL: SUBROUTINE PARTONS(x,Q2,XPDF)
void partons_(double &, double &, double *);
// CJK1LO: CJK1GRID(ISET,IOPT,X,XC,XB,Q2,XPDF,F2alfa)
void cjk1grid_(int &, int &, double &, double &, double &, double &, double *,
               double &);
// CJK2LO: SUBROUTINE CJK2GRID(ISET,IOPT,X,XC,XB,Q2,XPDF,F2alfa)
void cjk2grid_(int &, int &, double &, double &, double &, double &, double *,
               double &);
// CJKHO: SUBROUTINE CJKHOGRID(IOPT,X,Q2,XPDF)
void cjkhogrid_(int &, double &, double &, double *);
}

CJKph_Fortran_Interface::CJKph_Fortran_Interface(const ATOOLS::Flavour _bunch,
                                                 const std::string _set) {
  m_xmin = 1.e-5;
  m_xmax = 1.;
  m_q2min = 1;
  m_q2max = 2.e5;
  m_nf = 5;

  m_iset = 0;
  m_set = _set;
  m_path = "";
  if (_set == std::string("CJK1LO"))
    m_path = rpa->gen.Variable("SHERPA_SHARE_PATH") + "/CJK1Grid";
  else if (_set == std::string("CJK2LO"))
    m_path = rpa->gen.Variable("SHERPA_SHARE_PATH") + "/CJK2Grid";
  else if (_set == std::string("CJKHO"))
    m_path = rpa->gen.Variable("SHERPA_SHARE_PATH") + "/CJKHOGrid";
  else {
    msg_Out() << "Warning: Unknown option for CJK photon PDF. Will fall back "
                 "to the default LO parametrizations";
    m_set = "CJKLLO";
  }

  m_bunch = _bunch;
  m_d = m_u = m_s = m_c = m_b = m_g = 0.;

  // TODO: Check the code below
  for (int i = 1; i < m_nf + 1; i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
}

PDF_Base *CJKph_Fortran_Interface::GetCopy() {
  return new CJKph_Fortran_Interface(m_bunch, m_set);
}

void CJKph_Fortran_Interface::CalculateSpec(const double &_x,
                                            const double &_Q2) {
  double x = _x / m_rescale, Q2 = _Q2;

  int iopt = 2;

  double pdf[11];

  if (m_set == std::string("CJKLLO"))
    partons_(x, Q2, pdf);
  else {
    // Change directory to the path, which contains the respective .dat files
    // save current dir to later change it back
    char buffer[1024];
    char *err = getcwd(buffer, 1024);
    if (chdir(m_path.c_str()) != 0 || err == nullptr)
      msg_Error() << "Error in CTEQ6_Fortran_Interface.C " << std::endl
                  << "   path " << m_path << " not found " << std::endl;

    // TODO: Check the use of the Chi^2 in the fortran code
    if (m_set == std::string("CJK1LO")) {
      double chi_c, chi_b = 1.;
      double f_gamma = 0.;
      cjk1grid_(m_iset, iopt, x, chi_c, chi_b, Q2, pdf, f_gamma);
    } else if (m_set == std::string("CJK2LO")) {
      double chi_c, chi_b = 1.;
      double f_gamma = 0.;
      cjk2grid_(m_iset, iopt, x, chi_c, chi_b, Q2, pdf, f_gamma);
    } else if (m_set == std::string("CJKHO")) {
      cjkhogrid_(iopt, x, Q2, pdf);
    }
    if (chdir(buffer) != 0)
      msg_Error() << "Error in CTEQ6_Fortran_Interface.C " << std::endl
                  << "   path " << m_path << " not found." << std::endl;
  }

  m_g = pdf[5];
  m_d = pdf[6];
  m_u = pdf[7];
  m_s = pdf[8];
  m_c = pdf[9];
  m_b = pdf[10];
}

double CJKph_Fortran_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
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

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                          sqr(rpa->gen.Ecms()));

  return m_rescale * value;
}

double CJKph_Fortran_Interface::GetXPDF(const kf_code &kf, bool anti) {
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

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"),
                                          sqr(rpa->gen.Ecms()));

  return m_rescale * value;
}

DECLARE_PDF_GETTER(CJKph_Getter);

PDF_Base *CJKph_Getter::operator()(const Parameter_Type &args) const {
  if (!args.m_bunch.IsPhoton())
    return NULL;
  return new CJKph_Fortran_Interface(args.m_bunch, args.m_set);
}

void CJKph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
  str << "CJK photon PDF parametrizations, see "
         "https://www.fuw.edu.pl/~pjank/param.html \n"
      << "The different parametrizations are \n"
      << " - CJKLLO: Phys.Rev.D68:014010,2003 (hep-ph/0212160) \n"
      << " - CJK1LO and CJK2LO: Nucl.Phys.Proc.Suppl.126:28-37,2004 "
         "(hep-ph/0310029) and hep-ph/0404244\n"
      << " - CJKHO: hep-ph/0404063\n";
}

CJKph_Getter *p_get_cjk[4];

extern "C" void InitPDFLib() {
  p_get_cjk[0] = new CJKph_Getter("CJKLLO");
  p_get_cjk[1] = new CJKph_Getter("CJK1LO");
  p_get_cjk[2] = new CJKph_Getter("CJK2LO");
  p_get_cjk[3] = new CJKph_Getter("CJKHO");
}

extern "C" void ExitPDFLib() {
  for (int i(0); i < 4; ++i)
    delete p_get_cjk[i];
}
