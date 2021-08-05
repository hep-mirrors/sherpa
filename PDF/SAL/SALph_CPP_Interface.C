#include "SALph_CPP_Interface.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include <unistd.h>

#include "sal.h"

using namespace PDF;
using namespace ATOOLS;

SALph_CPP_Interface::SALph_CPP_Interface(const ATOOLS::Flavour _bunch) {
  m_xmin = 1.e-5;
  m_xmax = 0.9999;
  m_q2min = 2.;
  m_q2max = 8.e4;
  m_nf = 6;

  m_set = "SAL";
  m_bunch = _bunch;
  m_d = m_u = m_s = m_c = m_b = m_t = m_g = 0.;

  for (int i = 1; i < m_nf + 1; i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
}

PDF_Base *SALph_CPP_Interface::GetCopy() {
  return new SALph_CPP_Interface(m_bunch);
}

void SALph_CPP_Interface::CalculateSpec(const double &_x, const double &_Q2) {
  double x = _x / m_rescale, Q2 = _Q2;

  if (x < m_xmin || x > m_xmax)
    return;
  // indeces correspond to G,d,u,s,c,b,t
  double f[7];

  std::string path = rpa->gen.Variable("SHERPA_SHARE_PATH") + "/SALGrid";
  char buffer[1024];
  char *err = getcwd(buffer, 1024);
  if (chdir(path.c_str()) != 0 || err == nullptr)
    msg_Error() << "Error in SALph_Fortran_Interface.C " << std::endl
                << "   path " << path << " not found " << std::endl;
  SALPDF(x, Q2, f);
  if (chdir(buffer) != 0)
    msg_Error() << "Error in SALph_Fortran_Interface.C " << std::endl
                << "   path " << path << " not found." << std::endl;

  m_g = x * f[0];
  m_d = x * f[1];
  m_u = x * f[2];
  m_s = x * f[3];
  m_c = x * f[4];
  m_b = x * f[5];
  m_t = x * f[6];
}

double SALph_CPP_Interface::GetXPDF(const ATOOLS::Flavour &infl) {
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

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"), 0);

  return m_rescale * value;
}

double SALph_CPP_Interface::GetXPDF(const kf_code &kf, bool anti) {
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

  value *= MODEL::s_model->ScalarFunction(std::string("alpha_QED"), 0);

  return m_rescale * value;
}

DECLARE_PDF_GETTER(SALph_Getter);

PDF_Base *SALph_Getter::operator()(const Parameter_Type &args) const {
  if (!args.m_bunch.IsPhoton())
    return NULL;
  return new SALph_CPP_Interface(args.m_bunch);
}

void SALph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
  str << "SAL photon PDF, see Eur.Phys.J.C 45 (2006) 633-641 (hep-ph/0504003)";
}

SALph_Getter *p_get_sal;

extern "C" void InitPDFLib() { p_get_sal = new SALph_Getter("SAL"); }

extern "C" void ExitPDFLib() { delete p_get_sal; }
