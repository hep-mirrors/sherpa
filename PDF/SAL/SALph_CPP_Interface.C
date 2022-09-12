#include "SALph_CPP_Interface.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include <unistd.h>

#include "sal.h"

using namespace PDF;
using namespace ATOOLS;

SALph_CPP_Interface::SALph_CPP_Interface(const ATOOLS::Flavour _bunch)
    : Photon_PDF_Base(_bunch, "SAL", 6) {
  // This is acutally an upper limit, as the x_min is dependent on the flavour
  m_xmin = 1.e-5;
  m_xmax = 0.9999;
  m_q2min = 2.;
  m_q2max = 8.e4;

  rpa->gen.AddCitation(
      1, "The SAL photon PDF is published under \\cite{Slominski:2005bw}.");
}

PDF_Base *SALph_CPP_Interface::GetCopy() {
  return new SALph_CPP_Interface(m_bunch);
}

void SALph_CPP_Interface::CalculateSpec(const double &_x, const double &_Q2) {
  if (m_include_photon_in_photon)
    m_ph = GetPhotonCoefficient(_x, _Q2);

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

  // Adapt from SAL to Sherpa convention
  double alphaem = MODEL::s_model->ScalarFunction(std::string("alpha_QED"), 0);
  for (int i = 0; i < 7; ++i) {
    f[i] *= x * alphaem;
  }

  m_g = x * f[0];
  m_d = x * f[1];
  m_u = x * f[2];
  m_s = x * f[3];
  m_c = x * f[4];
  m_b = x * f[5];
  m_t = x * f[6];
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
