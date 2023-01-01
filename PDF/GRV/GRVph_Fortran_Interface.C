#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PDF/GRV/GRVph_Fortran_Interface.H"
#include <unistd.h>

using namespace PDF;
using namespace ATOOLS;

extern "C" {
void grvglo_(float &, float &, float &, float &, float &, float &, float &,
             float &);
void grvgho_(float &, float &, float &, float &, float &, float &, float &,
             float &);
}

GRVph_Fortran_Interface::GRVph_Fortran_Interface(const ATOOLS::Flavour _bunch,
                                                 const std::string _set)
    : Photon_PDF_Base(_bunch, _set, 5) {
  m_xmin = 1.e-5;
  m_xmax = 1.;
  m_q2min = .25;
  m_q2max = 1.e6;

  rpa->gen.AddCitation(1, "The GRV photon PDF is published under "
                          "\\cite{Gluck:1991jc} and \\cite{Gluck:1991ee}.");
}

PDF_Base *GRVph_Fortran_Interface::GetCopy() {
  return new GRVph_Fortran_Interface(m_bunch, m_set);
}

void GRVph_Fortran_Interface::CalculateSpec(const double &_x,
                                            const double &_Q2) {
  if (m_include_photon_in_photon)
    m_ph = GetPhotonCoefficient(_x, _Q2);

  float x = _x / m_rescale, Q2 = _Q2;
  float pdf[m_nf + 1];

  if (m_set == std::string("GRVLO"))
    grvglo_(x, Q2, pdf[0], pdf[1], pdf[2], pdf[3], pdf[4], pdf[5]);
  else if (m_set == std::string("GRVHO"))
    grvgho_(x, Q2, pdf[0], pdf[1], pdf[2], pdf[3], pdf[4], pdf[5]);
  else
    msg_Error() << "Error in GRVph_Fortran_Interface.C " << std::endl
                << "   path " << m_set << " not found." << std::endl;

  // Adapt from GRV to Sherpa convention
  double alphaem = MODEL::s_model->ScalarFunction(std::string("alpha_QED"), 0);
  for (int i = 0; i < m_nf + 1; ++i) {
    pdf[i] *= alphaem;
  }

  m_u = pdf[0];
  m_d = pdf[1];
  m_s = pdf[2];
  m_c = pdf[3];
  m_b = pdf[4];
  m_g = pdf[5];
}

DECLARE_PDF_GETTER(GRVph_Getter);

PDF_Base *GRVph_Getter::operator()(const Parameter_Type &args) const {
  if (!args.m_bunch.IsPhoton())
    return NULL;
  return new GRVph_Fortran_Interface(args.m_bunch, args.m_set);
}

void GRVph_Getter::PrintInfo(std::ostream &str, const size_t width) const {
  str << "GRV photon PDF library, see PRD45(1992)3986 and PRD46(1992)1973 \n"
      << "The two sets are \n"
      << " - GRVLO \n"
      << " - GRVHO \n";
}

GRVph_Getter *p_get_grv[2];

extern "C" void InitPDFLib() {
  p_get_grv[0] = new GRVph_Getter("GRVLO");
  p_get_grv[1] = new GRVph_Getter("GRVHO");
}

extern "C" void ExitPDFLib() {
  for (int i = 0; i < 2; ++i) {
    delete p_get_grv[i];
  }
}
