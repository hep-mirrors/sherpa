#include "LHAPDF_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "CXXFLAGS_PACKAGES.H"
#include <cstring>

using namespace PDF;
using namespace ATOOLS;

#ifdef LHAPDF__NATIVE__WRAPPER
#include "LHAPDF/LHAPDF.h"
#else
extern "C" {
  void   lhapdreset_();
  void   lhapdfinitset_(const char *, int len);
  void   lhapdfinitsetbyname_(const char *, int len);
  void   lhapdfinit_(int &);
  void   lhapdfevolve_(double &,double &,double *);
  double lhapdfalphas_(double &);
  void   lhapdfgetdesc_();
}
#endif

LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour _bunch,
						   const std::string _set,const int _member,
						   bool & initlhapdf) :
  m_set(_set), m_member(_member), m_anti(1)
{
  m_xmin=0.;
  m_xmax=1.;
  m_q2min=1.;
  m_q2max=1.e12;
  m_type="LHA["+m_set+"]";

  m_bunch = _bunch; 
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  if (!initlhapdf) {
    initlhapdf = true;
#ifdef LHAPDF__NATIVE__WRAPPER
    LHAPDF::initPDFByName(m_set, m_member);
#else
    std::string full = m_set;
    const char * help;
    help = full.c_str();
    lhapdfinitsetbyname_(help, strlen(help));
    lhapdfinit_(m_member);
#endif
  }

  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour((kf_code)(i)));
    m_partons.push_back(Flavour((kf_code)(i)).Bar());
  }
  m_partons.push_back(Flavour(kf_gluon));
  m_partons.push_back(Flavour(kf_jet));
  m_partons.push_back(Flavour(kf_quark));
  m_partons.push_back(Flavour(kf_quark).Bar());                               
}

PDF_Base * LHAPDF_Fortran_Interface::GetCopy() 
{
  bool init=0;
  return new LHAPDF_Fortran_Interface(m_bunch,m_set,m_member,init);
}


double LHAPDF_Fortran_Interface::AlphaSPDF(double scale2) {
  double scale = sqrt(scale2);
#ifdef LHAPDF__NATIVE__WRAPPER
  double as    = LHAPDF::alphasPDF(scale);
#else
  double as    = lhapdfalphas_(scale);
#endif
  return as;
}

void LHAPDF_Fortran_Interface::Output() {
  double scale = Flavour(kf_Z).Mass();
#ifdef LHAPDF__NATIVE__WRAPPER
  double as    = LHAPDF::alphasPDF(scale);
#else
  double as    = lhapdfalphas_(scale);
#endif
  msg_Out()<<" LHAPDF : "<<m_set<<" / "<<m_member<<std::endl
	   <<"          alpha_s(MZ) = "<<as<<std::endl;
#ifdef LHAPDF__NATIVE__WRAPPER
  LHAPDF::getDescription();
#else
  lhapdfgetdesc_();
#endif
}

void LHAPDF_Fortran_Interface::Calculate(double x,double Q2) {
  double Q = sqrt(Q2*m_fac_scale_factor);
#ifdef LHAPDF__NATIVE__WRAPPER
  m_fv=LHAPDF::xfx(x,Q);
#else
  lhapdfevolve_(x,Q,m_f);
#endif
}

double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) {
  int kfc = m_anti*int(infl);
  if (infl == Flavour(kf_gluon)) kfc=0;
  if (kfc<-6 || kfc>6) {
    msg_Out()<<"WARNING in LHAPDF_Fortran_Interface::GetXPDF("<<infl<<") not supported by this PDF!"<<std::endl;
    return 0.;
  }
#ifdef LHAPDF__NATIVE__WRAPPER
  return m_fv[6+kfc];
#else
  return m_f[6+kfc];
#endif
}

void LHAPDF_Fortran_Interface::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
