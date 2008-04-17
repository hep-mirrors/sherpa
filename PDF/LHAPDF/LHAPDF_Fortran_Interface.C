#include "LHAPDF_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace PDF;
using namespace ATOOLS;


extern "C" {
  void   lhapdreset_();
  void   lhapdfinitset_(const char *, int len);
  void   lhapdfinit_(int &);
  void   lhapdfevolve_(double &,double &,double *);
  double lhapdfalphas_(double &);
  void   lhapdfgetdesc_();
}

LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour _bunch,
						   const std::string _set,const int _member,
						   const std::string _path, bool & initlhapdf) :
  m_set(_set), m_member(_member), m_path(_path), m_anti(1)
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
    std::string full = m_path+std::string("/")+m_set;
    const char * help;
    help = full.c_str();
    lhapdfinitset_(help, strlen(help));
    lhapdfinit_(m_member);
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
  return new LHAPDF_Fortran_Interface(m_bunch,m_set,m_member,m_path,init);
}


double LHAPDF_Fortran_Interface::AlphaSPDF(double scale2) {
  double scale = sqrt(scale2);
  double as    = lhapdfalphas_(scale);
  return as;
}

void LHAPDF_Fortran_Interface::Output() {
  double scale = Flavour(kf_Z).Mass();
  msg_Out()<<" LHAPDF : "<<m_set<<" / "<<m_member<<std::endl
	   <<"          alpha_s(MZ) = "<<lhapdfalphas_(scale)<<std::endl;
  lhapdfgetdesc_();
}

void LHAPDF_Fortran_Interface::Calculate(double x,double z,double kp2,double Q2) {
  double Q = sqrt(Q2);
  lhapdfevolve_(x,Q,m_f);
}

double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) {
  if (infl == Flavour(kf_gluon)) return m_f[6];
  int kfc = m_anti*int(infl);
  if (kfc<-6 || kfc>6) {
    msg_Out()<<"WARNING in LHAPDF_Fortran_Interface::GetXPDF("<<infl<<") not supported by this PDF!"<<std::endl;
    return 0.;
  }
  return m_f[6+kfc];
}

void LHAPDF_Fortran_Interface::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
