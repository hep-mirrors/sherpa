#include "LHAPDF_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"

using namespace PDF;
using namespace ATOOLS;


extern "C" {
  void   lhapdreset_();
  void   lhapdfinitset_(const char *);
  void   lhapdfinit_(int &);
  void   lhapdfevolve_(double &,double &,double *);
  double lhapdfalphas_(double &);
}

LHAPDF_Fortran_Interface::LHAPDF_Fortran_Interface(const ATOOLS::Flavour _bunch,
						   const std::string _set,const int _member,
						   const std::string _path, bool & initlhapdf) :
  m_set(_set), m_member(_member), m_path(_path), m_anti(1)
{
  m_bunch = _bunch; 
  if (m_bunch==Flavour(kf::p_plus).Bar()) m_anti=-1;
  if (!initlhapdf) {
    initlhapdf = 1;
    std::string full = m_path+std::string("/")+m_set+std::string(".LHpdf");
    const char * help;
    help = full.c_str();
    lhapdfinitset_(help);
    lhapdfinit_(m_member);
  }

  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour(kf::code(i)));
    m_partons.push_back(Flavour(kf::code(i)).Bar());
  }
  m_partons.push_back(Flavour(kf::gluon));
  m_partons.push_back(Flavour(kf::jet));
  m_partons.push_back(Flavour(kf::quark));
  m_partons.push_back(Flavour(kf::quark).Bar());                               
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
  double scale = 91.2;
  msg.Out()<<" LHAPDF : "<<m_set<<" / "<<m_member<<std::endl
	   <<"          alpha_s(MZ) = "<<lhapdfalphas_(scale)<<std::endl;
}

void LHAPDF_Fortran_Interface::Calculate(const double _x, const double _Q2) {
  double x = _x, Q = sqrt(_Q2);
  lhapdfevolve_(x,Q,m_f);
}

double LHAPDF_Fortran_Interface::GetXPDF(const ATOOLS::Flavour & infl) {
  if (infl == Flavour(kf::gluon)) return m_f[6];
  int kfc = m_anti*int(infl);
  return m_f[6+kfc];
}

