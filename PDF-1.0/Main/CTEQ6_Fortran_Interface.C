#include "CTEQ6_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include <unistd.h> 

using namespace PDF;
using namespace ATOOLS;

extern "C" {
    void    ctq6initset_(int &);
    double  ctq6evolve_(int &,double &, double &);
}

CTEQ6_Fortran_Interface::CTEQ6_Fortran_Interface(const ATOOLS::Flavour _bunch,
						 const std::string _set,const int _member,
						 const std::string _path):
  m_set(_set), m_path(_path), m_member(_member), m_anti(1) 
{
  m_xmin=0.;
  m_xmax=1.;
  m_q2min=.5;
  m_q2max=1.e12;

  m_type=m_set;
  m_bunch = _bunch;
  if (m_bunch==Flavour(kf::p_plus).Bar()) m_anti=-1;
  int iset = 0;
  
  if (m_set==std::string("cteq6m"))  iset = 1;
  if (m_set==std::string("cteq6d"))  iset = 2;
  if (m_set==std::string("cteq6l"))  iset = 3;
  if (m_set==std::string("cteq6l1")) iset = 4;
  
  char buffer[1024];
  char * err = getcwd(buffer,1024);
  if (err==NULL) {
    msg.Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl;
  }
  int stat=chdir(m_path.c_str());
  ctq6initset_(iset);
  if (stat==0) {
    chdir(buffer);
  }
  else {
    msg.Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl
	       <<"   path "<<m_path<<" not found "<<std::endl;
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

PDF_Base *CTEQ6_Fortran_Interface::GetCopy()
{
  PDF_Base *copy = new CTEQ6_Fortran_Interface(m_bunch,m_set,m_member,m_path);
  m_copies.push_back(copy);
  return copy;
}



double CTEQ6_Fortran_Interface::AlphaSPDF(double scale2) 
{
  //  ** ALL fits are obtained by using the same coupling strength
  //   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
  //   which uses the LO running \alpha_s and its value determined from the fit.

  double asmz = 0.;
  if (m_set==std::string("cteq6m"))  asmz = 0.118;
  if (m_set==std::string("cteq6d"))  asmz = 0.118;
  if (m_set==std::string("cteq6l"))  asmz = 0.118;
  if (m_set==std::string("cteq6l1")) asmz = 0.130;

  return asmz;
}

void CTEQ6_Fortran_Interface::Output() {}

void CTEQ6_Fortran_Interface::Calculate(double x,double z,double kp2,double _Q2) 
{
  for (size_t i=0;i<11;++i) m_calculated[i]=false;
  m_x=x/m_rescale;
  m_Q=sqrt(_Q2);
}

double CTEQ6_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) 
{
  if ((m_x>m_xmax && m_rescale<1.) || m_rescale<0.) return 0.;
  int cteqindex;
  switch (infl.Kfcode()) {
  case ATOOLS::kf::gluon: cteqindex=0;                  break;
  case ATOOLS::kf::d:     cteqindex=m_anti*int(infl)*2; break;
  case ATOOLS::kf::u:     cteqindex=m_anti*int(infl)/2; break;
  default:                cteqindex=m_anti*int(infl);   break;
  }
  if (!m_calculated[5-cteqindex]) {
    m_f[5-cteqindex]=ctq6evolve_(cteqindex,m_x,m_Q)*m_x; 
    m_calculated[5-cteqindex]=true;
  }
  return m_rescale*m_f[5-cteqindex];     
}

void CTEQ6_Fortran_Interface::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
