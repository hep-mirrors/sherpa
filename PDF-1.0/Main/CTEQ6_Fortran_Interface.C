#include "CTEQ6_Fortran_Interface.H"
#include "Run_Parameter.H"
#include "Message.H"
#include <unistd.h> 

using namespace PDF;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;

extern "C" {
    void    ctq6initset_(int &);
    double  ctq6evolve_(int &,double &, double &);
}

CTEQ6_Fortran_Interface::CTEQ6_Fortran_Interface(const APHYTOOLS::Flavour _bunch,
						   const std::string _set,const int _member,
						   const std::string _path) :
  m_set(_set), m_member(_member), m_path(_path), m_anti(1) 
{

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

PDF_Base * CTEQ6_Fortran_Interface::GetCopy()
{
  return new CTEQ6_Fortran_Interface(m_bunch,m_set,m_member,m_path);
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

void CTEQ6_Fortran_Interface::Calculate(const double _x, const double _Q2) 
{
  double x = _x, Q = sqrt(_Q2);
  
  int j;
  for (int i=0;i<11;i++) {
    j       = 5-i;
    m_f[i]  = ctq6evolve_(j,x,Q)*x;
  }
}

double CTEQ6_Fortran_Interface::GetXPDF(const APHYTOOLS::Flavour & infl) 
{
  if (infl == Flavour(kf::gluon)) return m_f[5];
  if (infl.Kfcode()==2)   return m_f[5-m_anti*int(infl)/2];     // +/- 1
  if (infl.Kfcode()==1)   return m_f[5-m_anti*int(infl)*2];     // +/- 2
  if (infl.Kfcode()==3)   return m_f[5-m_anti*int(infl)];   // +/- 3
  if (infl.Kfcode()==4)   return m_f[5-m_anti*int(infl)];   // +/- 4
  if (infl.Kfcode()==5)   return m_f[5-m_anti*int(infl)];       // +/- 5
  return 0.;
}

