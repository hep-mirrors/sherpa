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
    msg.Tracking()<<"Try to initialize cteq6 PDF set."<<endl
		  <<"  Set = "<<m_set<<" v "<<m_member<<" for "<<m_bunch<<endl;

    int iset = 0;
    
    if (m_set==std::string("cteq6m"))  iset = 1;
    if (m_set==std::string("cteq6d"))  iset = 2;
    if (m_set==std::string("cteq6l"))  iset = 3;
    if (m_set==std::string("cteq6l1")) iset = 4;

    char buffer[1024];
    char * err = getcwd(buffer,1024);
    if (err==NULL) {
      cout<<" something wrong in CTEQ6_Fortran_Interface.C "<<endl;
    }
    int stat=chdir(m_path.c_str());
    ctq6initset_(iset);
    if (stat==0) {
      chdir(buffer);
    }
    else {
      cout<<" path "<<m_path<<" not found "<<endl;
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



double CTEQ6_Fortran_Interface::AlphaSPDF(double scale2) {}


void CTEQ6_Fortran_Interface::Output() {}

void CTEQ6_Fortran_Interface::Calculate(const double _x, const double _Q2) {
  
    //cout<<"In Calculate "<<_x<<" "<<_Q2<<endl;

    double x = _x, Q = sqrt(_Q2);

    //cout<<"In Calculate "<<x<<" "<<Q<<endl;
       
    int j;
    for (int i=0;i<11;i++) {
	j       = 5-i;
	m_f[i]  = ctq6evolve_(j,x,Q)*x;
	//cout<<"m_f["<<i<<"] = "<<m_f[i]<<endl;
    }
}

double CTEQ6_Fortran_Interface::GetXPDF(const APHYTOOLS::Flavour & infl) {

    if (infl == Flavour(kf::gluon)) return m_f[5];
    if (infl.Kfcode()==2)   return m_f[5-m_anti*int(infl)/2];     // +/- 1
    if (infl.Kfcode()==1)   return m_f[5-m_anti*int(infl)*2];     // +/- 2
    if (infl.Kfcode()==4)   return m_f[5-m_anti*int(infl)/4*3];   // +/- 3
    if (infl.Kfcode()==3)   return m_f[5-m_anti*int(infl)/3*4];   // +/- 4
    if (infl.Kfcode()==5)   return m_f[5-m_anti*int(infl)];       // +/- 5
  
}

