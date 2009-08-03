#include "PDF/CTEQ/CTEQ6_Fortran_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <unistd.h> 

using namespace PDF;
using namespace ATOOLS;

extern "C" {
  void    ctq6initset_(int &);
  double  ctq6evolve_(int &,double &, double &);
  void    errmsg_();
}

void errmsg_() {
  CTEQ6_Fortran_Interface::Error();
}


CTEQ6_Fortran_Interface::CTEQ6_Fortran_Interface(const ATOOLS::Flavour _bunch,
						 const std::string _set,const int _member,
						 const std::string _path):
  m_set(_set), m_path(_path), m_member(_member), m_anti(1) 
{
  m_xmin=1.e-6;
  m_xmax=1.;
  m_q2min=.5;
  m_q2max=1.e12;

  m_type=m_set;
  m_bunch = _bunch;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  int iset = 0;
  
  if (m_set==std::string("cteq6m"))  iset = 1;
  if (m_set==std::string("cteq6d"))  iset = 2;
  if (m_set==std::string("cteq6l"))  iset = 3;
  if (m_set==std::string("cteq6l1")) iset = 4;
  
  char buffer[1024];
  char * err = getcwd(buffer,1024);
  if (err==NULL) {
    msg_Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl;
  }
  int stat=chdir(m_path.c_str());
  ctq6initset_(iset);
  if (stat==0) {
    stat=chdir(buffer);
  }
  else {
    msg_Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl
	       <<"   path "<<m_path<<" not found "<<std::endl;
  }
  
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());                               
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

void CTEQ6_Fortran_Interface::CalculateSpec(double x,double _Q2) 
{
  for (size_t i=0;i<11;++i) m_calculated[i]=false;
  m_x=x/m_rescale;
  m_Q=sqrt(_Q2);
  if(m_Q*m_Q<m_q2min) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<m_Q
	       <<" < "<<sqrt(m_q2min)<<". Set Q -> "
	       <<sqrt(m_q2min)<<"."<<std::endl;
    m_Q=1.000001*sqrt(m_q2min);
  }
  if(m_Q*m_Q>m_q2max) {
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<m_Q
	       <<" > "<<sqrt(m_q2max)<<". Set Q -> "
	       <<sqrt(m_q2max)<<"."<<std::endl;
    m_Q=0.999999*sqrt(m_q2max);
  }
  if(m_x<m_xmin) {
    msg_Error()<<METHOD<<"("<<this<<"): x = "<<m_x<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Set x -> "<<m_xmin<<".\n";
    m_x=m_xmin;
  }
  if(m_x>m_xmax) {
    msg_Error()<<METHOD<<"("<<this<<"): x = "<<m_x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Set x -> "<<m_xmax<<".\n";
    m_x=m_xmax;
  }
}

double CTEQ6_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) 
{
  if ((m_x>m_xmax && m_rescale<1.) || m_rescale<0.) return 0.;
  int cteqindex;
  switch (infl.Kfcode()) {
  case kf_gluon: cteqindex=0;                  break;
  case kf_d:     cteqindex=m_anti*int(infl)*2; break;
  case kf_u:     cteqindex=m_anti*int(infl)/2; break;
  default:                cteqindex=m_anti*int(infl);   break;
  }
  if (!m_calculated[5-cteqindex]) {
    m_f[5-cteqindex]=ctq6evolve_(cteqindex,m_x,m_Q)*m_x; 
    m_calculated[5-cteqindex]=true;
  }
  return m_rescale*m_f[5-cteqindex];     
}

void CTEQ6_Fortran_Interface::Error()
{
  THROW(critical_error,"Cteq6Pdf called ERRORMSG ");
}

DECLARE_PDF_GETTER(CTEQ6_Getter);

PDF_Base *CTEQ6_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  int mode=args.p_read->GetValue<int>("PDF_SET_VERSION",1);
  return new CTEQ6_Fortran_Interface(args.m_bunch,m_key,mode,args.m_path);
}

void CTEQ6_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"CTEQ 6 fit, see hep-ph/0201195";
}

CTEQ6_Getter *p_get[3];

extern "C" void InitPDFLib(const std::string &path)
{
  p_get[0] = new CTEQ6_Getter("cteq6l1");
  p_get[1] = new CTEQ6_Getter("cteq6l");
  p_get[2] = new CTEQ6_Getter("cteq6m");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<3;++i) delete p_get[i];
}
