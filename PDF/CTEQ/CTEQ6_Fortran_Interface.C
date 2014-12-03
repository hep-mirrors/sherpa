#include "PDF/CTEQ/CTEQ6_Fortran_Interface.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <unistd.h> 

using namespace PDF;
using namespace ATOOLS;

extern "C" {
  void    setctq6_(int &);
  double  ctq6pdf_(int &,double &, double &);
  void    errmsg_();
}

void errmsg_() {
  CTEQ6_Fortran_Interface::Error();
}


CTEQ6_Fortran_Interface::CTEQ6_Fortran_Interface(const ATOOLS::Flavour _bunch,
						 const std::string _set,const int _member):
  m_set(_set), m_anti(1) 
{
  m_member=_member;
  m_xmin=1.e-6;
  m_xmax=1.;
  m_q2min=.5;
  m_q2max=1.e12;

  m_type=m_set;
  m_bunch = _bunch;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  int iset = 0;
  std::string path = rpa->gen.Variable("SHERPA_SHARE_PATH")+"/";
  
  if (m_set==std::string("cteq6.6m")) {
    iset = 400;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.118;
    path+="CTEQ66Grid";
    m_lhef_number=10550+m_member;
  }
  if (m_set==std::string("cteq6.6a1")) {
    iset = 460;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.125;
    path+="CTEQ66Grid";
  }
  if (m_set==std::string("cteq6.6a2")) {
    iset = 461;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.122;
    path+="CTEQ66Grid";
  }
  if (m_set==std::string("cteq6.6a3")) {
    iset = 462;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.114;
    path+="CTEQ66Grid";
  }
  if (m_set==std::string("cteq6.6a4")) {
    iset = 463;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.112;
    path+="CTEQ66Grid";
  }
  if (m_set==std::string("cteq6m")) {
    iset = 1;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.118;
    path+="CTEQ6Grid";
  }
  if (m_set==std::string("cteq6d")) {
    iset = 2;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.118;
    path+="CTEQ6Grid";
  }
  if (m_set==std::string("cteq6l")) {
    iset = 3;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.117981;
    path+="CTEQ6Grid";
  }
  if (m_set==std::string("cteq6l1")) {
    iset = 4;
    m_asinfo.m_order=0;
    m_asinfo.m_asmz=0.129783;
    path+="CTEQ6Grid";
  }
  if (iset==1 && m_member>0 && m_member<=40) {
    iset=100+m_member;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.118;
    path+="CTEQ6Grid";
  }
  if (iset==400 && m_member>0 && m_member<=44) {
    iset+=m_member;
    m_asinfo.m_order=1;
    m_asinfo.m_asmz=0.118;
    path+="CTEQ6Grid";
  }
  
  char buffer[1024];
  char * err = getcwd(buffer,1024);
  if (err==NULL) {
    msg_Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl;
  }
  int stat=chdir(path.c_str());
  msg_Info()<<METHOD<<"(): Init member "<<iset<<"."<<std::endl;
  setctq6_(iset);
  if (stat==0) {
    stat=chdir(buffer);
  }
  else {
    msg_Error()<<"Error in CTEQ6_Fortran_Interface.C "<<std::endl
	       <<"   path "<<path<<" not found "<<std::endl;
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
  PDF_Base *copy = new CTEQ6_Fortran_Interface(m_bunch,m_set,m_member);
  m_copies.push_back(copy);
  return copy;
}

void CTEQ6_Fortran_Interface::CalculateSpec(double x,double _Q2) 
{
  for (size_t i=0;i<11;++i) m_calculated[i]=false;
  m_x=x;
  m_Q=sqrt(_Q2);
}

double CTEQ6_Fortran_Interface::GetXPDF(const ATOOLS::Flavour infl) 
{
  if (m_x>m_xmax) return 0.;
  int cteqindex;
  switch (infl.Kfcode()) {
  case kf_gluon: cteqindex=0;                  break;
  case kf_d:     cteqindex=m_anti*int(infl)*2; break;
  case kf_u:     cteqindex=m_anti*int(infl)/2; break;
  default:                cteqindex=m_anti*int(infl);   break;
  }
  if (!m_calculated[5-cteqindex]) {
    m_f[5-cteqindex]=ctq6pdf_(cteqindex,m_x,m_Q)*m_x; 
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
  int mode=args.p_read->GetValue<int>("PDF_SET_VERSION",0);
  int ibeam=args.m_ibeam;
  mode=args.p_read->GetValue<int>("PDF_SET_VERSION_"+ToString(ibeam+1),mode);
  return new CTEQ6_Fortran_Interface(args.m_bunch,m_key,mode);
}

void CTEQ6_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"CTEQ 6 fit, see hep-ph/0201195"
     <<" / CTEQ 6.6 fit, see arXiv:0802.0007 [hep-ph]";
}

CTEQ6_Getter *p_get_cteq6[8];

extern "C" void InitPDFLib()
{
  p_get_cteq6[0] = new CTEQ6_Getter("cteq6l1");
  p_get_cteq6[1] = new CTEQ6_Getter("cteq6l");
  p_get_cteq6[2] = new CTEQ6_Getter("cteq6m");
  p_get_cteq6[3] = new CTEQ6_Getter("cteq6.6m");
  p_get_cteq6[4] = new CTEQ6_Getter("cteq6.6a1");
  p_get_cteq6[5] = new CTEQ6_Getter("cteq6.6a2");
  p_get_cteq6[6] = new CTEQ6_Getter("cteq6.6a3");
  p_get_cteq6[7] = new CTEQ6_Getter("cteq6.6a4");
}

extern "C" void ExitPDFLib()
{
  for (int i(0);i<8;++i) delete p_get_cteq6[i];
}
