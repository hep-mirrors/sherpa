#include "PDF_MRST99.H"

#include "Message.H"
#include "MyStrStream.H"

using namespace std;
using namespace PDF;
using namespace ATOOLS;
c_mrst * PDF_MRST99::p_proton = NULL;


PDF_MRST99::PDF_MRST99(const ATOOLS::Flavour _bunch,
		       const int _set,const std::string _path) :
  m_set(_set), m_path(_path)//, p_proton(NULL)
{
  if ((m_set<1)||(m_set>12)) {
    msg_Error()<<"Error in PDF_MRST99::PDF_MRST99 : Wrong set : "<<m_set<<std::endl
	       <<"    will continue with set 1."<<std::endl;
    m_set  = 1;
  }
  m_type=std::string("MRST99")+ATOOLS::ToString(m_set);
  m_bunch  = _bunch;
  m_anti   = 1;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti = -1;

  if (p_proton==NULL) p_proton = new c_mrst(m_path);

  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour((kf_code)(i)));
    m_partons.push_back(Flavour((kf_code)(i)).Bar());
  }
  m_partons.push_back(Flavour(kf_gluon));
  m_partons.push_back(Flavour(kf_jet));
  m_partons.push_back(Flavour(kf_quark));
  m_partons.push_back(Flavour(kf_quark).Bar());

  m_xmin=MRST99::xmin;
  m_xmax=MRST99::xmax;
  m_q2min=MRST99::qsqmin;
  m_q2max=MRST99::qsqmax;
}


PDF_Base *PDF_MRST99::GetCopy() 
{
  PDF_Base *copy = new PDF_MRST99(m_bunch,m_set,m_path);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MRST99::Output() {
  msg_Out()<<" internal PDF_MRST99 : "<<endl;
  switch(m_set) {
  case 1:  msg_Out()<<"          central gluon, central alpha_s = "<<0.1175<<endl;break;
  case 2:  msg_Out()<<"          higher gluon, central alpha_s = "<<0.1175<<endl;break;
  case 3:  msg_Out()<<"          lower gluon, central alpha_s = "<<0.1175<<endl;break;
  case 4:  msg_Out()<<"          lower alpha_s = "<<0.1125<<endl;break;
  case 5:  msg_Out()<<"          higher alpha_s = "<<0.1225<<endl;break;
  case 6:  msg_Out()<<"          quarks up, alpha_s = "<<0.1178<<endl;break;
  case 7:  msg_Out()<<"          quarks down, alpha_s = "<<0.1171<<endl;break;
  case 8:  msg_Out()<<"          strange up, central alpha_s = "<<0.1175<<endl;break;
  case 9:  msg_Out()<<"          strange down, central alpha_s = "<<0.1175<<endl;break;
  case 10: msg_Out()<<"          charm up, central alpha_s = "<<0.1175<<endl;break;
  case 11: msg_Out()<<"          charm down, central alpha_s = "<<0.1175<<endl;break;
  case 12: msg_Out()<<"          larger d/u, central alpha_s = "<<0.1175<<endl;break;
  default : msg_Out()<<"Error in  PDF_MRST99::Output() : Set = "<<m_set<<endl;
  }
}

void PDF_MRST99::Calculate(double x,double _Q2) 
{
  double Q2(_Q2*m_fac_scale_factor);
  m_overscaled=false;
  if (x/m_rescale>m_xmax || m_rescale<0.) {
    m_overscaled=true;
    return;
  }
  p_proton->mrst99(x/m_rescale,Q2,m_set);
  m_content = p_proton->cont;
}


double PDF_MRST99::GetXPDF(const ATOOLS::Flavour infl) 
{
  if (m_overscaled) return 0.;
  int kfc=m_anti*int(infl);
  switch (kfc) {
  case  kf_d : return m_rescale*(m_content.dnv + m_content.dsea);
  case -kf_d : return m_rescale*m_content.dsea; 
  case  kf_u : return m_rescale*(m_content.upv + m_content.usea);
  case -kf_u : return m_rescale*m_content.usea; 
  case  kf_s :
  case -kf_s : return m_rescale*m_content.str;
  case  kf_c : 
  case -kf_c : return m_rescale*m_content.chm;
  case  kf_b : 
  case -kf_b : return m_rescale*m_content.bot;
  case kf_gluon : 
  // pseudo anti gluon for anti-proton
  case -kf_gluon :return m_rescale*m_content.glu; 
  default: return 0.;
  }
}

void PDF_MRST99::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
