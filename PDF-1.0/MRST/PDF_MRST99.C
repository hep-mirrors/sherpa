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
    msg.Error()<<"Error in PDF_MRST99::PDF_MRST99 : Wrong set : "<<m_set<<std::endl
	       <<"    will continue with set 1."<<std::endl;
    m_set  = 1;
  }
  m_type=std::string("MRST99")+ATOOLS::ToString(m_set);
  m_bunch  = _bunch;
  m_anti   = 1;
  if (m_bunch==Flavour(kf::p_plus).Bar()) m_anti = -1;

  if (p_proton==NULL) p_proton = new c_mrst(m_path);

  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour(kf::code(i)));
    m_partons.push_back(Flavour(kf::code(i)).Bar());
  }
  m_partons.push_back(Flavour(kf::gluon));
  m_partons.push_back(Flavour(kf::jet));
  m_partons.push_back(Flavour(kf::quark));
  m_partons.push_back(Flavour(kf::quark).Bar());

  m_xmin=MRST99::xmin;
  m_xmax=MRST99::xmax;
  m_q2min=MRST99::qsqmin;
  m_q2max=MRST99::qsqmax;
};


PDF_Base *PDF_MRST99::GetCopy() 
{
  PDF_Base *copy = new PDF_MRST99(m_bunch,m_set,m_path);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MRST99::Output() {
  msg.Out()<<" internal PDF_MRST99 : "<<endl;
  switch(m_set) {
  case 1:  msg.Out()<<"          central gluon, central alpha_s = "<<0.1175<<endl;break;
  case 2:  msg.Out()<<"          higher gluon, central alpha_s = "<<0.1175<<endl;break;
  case 3:  msg.Out()<<"          lower gluon, central alpha_s = "<<0.1175<<endl;break;
  case 4:  msg.Out()<<"          lower alpha_s = "<<0.1125<<endl;break;
  case 5:  msg.Out()<<"          higher alpha_s = "<<0.1225<<endl;break;
  case 6:  msg.Out()<<"          quarks up, alpha_s = "<<0.1178<<endl;break;
  case 7:  msg.Out()<<"          quarks down, alpha_s = "<<0.1171<<endl;break;
  case 8:  msg.Out()<<"          strange up, central alpha_s = "<<0.1175<<endl;break;
  case 9:  msg.Out()<<"          strange down, central alpha_s = "<<0.1175<<endl;break;
  case 10: msg.Out()<<"          charm up, central alpha_s = "<<0.1175<<endl;break;
  case 11: msg.Out()<<"          charm down, central alpha_s = "<<0.1175<<endl;break;
  case 12: msg.Out()<<"          larger d/u, central alpha_s = "<<0.1175<<endl;break;
  default : msg.Out()<<"Error in  PDF_MRST99::Output() : Set = "<<m_set<<endl;
  }
}

void PDF_MRST99::Calculate(double x,double z,double kp2,double Q2) 
{
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
  case  ATOOLS::kf::d : return m_rescale*(m_content.dnv + m_content.dsea);
  case -ATOOLS::kf::d : return m_rescale*m_content.dsea; 
  case  ATOOLS::kf::u : return m_rescale*(m_content.upv + m_content.usea);
  case -ATOOLS::kf::u : return m_rescale*m_content.usea; 
  case  ATOOLS::kf::s :
  case -ATOOLS::kf::s : return m_rescale*m_content.str;
  case  ATOOLS::kf::c : 
  case -ATOOLS::kf::c : return m_rescale*m_content.chm;
  case  ATOOLS::kf::b : 
  case -ATOOLS::kf::b : return m_rescale*m_content.bot;
  case ATOOLS::kf::gluon : 
  // pseudo anti gluon for anti-proton
  case -ATOOLS::kf::gluon :return m_rescale*m_content.glu; 
  default: return 0.;
  }
}

void PDF_MRST99::AssignKeys(ATOOLS::Integration_Info *const info)
{
}
