#include "PDF_MRST99.H"
#include "Message.H"

using namespace std;
using namespace PDF;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
c_mrst * PDF_MRST99::p_proton = NULL;


PDF_MRST99::PDF_MRST99(const APHYTOOLS::Flavour _bunch,
		       const int _set,const std::string _path) :
  m_set(_set), m_path(_path), p_proton(NULL)
{
  if ((m_set<1)||(m_set>12)) {
    std::cout <<"Error in PDF_MRST99::PDF_MRST99 : Wrong set : "<<m_set<<std::endl
	      <<"    will continue with set 1."<<std::endl;
    m_set  = 1;
  }
  m_bunch  = _bunch;
  m_anti   = 1;
  if (m_bunch==Flavour(kf::p_plus).Bar()) m_anti = -1;

  msg.Debugging()<<"Initialise MRST99/"<<m_set<<" from "<<m_path<<endl;
  if (p_proton==NULL) p_proton = new c_mrst(m_path);

  for (int i=1;i<6;i++) {
    m_partons.push_back(Flavour(kf::code(i)));
    m_partons.push_back(Flavour(kf::code(i)).Bar());
  }
  m_partons.push_back(Flavour(kf::gluon));
  m_partons.push_back(Flavour(kf::jet));
  m_partons.push_back(Flavour(kf::quark));
  m_partons.push_back(Flavour(kf::quark).Bar());
};


PDF_Base * PDF_MRST99::GetCopy() {
  return new PDF_MRST99(m_bunch,m_set,m_path);
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

void PDF_MRST99::Calculate(const double x, const double Q2) {
  p_proton->mrst99(x,Q2,m_set);
  m_content = p_proton->cont;
}


double PDF_MRST99::GetXPDF(const APHYTOOLS::Flavour & infl) {
  int kfc = m_anti*infl.Kfcode();
  switch (kfc) {
  case  APHYTOOLS::kf::d : return (m_content.dnv + m_content.dsea);
  case -APHYTOOLS::kf::d : return m_content.dsea; 
  case  APHYTOOLS::kf::u : return (m_content.upv + m_content.usea);
  case -APHYTOOLS::kf::u : return m_content.usea; 
  case  APHYTOOLS::kf::s :
  case -APHYTOOLS::kf::s : return m_content.str;
  case  APHYTOOLS::kf::c : 
  case -APHYTOOLS::kf::c : return m_content.chm;
  case  APHYTOOLS::kf::b : 
  case -APHYTOOLS::kf::b : return m_content.bot;
  case APHYTOOLS::kf::gluon : 
  case -APHYTOOLS::kf::gluon :return m_content.glu; // pseudo anti gluon for anti-proton
  default: return 0.;
  }
}
