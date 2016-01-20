#include "PDF/MRST/PDF_MRST99.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace std;
using namespace PDF;
using namespace ATOOLS;
c_mrst * PDF_MRST99::p_proton = NULL;


PDF_MRST99::PDF_MRST99(const ATOOLS::Flavour _bunch,
                       const int _member) :
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/MRST99Grid")
{
  m_set="MRST99";
  // default set is member=1
  m_member=(_member==0?1:_member);
  if ((m_member<1)||(m_member>12)) {
    msg_Error()<<"Error in "<<METHOD<<" : Wrong member : "<<m_member<<std::endl
	       <<"    will continue with set 1."<<std::endl;
    m_member  = 1;
  }
  m_type="MRST99["+ATOOLS::ToString(m_member)+"]";
  m_bunch  = _bunch;
  m_anti   = 1;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti = -1;

  if (p_proton==NULL) p_proton = new c_mrst(m_path);

  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());

  m_xmin=MRST99::xmin;
  m_xmax=MRST99::xmax;
  m_q2min=MRST99::qsqmin;
  m_q2max=MRST99::qsqmax;
  m_nf=5;
}


PDF_Base *PDF_MRST99::GetCopy() 
{
  PDF_Base *copy = new PDF_MRST99(m_bunch,m_member);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MRST99::CalculateSpec(const double& x,const double& Q2)
{
  m_overscaled=false;
  if (x/m_rescale>m_xmax || m_rescale<0.) {
    m_overscaled=true;
    return;
  }
  p_proton->mrst99(x/m_rescale,Q2,m_member);
  m_content = p_proton->cont;
}


double PDF_MRST99::GetXPDF(const ATOOLS::Flavour& infl)
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

double PDF_MRST99::GetXPDF(const kf_code& kf, bool anti)
{
  if (m_overscaled) return 0.;
  int kfc=m_anti*(anti?-kf:kf);
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

DECLARE_PDF_GETTER(MRST99_Getter);

PDF_Base *MRST99_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new PDF_MRST99(args.m_bunch,args.m_member);
}

void MRST99_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"MRST 1999 fit\n"
     <<std::string(width+4,' ')<<"see hep-ph/9907231";
}

MRST99_Getter *p_get_mrst99;

extern "C" void InitPDFLib()
{
  p_get_mrst99 = new MRST99_Getter("MRST99");
}

extern "C" void ExitPDFLib()
{
  delete p_get_mrst99;
}
