#include "PDF/MRST/PDF_MRST01LO.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PDF;
using namespace ATOOLS;

extern "C" {
  void mrstlo_(double *x,double *q2,int *mode,
	       double *upv,double *dnv,double *usea,
	       double *dsea,double *str,double *chm,
	       double *bot,double *glu);
}

void mrstlo(double x,double q2,int mode,
	    double &upv,double &dnv,double &usea,
	    double &dsea,double &str,double &chm,
	    double &bot,double &glu)
{
  mrstlo_(&x,&q2,&mode,&upv,&dnv,&usea,&dsea,&str,&chm,&bot,&glu);
}

extern "C" {
extern struct {
  char mfile[128];
} mrinput_;
#define input mrinput_
}

inline void MakeFortranString(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
}

PDF_MRST01LO::PDF_MRST01LO(const ATOOLS::Flavour bunch):
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/MRST01Grid"),
  m_anti(1),
  m_mode(1)
{
  m_type=std::string("MRST01LO");
  m_bunch=bunch;
  if (m_bunch==Flavour(kf_p_plus).Bar()) m_anti=-1;
  for (int i=1;i<6;i++) {
    m_partons.insert(Flavour((kf_code)(i)));
    m_partons.insert(Flavour((kf_code)(i)).Bar());
  }
  m_partons.insert(Flavour(kf_gluon));
  m_partons.insert(Flavour(kf_jet));
  m_partons.insert(Flavour(kf_quark));
  m_partons.insert(Flavour(kf_quark).Bar());
  m_xmin=1.e-5;
  m_xmax=1.;
  m_q2min=1.25;
  m_q2max=1.e7;
  m_nf=5;
  MakeFortranString(input.mfile,m_path+std::string("/lo2002.dat"),128);
}


PDF_Base *PDF_MRST01LO::GetCopy() 
{
  PDF_Base *copy = new PDF_MRST01LO(m_bunch);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MRST01LO::CalculateSpec(const double& x,const double& Q2)
{
  m_overscaled=false;
  if (x/m_rescale>m_xmax || m_rescale<0.) {
    m_overscaled=true;
    return;
  }
  mrstlo(x/m_rescale,Q2,m_mode,p_xpdfv[1],p_xpdfv[0],p_xpdf[1],
	 p_xpdf[0],p_xpdf[2],p_xpdf[3],p_xpdf[4],p_xpdf[5]);
}


double PDF_MRST01LO::GetXPDF(const ATOOLS::Flavour& infl)
{
  if (m_overscaled) return 0.;
  int kfc=m_anti*int(infl);
  switch (kfc) {
  case  kf_d : return m_rescale*(p_xpdfv[0]+p_xpdf[0]);
  case -kf_d : return m_rescale*p_xpdf[0]; 
  case  kf_u : return m_rescale*(p_xpdfv[1]+p_xpdf[1]);
  case -kf_u : return m_rescale*p_xpdf[1]; 
  case  kf_s :
  case -kf_s : return m_rescale*p_xpdf[2];
  case  kf_c : 
  case -kf_c : return m_rescale*p_xpdf[3];
  case  kf_b : 
  case -kf_b : return m_rescale*p_xpdf[4];
  case  kf_gluon : 
  case -kf_gluon :return m_rescale*p_xpdf[5]; 
  default: return 0.;
  }
}

double PDF_MRST01LO::GetXPDF(const kf_code& kf, bool anti)
{
  if (m_overscaled) return 0.;
  int kfc=m_anti*(anti?-kf:kf);
  switch (kfc) {
  case  kf_d : return m_rescale*(p_xpdfv[0]+p_xpdf[0]);
  case -kf_d : return m_rescale*p_xpdf[0];
  case  kf_u : return m_rescale*(p_xpdfv[1]+p_xpdf[1]);
  case -kf_u : return m_rescale*p_xpdf[1];
  case  kf_s :
  case -kf_s : return m_rescale*p_xpdf[2];
  case  kf_c :
  case -kf_c : return m_rescale*p_xpdf[3];
  case  kf_b :
  case -kf_b : return m_rescale*p_xpdf[4];
  case  kf_gluon :
  case -kf_gluon :return m_rescale*p_xpdf[5];
  default: return 0.;
  }
}

DECLARE_PDF_GETTER(MRST01LO_Getter);

PDF_Base *MRST01LO_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new PDF_MRST01LO(args.m_bunch);
}

void MRST01LO_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"MRST 2001 LO fit\n"
     <<std::string(width+4,' ')<<"see hep-ph/0201127";
}

MRST01LO_Getter *p_get_mrst01lo;

extern "C" void InitPDFLib()
{
  p_get_mrst01lo = new MRST01LO_Getter("MRST01LO");
}

extern "C" void ExitPDFLib()
{
  delete p_get_mrst01lo;
}
