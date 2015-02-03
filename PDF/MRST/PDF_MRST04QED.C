#include "PDF/MRST/PDF_MRST04QED.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PDF;
using namespace ATOOLS;

extern "C" {
  void mrstqed_(double *x,double *q2,int *mode,
		double *upv,double *dnv,double *usea,
		double *dsea,double *str,double *chm,
		double *bot,double *glu,double *phot);
}

void mrstqed(double x,double q2,int mode,
	     double &upv,double &dnv,double &usea,
	     double &dsea,double &str,double &chm,
	     double &bot,double &glu,double &phot)
{
  mrstqed_(&x,&q2,&mode,&upv,&dnv,&usea,&dsea,&str,&chm,&bot,&glu,&phot);
}

extern "C" {
extern struct {
  char nfile[128], pfile[128];
} mrinput_;
#define input mrinput_
}

inline void MakeFortranString(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) output[j]=(char)input[j];
}

PDF_MRST04QED::PDF_MRST04QED(const ATOOLS::Flavour bunch):
  m_path(rpa->gen.Variable("SHERPA_SHARE_PATH")+"/MRST04Grid"),
  m_anti(1),
  m_mode(1)
{
  m_set="MRST04QED";
  m_type=m_set;
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
  m_partons.insert(Flavour(kf_photon));
  m_partons.insert(Flavour(kf_resummed));
  m_xmin=1.e-5;
  m_xmax=1.;
  m_q2min=1.25;
  m_q2max=1.e7;
  MakeFortranString(input.nfile,m_path+std::string("/qed6-10gridn.dat"),128);
  MakeFortranString(input.pfile,m_path+std::string("/qed6-10gridp.dat"),128);
}


PDF_Base *PDF_MRST04QED::GetCopy() 
{
  PDF_Base *copy = new PDF_MRST04QED(m_bunch);
  m_copies.push_back(copy);
  return copy;
}

void PDF_MRST04QED::CalculateSpec(const double& x,const double& Q2)
{
  m_overscaled=false;
  double xx(x);
  if(xx<m_xmin) xx=m_xmin;
  if (xx>m_xmax) {
    m_overscaled=true;
    return;
  }
  mrstqed(xx,Q2,m_mode,p_xpdfv[1],p_xpdfv[0],p_xpdf[1],
	  p_xpdf[0],p_xpdf[2],p_xpdf[3],p_xpdf[4],p_xpdf[5],p_xpdf[6]);
}


double PDF_MRST04QED::GetXPDF(const ATOOLS::Flavour& infl)
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
  case -kf_gluon : return m_rescale*p_xpdf[5]; 
  case  kf_photon :
  case -kf_photon : return m_rescale*p_xpdf[6]; 
  default: return 0.;
  }
}

double PDF_MRST04QED::GetXPDF(const kf_code& kf, bool anti)
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
  case -kf_gluon : return m_rescale*p_xpdf[5];
  case  kf_photon :
  case -kf_photon : return m_rescale*p_xpdf[6];
  default: return 0.;
  }
}

DECLARE_PDF_GETTER(MRST04QED_Getter);

PDF_Base *MRST04QED_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.IsHadron()) return NULL;
  return new PDF_MRST04QED(args.m_bunch);
}

void MRST04QED_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"MRST 2004 fit including O(alpha) contributions\n"
     <<std::string(width+4,' ')<<"see hep-ph/0411040";
}

MRST04QED_Getter *p_get_mrst04qed;

extern "C" void InitPDFLib()
{
  p_get_mrst04qed = new MRST04QED_Getter("MRST04QED");
}

extern "C" void ExitPDFLib()
{
  delete p_get_mrst04qed;
}
