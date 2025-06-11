#include "PDF/Electron/PDF_Electron.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"
//#include <iostream>


using namespace MODEL;
using namespace ATOOLS;
using namespace PDF;

PDF_Electron::PDF_Electron(const Flavour _bunch,const int _izetta,const int _order) : 
  m_izetta(_izetta), m_order(_order)
{
  Settings& s = Settings::GetMainSettings();

  m_eps = s["PDFE_EPS"].Get<double>();
  m_delta = s["PDFE_DELTA"].Get<double>();
  m_xmin=1.e-6;
  m_xmax = 1.-m_eps;
  m_q2min=1.e-12; // electron mass handled separately
  m_q2max=1.e14;

  s["PDFE_RESCALE_MODE"].SetDefault(1);
  m_rescale_mode = s["PDFE_RESCALE_MODE"].Get<int>();

  s["PDFE_RUNNING_ALPHA"].SetDefault(false);
  m_runningalpha = s["PDFE_RUNNING_ALPHA"].Get<bool>();

  m_set    = "PDFE";
  m_bunch  = _bunch;
  m_partons.insert(m_bunch);
  m_type   = std::string("PDF_")+std::string(m_bunch.IDName());
  
  m_mass   = m_bunch.Mass(true);

  m_init=false;
}

double PDF_Electron::GetXPDF(const ATOOLS::Flavour& fl)
{
  if (fl==m_bunch) return m_xpdf;
  return 0.;
}

double PDF_Electron::GetXPDF(const kf_code& kf, bool anti)
{
  if (kf==m_bunch.Kfcode() && anti==m_bunch.IsAnti()) return m_xpdf;
  return 0.;
}

PDF_Base * PDF_Electron::GetCopy() { return new PDF_Electron(m_bunch,m_izetta,m_order); }

void PDF_Electron::Calculate(double x,double Q2)
{
  if(Q2<m_q2min) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" < "<<sqrt(m_q2min)<<". Return.\n";
    lasterr=Q2;
    return;
  }
  if(Q2>m_q2max) {
    static double lasterr(-1.0);
    if (Q2!=lasterr)
    msg_Error()<<METHOD<<"(): Q-range violation Q = "<<sqrt(Q2)
	       <<" > "<<sqrt(m_q2max)<<". Return.\n";
    lasterr=Q2;
    return;
  }
  double xR = x*(m_rescX?m_rescale:1.);
  if(xR<m_xmin*m_rescale) {
    static double lasterr(-1.0);
    if (xR!=lasterr)
    msg_Error()<<METHOD<<"(): x = "<<xR<<" ("<<m_rescale
	       <<") < "<<m_xmin<<". Return.\n";
    lasterr=xR;
    return;
  }
  if(xR>m_xmax*m_rescale) {
    static double lasterr(-1.0);
    if (xR!=lasterr)
    msg_Tracking()<<METHOD<<"(): x = "<<x<<" ("<<m_rescale
	       <<") > "<<m_xmax<<". Return.\n";
    lasterr=xR;
    return;
  }
  return CalculateSpec(xR,Q2);
}

void PDF_Electron::CalculateSpec(const double& x, const double& Q2)
{
  if (!m_init) {
    m_alpha  = m_runningalpha?(*aqed)(sqr(rpa->gen.Ecms())):aqed->AqedThomson();
    double L = log(sqr(rpa->gen.Ecms()/m_bunch.Mass(true)));
    m_exponent = m_beta = (*aqed)(sqr(m_bunch.Mass(true)))/M_PI*(L-1.);
    m_init = true;
  }

  m_xpdf  = 0.;
  m_alpha = m_runningalpha?(*aqed)(Q2):aqed->AqedThomson();
  if (x>m_xmax) return;

  double L       = 2.*log(sqrt(Q2)/m_mass);
  double beta_e  = 2.*m_alpha/M_PI*(L-1.);
  double eta     = 2.*m_alpha/M_PI*L;

  if (beta_e <= 0 || eta <= 0) {
    msg_Error()<<"Warning: freezing out electron structure function at m_e\n";
    // set Q^2 to approximately e m^2
    double betamin = 1.e-6;
    L = 1 + betamin*M_PI/(2.*m_alpha);
    beta_e = betamin;
    eta = 2.*m_alpha/M_PI + betamin;
  }

  double betaS,betaH;
  switch(m_izetta) {
  case 0:
    m_beta = beta_e;
    betaS  = betaH = eta;
    break;
  case 1:
    m_beta = betaS = beta_e;
    betaH  = eta;
    break;
  case 2:
    // default
    m_beta = betaS = betaH = beta_e;
  default:
    THROW(fatal_error, "Undefined scheme for electron structure function!");
  }
  double gamma   = ::exp(Gammln(1.+m_beta/2.)); 

  // Produces collinear bremsstrahlung in exponentiated LLA,
  double S=0.,h0=0.,h1=0.,h2=0.;
  S  = exp(-.5*GAMMA_E*m_beta+.375*betaS)/gamma*m_beta/2.;
  
  // default order is 1
  if (m_order>=1) {
    h0 = -.25*(1.+x)*betaH;
  }
  if (m_order>=2) {  
    h1   = -1./32.*betaH*betaH*
      ((1.+3.*x*x)/(1.-x)*log(x)+4.*(1.+x)*log(1.-x)+5.+x); 
  }
  if (m_order>=3) {
    int i = 1;
    double Li = 0.;
    double del= 1.;
    double del2=1.;
    for (i=1;del2>rpa->gen.Accu();++i) {  
      del *=x;
      del2 = del/double(i*i);
      Li  += del2;
    }
    h2 = -1./384.*betaH*betaH*betaH*
      ((1.+x)*(6.*Li+12.*sqr(log(1.-x))-3.*M_PI*M_PI)+
       1./(1.-x)*(1.5*(1.+8.*x+3.*x*x)*log(x)+6.*(x+5.)*(1.-x)*log(1.-x)+
       12.*(1.+x*x)*log(x)*log(1.-x)-(.5+3.5*x*x)*sqr(log(x))+
      .25*(39.-24.*x-15.*x*x)));                
  } 

  m_xpdf = x * (S*pow(1.-x,m_beta/2.-1.)+(h0+h1+h2));

  // Rescale the region 1-delta<x<1-epsilon to contain the missing XS contributions
  // from 1-epsilon<x<1
  if (m_rescale_mode==1) {
    // lambda rescaling (default)
    if (x>1.-m_delta) {
      m_xpdf *= pow(m_delta/m_eps,m_beta/2)/(pow(m_delta/m_eps,m_beta/2)-1.);
    }
  }
  else if (m_rescale_mode==2) {
    // linear rescaling, used for QED shower as it is continuous at x=1-delta
    if (x>1.-m_delta) {
      double a = 1./(m_beta/(m_beta+2)*m_eps-m_delta+m_delta*(1.-m_beta/(m_beta+2))*pow(m_delta/m_eps,m_beta/2));
      double b = 1.-a*(1-m_delta);
      m_xpdf *= a*x+b;
    }
  }
}

DECLARE_PDF_GETTER(PDFE_Getter);

PDF_Base *PDFE_Getter::operator()
  (const Parameter_Type &args) const
{
  if (args.m_bunch.Kfcode()!=kf_e && args.m_bunch.Kfcode()!=kf_mu) return NULL;
  return new PDF_Electron(args.m_bunch,args.m_scheme,args.m_order);
}

void PDFE_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"electron PDF";
}

PDFE_Getter *p_get_pdfe;

extern "C" void InitPDFLib()
{
  p_get_pdfe = new PDFE_Getter("PDFe");
}

extern "C" void ExitPDFLib()
{
  delete p_get_pdfe;
}
