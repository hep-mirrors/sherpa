#include "PDF/Electron/PDF_Electron.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/MathTools.H"
//#include <iostream>


using namespace MODEL;
using namespace ATOOLS;
using namespace PDF;

PDF_Electron::PDF_Electron(const Flavour _bunch,const int _izetta,const int _order) : 
  m_izetta(_izetta), m_order(_order)
{
  m_xmin=1.e-6;
  m_xmax=.999999;
  m_q2min=0.25;
  m_q2max=1.e14;

  m_bunch  = _bunch;
  m_partons.insert(m_bunch);
  m_type   = std::string("PDF_")+std::string(m_bunch.IDName());
  
  m_mass   = m_bunch.Mass(true);

  m_init=false;
}

double PDF_Electron::GetXPDF(const ATOOLS::Flavour  fl) {
  if (fl==m_bunch) return m_xpdf;
  return 0.;
}

PDF_Base * PDF_Electron::GetCopy() { return new PDF_Electron(m_bunch,m_order,m_izetta); }

void PDF_Electron::CalculateSpec(double x,double Q2)
{
  if (!m_init) {
    m_alpha  = (*aqed)(sqr(rpa.gen.Ecms()));
    double L = log(sqr(rpa.gen.Ecms()/m_bunch.Mass(true)));
    m_exponent = m_beta = (*aqed)(sqr(m_bunch.Mass(true)))/M_PI*(L-1.);
    m_init = true;
  }

  m_xpdf  = 0.;
  m_alpha = (*aqed)(Q2);
  if (x>=0.999999) return;

  double L       = 2.*log(sqrt(Q2)/m_mass);
  double beta_e  = 2.*m_alpha/M_PI*(L-1.);
  double eta     = 2.*m_alpha/M_PI*L;

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
  default:
    m_beta = betaS = betaH = beta_e; 
  }
  double gamma   = ::exp(Gammln(1.+m_beta/2.)); 

  // Produces collinear bremsstrahlung in exponentiated LLA,
  double S=0.,h0=0.,h1=0.,h2=0.;
  S  = exp(-.5*GAMMA_E*m_beta+.375*betaS)/gamma*m_beta/2.;
  h0 = -.25*(1.+x)*betaH;

  if (m_order>=2) {  
    h1   = -1./32.*betaH*betaH*
      ((1.+3.*x*x)/(1.-x)*log(x)+4.*(1.+x)*log(1.-x)+5.+x); 
  }
  if (m_order==3) {
    int i = 1;
    double Li = 0.;
    double del= 1.;
    double del2=1.;
    for (i=1;del2>rpa.gen.Accu();++i) {  
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
  if (x>0.9999) m_xpdf *= pow(100.,m_beta/2)/(pow(100.,m_beta/2)-1.);
}

DECLARE_PDF_GETTER(PDFE_Getter);

PDF_Base *PDFE_Getter::operator()
  (const Parameter_Type &args) const
{
  if (!args.m_bunch.Kfcode()==kf_e) return NULL;
  int izetta=args.p_read->GetValue<int>("ISR_E_ORDER",1);
  int order=args.p_read->GetValue<int>("ISR_E_SCHEME",2);
  return new PDF_Electron(args.m_bunch,izetta,order);
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
