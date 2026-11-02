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

std::ostream & PDF::operator<<(std::ostream & s,
                               const PDF::alpha::scheme & as)
{
  if (as==alpha::thomson)            s<<"Thomson: alpha(0)";
  else if (as==alpha::shortdistance) s<<"Short-Distance: alpha_sd";
  else if (as==alpha::running)       s<<"Running: alpha(Q2)";
  else                               s<<"Undefined";
  return s;
}

std::istream &PDF::operator>>(std::istream &s,
                              PDF::alpha::scheme &as)
{
  std::string tag;
  getline(s,tag);
  as=alpha::running;
  if (tag.find("Thomson")!=std::string::npos)
    as=alpha::thomson;
  if (tag.find("alpha(0)")!=std::string::npos)
    as=alpha::thomson;
  if (tag.find("Short-Distance")!=std::string::npos)
    as=alpha::shortdistance;
  if (tag.find("alpha_sd")!=std::string::npos)
    as=alpha::shortdistance;
  if (tag.find("Running")!=std::string::npos)
    as=alpha::running;
  if (tag.find("alpha(Q2)")!=std::string::npos)
    as=alpha::running;
  return s;
}

std::ostream & PDF::operator<<(std::ostream & s,
                               const PDF::overweight::scheme & os)
{
  if (os==overweight::none)               s<<"None";
  else if (os==overweight::constant)      s<<"Constant";
  else if (os==overweight::linear)        s<<"Linear";
  else if (os==overweight::quadratic)     s<<"Quadratic";
  else if (os==overweight::sigmoid)       s<<"Sigmoid";
  return s;
}

std::istream &PDF::operator>>(std::istream &s,
                              PDF::overweight::scheme &os)
{
  std::string tag;
  getline(s,tag);
  os=overweight::none;
  if (tag.find("Constant")!=std::string::npos)
    os=overweight::constant;
  if (tag.find("Linear")!=std::string::npos)
    os=overweight::linear;
  if (tag.find("Quadratic")!=std::string::npos)
    os=overweight::quadratic;
  if (tag.find("Sigmoid")!=std::string::npos)
    os=overweight::sigmoid;
  return s;
}

PDF_Electron::PDF_Electron(const Flavour _bunch,const int _izetta,const int _order) : 
  m_mass(0.), m_xpdf(0.),
  m_izetta(_izetta), m_order(_order), m_init(0),
  m_ascheme(alpha::running), m_owscheme(overweight::none)
{
  Settings& s = Settings::GetMainSettings();
  m_eps = s["PDFE_EPS"].Get<double>();
  m_delta = s["PDFE_DELTA"].Get<double>();
  m_xmax = 1.-m_eps;
  m_xmin=1.e-6;
  m_q2min=1e-12;
  m_q2max=1.e14;

  m_owscheme=s["PDFE_OVERWEIGHT_SCHEME"].SetDefault(overweight::linear)
	                                .Get<overweight::scheme>();
  msg_Info()<<METHOD<<"(): overweight scheme: "<<m_owscheme
	                <<", delta = "<<m_delta<<" , eps = "<<m_eps<<std::endl;

  m_ascheme = s["PDFE_ALPHA_SCHEME"].SetDefault(alpha::thomson)
	                            .Get<alpha::scheme>();
  msg_Info()<<METHOD<<"(): alpha scheme: "<<m_ascheme<<std::endl;

  m_set    = "PDFE";
  m_bunch  = _bunch;
  m_partons.insert(m_bunch);
  m_type   = std::string("PDF_")+std::string(m_bunch.IDName());
  
  m_mass   = m_bunch.Mass(true);

  m_init=false;
}

double PDF_Electron::Alpha(const double &Q2) const
{
  switch (m_ascheme) {
  case alpha::thomson:       return aqed->AqedThomson();
  case alpha::shortdistance: return aqed->Default();
  case alpha::running:       return (*aqed)(Q2);
  default:
    THROW(not_implemented,"Specified alpha-scheme for ePDF not implemented.");
  }
  return 0.;
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

double PDF_Electron::GetIntPDF(const ATOOLS::Flavour& fl,
                               const double& eps, const double& Q2)
{
  return GetIntPDF(fl.Kfcode(),fl.IsAnti(),eps,Q2);
}

double PDF_Electron::GetIntPDF(const kf_code& kf, bool anti,
	                       const double& eps, const double& Q2)
{
  // with the above structure function in the limit z->1, 
  // f(z) = beta/2 (1-z)^{beta/2-1} ,
  // the general form \int_{eta}^{eta/(1-eps)} dx f_a(eta/x)
  // returns for a the electron 
  // \int_{eta}^{eta/(1-eps)} dx/eta beta/2 (1-eta/x)^{beta/2-1}
  // -> eps^{beta/2}  for eps->0
  // take values of beta = beta(Q2) from already executed CalculateSpec function
  if (kf==m_bunch.Kfcode() && anti==m_bunch.IsAnti()) {
    double L(log(Q2/sqr(m_mass)));
    double beta(2.*Alpha(Q2)/M_PI*(L-1.));
    return pow(eps,beta/2.);
  }
  else return 0.;
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
    // this calculates the largest possible beta: alpha(s)/pi * [log(S/me2)-1]
    double L = log(sqr(rpa->gen.Ecms()/m_mass));
    m_exponent = (*aqed)(sqr(rpa->gen.Ecms()))/M_PI*(L-1.);
    m_init = true;
  }

  double alpha = Alpha(Q2);
  m_xpdf  = 0.;
  if (x>m_xmax) return;

  double L       = log(Q2/sqr(m_mass));
  double beta    = 2.*alpha/M_PI*(L-1.);
  double eta     = 2.*alpha/M_PI*L;

  if (beta <= 0 || eta <= 0) {
    msg_Error()<<"Warning: freezing out electron structure function at m_e\n";
    // set Q^2 to approximately e m^2
    double eps = 1.e-6;
    L = 1. + eps*M_PI/(2.*alpha);
    beta = eps;
    eta = 2.*alpha/M_PI + eps;
  }

  double betaS,betaH;
  switch(m_izetta) {
  case 0:
    betaS = betaH = eta;
    break;
  case 1:
    betaS = beta;
    betaH = eta;
    break;
  case 2:
    // default
    betaS = betaH = beta; 
    break;
  default:
    THROW(fatal_error, "Undefined scheme for electron structure function!");
  }
  // Gamma[1+beta/2]
  double gamma   = ::exp(Gammln(1.+beta/2.)); 

  // Produces collinear bremsstrahlung in exponentiated LLA,
  double S=0.,h0=0.,h1=0.,h2=0.;
  S  = beta/2.*exp(-.5*GAMMA_E*beta+.375*betaS)/gamma;

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

  m_xpdf = x * (S*pow(1.-x,beta/2.-1.)+(h0+h1+h2));  
  
  // Rescale the region 1-delta<x<1-epsilon to contain the missing XS contributions
  // from 1-epsilon<x<1
  if (x>1.-m_delta) {
    if (m_owscheme==overweight::constant) {
      // rescaling by a constant
      m_xpdf *= pow(m_delta/m_eps,beta/2.)/(pow(m_delta/m_eps,beta/2.)-1.);
    }
    else if (m_owscheme==overweight::linear) {
      // rescaling with a linear function which is continuous
      double a = 1./( beta/(beta+2.)*m_eps
		     -m_delta
		     +m_delta*(1.-beta/(beta+2.))*pow(m_delta/m_eps,beta/2.));
      double b = 1.-a*(1.-m_delta);
      m_xpdf *= a*x+b;
    }
    else if (m_owscheme==overweight::quadratic) {
      // rescaling with quadratic function which is differentiable
      double a = -(2.+beta)*(4.+beta)
                  /( (8.+6.*beta+sqr(beta))*sqr(m_delta)
                     -(8.*beta+2.*sqr(beta))*m_delta*m_eps
                     +(2.*beta+sqr(beta))*sqr(m_eps)
                     -8.*sqr(m_delta)*pow(m_delta/m_eps,beta/2.) );
      double b = -2.*a*(1.-m_delta);
      double c = 1.+a*sqr(1.-m_delta);
      m_xpdf *= a*x*x+b*x+c;
    }
    else if (m_owscheme==overweight::sigmoid) {
      // rescaleing with sigmoid function which is infinitely differentiable
      THROW(not_implemented,"Sigmoid rescaling not implemented yet.");
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
