#include "APACIC++/Showers/QCD_Splitting_Functions.H"
#include "APACIC++/Showers/Sudakov_Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"

#include "ATOOLS/Org/Message.H"

using namespace APACIC;
using namespace std;

// quark to quark + gluon  splitting function
q_qg::q_qg(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour quarkflavour,double fmed) : 
  Splitting_Function(ms), p_tools(0), m_kfactor(1.0), m_fmed(fmed) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon);
  m_alpha    = 1.;
}

q_qg::q_qg(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed) :
  Splitting_Function(ms), p_tools(_tools), m_kfactor(1.0), m_fmed(fmed)  
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon);
  m_alpha    = p_tools->GetASmax();
  if (s_kfactorscheme==1) m_kfactor  = 1.0+m_alpha/(2.0*M_PI)*s_kappa;
}

double q_qg::operator()(double z) 
{ 
  return m_kfactor*s_CF*(2.*(1.+m_fmed)/(1.-z)-(1.+z));
}    

double q_qg::GetZ() 
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());   
}

double q_qg::GetPhi(double z)
{
  return ATOOLS::ran.Get()*2.*M_PI;
}

const ATOOLS::Simple_Polarisation_Info q_qg::GetPolB(double z, double phi) 
{
  return ATOOLS::Simple_Polarisation_Info();
}

const ATOOLS::Simple_Polarisation_Info q_qg::GetPolC(double z, double phi, double phi_b) 
{
  double phi_c,f;
  double deno = 1./ATOOLS::sqr(1+z);
  double nume = 1 + z*z;
  do {
    phi_c=ATOOLS::ran.Get()*2.*M_PI;
    f=(nume + 2.*z*cos(2.*phi_c))*deno;
  }
  while (f<ATOOLS::ran.Get());
  return ATOOLS::Simple_Polarisation_Info(1,phi_c);
}


double q_qg::GetCoupling()         { return m_alpha;}
double q_qg::GetCoupling(double t) { return p_tools->AlphaS(t);}

double q_qg::GetWeight(double z,double pt2,bool massterm) 
{ 
  double weight(1.-(1.-z*z)/(2.*(1.+m_fmed)));
  if (massterm) {
    double mu(z*(1.-z)*ATOOLS::sqr(p_ms->Mass(m_flavs[0]))/
	      (pt2+ATOOLS::sqr((1.-z)*p_ms->Mass(m_flavs[0]))));
    weight *= 1.-(1.-z)*2.*mu/(2.*(1.+m_fmed)-(1.-z*z));
  }
  if (s_kfactorscheme==1) {
    double kappa(s_kappaG-p_tools->Nf(pt2)*s_kappaQ);
    weight *= (1.+kappa*p_tools->AlphaS(pt2)/(2.*M_PI))/m_kfactor;
  } 
  return weight;
}
                   
double q_qg::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*s_CF*m_alpha*m_kfactor*(1.+m_fmed)*log((1.-m_zmin)/(1.-m_zmax));         
}

double q_qg::Integral(double zmin, double zmax) 
{
  return s_CF*(-0.5*((1.+zmax)*(1.+zmax)-(1.+zmin)*(1.+zmin))-2.*log((1.-zmax)/(1.-zmin)));
}

// gluon to gluon + gluon splitting function (needed twice for initial state shower)
g_gg::g_gg(ATOOLS::Mass_Selector *&ms,double fmed) : 
  Splitting_Function(ms), p_tools(0), m_kfactor(1.0), m_fmed(fmed)  
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon); 
  m_alpha    = 1.;
}

g_gg::g_gg(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * _tools,double fmed) : 
  Splitting_Function(ms), p_tools(_tools), m_kfactor(1.0), m_fmed(fmed)  
{ 
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon); 
  m_alpha    = p_tools->GetASmax();
  if (s_kfactorscheme==1)
    m_kfactor  = 1.0+m_alpha/(2.0*M_PI)*s_kappa;
}

double g_gg::operator()(double z) 
{
  return m_kfactor*s_CA*( (1.+m_fmed)/(z*(1.-z)) + z*(1.-z)-2 );
}

double g_gg::GetZ() 
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg::GetPhi(double z)
{
  //  std::cout<<"g_gg::GetPhi("<<z<<") --- "<<std::endl;

  double phi,f;
  double nume1 = ATOOLS::sqr(1-z) + ATOOLS::sqr(z);
  double nume2 = ATOOLS::sqr(z*(1-z));
  double deno  = 1./(nume1 + 2.*nume2);
  do {
    phi=ATOOLS::ran.Get()*2.*M_PI;
    f=(nume1 + nume2*cos(2.*phi))*deno;
  } 
  while (f<ATOOLS::ran.Get());
  return phi;
}

const ATOOLS::Simple_Polarisation_Info g_gg::GetPolB(double z, double phi) 
{
  double phi_b,f;
  double c2phi = ATOOLS::sqr(cos(phi));
  double s2phi = 1. - c2phi;

  double nume1 = (1.-z)/z;
  double nume2 = z/(1.-z);
  double nume3 = z*(1.-z);
  double deno = 1./(nume1 + nume2 + nume3*c2phi);
  do {
    phi_b=ATOOLS::ran.Get()*2.*M_PI;
    double c2phi_b = ATOOLS::sqr(cos(phi_b));
    double s2phi_b = 1. - c2phi_b;

    f=(nume1*c2phi_b + nume2*(s2phi*s2phi_b + c2phi*c2phi_b) + nume3*c2phi)*deno;
  }
  while (f<ATOOLS::ran.Get());
  return ATOOLS::Simple_Polarisation_Info(1,phi_b);
}

const ATOOLS::Simple_Polarisation_Info g_gg::GetPolC(double z, double phi, double phi_b) 
{
  double phi_c,f;
  double c2phi = ATOOLS::sqr(cos(phi));
  double s2phi = 1. - c2phi;
  double c2phi_b = ATOOLS::sqr(cos(phi_b));
  double s2phi_b = 1. - c2phi_b;

  double nume1 = (1.-z)/z;
  double nume2 = z/(1.-z);
  double nume3 = z*(1.-z);
  double deno; // = 1./(nume1 + nume2 + nume3);
  deno = 1./(nume1*c2phi_b + nume2*(s2phi*s2phi_b + c2phi*c2phi_b) + nume3*c2phi);
  do {
    phi_c=ATOOLS::ran.Get()*2.*M_PI;
    double c2phi_c = ATOOLS::sqr(cos(phi_c));
    double s2phi_c = 1. - c2phi_c;

    f=(nume1*c2phi_b*(s2phi*s2phi_c + c2phi*c2phi_c) +
       nume2*c2phi_c*(s2phi*s2phi_b + c2phi*c2phi_b) + 
       nume3*c2phi*(s2phi_b*s2phi_c + c2phi_b*c2phi_c))*deno;
  }
  while (f<ATOOLS::ran.Get());
  return ATOOLS::Simple_Polarisation_Info(1,phi_c);
}


double g_gg::GetCoupling()         { return m_alpha;}
double g_gg::GetCoupling(double t) { return p_tools->AlphaS(t);}
double g_gg::GetWeight(double z,double pt2,bool masses)   
{ 
  double weight((ATOOLS::sqr(1.-z*(1.-z))+m_fmed)/(1.+m_fmed));
  if (s_kfactorscheme==1) {
    double kappa(s_kappaG-p_tools->Nf(pt2)*s_kappaQ);
    weight *= (1.+kappa*p_tools->AlphaS(pt2)/(2.*M_PI))/m_kfactor;
  } 
  return weight;
}
    
double g_gg::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return m_kfactor*s_CA*m_alpha*log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax)))*(1.+m_fmed);                    
} 

double g_gg::Integral(double zmin, double zmax) 
{
  return s_CA*(log(zmax*(1.-zmin)/(zmin*(1.-zmax)))-
	       zmax*(zmax*(zmax/3.-0.5)+2.)+zmin*(zmin*(zmin/3.-0.5)+2.));
}

//! gluon to quark + anti-quark splitting function
g_qq::g_qq(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour quarkflavour,double fmed) : 
  Splitting_Function(ms), p_tools(0), m_fmed(fmed)  
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = quarkflavour.Bar(); 
  m_alpha    = 1.;
}

g_qq::g_qq(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed) :
  Splitting_Function(ms), p_tools (_tools), m_fmed(fmed)  
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = quarkflavour.Bar(); 
  m_alpha    = p_tools->GetASmax();
}

double g_qq::operator()(double z) 
{
  return s_TR*(z*z + (1.-z)*(1.-z));
}

double g_qq::GetZ()      
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran.Get();                       
}

double g_qq::GetPhi(double z)
{
  double phi,f;
  double nume1 = z*z + ATOOLS::sqr(1-z);
  double nume2 = -2.*z*(1-z);
  do {
    phi = ATOOLS::ran.Get()*2.*M_PI;
    f   = nume1 + nume2*cos(2.*phi);
  }
  while (f<ATOOLS::ran.Get());
  return phi;
}

const ATOOLS::Simple_Polarisation_Info g_qq::GetPolB(double z, double phi) 
{
  return ATOOLS::Simple_Polarisation_Info();
}

const ATOOLS::Simple_Polarisation_Info g_qq::GetPolC(double z, double phi, double phi_b) 
{
  return ATOOLS::Simple_Polarisation_Info();
}

double g_qq::GetCoupling()         { return m_alpha; }
double g_qq::GetCoupling(double t) { return p_tools->AlphaS(t); }
double g_qq::GetWeight(double z,double pt2,bool massterm) 
{ 
  double weight((*this)(z)/s_TR);
  if (massterm) {
    double m2(ATOOLS::sqr(p_ms->Mass(m_flavs[1])));
    double mu(2.*m2*z*(1.-z)/(pt2+m2));
    weight *= 1.+mu/(1.-2.*z*(1.-z));
  }
  return weight;
}
                 
double g_qq::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return s_TR*m_alpha*(m_zmax-m_zmin);                                             
}

double g_qq::Integral(double zmin, double zmax) 
{
  return s_TR*(zmax*(zmax*(2.0*zmax/3.0-1.0)+1.0)+
 	       -zmin*(zmin*(2.0*zmin/3.0-1.0)+1.0));
}

// quark to qluon + quark splitting function (only used in Initial State Shower)
q_gq::q_gq(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour quarkflavour,double fmed) :
  Splitting_Function(ms), m_kfactor(1.0), m_fmed(fmed) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = quarkflavour; 
  m_alpha    = 1.;
}

q_gq::q_gq(ATOOLS::Mass_Selector *&ms,
	   ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools,double fmed) : 
  Splitting_Function(ms), p_tools(_tools), m_kfactor(1.0), m_fmed(fmed)  
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = quarkflavour; 
  m_alpha    = p_tools->GetASmax();
  if (s_kfactorscheme==1) m_kfactor  = 1.0+m_alpha/(2.0*M_PI)*s_kappa;
}

double q_gq::operator()(double z) {return s_CF*(2.*(1.+m_fmed))/z-(2.-z);} 
double q_gq::GetZ()      
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}
double q_gq::GetPhi(double z)
{
  //  std::cout<<"q_gq::GetPhi --- uniform"<<std::endl;

  return ATOOLS::ran.Get()*2.*M_PI;
}

const ATOOLS::Simple_Polarisation_Info q_gq::GetPolB(double z, double phi) 
{
  double phi_b,f;
  double deno = 1./ATOOLS::sqr(2-z);
  double nume = 1. + ATOOLS::sqr(1-z);
  do {
    phi_b=ATOOLS::ran.Get()*2.*M_PI;
    f=(nume + 2.*(1.-z)*cos(2.*phi_b))*deno;
  }
  while (f<ATOOLS::ran.Get());
  return ATOOLS::Simple_Polarisation_Info(1,phi_b);
}

const ATOOLS::Simple_Polarisation_Info q_gq::GetPolC(double z, double phi, double phi_b) 
{
  return ATOOLS::Simple_Polarisation_Info();
}


double q_gq::GetCoupling()         { return m_alpha;}
double q_gq::GetCoupling(double t) { return p_tools->AlphaS(t);}
double q_gq::GetWeight(double z,double pt2,bool massterm) 
{ 
  double weight(1.-(2.-z)*z/(2.*(1.+m_fmed)));
  if (massterm) {
    double mu(z*(1.-z)*ATOOLS::sqr(p_ms->Mass(m_flavs[0]))/
	      (pt2+ATOOLS::sqr(z*p_ms->Mass(m_flavs[0]))));
    weight *= 1.-z*mu/(2.*(1.+m_fmed)-z*(2.-z));
  }
  if (s_kfactorscheme==1) {
    double kappa(s_kappaG-p_tools->Nf(pt2)*s_kappaQ);
    weight *= (1.+kappa*p_tools->AlphaS(pt2)/(2.*M_PI))/m_kfactor;
  } 
  return weight;
}
                   
double q_gq::CrudeInt(double _zmin, double _zmax) {
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*(1.+m_fmed)*m_kfactor*s_CF*m_alpha*log(m_zmax/m_zmin);
}

double q_gq::Integral(double zmin, double zmax) 
{
  return s_CF*(0.5*(zmax*zmax-zmin*zmin)-2.*(zmax-zmin)+2.0*log(zmax/zmin));
}
