#include "QCD_Splitting_Functions.H"
#include "Sudakov_Tools.H"
#include "MathTools.H"
#include "Random.H"

#include "Message.H"

using namespace APACIC;
using namespace std;

// quark to quark + gluon  splitting function
q_qg::q_qg(ATOOLS::Flavour quarkflavour) : p_tools(0) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(ATOOLS::kf::gluon);
  m_alpha    = 1.;
}

q_qg::q_qg(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
  p_tools(_tools) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(ATOOLS::kf::gluon);
  m_alpha    = p_tools->GetASmax();
}

double q_qg::operator()(double z) {return s_CF*(1.+z*z)/(1.-z);};             

double q_qg::GetZ()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());   
}

double q_qg::GetCoupling()         { return m_alpha;}
double q_qg::GetCoupling(double t) { return p_tools->AlphaS(t);}

double q_qg::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+z*z);
  return ( (1.+z*z)/2. - 
	   z*ATOOLS::sqr((1.-z)*m_flavs[0].PSMass())/
	   (pt2+ATOOLS::sqr((1.-z)*m_flavs[0].PSMass())) );
}
                   
double q_qg::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*s_CF*m_alpha*log((1.-m_zmin)/(1.-m_zmax));                              
}

double q_qg::Integral(double zmin, double zmax) 
{
  return s_CF*(-0.5*((1.+zmax)*(1.+zmax)-(1.+zmin)*(1.+zmin))-2.*log((1.-zmax)/(1.-zmin)));
}

// gluon to gluon + gluon splitting function (needed twice for initial state shower)
g_gg::g_gg() : p_tools(0) 
{
  m_flavs[0] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[1] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[2] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_alpha    = 1.;
}

g_gg::g_gg(Sudakov_Tools * _tools) : p_tools(_tools) 
{ 
  m_flavs[0] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[1] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[2] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_alpha    = p_tools->GetASmax();
}

double g_gg::operator()(double z) 
{
  return s_CA*ATOOLS::sqr(1.-z*(1.-z))/(z*(1.-z));
}

double g_gg::GetZ() 
{
  return 1./(1. + ((1.-m_zmin)/m_zmin) *
	     pow( m_zmin*(1.-m_zmax)/((1.-m_zmin)*m_zmax), ATOOLS::ran.Get()));
}

double g_gg::GetCoupling()         { return m_alpha;}
double g_gg::GetCoupling(double t) { return p_tools->AlphaS(t);}
double g_gg::GetWeight(double z,double pt2,bool masses)   
{ 
  return ATOOLS::sqr(1.-z*(1.-z));
}
    
double g_gg::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return s_CA*m_alpha*log((1.-m_zmin)*m_zmax/(m_zmin*(1.-m_zmax)));                    
} 

double g_gg::Integral(double zmin, double zmax) 
{
  return s_CA*(log(zmax*(1.-zmin)/(zmin*(1.-zmax)))-
	       zmax*(zmax*(zmax/3.-0.5)+2.)+zmin*(zmin*(zmin/3.-0.5)+2.));
}

//! gluon to quark + anti-quark splitting function
g_qq::g_qq(ATOOLS::Flavour quarkflavour): p_tools(0) 
{
  m_flavs[0] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[1] = quarkflavour; 
  m_flavs[2] = quarkflavour.Bar(); 
  m_alpha    = 1.;
}

g_qq::g_qq(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
  p_tools (_tools) 
{
  m_flavs[0] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
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

double g_qq::GetCoupling()         { return m_alpha; }
double g_qq::GetCoupling(double t) { return p_tools->AlphaS(t); }
double g_qq::GetWeight(double z,double pt2,bool masses) 
{ 
  if (!masses) return (*this)(z)/s_TR;
  return (1. - 2.*z*(1.-z)*(1.- pt2/(pt2+ATOOLS::sqr(m_flavs[1].PSMass())))); // /s_TR;
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
q_gq::q_gq(ATOOLS::Flavour quarkflavour) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[2] = quarkflavour; 
  m_alpha    = 1.;
}

q_gq::q_gq(ATOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) : 
  p_tools(_tools) 
{
  m_flavs[0] = quarkflavour; 
  m_flavs[1] = ATOOLS::Flavour(ATOOLS::kf::gluon); 
  m_flavs[2] = quarkflavour; 
  m_alpha    = p_tools->GetASmax();
}

double q_gq::operator()(double z) {return s_CF*(1.+(1.-z)*(1.-z))/z;} 
double q_gq::GetZ()      
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double q_gq::GetCoupling()         { return m_alpha;}
double q_gq::GetCoupling(double t) { return p_tools->AlphaS(t);}
double q_gq::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+ATOOLS::sqr(1.-z));
  return ( (1.+ATOOLS::sqr(1.-z))/2. - 
	   (1.-z)*ATOOLS::sqr(z*m_flavs[0].PSMass())/
	   (pt2+ATOOLS::sqr(z*m_flavs[0].PSMass())) );
}
                   
double q_gq::CrudeInt(double _zmin, double _zmax) {
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*s_CF*m_alpha*log(m_zmax/m_zmin);
}

double q_gq::Integral(double zmin, double zmax) 
{
  return s_CF*(0.5*(zmax*zmax-zmin*zmin)-2.*(zmax-zmin)+2.0*log(zmax/zmin));
}
