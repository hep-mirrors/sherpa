#include "APACIC++/Showers/SUSY_QCD_Splitting_Functions.H"
#include "APACIC++/Showers/Sudakov_Tools.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"

#include "ATOOLS/Org/Message.H"

using namespace APACIC;
using namespace std;

// squark to squark + gluon  splitting function
SQuark__SQuark_Gluon::SQuark__SQuark_Gluon
(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour) : 
  Splitting_Function(ms), p_tools(0) 
{
  m_flavs[0] = squarkflavour; 
  m_flavs[1] = squarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon);
  m_alpha    = 1.;
}

SQuark__SQuark_Gluon::SQuark__SQuark_Gluon
(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour,Sudakov_Tools *tools) :
  Splitting_Function(ms), p_tools(tools) 
{
  m_flavs[0] = squarkflavour; 
  m_flavs[1] = squarkflavour; 
  m_flavs[2] = ATOOLS::Flavour(kf_gluon);
  m_alpha    = p_tools->GetASmax();
}

double SQuark__SQuark_Gluon::operator()(double z) {return s_CF*(2.*z)/(1.-z);}

double SQuark__SQuark_Gluon::GetZ()
{
  return 1.-(1.-m_zmin)*pow((1.-m_zmax)/(1.-m_zmin),ATOOLS::ran.Get());   
}

double SQuark__SQuark_Gluon::GetCoupling()         { return m_alpha;}
double SQuark__SQuark_Gluon::GetCoupling(double t) { return p_tools->AlphaS(t);}

double SQuark__SQuark_Gluon::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return z;
  return ( z - 
	   z*ATOOLS::sqr((1.-z)*p_ms->Mass(m_flavs[0]))/
	   (pt2+ATOOLS::sqr((1.-z)*p_ms->Mass(m_flavs[0]))));
}
                   
double SQuark__SQuark_Gluon::CrudeInt(double zmin, double zmax) 
{
  m_zmin = zmin;
  m_zmax = zmax;
  return 2.*s_CF*m_alpha*log((1.-m_zmin)/(1.-m_zmax));                              
}

double SQuark__SQuark_Gluon::Integral(double zmin, double zmax) 
{
  return 0.;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gluon__Gluino_Gluino::Gluon__Gluino_Gluino(ATOOLS::Mass_Selector *&ms):
  Splitting_Function(ms), p_tools(0) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = ATOOLS::Flavour(kf_Gluino); 
  m_flavs[2] = ATOOLS::Flavour(kf_Gluino); 
  m_alpha    = 1.;
}

Gluon__Gluino_Gluino::Gluon__Gluino_Gluino
(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * _tools) :
  Splitting_Function(ms), p_tools (_tools) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = ATOOLS::Flavour(kf_Gluino); 
  m_flavs[2] = ATOOLS::Flavour(kf_Gluino); 
  m_alpha    = p_tools->GetASmax();
}

double Gluon__Gluino_Gluino::operator()(double z) 
{
  return s_CA*(z*z + (1.-z)*(1.-z));
}

double Gluon__Gluino_Gluino::GetZ()      
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran.Get();                       
}

double Gluon__Gluino_Gluino::GetCoupling()         { return m_alpha; }
double Gluon__Gluino_Gluino::GetCoupling(double t) { return p_tools->AlphaS(t); }
double Gluon__Gluino_Gluino::GetWeight(double z,double pt2,bool masses) 
{ 
  if (!masses) return (*this)(z)/s_CA;
  return (1. - 2.*z*(1.-z)*(1.- pt2/(pt2+ATOOLS::sqr(p_ms->Mass(m_flavs[1]))))); // /s_CA;
}
                 
double Gluon__Gluino_Gluino::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return s_CA*m_alpha*(m_zmax-m_zmin);                                             
}

double Gluon__Gluino_Gluino::Integral(double zmin, double zmax) 
{
  return 0.;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gluino__Gluino_Gluon::Gluino__Gluino_Gluon(ATOOLS::Mass_Selector *&ms):
  Splitting_Function(ms)
{
  m_flavs[0] = ATOOLS::Flavour(kf_Gluino); 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = ATOOLS::Flavour(kf_Gluino); 
  m_alpha    = 1.;
}

Gluino__Gluino_Gluon::Gluino__Gluino_Gluon
(ATOOLS::Mass_Selector *&ms,Sudakov_Tools * _tools) : 
  Splitting_Function(ms), p_tools(_tools) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_Gluino); 
  m_flavs[1] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[2] = ATOOLS::Flavour(kf_Gluino); 
  m_alpha    = p_tools->GetASmax();
}

double Gluino__Gluino_Gluon::operator()(double z) {return s_CA*(1.+(1.-z)*(1.-z))/z;} 
double Gluino__Gluino_Gluon::GetZ()      
{
  return m_zmin*pow(m_zmax/m_zmin,ATOOLS::ran.Get());
}

double Gluino__Gluino_Gluon::GetCoupling()         { return m_alpha;}
double Gluino__Gluino_Gluon::GetCoupling(double t) { return p_tools->AlphaS(t);}
double Gluino__Gluino_Gluon::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+ATOOLS::sqr(1.-z));
  return ( (1.+ATOOLS::sqr(1.-z))/2. - 
	   (1.-z)*ATOOLS::sqr(z*p_ms->Mass(m_flavs[0]))/
	   (pt2+ATOOLS::sqr(z*p_ms->Mass(m_flavs[0]))) );
}
                   
double Gluino__Gluino_Gluon::CrudeInt(double _zmin, double _zmax) {
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 2.*s_CA*m_alpha*log(m_zmax/m_zmin);
}

double Gluino__Gluino_Gluon::Integral(double zmin, double zmax) 
{
  return 0.;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gluon__SQuark_SQuark::Gluon__SQuark_SQuark
(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour): 
  Splitting_Function(ms), p_tools(0) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = squarkflavour; 
  m_flavs[2] = squarkflavour.Bar(); 
  m_alpha    = 1.;
}

Gluon__SQuark_SQuark::Gluon__SQuark_SQuark
(ATOOLS::Mass_Selector *&ms,ATOOLS::Flavour squarkflavour,Sudakov_Tools * _tools) :
  Splitting_Function(ms), p_tools (_tools) 
{
  m_flavs[0] = ATOOLS::Flavour(kf_gluon); 
  m_flavs[1] = squarkflavour; 
  m_flavs[2] = squarkflavour.Bar(); 
  m_alpha    = p_tools->GetASmax();
}

double Gluon__SQuark_SQuark::operator()(double z) 
{
  return s_TR*(z*(1.-z));
}

double Gluon__SQuark_SQuark::GetZ()      
{
  return m_zmin+(m_zmax-m_zmin)*ATOOLS::ran.Get();                       
}

double Gluon__SQuark_SQuark::GetCoupling()         { return m_alpha; }
double Gluon__SQuark_SQuark::GetCoupling(double t) { return p_tools->AlphaS(t); }
double Gluon__SQuark_SQuark::GetWeight(double z,double pt2,bool masses) 
{ 
  if (!masses) return (*this)(z)/(s_TR*.25);
  return (4.*z*(1.-z)*(1.- pt2/(pt2+ATOOLS::sqr(p_ms->Mass(m_flavs[1]))))); // /(s_TR*.25);
}
                 
double Gluon__SQuark_SQuark::CrudeInt(double _zmin, double _zmax) 
{
  m_zmin = _zmin;
  m_zmax = _zmax;
  return 0.25*s_TR*m_alpha*(m_zmax-m_zmin);
}

double Gluon__SQuark_SQuark::Integral(double zmin, double zmax) 
{
  return 0.;	    
}
