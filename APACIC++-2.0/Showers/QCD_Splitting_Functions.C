#include "QCD_Splitting_Functions.H"

using namespace APACIC;

// quark to quark + gluon  splitting function
q_qg::q_qg(APHYTOOLS::Flavour quarkflavour) : tools(0) 
{
  flavs[0] = quarkflavour; 
  flavs[1] = quarkflavour; 
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  alpha    = 1.;
}

q_qg::q_qg(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
  tools(_tools) 
{
  flavs[0] = quarkflavour; 
  flavs[1] = quarkflavour; 
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
  alpha    = tools->GetASmax();
  cout<<" q_qg: alpha ="<<alpha<<endl;
}

double q_qg::operator()(double z) {return CF*(1.+z*z)/(1.-z);};             

double q_qg::GetZ()
{
  return 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),AMATOOLS::ran.Get());   
}

double q_qg::GetCoupling()         { return alpha;}
double q_qg::GetCoupling(double t) { return tools->AlphaS(t);}

double q_qg::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+z*z);
  return ( (1.+z*z)/2. - 
	   z*AMATOOLS::sqr((1.-z)*flavs[0].PSMass())/
	   (pt2+AMATOOLS::sqr((1.-z)*flavs[0].PSMass())) );
}
                   
double q_qg::CrudeInt(double _zmin, double _zmax) 
{
  zmin = _zmin;
  zmax = _zmax;
  return 2.*CF*alpha*log((1.-zmin)/(1.-zmax));                              
}

// gluon to gluon + gluon splitting function (needed twice for initial state shower)
g_gg::g_gg() : tools(0) 
{
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  alpha    = 1.;
}

g_gg::g_gg(Sudakov_Tools * _tools) : tools(_tools) 
{ 
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  alpha    = tools->GetASmax();
  cout<<" q_gg: alpha ="<<alpha<<endl;

}

double g_gg::operator()(double z) 
{
  return CA*AMATOOLS::sqr(1.-z*(1.-z))/(z*(1.-z));
}

double g_gg::GetZ() 
{
  return 1./(1. + ((1.-zmin)/zmin) *
	     pow( zmin*(1.-zmax)/((1.-zmin)*zmax), AMATOOLS::ran.Get()));
}

double g_gg::GetCoupling()         { return alpha;}
double g_gg::GetCoupling(double t) { return tools->AlphaS(t);}
double g_gg::GetWeight(double z,double pt2,bool masses)   
{ 
  return AMATOOLS::sqr(1.-z*(1.-z));
}
    
double g_gg::CrudeInt(double _zmin, double _zmax) 
{
  zmin = _zmin;
  zmax = _zmax;
  return CA*alpha*log((1.-zmin)*zmax/(zmin*(1.-zmax)));                    
} 

//! gluon to quark + anti-quark splitting function
g_qq::g_qq(APHYTOOLS::Flavour quarkflavour): tools(0) 
{
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[1] = quarkflavour; 
  flavs[2] = quarkflavour.Bar(); 
  alpha    = 1.;
}

g_qq::g_qq(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
  tools (_tools) 
{
  flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[1] = quarkflavour; 
  flavs[2] = quarkflavour.Bar(); 
  alpha    = tools->GetASmax();
}

double g_qq::operator()(double z) 
{
  return TR*(z*z + (1.-z)*(1.-z));
}

double g_qq::GetZ()      
{
  return zmin+(zmax-zmin)*AMATOOLS::ran.Get();                       
}

double g_qq::GetCoupling()         { return alpha; }
double g_qq::GetCoupling(double t) { return tools->AlphaS(t); }
double g_qq::GetWeight(double z,double pt2,bool masses) 
{ 
  if (!masses) return (*this)(z)/TR;
  return (1. - 2.*z*(1.-z)*(1.- pt2/(pt2+AMATOOLS::sqr(flavs[1].PSMass())))); // /TR;
}
                 
double g_qq::CrudeInt(double _zmin, double _zmax) 
{
  zmin = _zmin;
  zmax = _zmax;
  return TR*alpha*(zmax-zmin);                                             
}

// quark to qluon + quark splitting function (only used in Initial State Shower)
q_gq::q_gq(APHYTOOLS::Flavour quarkflavour) 
{
  flavs[0] = quarkflavour; 
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[2] = quarkflavour; 
  alpha    = 1.;
}

q_gq::q_gq(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) : 
  tools(_tools) 
{
  flavs[0] = quarkflavour; 
  flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
  flavs[2] = quarkflavour; 
  alpha    = tools->GetASmax();
}

double q_gq::operator()(double z) {return CF*(1.+(1.-z)*(1.-z))/z;} 
double q_gq::GetZ()      
{
  return zmin*pow(zmax/zmin,AMATOOLS::ran.Get());
}

double q_gq::GetCoupling()         { return alpha;}
double q_gq::GetCoupling(double t) { return tools->AlphaS(t);}
double q_gq::GetWeight(double z,double pt2,bool massterm) 
{ 
  if (!massterm) return 0.5*(1.+AMATOOLS::sqr(1.-z));
  return ( (1.+AMATOOLS::sqr(1.-z))/2. - 
	   (1.-z)*AMATOOLS::sqr(z*flavs[0].PSMass())/
	   (pt2+AMATOOLS::sqr(z*flavs[0].PSMass())) );
}
                   
double q_gq::CrudeInt(double _zmin, double _zmax) {
  zmin = _zmin;
  zmax = _zmax;
  return 2.*CF*alpha*log(zmax/zmin);
}
