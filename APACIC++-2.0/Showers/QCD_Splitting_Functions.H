#ifndef QCD_Splitting_Function_H
#define QCD_Splitting_Function_H
//----------------------------------------------------------------------
//               QCD_Splitting_Functions.H
//----------------------------------------------------------------------
#include "Splitting_Function.H"
#include "MathTools.H"
#include "Random.H"

namespace APACIC {

  //! quark to quark + gluon  splitting function
  class q_qg : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
  public:
    q_qg(APHYTOOLS::Flavour quarkflavour) : tools(0) {
      flavs[0] = quarkflavour; 
      flavs[1] = quarkflavour; 
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
      alpha    = 1.;
    };
    q_qg(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
      tools(_tools) {
      flavs[0] = quarkflavour; 
      flavs[1] = quarkflavour; 
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon);
      alpha    = tools->GetASmax();
      cout<<" q_qg: alpha ="<<alpha<<endl;
    };
    double operator()(double z) {return CF*(1.+z*z)/(1.-z);};             
    double GetZ()      {
      return 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),AMATOOLS::Ran.get());   
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->AlphaS(t);};
    double GetWeight(double z,double pt2,bool massterm) { 
      if (!massterm) return 0.5*(1.+z*z);
      return ( (1.+z*z)/2. - 
	       z*AMATOOLS::sqr((1.-z)*flavs[0].PSmass())/
	       (pt2+AMATOOLS::sqr((1.-z)*flavs[0].PSmass())) );
    };                   
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return 2.*CF*alpha*log((1.-zmin)/(1.-zmax));                              
    };
  };

  //! gluon to gluon + gluon splitting function (needed twice for initial state shower)
  class g_gg : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
  public:
    g_gg() : tools(0) {
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      alpha    = 1.;
    };
    g_gg(Sudakov_Tools * _tools) : tools(_tools) { 
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      alpha    = tools->GetASmax();
      cout<<" q_gg: alpha ="<<alpha<<endl;

    };
    double operator()(double z) {
      return CA*AMATOOLS::sqr(1.-z*(1.-z))/(z*(1.-z));};                 
    double GetZ() {
      return 1./(1. + ((1.-zmin)/zmin) *
      		 pow( zmin*(1.-zmax)/((1.-zmin)*zmax), AMATOOLS::Ran.get()));
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->AlphaS(t);};
    double GetWeight(double z,double pt2,bool masses)   { 
      return AMATOOLS::sqr(1.-z*(1.-z));
    };    
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return CA*alpha*log((1.-zmin)*zmax/(zmin*(1.-zmax)));                    
    } 
  };

  //! gluon to quark + anti-quark splitting function
  class g_qq : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
  public:
    g_qq(APHYTOOLS::Flavour quarkflavour): tools(0) {
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[1] = quarkflavour; 
      flavs[2] = quarkflavour.bar(); 
      alpha    = 1.;
    };
    g_qq(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) :
      tools (_tools) {
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[1] = quarkflavour; 
      flavs[2] = quarkflavour.bar(); 
      alpha    = tools->GetASmax();
    };
    double operator()(double z) {
      return TR*(z*z + (1.-z)*(1.-z));};                                 
    double GetZ()      {
      return zmin+(zmax-zmin)*AMATOOLS::Ran.get();                       
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->AlphaS(t);};
    double GetWeight(double z,double pt2,bool masses) { 
      if (!masses) return (*this)(z)/TR;
      return (1. - 2.*z*(1.-z)*(1.- pt2/(pt2+AMATOOLS::sqr(flavs[1].PSmass())))); // /TR;
    };                 
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return TR*alpha*(zmax-zmin);                                             
    };
  };


  //! quark to qluon + quark splitting function (only used in Initial State Shower)
  class q_gq : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
  public:
    q_gq(APHYTOOLS::Flavour quarkflavour) {
      flavs[0] = quarkflavour; 
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[2] = quarkflavour; 
      alpha    = 1.;
    };
    q_gq(APHYTOOLS::Flavour quarkflavour,Sudakov_Tools * _tools) : 
      tools(_tools) {
      flavs[0] = quarkflavour; 
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::gluon); 
      flavs[2] = quarkflavour; 
      alpha    = tools->GetASmax();
    };
    double operator()(double z) {return CF*(1.+(1.-z)*(1.-z))/z;}; 
    double GetZ()      {
      return zmin*pow(zmax/zmin,AMATOOLS::Ran.get());
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->AlphaS(t);};
    double GetWeight(double z,double pt2,bool massterm) { 
      if (!massterm) return 0.5*(1.+AMATOOLS::sqr(1.-z));
      return ( (1.+AMATOOLS::sqr(1.-z))/2. - 
	       (1.-z)*AMATOOLS::sqr(z*flavs[0].PSmass())/
	       (pt2+AMATOOLS::sqr(z*flavs[0].PSmass())) );
    };                   
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return 2.*CF*alpha*log(zmax/zmin);
    };
  };

}
#endif








