#ifndef QED_Splitting_Function_H
#define QED_Splitting_Function_H
//----------------------------------------------------------------------
//               QED_Splitting_Functions.H
//----------------------------------------------------------------------
#include "Splitting_Function.H"
#include "MathTools.H"
#include "Random.H"

namespace APACIC {
  //! fermion to fermion + photon  splitting function
  class f_fp : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
    double          qsqr;
  public:
    f_fp(APHYTOOLS::Flavour fermionflavour) {
      flavs[0] = fermionflavour; // a
      flavs[1] = fermionflavour; // b
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = 1.;
    };
    f_fp(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
      tools (_tools) {
      flavs[0] = fermionflavour; // a
      flavs[1] = fermionflavour; // b
      flavs[2] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = tools->GetAQEDmax();
    };
    double operator()(double z) {return qsqr*(1.+z*z)/(1.-z);}; 
    double GetZ()      {
      return 1.-(1.-zmin)*pow((1.-zmax)/(1.-zmin),AMATOOLS::Ran.get());
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->Alpha(t);};
    double GetWeight(double z,double pt2,bool massterm) { 
      if (!massterm) return 0.5*(1.+z*z);
      return ( (1.+z*z)/2. - 
	       z*AMATOOLS::sqr((1.-z)*flavs[0].PSmass())/
	       (pt2+AMATOOLS::sqr((1.-z)*flavs[0].PSmass())) );
    };                   
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return 2.*qsqr*alpha*log((1.-zmin)/(1.-zmax));
    };
  };

  //! fermion to photon + fermion splitting function (only used in Initial State Shower)
  class f_pf : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
    double          qsqr;
  public:
    f_pf(APHYTOOLS::Flavour fermionflavour) {
      flavs[0] = fermionflavour; // a
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // b
      flavs[2] = fermionflavour; // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = 1.;
    };
    f_pf(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
      tools (_tools) {
      flavs[0] = fermionflavour; // a
      flavs[1] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // b
      flavs[2] = fermionflavour; // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = tools->GetAQEDmax();
    };
    double operator()(double z) {return qsqr*(1.+(1.-z)*(1.-z))/z;}; 
    double GetZ()      {
      return zmin*pow(zmax/zmin,AMATOOLS::Ran.get());
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->Alpha(t);};
    double GetWeight(double z,double pt2,bool massterm) { 
      if (!massterm) return 0.5*(1.+AMATOOLS::sqr(1.-z));
      return ( (1.+AMATOOLS::sqr(1.-z))/2. - 
	       (1.-z)*AMATOOLS::sqr(z*flavs[0].PSmass())/
	       (pt2+AMATOOLS::sqr(z*flavs[0].PSmass())) );
    };                   
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return 2.*qsqr*alpha*log(zmax/zmin);
    };
  };

  //! photon to fermion + anti-fermion splitting function
  class p_ff : public Splitting_Function {
    Sudakov_Tools * tools;
    double          alpha;
    double          qsqr;
  public:
    p_ff(APHYTOOLS::Flavour fermionflavour) {
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // a
      flavs[1] = fermionflavour; // b
      flavs[2] = fermionflavour.bar(); // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = 1.;
    };
    p_ff(APHYTOOLS::Flavour fermionflavour,Sudakov_Tools * _tools) :
      tools (_tools) {
      flavs[0] = APHYTOOLS::Flavour(APHYTOOLS::kf::photon); // a
      flavs[1] = fermionflavour; // b
      flavs[2] = fermionflavour.bar(); // c
      qsqr     = AMATOOLS::sqr(fermionflavour.charge());
      alpha    = tools->GetAQEDmax();
    };
    double operator()(double z) {
      return qsqr*(z*z+ AMATOOLS::sqr(1-z));};
    double GetZ()      {
      return zmin+(zmax-zmin)*AMATOOLS::Ran.get();
    };
    double GetCoupling()         { return alpha;};
    double GetCoupling(double t) { return tools->Alpha(t);};
    double GetWeight(double z,double pt2,bool masses) { 
      if (masses) return (*this)(z)/qsqr;
      return (1. - 2.*z*(1.-z)*pt2/(pt2+AMATOOLS::sqr(flavs[1].PSmass())));
    };                 
    double CrudeInt(double _zmin, double _zmax) {
      zmin = _zmin;
      zmax = _zmax;
      return (zmax-zmin)*alpha*qsqr;
    };
  };
}
#endif

