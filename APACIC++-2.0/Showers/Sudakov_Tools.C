#include "Sudakov_Tools.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include "Message.H"

using namespace APACIC;
using namespace MODEL;
using namespace ATOOLS;


Sudakov_Tools::Sudakov_Tools(MODEL::Model_Base * _model) {
  scalefac = 1.;
  p_as     = _model->GetScalarFunction(std::string("alpha_S"));
  p_aqed   = _model->GetScalarFunction(std::string("alpha_QED"));
  FixLambda2(sqr((Flavour(kf::Z)).Mass()));
  if (msg.LevelIsDebugging()) Output();
}

Sudakov_Tools::Sudakov_Tools(int _scheme,MODEL::Model_Base * _model,double tmin, double tmax) {
  p_as   = _model->GetScalarFunction(std::string("alpha_S"));
  p_aqed = _model->GetScalarFunction(std::string("alpha_QED"));
  scheme = _scheme;
  if (scheme>0) {
    alphaQEDmax = (*p_aqed)(tmax);    
    if (.25*tmin<static_cast<Running_AlphaS*>(p_as)->CutQ2()) {
      double cutq2 = static_cast<Running_AlphaS*>(p_as)->CutQ2();
      alphaSmax   = AlphaS(cutq2); 
    }
    else {
      //    alphaSmax   = AlphaS(tmin);
      alphaSmax   = AlphaS(0.25*tmin);      
    }
    FixLambda2(sqr((Flavour(kf::Z)).Mass())); 
    Setscalefac(tmin);   
  }
  else {
    alphaQEDmax     = 1./128.;     
    alphaSmax       = 0.2;         
    beta0 = lambda2 = 0.;          
    scalefac        = 1.;          
  }
}

double Sudakov_Tools::CrudeAlphaS(double t){
  if (t<0.) t = -t;
  return scalefac/(beta0*log(t/lambda2));
};

double Sudakov_Tools::AlphaS(double t){
  if (t<0.) t = -t;

  return (*p_as)(t);

  // const double b   =0.6100939485; 
  const double lam2 = sqr(0.29); 
  double thr[7];
  int nf =5;
  thr[0]=thr[1]=thr[2]=thr[3]=0.;
  thr[4]=sqr(Flavour(kf::c).PSMass());
  thr[5]=sqr(Flavour(kf::b).PSMass());
  thr[6]=sqr(Flavour(kf::t).PSMass());
  while (t<thr[nf]) {
    --nf;
  }
  double b_eff=1./(12.* M_PI) * (33. - 2.*nf);
  double alp= 1./(b_eff*log(t/lam2));

  return alp;
}

double Sudakov_Tools::Alpha(double t){
  if (t<0.) t = -t;
  return (*p_aqed)(t);
}

void Sudakov_Tools::FixLambda2(double t) { 
  beta0   = as->Beta0(t)/M_PI;  
  lambda2 = t*exp(-1./(beta0*AlphaS(t)));
}

void Sudakov_Tools::Setscalefac(double t0) {
  if (t0<0.) t0=-t0;
  scalefac = 1.; 
  scalefac = AlphaS(t0)/CrudeAlphaS(t0);
}

void Sudakov_Tools::Output() {
  msg_Debugging()<<"Initialise Sudakov-Tools with scheme : "<<scheme<<std::endl
		 <<"beta0      = "<<beta0<<std::endl
		 <<"lambda2    = "<<lambda2<<std::endl	
		 <<"alphaS(MZ) = "
		 <<CrudeAlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (estimated)"<<std::endl
		 <<"alphaS(MZ) = "
		 <<AlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (exact)"<<std::endl;
  msg_Debugging()<<" scalefac="<<scalefac<<std::endl;
}
