#include "Sudakov_Tools.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"
#include "Run_Parameter.H"
#include "Function_Base.H"

#include "Message.H"

using namespace APACIC;
using namespace MODEL;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;


Sudakov_Tools::Sudakov_Tools() {
  scalefac = 1.;
  FixLambda2(sqr((Flavour(kf::Z)).Mass()));
  if (rpa.gen.Debugging()) Output();
}

Sudakov_Tools::Sudakov_Tools(int _scheme,double tmin, double tmax) {
  scheme = _scheme;
  if (scheme>0) {
    alphaQEDmax = (*aqed)(tmax);    // max alpha_S 
    alphaSmax   = AlphaS(tmin);      // max alpha_S 
    FixLambda2(sqr((Flavour(kf::Z)).Mass()));                   
    // determine Lambda2 and beta0 at MZ;
    Setscalefac(tmin);              // fix scaling appropriate for q02/4
    //    Setscalefac(sqr((Flavour(kf::Z)).Mass()));              // fix scaling only for pythia alphaS !!!!!!
  }
  else {
    alphaQEDmax     = 1./128.;      // max alpha_QED has to be read in from
    alphaSmax       = 0.2;          // max alpha_S has to be read in from
    beta0 = lambda2 = 0.;           // parameter-file ...
    scalefac        = 1.;           // won't be used ....
  }
  if (rpa.gen.Debugging()) { 
    Output();
    cout<<" tmin= "<< tmin<<endl;
    cout<<" alpha_max="<< alphaSmax <<endl;
    cout<<" Checking alphaS "<<endl;
    double q2_max=sqr(91.2);
    double q2_min=sqr(.912);
    int    n =10;
    for (int i=0;i<=n;++i) {
      double q2=q2_min*pow((q2_max/q2_min),double(i)/double(n));
      cout<<" "<<q2<<" \t"<<CrudeAlphaS(q2)<<" \t"<<AlphaS(q2)<<" \t"<<(*as)(q2)<<endl;
    }
    //    exit(0);
  }
}

double Sudakov_Tools::CrudeAlphaS(double t){
  if (t<0.) t = -t;
  return scalefac/(beta0*log(t/lambda2));
};

double Sudakov_Tools::AlphaS(double t){
  if (t<0.) t = -t;

  // exact (LO) alphaS
  return (*as)(t);

  const double b   =0.6100939485; // 1/(12 Pi) * (33 - 2*5)
  const double lam2=sqr(0.29); 
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


  if (t<=0.25) {
    cout<<" small q2 in alphas ="<<t<<endl;
  }

  return alp;

  // effective (LO) alphaS  (according to spacelike pythia shower)
  /*
  const double b   =0.6100939485; // 1/(12 Pi) * (33 - 2*5)
  const double lam2=0.02106116817; //0.16;
  int nf =5;
  double lam2_eff=lam2;
//   cout<<"as b="<<b<<endl;
//   cout<<"as lam2="<<lam2<<endl;
//   cout<<"as t="<<t<<endl;
  double thr[7];
  thr[0]=thr[1]=thr[2]=thr[3]=0.;
  thr[4]=sqr(Flavour(kf::c).PSMass());
  thr[5]=sqr(Flavour(kf::b).PSMass());
  thr[6]=sqr(Flavour(kf::t).PSMass());
  while (t<thr[nf]) {
    //    cout<<thr[nf]<<endl;
    --nf;
    //    cout<<"lam was = "<<lam2_eff<<endl;
    lam2_eff=lam2_eff*pow(thr[nf+1]/lam2_eff,2./(33. - 2.* nf));
    //    cout<<"lam is = "<<lam2_eff<<endl;
  }
  double b_eff=1./(12.* M_PI) * (33. - 2.*nf);
  return 1./(b_eff*log(t/lam2_eff));
  */
};

double Sudakov_Tools::Alpha(double t){
  if (t<0.) t = -t;
  return (*aqed)(t);
};

void Sudakov_Tools::FixLambda2(double t) { 
  beta0   = as->Beta0(t)/M_PI;  // additional factor M_PI (since beta0 is "b" on Weber page 29)
  lambda2 = t*exp(-1./(beta0*AlphaS(t)));
};

void Sudakov_Tools::Setscalefac(double t0) {
  if (t0<0.) t0=-t0;
  scalefac = 1.; // will be used in CrudeAlphaS ...
  scalefac = AlphaS(t0)/CrudeAlphaS(t0);
};

void Sudakov_Tools::Output() {
  msg.Debugging()<<"Initialise Sudakov-Tools with scheme : "<<scheme<<std::endl
		 <<"beta0      = "<<beta0<<std::endl
		 <<"lambda2    = "<<lambda2<<std::endl	
		 <<"alphaS(MZ) = "
		 <<CrudeAlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (estimated)"<<std::endl
		 <<"alphaS(MZ) = "
		 <<AlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (exact)"<<std::endl;
  msg.Debugging()<<" scalefac="<<scalefac<<endl;
}
