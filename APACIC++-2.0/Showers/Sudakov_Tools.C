#include "Sudakov_Tools.H"
#include "Run_Parameter.H"
#include "Running_AlphaS.H"
#include "Running_AlphaQED.H"

#include "Message.H"

using namespace APACIC;
using namespace MODEL;
using namespace ATOOLS;


Sudakov_Tools::Sudakov_Tools(MODEL::Model_Base * model) {
  m_scalefac = 1.;
  p_as     = model->GetScalarFunction(std::string("alpha_S"));
  p_aqed   = model->GetScalarFunction(std::string("alpha_QED"));
  m_renormalization_scale_factor = rpa.gen.RenormalizationScaleFactor();
  //  std::cout<<"Sudakov_Tools::Sudakov_Tools("<<m_renormalization_scale_factor<<std::endl;
  FixLambda2(sqr((Flavour(kf::Z)).Mass()));
  //  if (msg.LevelIsDebugging()) 
  Output();
}

Sudakov_Tools::Sudakov_Tools(int _m_scheme,MODEL::Model_Base * model,double tmin, double tmax) {
  p_as   = model->GetScalarFunction(std::string("alpha_S"));
  p_aqed = model->GetScalarFunction(std::string("alpha_QED"));
  m_renormalization_scale_factor = rpa.gen.RenormalizationScaleFactor();
  m_scheme = _m_scheme;

  //  tmin*=m_renormalization_scale_factor;
  //  tmax*=m_renormalization_scale_factor;

  if (m_scheme>0) {
    m_alphaQEDmax = (*p_aqed)(tmax*m_renormalization_scale_factor);    
    if (.25*tmin*m_renormalization_scale_factor<static_cast<Running_AlphaS*>(p_as)->CutQ2()) {
      double cutq2 = static_cast<Running_AlphaS*>(p_as)->CutQ2();
      m_alphaSmax   = AlphaS(cutq2); 
    }
    else {
      //    m_alphaSmax   = AlphaS(tmin);
      m_alphaSmax   = AlphaS(0.25*tmin);      
    }
    FixLambda2(sqr((Flavour(kf::Z)).Mass())); 
    Setscalefac(tmin);   
  }
  else {
    m_alphaQEDmax       = 1./128.;     
    m_alphaSmax         = 0.2;         
    m_beta0 = m_lambda2 = 0.;          
    m_scalefac          = 1.;          
  }
  //  std::cout<<"Sudakov_Tools::Sudakov_Tools("<<m_renormalization_scale_factor<<std::endl;
  //if (msg.LevelIsDebugging()) 
  Output();
}

double Sudakov_Tools::CrudeAlphaS(double t){
  t*=m_renormalization_scale_factor;
  if (t<0.) t = -t;
  return m_scalefac/(m_beta0*log(t/m_lambda2));
};

double Sudakov_Tools::AlphaS(double t){
  //  std::cout<<"Sudakov_Tools::AlphaS("<<t<<")*"<<m_renormalization_scale_factor<<"\n";
  t*=m_renormalization_scale_factor;
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

double Sudakov_Tools::Nf(double t){
  //  std::cout<<"Sudakov_Tools::Nf("<<t<<")*"<<m_renormalization_scale_factor<<"\n";
  t*=m_renormalization_scale_factor;
  if (t<0.) t = -t;

  MODEL::Running_AlphaS* as=static_cast<MODEL::Running_AlphaS*>(p_as);
  return as->Nf(t);
}

double Sudakov_Tools::Alpha(double t){
  t*=m_renormalization_scale_factor;
  if (t<0.) t = -t;
  return (*p_aqed)(t);
}

void Sudakov_Tools::FixLambda2(double t) { 
  m_beta0   = as->Beta0(t)/M_PI;  
  m_lambda2 = t*exp(-1./(m_beta0*AlphaS(t)));
}

void Sudakov_Tools::Setscalefac(double t0) {
  if (t0<0.) t0=-t0;
  m_scalefac = 1.; 
  m_scalefac = AlphaS(t0)/CrudeAlphaS(t0);
}

void Sudakov_Tools::Output() {
  msg.Out()<<"Initialise Sudakov-Tools with scheme : "<<m_scheme<<std::endl
		 <<"beta0      = "<<m_beta0<<std::endl
		 <<"lambda2    = "<<m_lambda2<<std::endl	
		 <<"alphaS(MZ) = "
		 <<CrudeAlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (estimated)"<<std::endl
		 <<"alphaS(MZ) = "
		 <<AlphaS(sqr((Flavour(kf::Z)).Mass()))
		 <<"  (exact)"<<std::endl;
  msg.Out()<<" scalefac="<<m_scalefac<<std::endl;
}
