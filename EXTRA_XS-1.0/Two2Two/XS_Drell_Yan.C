#include "XS_Drell_Yan.H"
#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


/* 
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha. Hence one Pi remains in the game.
*/

XS_qqbar_pg::XS_qqbar_pg(int _nin,int _nout,APHYTOOLS::Flavour * _fl)
  : Single_XS(_nin,_nout,_fl) 
{ 
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0; 

  barred = _fl[0].isanti();

  colours[0][barred] = colours[2][barred] = 500; 
  colours[1][1-barred] = colours[2][1-barred] = 501;
  name=" q qbar -> g  photon ";
} 

double XS_qqbar_pg::operator()(double s,double t,double u) {  
  if (s<thres) return 0.;
  return 8. * (t*t + u*u + 2. * s * ( s + t + u)) / ( 9. * t*u);
} 

bool XS_qqbar_pg::SetColours(double s,double t,double u) { 
  //  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}

//======================================================================

XS_qg_qp::XS_qg_qp (int _nin,int _nout,APHYTOOLS::Flavour * _fl)
  : Single_XS(_nin,_nout,_fl) 
{ 
  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0; 

  barred = _fl[0].isanti();

  colours[0][barred] = colours[1][1-barred] = 500; 
  colours[1][barred] = colours[2][barred] = 501; 
  name=" q g -> q  photon ";
}
 
double XS_qg_qp::operator()(double s,double t,double u) { 
  if (s<thres) return 0.;
  return (-1.) * (t*t + u*u + 2. * s * ( s + t + u)) / ( 3. * s*u);
} 

bool XS_qg_qp::SetColours(double s,double t,double u) { 
  //  scale = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}

//======================================================================

XS_ee_ffbar::XS_ee_ffbar(int _nin,int _nout,APHYTOOLS::Flavour * _fl) 
  : Single_XS(_nin,_nout,_fl) 
{
  msg.Debugging()<<"In XS_ee_ffbar."<<std::endl;

  MZ2    = sqr(APHYTOOLS::Flavour(APHYTOOLS::kf::Z).mass());
  GZ2    = sqr(APHYTOOLS::Flavour(APHYTOOLS::kf::Z).width());
 
  alpha  = aqed->Aqed((sqr(rpa.gen.Ecms())));
  sin2tw = rpa.consts.Sin2TW();
  if (APHYTOOLS::Flavour(APHYTOOLS::kf::Z).ison()) 
    kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  else
    kappa  = 0.;

  mass     = _fl[2].mass();
  qe       = _fl[0].charge();
  qf       = _fl[2].charge();
  ae       = _fl[0].isoweak();      
  af       = _fl[2].isoweak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  colfac   = 1.;

  kswitch  = 2;  
  fac      = 2./(3.*M_PI);
  fin      = 2.*M_PI/9. - 7./(3.*M_PI) + 9./(3.*M_PI);

  for (short int i=0;i<4;i++) colours[i][0] = colours[i][1] = 0;
  if (_fl[2].isquark()) {
    barred = _fl[2].isanti();
    colours[2][barred] = colours[3][1-barred] = 500;
    colfac = 3.;
  }

  if (_fl[0].isquark())  {
    barred = _fl[0].isanti();
    colours[0][barred] = colours[1][1-barred] = 500;
    colfac  = 1./3.;
    kswitch = 1;
  }
}

double XS_ee_ffbar::operator()(double s,double t,double u) {
  if (s<thres) return 0.;
  chi1  = kappa * s * (s-MZ2)/(sqr(s-MZ2) + GZ2*MZ2);
  chi2  = sqr(kappa * s)/(sqr(s-MZ2) + GZ2*MZ2);

  term1 = (1+sqr(1.+2.*t/s)) * (sqr(qf*qe) + 2.*(qf*qe*vf*ve) * chi1 +
				(ae*ae+ve*ve) * (af*af+vf*vf) * chi2);
  term2 = (1.+2.*t/s) * (4. * qe*qf*ae*af * chi1 + 8. * ae*ve*af*vf * chi2);

  // Divide by two ???? no !
  return sqr(4.*M_PI*alpha) * colfac * (term1+term2); 
}

bool XS_ee_ffbar::SetColours(double s,double t,double u) { 
  //  scale = s;
  return 1; 
}


double XS_ee_ffbar::KFactor(double s) {
  double alphaS = (*as).AlphaS(s);
  if (kswitch==1) return exp(fac*alphaS) * (1. + fin*alphaS);
  if (kswitch==2) return 1. + alphaS/M_PI + (1.986 - 0.115*5.) * sqr(alphaS/M_PI);
  return 1.;
}
