#include "XS_Drell_Yan.H"
#include "Run_Parameter.H"
#include "Running_AlphaQED.H"
#include "Running_AlphaS.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;


/* 
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha. Hence one Pi remains in the game.
*/

XS_qqbar_pg::XS_qqbar_pg(const size_t nin,const size_t nout,const ATOOLS::Flavour *fl)
  : Single_XS(nin,nout,fl) 
{ 
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0; 

  barred = fl[0].IsAnti();

  p_colours[0][barred]   = p_colours[2][barred]   = 500; 
  p_colours[1][1-barred] = p_colours[2][1-barred] = 501;
  m_name=" q qbar -> g  photon ";
} 

double XS_qqbar_pg::operator()(double s,double t,double u) {  
  if (s<m_threshold) return 0.;
  return 8. * (t*t + u*u + 2. * s * ( s + t + u)) / ( 9. * t*u);
} 

bool XS_qqbar_pg::SetColours(double s,double t,double u) { 
  m_scale[PHASIC::stp::as] = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}

//======================================================================

XS_qg_qp::XS_qg_qp (const size_t nin,const size_t nout,const ATOOLS::Flavour *fl)
  : Single_XS(nin,nout,fl) 
{ 
  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0; 

  barred = fl[0].IsAnti();

  p_colours[0][barred] = p_colours[1][1-barred] = 500; 
  p_colours[1][barred] = p_colours[2][barred] = 501; 
  m_name=" q g -> q  photon ";
}
 
double XS_qg_qp::operator()(double s,double t,double u) { 
  if (s<m_threshold) return 0.;
  return (-1.) * (t*t + u*u + 2. * s * ( s + t + u)) / ( 3. * s*u);
} 

bool XS_qg_qp::SetColours(double s,double t,double u) { 
  m_scale[PHASIC::stp::as] = (2.*s*t*u)/(s*s+t*t+u*u);
  return 1; 
}

//======================================================================

template <> 
Single_XS *Single_XS::GetProcess<XS_ee_ffbar>(const size_t nin,const size_t nout,
					      const ATOOLS::Flavour *flavours,
					      const size_t nqed, const size_t nqcd)
{
  if ((flavours[2].IsLepton() && flavours[3]==flavours[2].Bar() && flavours[0].IsQuark() && 
       flavours[1]==flavours[0].Bar()) ||   
      (flavours[0].IsLepton() && flavours[1]==flavours[0].Bar() && flavours[2].IsQuark() && 
       flavours[3]==flavours[2].Bar())) { 
    if (nqcd==0 && nqed==2) {
      return new XS_ee_ffbar(nin,nout,flavours);  
    }
  }
  return NULL;
}

XS_ee_ffbar::XS_ee_ffbar(const size_t nin,const size_t nout,const ATOOLS::Flavour *fl) 
  : Single_XS(nin,nout,fl) 
{
  MZ2    = sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass());
  GZ2    = sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Width());
 
  alpha  = aqed->Aqed((sqr(rpa.gen.Ecms())));
  sin2tw = rpa.gen.ScalarConstant(string("sin2_thetaW"));
  if (ATOOLS::Flavour(ATOOLS::kf::Z).IsOn()) 
    kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  else
    kappa  = 0.;

  mass     = fl[2].Mass();
  qe       = fl[0].Charge();
  qf       = fl[2].Charge();
  ae       = fl[0].IsoWeak();      
  af       = fl[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  colfac   = 1.;

  kswitch  = 0;  
  fac      = 2./(3.*M_PI);
  fin      = 2.*M_PI/9. - 7./(3.*M_PI) + 9./(3.*M_PI);

  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (fl[2].IsQuark()) {
    barred = fl[2].IsAnti();
    p_colours[2][barred] = p_colours[3][1-barred] = 500;
    colfac = 3.;
  }

  if (fl[0].IsQuark())  {
    barred = fl[0].IsAnti();
    p_colours[0][barred] = p_colours[1][1-barred] = 500;
    colfac  = 1./3.;
    kswitch = 1;
  }
}

double XS_ee_ffbar::operator()(double s,double t,double u) {
  if (s<m_threshold) return 0.;
  chi1  = kappa * s * (s-MZ2)/(sqr(s-MZ2) + GZ2*MZ2);
  chi2  = sqr(kappa * s)/(sqr(s-MZ2) + GZ2*MZ2);

  term1 = (1+sqr(1.+2.*t/s)) * (sqr(qf*qe) + 2.*(qf*qe*vf*ve) * chi1 +
				(ae*ae+ve*ve) * (af*af+vf*vf) * chi2);
  term2 = (1.+2.*t/s) * (4. * qe*qf*ae*af * chi1 + 8. * ae*ve*af*vf * chi2);

  // Divide by two ????
  return sqr(4.*M_PI*alpha) * colfac * (term1+term2); 
}

bool XS_ee_ffbar::SetColours(double s,double t,double u) { 
  m_scale[PHASIC::stp::as] = s;
  return 1; 
}


double XS_ee_ffbar::KFactor(double s) {
  double alphaS = as->AlphaS(s);
  if (kswitch==1) return exp(fac*alphaS) * (1. + fin*alphaS);
  if (kswitch==2) return 1. + alphaS/M_PI + (1.986 - 0.115*5.) * sqr(alphaS/M_PI);
  return 1.;
}
