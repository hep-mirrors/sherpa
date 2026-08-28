#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Scoped_Settings.H"


#include <iostream>


using namespace MODEL;
using namespace ATOOLS;
using namespace std;



namespace MODEL {
  Running_AlphaQED * aqed  = 0 ;

  const double Running_AlphaQED::m_A[4]={0.0,0.0,0.00165,0.00221};
  const double Running_AlphaQED::m_B[4]={0.00835,0.00238,0.00299,0.00293};
  const double Running_AlphaQED::m_C[4]={1.0,3.927,1.0,1.0};
 #ifdef USING__HADALPHAQED
    extern "C"{
      // hadr5x declares its arguments `real` (hadr5x19.f:54), but the file is
      // built with -fdefault-real-8, so `real` is 8 bytes and double* is the
      // correct prototype. (dhadr5x_qedc19 is the explicit real*8 wrapper and
      // would work too; it is redundant under these flags.) If those Fortran
      // flags ever change, this breaks silently - Tests/Unit/
      // test_MODEL_Running_AlphaQED.C pins the value that would catch it.
      void hadr5x_(double *e, double *st2,
                  double *der, double *errdersta,
                  double *errdersys, double *deg,
                  double *errdegsta, double *errdegsys);
    }
    #endif

}



Running_AlphaQED::Running_AlphaQED(const double _alpha0) :
  m_alpha0(_alpha0)
{
  m_type = std::string("Running Coupling");
  m_name  = "Alpha_QED";
  m_defval = _alpha0;
  Scoped_Settings s{ Settings::GetMainSettings()[m_name] };
  m_mode     = s["VPMODE"].SetDefault(vpmode::off).Get<vpmode::code>();
}


double Running_AlphaQED::operator()(double t)
{
  double Q2    = t;
  if (t<0.) Q2 = -t; 

  int i = 3;
  if (Q2<0.3)        i=0;
  else if (Q2<3.0)   i=1;
  else if (Q2<100.0) i=2;

  if(IsZero(t) || m_mode==vpmode::off) return m_alpha0;

  double sig_lep_gg = m_alpha0/(3.*M_PI) * 
    (PiGamma(Flavour(kf_e),Q2)+PiGamma(Flavour(kf_mu),Q2)+PiGamma(Flavour(kf_tau),Q2));
  double sig_ha_gg  = m_A[i] + m_B[i]*log(1+m_C[i]*Q2);
  double sig_top_gg = m_alpha0/(3.*M_PI) * 3. * (PiGamma(Flavour(kf_t),Q2));
  double sigma_gg   = sig_lep_gg+sig_ha_gg+sig_top_gg;

  #ifdef USING__HADALPHAQED
    double delta_r(0.),errdersta,errdersys,deg,errdegsta,errdegsys;
    double sin2 = 0.23153;// MODEL::m_model->ComplexConstant("csin2_thetaW").real();
    if(m_mode!=vpmode::lp){
      // The routine takes E = signed ENERGY in GeV:
      //   "E energy ( momentum transfer ): E>0 timelike , E<0 spacelike"
      //   (hadr5x19.f:40)
      // NOT Q^2 in GeV^2. Passing Q^2 straight through evaluated
      // Delta_alpha_had at E = 8315 GeV instead of 91.19 at the Z pole (giving
      // 0.0542 where the measured value is 0.02766). The argument of this
      // operator is the signed Lorentz invariant - s>0 timelike, t<0 spacelike -
      // so the energy is sign(t)*sqrt(|t|).
      double e = (t<0. ? -sqrt(Q2) : sqrt(Q2));
      hadr5x_(&e, &sin2, &delta_r, &errdersta,
              &errdersys, &deg, &errdegsta, &errdegsys);
    }
    switch(m_mode){
      case vpmode::full:
        // delta_r is the 5-flavour HADRONIC piece alone ("the routine returns
        // the hadronic contribution of 5 flavors (u,d,s,c,b)", hadr5x19.f:43),
        // so the leptonic and top parts have to be added explicitly. They used
        // to be dropped here, which made `full` identical to `hp` and silently
        // discarded the leptonic contribution - at M_Z the LARGER half
        // (0.0314 leptonic vs 0.0277 hadronic).
        sigma_gg = sig_lep_gg + delta_r + sig_top_gg;
        break;
      case vpmode::lp:
        sigma_gg = sig_lep_gg;
        break;
      case vpmode::hp:
        sigma_gg = delta_r;
        break;
      case vpmode::off:
        break;   // unreachable, handled by the early return above
      }
  #endif

  return m_alpha0/(1.-sigma_gg);
}  

double Running_AlphaQED::PiGamma(const Flavour & fl,double scale) {
  double mass2  = sqr(fl.Mass(true)); // onshell mass
  if(mass2==0.) THROW(fatal_error, "Cannot evolve QED coupling with zero fermion masses");    
  double mqs    = mass2/scale;
  if (scale==0.) return 0.;
  if (4.*mqs<1.e-3) return (-5./3.-log(mqs));
  else if (4.*mqs<=1.) {
    double beta = sqrt(1.-4.*mqs);
    return (1./3.-(1.+2.*mqs)*(2.+beta*log((1.-beta)/(1.+beta))));
  } 
  else { 
    return 0.; 
  }
}

void Running_AlphaQED::PrintSummary()
{
  msg_Info()<<"Set \\alpha according to EW scheme"
            <<"\n  1/\\alpha(0)   = "<<1./m_alpha0
            <<"\n  1/\\alpha(def) = "<<1./m_defval<<"\n";
}

std::istream &MODEL::operator>>(std::istream &str,vpmode::code &vp)
{
  std::string tag;
  str>>tag;
  vp=vpmode::off;
  if      (tag.find("None")!=std::string::npos) vp=vpmode::off;
  else if (tag.find("Full")!=std::string::npos) vp=vpmode::full;
  else if (tag.find("1")!=std::string::npos)    vp=vpmode::full;
  else if (tag.find("HP")!=std::string::npos)   vp=vpmode::hp;
  else if (tag.find("2")!=std::string::npos)    vp=vpmode::hp;
  else if (tag.find("LP")!=std::string::npos)   vp=vpmode::lp;
  else if (tag.find("3")!=std::string::npos)    vp=vpmode::lp;
  return str;
}