#include "HADRON_RESCATTERING/XSecs/MesonMeson.H"
#include "HADRON_RESCATTERING/XSecs/PiPi.H"

#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

MesonMeson::MesonMeson() :
  m_test(true)
{
  InitialiseResonances();
  if (m_test) { Tests(); exit(1); }
}

MesonMeson::~MesonMeson() {}

void MesonMeson::InitialiseResonances() {
  Flavour pizero  = Flavour(kf_pi);
  Flavour piplus  = Flavour(kf_pi_plus);
  Flavour piminus = Flavour(kf_pi_plus).Bar();
  Flavour Kzero   = Flavour(kf_K);
  Flavour KS      = Flavour(kf_K_S);
  Flavour KL      = Flavour(kf_K_L);
  Flavour Kplus   = Flavour(kf_K_plus);
  // pi^+ pi^- resonances: rho(770), f_2(1270), f_0(980), f_0(600)
  m_resonances[HR_Res::pip_pim].push_back(HR_Resonance(Flavour(kf_rho_770),
						       piplus,piminus,1./2.,1.) );
  m_resonances[HR_Res::pip_pim].push_back(HR_Resonance(Flavour(kf_f_2_1270),
						       piplus,piminus,1./3.,0.84) );
  m_resonances[HR_Res::pip_pim].push_back(HR_Resonance(Flavour(kf_f_0_980),
						       piplus,piminus,1./3.,0.781) );
  m_resonances[HR_Res::pip_pim].push_back(HR_Resonance(Flavour(kf_f_0_600),
						       piplus,piminus,1./3.,1.,
						       0.559, 0.370) );
  // pi^+ pi^0 resonances: rho^+(770)
  m_resonances[HR_Res::pip_pi0].push_back(HR_Resonance(Flavour(kf_rho_770_plus),
						       piplus,pizero,1./2.,1.) );
  // pi^0 pi^0 resonances: f_0(980), f_0(600)
  // Note: CG factor of 1/3 obtains extra symmetrisation of 1/2
  m_resonances[HR_Res::pi0_pi0].push_back(HR_Resonance(Flavour(kf_f_0_980),
						       pizero,pizero,1./6.,0.781) );
  m_resonances[HR_Res::pi0_pi0].push_back(HR_Resonance(Flavour(kf_f_0_600),
						       pizero,pizero,1./6.,1.,
						       0.559, 0.370) );
  // K^+ pi^- resonances: K^*(892), K_0^*(1430), K_2^*(1430), K^*(1680)
  m_resonances[HR_Res::Kp_pim].push_back(HR_Resonance(Flavour(kf_K_star_892),
						      Kplus,piminus,2./3.,0.999) );
  m_resonances[HR_Res::Kp_pim].push_back(HR_Resonance(Flavour(kf_K_0_star_1430),
						      Kplus,piminus,2./3.,0.93) );
  m_resonances[HR_Res::Kp_pim].push_back(HR_Resonance(Flavour(kf_K_2_star_1430),
						      Kplus,piminus,2./3.,0.499) );
  m_resonances[HR_Res::Kp_pim].push_back(HR_Resonance(Flavour(kf_K_star_1680),
						      Kplus,piminus,2./3.,0.93) );
  // K^+ pi^0 resonances: K^{*+}(892), K_0^{*+}(1430), K_2^{*+}(1430), K^{*+}(1680)
  m_resonances[HR_Res::Kp_pi0].push_back(HR_Resonance(Flavour(kf_K_star_892_plus),
						      Kplus,piminus,1./3.,0.999) );
  m_resonances[HR_Res::Kp_pi0].push_back(HR_Resonance(Flavour(kf_K_0_star_1430_plus),
						      Kplus,piminus,1./3.,0.93) );
  m_resonances[HR_Res::Kp_pi0].push_back(HR_Resonance(Flavour(kf_K_2_star_1430_plus),
						      Kplus,piminus,1./3.,0.499) );
  m_resonances[HR_Res::Kp_pi0].push_back(HR_Resonance(Flavour(kf_K_star_1680_plus),
						      Kplus,piminus,1./3.,0.387) );
  // K^0 pi^+ resonances: K^{*+}(892), K_0^{*+}(1430), K_2^{*+}(1430), K^{*+}(1680)
  m_resonances[HR_Res::K0_pip].push_back(HR_Resonance(Flavour(kf_K_star_892_plus),
						      Kplus,piminus,2./3.,0.999) );
  m_resonances[HR_Res::K0_pip].push_back(HR_Resonance(Flavour(kf_K_0_star_1430_plus),
						      Kplus,piminus,2./3.,0.93) );
  m_resonances[HR_Res::K0_pip].push_back(HR_Resonance(Flavour(kf_K_2_star_1430_plus),
						      Kplus,piminus,2./3.,0.499) );
  m_resonances[HR_Res::K0_pip].push_back(HR_Resonance(Flavour(kf_K_star_1680_plus),
						      Kplus,piminus,2./3.,0.387) );
  // K^0 pi^0 resonances: none?  Will need to check this!
}

double MesonMeson::XStot(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
			 const double & s) {
  if (!(A.IsMeson() && B.IsMeson())) return 0.;
  HR_Res psps    = HR_Res::none;
  // pi^+ pi^- scattering
  if ( A.Kfcode()==211 && B.Kfcode()==211 && 
       (( A.IsAnti() && !B.IsAnti()) ||
	(!A.IsAnti() &&  B.IsAnti())) )            psps = HR_Res::pip_pim;
  // pi^+/- pi^0 scattering
  else if ( (A.Kfcode()==211 && B.Kfcode()==111) ||
	    (A.Kfcode()==111 && B.Kfcode()==211) ) psps = HR_Res::pip_pi0;
  // pi^0 pi^0 scattering
  else if ( (A.Kfcode()==111 && B.Kfcode()==111) ) psps = HR_Res::pip_pi0;
  // K^+/- pi^-/+ scattering
  else if ( A.Kfcode()==321 && B.Kfcode()==211 && 
	    (( A.IsAnti() && !B.IsAnti()) ||
	     (!A.IsAnti() &&  B.IsAnti())) )       psps = HR_Res::Kp_pim;
  else if ( A.Kfcode()==321 && B.Kfcode()==111 )   psps = HR_Res::Kp_pi0;
  else if ( ((A.Kfcode()==311 ||
	      A.Kfcode()==310 || A.Kfcode()==130) &&
	     B.Kfcode()==211) ||
	    ((B.Kfcode()==311 ||
	      B.Kfcode()==310 || B.Kfcode()==130) &&
	     A.Kfcode()==211) )                    psps = HR_Res::K0_pip;
  map<HR_Res,HR_resonances>::iterator lrit = m_resonances.find(psps);
  double xstot = 5.;
  msg_Out()<<METHOD<<" ["<<A<<", "<<B<<", E = "<<sqrt(s)<<"] : "<<int(psps)<<"\n";
  if (lrit!=m_resonances.end()) {
    for (HR_resonances::iterator rit=lrit->second.begin();
	 rit!=lrit->second.end();rit++) xstot += rit->XStot_R(s);
  }
  return xstot;
}

double MesonMeson::XSel(const ATOOLS::Flavour & A,const ATOOLS::Flavour & B,
			const double & s) {
  return 0.;
}

void MesonMeson::Tests() {
  PiPi pipi;
}

