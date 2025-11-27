#include "HADRON_RESCATTERING/XSecs/NNfits.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;


NNfits::NNfits() :
  m_mp(Flavour(kf_p_plus).HadMass()), m_mp2(sqr(m_mp)),
  m_mn(Flavour(kf_n).HadMass()), m_mn2(sqr(m_mn)),
  m_bp(2.3), m_bn(2.3), m_alphaP(0.25), m_s0(4.)
{}

//////////////////////////////////////////////////////////////////////////////
//
// Neutron-proton total and elastic cross section and B parameter
//
// Total cross section:
// - for E>5 GeV: HPR1R2 parametrization
//                PDG, Chin.Phys. C40 (2016) 591
// - for lower energies: GiBUU parametrization, inherited from
//                       Cugnon, Hote, Vandermeulen
//                       Nucl.Inst.Meth. B111 (1996) 215, Eqs (1,2)
//                       https://inspirehep.net/literature/430825
//
// Elastic cross section:
// - may want to use some HPR1R2 parametrization from
//   PDG, Chin.Phys. C40 (2016) 591
//   but will need to refine and understand the inputs.
// - for E>5GeV: CERN/HERA parametrization 
//               TODO: Add reference
// - for lower energies: GiBUU parametrization, inherited from
//                       Cugnon, Hote, Vandermeulen
//                       Nucl.Inst.Meth. B111 (1996) 215, Eq (5)
//                       https://inspirehep.net/literature/430825
//                       supplemtend with Eqs. (B1) and (B5) from
//                       O.Buss et al., 1106.1344
//
// B parameter:
// - B = 2 bA + 2bB + 2alpha' log(s/s0) 
// - could alternatively use Cugnon, Hote, Vandermeulen
//   Nucl.Inst.Meth. B111 (1996) 215, Eq (8)
//
//////////////////////////////////////////////////////////////////////////////

double NNfits::pptot(const double & s) {
  if (s<4.*m_mn2) return 0.;
  if (sqrt(s)>=5.) return m_hpr1r2.xs_tot(hpr1r2::pp,s);
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>=1.5)     return 41.+60.*(plab-0.9)*exp(-1.2*plab);
  if (plab>=0.8)     return 23.5+24.6/(1.+exp(-10.*(plab-1.2)));
  if (plab>=0.4)     return 23.5+1000.*pow(plab-0.7,4.);
  if (plab>=0.1)     return 34.*pow(plab/0.4,-2.104);
  return 0.;
}

double NNfits::ppel(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (sqrt(s)>2.6)   return ( 11.9 + 26.9*pow(plab,-1.21) +
			      0.169*sqr(log(plab)) - 1.85*log(plab) );
  if (plab>=2.0)     return 77./(plab+1.5);
  if (plab>=0.8)     return 1250./(plab+50.)-4.*sqr(plab-1.3);
  if (plab>=0.435)   return 23.5+1000.*pow(plab-0.7,4.);
  if (plab>=0.1)     return 5.12*m_mp/(s-4.*m_mp2)+1.67;
  return 0.;
}

double NNfits::ppB(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab  = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>5.) return 4.*m_bp+m_alphaP*log(s/m_s0);
  if (plab>2.) return 5.334+0.67*(plab-2.);
  double plab8 = pow(plab,8.); 
  return 5.5*plab8/(7.7+plab8);
}

//////////////////////////////////////////////////////////////////////////////
//
// Neutron-proton total and elastic cross section and B parameter
//
// Total cross section:
// - for E>5 GeV: HPR1R2 parametrization
//                PDG, Chin.Phys. C40 (2016) 591
// - for lower energies: GiBUU parametrization, inherited from
//                       Cugnon, Hote, Vandermeulen
//                       Nucl.Inst.Meth. B111 (1996) 215, Eqs (3,4)
//                       https://inspirehep.net/literature/430825
//
// Elastic cross section:
// - may want to use some HPR1R2 parametrization from
//   PDG, Chin.Phys. C40 (2016) 591
//   but will need to refine and understand the inputs.
// - for E>5GeV: CERN/HERA parametrization 
//               TODO: Add reference
// - for lower energies: GiBUU parametrization, inherited from
//                       Cugnon, Hote, Vandermeulen
//                       Nucl.Inst.Meth. B111 (1996) 215, Eq (6)
//                       https://inspirehep.net/literature/430825
//                       supplemtend with Eqs. (B1) and (B5) from
//                       O.Buss et al., 1106.1344
//
// B parameter:
// - B = 2 bA + 2bB + 2alpha' log(s/s0) 
// - could alternatively use Cugnon, Hote, Vandermeulen
//   Nucl.Inst.Meth. B111 (1996) 215, Eq (8)
//
//////////////////////////////////////////////////////////////////////////////

double NNfits::pntot(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s-m_mp2-m_mn2)-4.*m_mn2*m_mp2)/(m_mn+m_mp);
  if (sqrt(s)>=5.) return m_hpr1r2.xs_tot(hpr1r2::pn,s);
  if (plab>2.)     return 42.;
  if (plab>1.)     return 24.2+8.9*plab;
  if (plab>0.4)    return 33.+196.*pow(dabs(plab-0.95),2.5);
  if (plab>0.05)   return 6.3555*pow(plab,-3.2481)*exp(-0.377*sqr(log(plab)));
  return 0.;
}

double NNfits::pnel(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s-m_mp2-m_mn2)-4.*m_mn2*m_mp2)/(m_mn+m_mp);
  if (sqrt(s)>2.6)   return ( 11.9 + 26.9*pow(plab,-1.21) +
			      0.169*sqr(log(plab)) - 1.85*log(plab) );
  if (plab>=2.0)     return 77./(plab+1.5);
  if (plab>=0.8)     return 31./sqrt(plab);
  if (plab>=0.525)   return 33.+196.*pow(dabs(plab-0.95),2.5);
  if (plab>=0.1)     return 17.05*m_mp/(s-4.*m_mp2)-6.83;
  return 0.;
}

double NNfits::pnB(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab  = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>1.6)   return ppB(s);
  if (plab>0.6)   return -1.63*plab+7.16;
  if (plab>0.225) return 16.53*(plab-0.225);
  return 0.;
}


//////////////////////////////////////////////////////////////////////////////
//
// Neutron-proton total and elastic cross section and B parameter
//
// Total cross section:
// - for E>6.5 GeV: HPR1R2 parametrization
//                  PDG, Chin.Phys. C40 (2016) 591
// - for lower energies: URQMD parametrization
//   transition point chosen with help of
//   Sjostrand & Utheim, 2005.05658
//
// Elastic cross section:
// - for E>5GeV: CERN/HERA parametrization 
//               TODO: Add reference
// - for lower energies: URQMD parametrization according to
//                       Sjostrand & Utheim, 2005.05658
// - could use alternative GiBUU parametrisation from
//   O.Buss et al., 1106.1344, Eq (B21)
//   modified parameters of original
//   Cugnon, Vandermeulen, Ann.Phys. (France) 14 (1989) 49-88
//   https://orbi.uliege.be/bitstream/2268/210870/1/94.pdf, Eq (2.5c).
//
// B parameter:
// - B = 2 bA + 2bB + 2alpha' log(s/s0)
//   identical to pp scattering.
//
//////////////////////////////////////////////////////////////////////////////

double NNfits::ppbartot(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  //if (plab>=6.5) return m_hpr1r2.xs_tot(hpr1r2::ppbar,s);
  if (plab>=5) return 38.4+77.6*pow(plab, -0.64)+0.26*sqr(log(plab))-1.2*log(plab);
  if (plab>0.3) return 75.+43.1/plab+2.6/sqr(plab)-3.9*plab;
  return 271.6*exp(-1.1*sqr(plab));
}

double NNfits::ppbarel(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>5.) return (10.2 + 52.7*pow(plab,-1.16) +
			  0.125*sqr(log(plab)) - 1.28*log(plab) );
  if (plab>0.3)   return 31.6+18.3/plab-1.1/sqr(plab)-3.8*plab;
  return 78.6;
  // alternative - may want to check that one:
  // return 40.*pow(plab,-0.56)+5.8*exp(-sqr(plab-1.85));  
}

double NNfits::ppbarB(const double & s) {
  return ppB(s);
}

//////////////////////////////////////////////////////////////////////////////
//
// Baryon-antibaryon annihilation cross sections taken from
// Cugnon, Vandermeulen, Ann.Phys. (France) 14 (1989) 49-88
// https://orbi.uliege.be/bitstream/2268/210870/1/94.pdf, Eq (2.5a)
// and O.Buss et al., 1106.1344, Eqs (B19-B21)
//
// Charge exchange cross section taken from
// Cugnon, Vandermeulen, Ann.Phys. (France) 14 (1989) 49-88
// https://orbi.uliege.be/bitstream/2268/210870/1/94.pdf, Eq (2.5b),
// and B parameter from 
// 
//////////////////////////////////////////////////////////////////////////////

double NNfits::ppbarAnnihil(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>6.34) return 38./sqrt(plab)+24.*pow(plab,-1.1);
  if (plab>0.51) return 88.*pow(plab,-0.4)-24.2;
  if (plab>0.3)  return 51.52*pow(plab,-0.85)+0.0034*pow(plab,-2.94);
  return ppbartot(s);
}

double NNfits::pnbarAnnihil(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>0.382) return 41.4+29./plab;
  return ppbarAnnihil(s);
}

double NNfits::ppbarCEX(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  if (plab>0.5)   return 7.1*pow(plab,-0.9);
  if (plab>0.001) return 10.9*(plab-0.1)*pow(plab,-1.6);
  return 0.;
}

double NNfits::ppbarCEXB(const double & s) {
  if (s<4.*m_mn2) return 0.;
  double plab = sqrt(sqr(s)-4.*m_mp2*s)/(2.*m_mp);
  return 11.*exp(-0.23*plab)+8./(1.+254.*pow(plab,-2.2));
}

