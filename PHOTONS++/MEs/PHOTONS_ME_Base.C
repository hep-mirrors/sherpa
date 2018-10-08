#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Data_Reader.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHOTONS::PHOTONS_ME_Base
#define PARAMETER_TYPE PHOTONS::Particle_Vector_Vector
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

PHOTONS_ME_Base::PHOTONS_ME_Base(const Particle_Vector_Vector& pvv) :
  m_alpha(Photons::s_alpha),
  m_e(sqrt(4.*M_PI*m_alpha)),
  m_GF(1.16639e-5),
  m_sqrt2(1.41421356237),
  m_i(Complex(0.,1.)),
  p_boost(NULL), p_rot(NULL),
  m_pvv_zero(pvv)
{
  Data_Reader reader(" ",";","#","=");
  // NLO EW corrections: 1 - all QED and weak corrections, 0 - only QED corrections
  m_ew = reader.GetValue<int>("YFS_EW_CORRECTIONS",1);
  // EW scheme: 1 = alpha(MZ) (default), 2 = alpha(0), 3 = Gmu
  m_ew_scheme = Photons::s_ew_scheme;
  // Calculate MEs in limit ml^2 << s (only available for some MEs)
  m_limit = reader.GetValue<int>("YFS_SMALL_MASS_LIMIT",0);
  // Whether to calculate NNLO QED corrections (only available for Z and Higgs decays):
  // 0 - no, i.e. only NLO, 1 - yes
  m_nnlo_qed = reader.GetValue<int>("YFS_NNLO_QED",0); 
#ifndef USING__YFS_NNLO
  if (m_nnlo_qed == 1) {
    msg_Error() << METHOD << " YFS NNLO corrections not compiled properly. Please reconfigure "
		<< "with option '--enable-yfs-nnlo'. Will calculate using NLO corrections for the "
		<< "being.\n";
  }
#endif
  // cutoff for stability check
  m_dev = reader.GetValue<double>("YFS_RV_STABILITY",0.1);
  // Factor to scale the dimensionful quantities in RV stability check
  m_xi = reader.GetValue<double>("YFS_RV_SCALING",0.1);  

  // Electroweak parameters.
  double  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  double  MH  = Flavour(kf_h0).Mass();
  double  mt  = Flavour(kf_t).Mass();
  MW2 = pow(MW,2.);
  MZ2 = pow(MZ,2.);
  MH2 = pow(MH,2.);
  mt2 = pow(mt,2.);
  double  GH  = 0.;//Flavour(kf_h0).Width();
  double  GW  = 0.;//*/Flavour(kf_Wplus).Width();
  double  GZ  = 0.;//*/Flavour(kf_Z).Width();
  Complex I   = Complex(0.,1.);
  Complex sw2 = 1.-sqr(MW/MZ);
  Complex cw2 = 1.-sw2;
  m_cW2=MW2/MZ2;
  m_sW2=1.-m_cW2;
  muW2 = MW*(MW-I*GW);
  muZ2 = MZ*(MZ-I*GZ);
  muH2 = MH*(MH-I*GH);
  cw2=muW2/muZ2;
  sw2=1.-cw2;
  m_sW = sqrt(std::abs(m_sW2));
  m_cW = sqrt(std::abs(m_cW2));
}

PHOTONS_ME_Base::~PHOTONS_ME_Base() {
  if (p_boost) delete p_boost;
  if (p_rot) delete p_rot;
}

PHOTONS_ME_Base * PHOTONS_ME_Base::GetIRsubtractedME
(const Particle_Vector_Vector& pvv)
{
  PHOTONS_ME_Getter::Getter_List glist(PHOTONS_ME_Getter::GetGetters());
  for (PHOTONS_ME_Getter::Getter_List::const_iterator git(glist.begin());
       git!=glist.end();++git) {
    PHOTONS_ME_Base * pme = (*git)->GetObject(pvv);
    if (pme && pme->Name() != "Collinear_Approximation_FF" && pme->Name() != "Collinear_Approximation_FI") return pme;
  }
  return NULL;
}

PHOTONS_ME_Base * PHOTONS_ME_Base::GetIRsubtractedME
(const std::string& tag, const Particle_Vector_Vector& pvv)
{
  PHOTONS_ME_Base * pme = PHOTONS_ME_Getter::GetObject(tag, pvv);
  if (!pme && (tag == "Collinear_Approximation_FF" || tag == "Collinear_Approximation_FI")) return NULL;
  if (!pme) THROW(fatal_error, "Did not find IR subtracted ME "+tag);
  return pme;
}

void PHOTONS_ME_Base::Print_EW_Parameters() {
  msg_Out() << "\nPHOTONS_ME_Base EW Parameters\n" 
	    << "Masses\n"
	    << "MZ    " << sqrt(MZ2) << "\n"
	    << "MZ2   " << MZ2 << "\n"
	    << "MW    " << sqrt(MW2) << "\n"
	    << "MW2   " << MW2 << "\n"
	    << "MH    " << sqrt(MH2) << "\n"
	    << "MH2   " << MH2 << "\n"
	    << "mt    " << sqrt(mt2) << "\n"
	    << "mt2   " << mt2 << "\n"
	    << "mtau  " << Flavour(kf_tau).Mass() << "\n"
	    << "mtau2 " << sqr(Flavour(kf_tau).Mass()) << "\n\n"
	    << "Weinberg Angle\n"
	    << "sw    " << m_sW << "\n"
	    << "sw2   " << m_sW2 << "\n"
	    << "cw    " << m_cW << "\n"
	    << "cw2   " << m_cW2 << "\n\n"
	    << "Couplings\n"
	    << "e     " << m_e << "\n"
	    << "gw    " << m_e/m_sW << "\n"
	    << "alpha " << m_alpha << "\n"
	    << "GF    " << /*m_GF*/m_sqrt2*sqr(m_e/m_sW)/(8.*MW2) << "\n"
	    << "vev   " << 2.*sqrt(MW2)*m_sW/m_e << "\n\n"; 
}
