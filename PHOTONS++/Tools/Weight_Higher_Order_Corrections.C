#include "PHOTONS++/Tools/Weight_Higher_Order_Corrections.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Main/Photons.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

Weight_Higher_Order_Corrections::Weight_Higher_Order_Corrections
(const Particle_Vector_Vector& pvv_old, const Particle_Vector_Vector& pvv_new,
 Dipole_Type::code dtype) : m_n(pvv_new[4].size()), p_pme(NULL), p_app(NULL) {
  DEBUG_FUNC(PHOTONS::Photons::s_useme);
  if (PHOTONS::Photons::s_useme)
    p_pme = PHOTONS_ME_Base::GetIRsubtractedME(pvv_old);
  if (p_pme) {
    msg_Debugging()<<"ME -> "<<p_pme->Name()<<std::endl;
    p_pme->FillMomentumArrays(pvv_new);
    // Reject calculation of higher order corrections entirely if momenta are not corrected properly
    // for some case in momentum reconstruction.
    // In that case return only the soft approximation.
    if (p_pme->RejectCalc()) {
      m_weight = 1.;
      m_maxweight = 1.;
    }
    else {
      CalculateWeightAndMaxWithME();
    }
  }
  else {
    msg_Debugging()<<"ME -> none"<<std::endl;
    m_dtype       = dtype;
    m_newdipole   = pvv_new[2];
    m_olddipole   = pvv_old[2];
    m_softphotons = pvv_new[4];
    // Get collinear approximation
    if (m_dtype == Dipole_Type::ff) {
      m_M = pvv_old[1][0]->FinalMass();
      p_app = PHOTONS_ME_Base::GetIRsubtractedME("Collinear_Approximation_FF",pvv_old);
    }
    if (m_dtype == Dipole_Type::fi) {
      m_M = pvv_old[0][0]->FinalMass();
      p_app = PHOTONS_ME_Base::GetIRsubtractedME("Collinear_Approximation_FI",pvv_old);
    }
    if (p_app) {
      msg_Debugging()<<"ME -> "<<p_app->Name()<<std::endl;
      p_app->FillMomentumArrays(pvv_new);
      // Reject calculation of higher order corrections entirely if momenta are not corrected properly
      // for some case in momentum reconstruction.
      // In that case return only the soft approximation.
      if (p_pme->RejectCalc()) {
	m_weight = 1.;
	m_maxweight = 1.;
      }
      else {
	CalculateWeightAndMaxWithME();
      }
    }
    else {
      msg_Debugging() << METHOD << ": Could not find collinear approximation in this case. Weird!\n";
    }
  }
}

Weight_Higher_Order_Corrections::~Weight_Higher_Order_Corrections() {
  if (p_pme) delete p_pme;
  if (p_app) delete p_app;
}

void Weight_Higher_Order_Corrections::CalculateWeight() {
  m_weight = 1.;
}

void Weight_Higher_Order_Corrections::CalculateMax() {
  m_maxweight = 1.;
}

void Weight_Higher_Order_Corrections::CalculateWeightAndMaxWithME() {
  DEBUG_FUNC("");
  double B_0_0 = p_pme->GetBeta_0_0();
  if (IsNan(B_0_0)) {
    msg_Out() << METHOD << ": Tree level ME is nan." <<
      "ME responsible for this is " << p_pme->Name() << ". Set ME weight to 1." << std::endl;
    m_weight = 1.;
    m_maxweight = 1.;
    return;
  }
  double B_0_1 = p_pme->GetBeta_0_1();
  double B_0_2 = p_pme->GetBeta_0_2();
  double Sum_1_1(0.), Sum_1_2(0.);
  for (unsigned int i=0; i<m_n; i++) {
    /// same safeguard as below
    Sum_1_1 = Sum_1_1 + p_pme->GetBeta_1_1(i)/p_pme->SmodFull(i);
    Sum_1_2 = Sum_1_2 + p_pme->GetBeta_1_2(i)/p_pme->SmodFull(i);
  }
  double Sum_2_2(0.);
  // Choose unconstrained pairs as MEs implemented this way
  for (unsigned int j=0; j<m_n; j++) {
    for (unsigned int i=0; i<j; i++) {
      double beta = p_pme->GetBeta_2_2(i,j);
      if (beta != 0.) {
	Sum_2_2 = Sum_2_2 + beta
	  /(p_pme->SmodFull(i)*p_pme->SmodFull(j));
      }
    }
  }
  if (msg_LevelIsDebugging()) {
    std::string name(p_pme->Name());
    std::string space(name.size(),' ');
    msg_Out()<<name <<"  \\beta_0^0 = "<<B_0_0<<std::endl
             <<space<<"  \\beta_0^1 = "<<B_0_1<<std::endl
             <<space<<"  \\beta_0^2 = "<<B_0_2<<std::endl
             <<space<<"  \\sum_i[\\beta_1^1(i)/S(i)] = "<<Sum_1_1<<std::endl
             <<space<<"  \\sum_i[\\beta_1^2(i)/S(i)] = "<<Sum_1_2<<std::endl
             <<space<<"  \\sum_i\\sum_j[\\beta_2^2(i,j)/S(i)S(j)] = "<<Sum_2_2
             <<std::endl;
  }
  m_weight    = 1. + 1./B_0_0*( B_0_1 + B_0_2 + Sum_1_1 + Sum_1_2 + Sum_2_2 );
  m_maxweight = 1. + 1./B_0_0*( B_0_1 + B_0_2 );
}

