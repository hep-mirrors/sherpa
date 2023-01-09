#include "PHOTONS++/MEs/Z_Decay_RV_Diagrams.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)
#define D_0(A,B,C,D,E,F,G,H,I,J,M) Master_Box(A,B,C,D,E,F,G,H,I,J,M)


using namespace ATOOLS;
using namespace METOOLS;
using namespace PHOTONS;

Z_Decay_RV_Diagrams::Z_Decay_RV_Diagrams
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  m_m(m), m_s(s), m_mu2(mu2), m_p1(p1), m_p2(p2), m_pP(pP),
  m_cL(cL), m_cR(cR)
{
  m_p12 = m_p1*m_p1;
  m_p22 = m_p2*m_p2;
  m_pP2 = m_pP*m_pP;
  m_s12 = (m_p1+m_p2)*(m_p1+m_p2);
  m_s1k = (m_p1+m_pP)*(m_p1+m_pP);
  m_s2k = (m_p2+m_pP)*(m_p2+m_pP);
  Init_Powers();
  Zero = DivArrC(0.,0.,0.,0.,0.,0.);
  One = DivArrC(0.,0.,0.,1.,0.,0.);
}

Z_Decay_RV_Diagrams::~Z_Decay_RV_Diagrams()
{
}

// Initialize abbreviations
void Z_Decay_RV_Diagrams::Init_Powers() 
{
  m_m2 = m_m*m_m;
  m_m3 = m_m2*m_m;
  m_m4 = m_m3*m_m;
  m_x = m_m2/m_s;
  m_z12 = m_s12/m_s;
  m_z1k = m_s1k/m_s;
  m_z2k = m_s2k/m_s;
  m_s_2 = m_s*m_s;
  m_s_3 = m_s_2*m_s;
  m_x_2 = m_x*m_x;
  m_x_3 = m_x_2*m_x;
  m_x_4 = m_x_3*m_x;
  m_x_5 = m_x_4*m_x;
  m_x_6 = m_x_5*m_x;
  m_x_7 = m_x_6*m_x;
  m_x_8 = m_x_7*m_x;
  m_x_9 = m_x_8*m_x;
  m_x_10 = m_x_9*m_x;
  m_x_11 = m_x_10*m_x;
  m_z12_2 = m_z12*m_z12;
  m_z12_3 = m_z12_2*m_z12;
  m_z12_4 = m_z12_3*m_z12;
  m_z12_5 = m_z12_4*m_z12;
  m_z12_6 = m_z12_5*m_z12;
  m_z12_7 = m_z12_6*m_z12;
  m_z12_8 = m_z12_7*m_z12;
  m_z12_9 = m_z12_8*m_z12;
  m_z1k_2 = m_z1k*m_z1k;
  m_z1k_3 = m_z1k_2*m_z1k;
  m_z1k_4 = m_z1k_3*m_z1k;
  m_z1k_5 = m_z1k_4*m_z1k;
  m_z1k_6 = m_z1k_5*m_z1k;
  m_z1k_7 = m_z1k_6*m_z1k;
  m_z1k_8 = m_z1k_7*m_z1k;
  m_z1k_9 = m_z1k_8*m_z1k;
  m_z2k_2 = m_z2k*m_z2k;
  m_z2k_3 = m_z2k_2*m_z2k;
  m_z2k_4 = m_z2k_3*m_z2k;
  m_z2k_5 = m_z2k_4*m_z2k;
  m_z2k_6 = m_z2k_5*m_z2k;
  m_z2k_7 = m_z2k_6*m_z2k;
  m_z2k_8 = m_z2k_7*m_z2k;
  m_z2k_9 = m_z2k_8*m_z2k;
  return;
}
