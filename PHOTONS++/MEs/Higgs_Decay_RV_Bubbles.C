#include "PHOTONS++/MEs/Higgs_Decay_RV_Diagrams.H"
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

Higgs_Decay_RV_Bubbles::Higgs_Decay_RV_Bubbles
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Higgs_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Bubble1Abb = new Complex[28];
  Bubble2Abb = new Complex[28];
  Init_Coefficients();
}

Higgs_Decay_RV_Bubbles::~Higgs_Decay_RV_Bubbles()
{
  delete [] Bubble1Abb;
  delete [] Bubble2Abb;
}

void Higgs_Decay_RV_Bubbles::Init_Coefficients() 
{
  Init_Bubble_1_Coefficients();
  Init_Bubble_2_Coefficients();
  return;
}

void Higgs_Decay_RV_Bubbles::Init_Bubble_1_Coefficients() 
{
  Bubble1Abb[0]=m_m2 - m_s1k;

  Bubble1Abb[1]=m_m2 + m_s1k;

  Bubble1Abb[2]=m_m4 - 6.*m_m2*m_s1k + m_s1k_2;

  Bubble1Abb[6]=-3.*m_m2;

  Bubble1Abb[7]=m_s1k;

  Bubble1Abb[8]=m_m2/m_s1k;

  Bubble1Abb[9]=-1./m_s1k;

  Bubble1Abb[13]=m_m2 - 3.*m_s1k;

  Bubble1Abb[16]=m_m;

  Bubble1Abb[17]=(-m_m3)/m_s1k;

  Bubble1Abb[18]=m_m/m_s1k;

  Bubble1Abb[21]=-3.*m_m2 + m_s1k;

  Bubble1Abb[26]=(2.*m_m2)/m_s1k;

  Bubble1Abb[27]=-2./m_s1k;

  return;
}


void Higgs_Decay_RV_Bubbles::Init_Bubble_2_Coefficients() 
{
  Bubble2Abb[0]=m_m2 - m_s2k;

  Bubble2Abb[1]=m_m2 + m_s2k;

  Bubble2Abb[2]=m_m4 - 6.*m_m2*m_s2k + m_s2k_2;

  Bubble2Abb[6]=-3.*m_m2;

  Bubble2Abb[7]=m_s2k;

  Bubble2Abb[8]=m_m2/m_s2k;

  Bubble2Abb[9]=-1./m_s2k;

  Bubble2Abb[13]=m_m2 - 3.*m_s2k;

  Bubble2Abb[16]=m_m;

  Bubble2Abb[17]=(-m_m3)/m_s2k;

  Bubble2Abb[18]=m_m/m_s2k;

  Bubble2Abb[21]=3.*m_m2 - m_s2k;

  Bubble2Abb[26]=(-2.*m_m2)/m_s2k;

  Bubble2Abb[27]=2./m_s2k;

  return;
}



DivArrC Higgs_Decay_RV_Bubbles::RV_Bubble_Insertion_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Box corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(2.*Bubble1Abb[21]*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[1]*Bubble1Abb[26]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[27]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2);
    }
    else if (LR == 1) {
      return One*(2.*Bubble1Abb[21]*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[1]*Bubble1Abb[26]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[27]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2);
    }
  }
  // ubar1 \slashed{k} P_i v2
  else if (ME == 2) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsP} P_i v2
  else if (ME == 3) {
    if (LR == 0) {
      return One*Bubble1Abb[16]/Bubble1Abb[0]

	+(Bubble1Abb[17]*B_0(0.,0.,m_m2,m_mu2))/Bubble1Abb[0]

	+(Bubble1Abb[13]*Bubble1Abb[18]*B_0(m_s1k,0.,m_m2,m_mu2))/Bubble1Abb[0];
    }
    else if (LR == 1) {
      return One*Bubble1Abb[16]/Bubble1Abb[0]

	+(Bubble1Abb[17]*B_0(0.,0.,m_m2,m_mu2))/Bubble1Abb[0]

	+(Bubble1Abb[13]*Bubble1Abb[18]*B_0(m_s1k,0.,m_m2,m_mu2))/Bubble1Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(Bubble1Abb[6]/pow(Bubble1Abb[0],2)+Bubble1Abb[7]/pow(Bubble1Abb[0],2))

	+(Bubble1Abb[1]*Bubble1Abb[8]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[9]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble1Abb[6]/pow(Bubble1Abb[0],2)+Bubble1Abb[7]/pow(Bubble1Abb[0],2))

	+(Bubble1Abb[1]*Bubble1Abb[8]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[9]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}


DivArrC Higgs_Decay_RV_Bubbles::RV_Bubble_Insertion_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Box corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // PRINT_VAR("\n\nBox 1 Finite " << m_s << " , " << m_s12 <<  " , " << m_s2k << " , " << m_m2 << " , " << m_mu2 
  // 	    << "\n" << B_0(0.,0.,m_m2,m_mu2).Finite()
  // 	    << "\n" << B_0(m_s2k,0.,m_m2,m_mu2).Finite());
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(2.*Bubble2Abb[21]*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[1]*Bubble2Abb[26]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[2]*Bubble2Abb[27]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2);
    }
    else if (LR == 1) {
      return One*(2.*Bubble2Abb[21]*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[1]*Bubble2Abb[26]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[2]*Bubble2Abb[27]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2);
    }
  }
  // ubar1 \slashed{k} P_i v2
  else if (ME == 2) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsP} P_i v2
  else if (ME == 3) {
    if (LR == 0) {
      return One*Bubble2Abb[16]/Bubble2Abb[0]

	+(Bubble2Abb[17]*B_0(0.,0.,m_m2,m_mu2))/Bubble2Abb[0]

	+(Bubble2Abb[13]*Bubble2Abb[18]*B_0(m_s2k,0.,m_m2,m_mu2))/Bubble2Abb[0];
    }
    else if (LR == 1) {
      return One*Bubble2Abb[16]/Bubble2Abb[0]

	+(Bubble2Abb[17]*B_0(0.,0.,m_m2,m_mu2))/Bubble2Abb[0]

	+(Bubble2Abb[13]*Bubble2Abb[18]*B_0(m_s2k,0.,m_m2,m_mu2))/Bubble2Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(Bubble2Abb[6]/pow(Bubble2Abb[0],2)+Bubble2Abb[7]/pow(Bubble2Abb[0],2))

	+(Bubble2Abb[1]*Bubble2Abb[8]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[2]*Bubble2Abb[9]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble2Abb[6]/pow(Bubble2Abb[0],2)+Bubble2Abb[7]/pow(Bubble2Abb[0],2))

	+(Bubble2Abb[1]*Bubble2Abb[8]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[2]*Bubble2Abb[9]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}


