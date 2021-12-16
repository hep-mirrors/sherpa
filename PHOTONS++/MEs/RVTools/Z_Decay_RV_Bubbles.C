#include "PHOTONS++/MEs/RVTools/Z_Decay_RV_Diagrams.H"
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

// Calculation of diagrams involving internal bubbles on the fermion legs.
Z_Decay_RV_Bubbles::Z_Decay_RV_Bubbles
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Z_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Bubble1Abb = new Complex[50];
  Bubble2Abb = new Complex[50];
  Init_Coefficients();
}

Z_Decay_RV_Bubbles::~Z_Decay_RV_Bubbles()
{
  delete [] Bubble1Abb;
  delete [] Bubble2Abb;
}

void Z_Decay_RV_Bubbles::Init_Coefficients() 
{
  Init_Bubble_1_Coefficients();
  Init_Bubble_2_Coefficients();
  return;
}

// Set up coefficients multiplying master integrals for emission from leg 1
void Z_Decay_RV_Bubbles::Init_Bubble_1_Coefficients() 
{
  Bubble1Abb[0]=m_x - m_z1k;

  Bubble1Abb[1]=-3.*m_x + m_z1k;

  Bubble1Abb[2]=m_x + m_z1k;

  Bubble1Abb[3]=-6.*m_x + m_z1k;

  Bubble1Abb[4]=m_m4 + Bubble1Abb[3]*m_s_2*m_z1k;

  Bubble1Abb[9]=m_cL/m_s;

  Bubble1Abb[10]=m_s*m_x;

  Bubble1Abb[11]=(m_cL*m_x)/(m_s*m_z1k);

  Bubble1Abb[12]=m_s*m_z1k;

  Bubble1Abb[13]=(-m_cL)/(m_s_3*m_z1k);

  Bubble1Abb[18]=m_cR/m_s;

  Bubble1Abb[19]=(m_cR*m_x)/(m_s*m_z1k);

  Bubble1Abb[20]=(-m_cR)/(m_s_3*m_z1k);

  Bubble1Abb[21]=m_s*m_x - m_s*m_z1k;

  Bubble1Abb[22]=-m_x + m_z1k;

  Bubble1Abb[23]=m_x - 3.*m_z1k;

  Bubble1Abb[27]=m_cL*m_m;

  Bubble1Abb[28]=(m_cL*m_m3)/(m_s_2*m_z1k);

  Bubble1Abb[29]=(m_cL*m_m)/m_z1k;

  Bubble1Abb[33]=m_cR*m_m;

  Bubble1Abb[34]=(m_cR*m_m3)/(m_s_2*m_z1k);

  Bubble1Abb[35]=(m_cR*m_m)/m_z1k;

  Bubble1Abb[40]=(2.*m_cL)/m_s;

  Bubble1Abb[41]=(2.*m_cL*m_x)/(m_s*m_z1k);

  Bubble1Abb[42]=(-2.*m_cL)/(m_s_3*m_z1k);

  Bubble1Abb[47]=(2.*m_cR)/m_s;

  Bubble1Abb[48]=(2.*m_cR*m_x)/(m_s*m_z1k);

  Bubble1Abb[49]=(-2.*m_cR)/(m_s_3*m_z1k);
}


// Set up coefficients multiplying master integrals for emission from leg 2
void Z_Decay_RV_Bubbles::Init_Bubble_2_Coefficients() 
{
  Bubble2Abb[0]=m_x - m_z2k;

  Bubble2Abb[1]=3.*m_x - m_z2k;

  Bubble2Abb[2]=m_x + m_z2k;

  Bubble2Abb[3]=-6.*m_x + m_z2k;

  Bubble2Abb[4]=m_m4 + Bubble2Abb[3]*m_s_2*m_z2k;

  Bubble2Abb[9]=m_cL/m_s;

  Bubble2Abb[10]=m_s*m_x;

  Bubble2Abb[11]=(-m_cL*m_x)/(m_s*m_z2k);

  Bubble2Abb[12]=m_s*m_z2k;

  Bubble2Abb[13]=m_cL/(m_s_3*m_z2k);

  Bubble2Abb[18]=m_cR/m_s;

  Bubble2Abb[19]=(-m_cR*m_x)/(m_s*m_z2k);

  Bubble2Abb[20]=m_cR/(m_s_3*m_z2k);

  Bubble2Abb[21]=m_s*m_x - m_s*m_z2k;

  Bubble2Abb[22]=-m_x + m_z2k;

  Bubble2Abb[23]=m_x - 3.*m_z2k;

  Bubble2Abb[27]=m_cR*m_m;

  Bubble2Abb[28]=(m_cR*m_m3)/(m_s_2*m_z2k);

  Bubble2Abb[29]=(m_cR*m_m)/m_z2k;

  Bubble2Abb[33]=m_cL*m_m;

  Bubble2Abb[34]=(m_cL*m_m3)/(m_s_2*m_z2k);

  Bubble2Abb[35]=(m_cL*m_m)/m_z2k;

  Bubble2Abb[40]=(2.*m_cL)/m_s;

  Bubble2Abb[41]=(-2.*m_cL*m_x)/(m_s*m_z2k);

  Bubble2Abb[42]=(2.*m_cL)/(m_s_3*m_z2k);

  Bubble2Abb[47]=(2.*m_cR)/m_s;

  Bubble2Abb[48]=(-2.*m_cR*m_x)/(m_s*m_z2k);

  Bubble2Abb[49]=(2.*m_cR)/(m_s_3*m_z2k);
}



DivArrC Z_Decay_RV_Bubbles::RV_Bubble_Insertion_1(const int& ME, const int& LR,
						  const Vec4C& epsV, const Vec4C& epsP)
{
  // Bubble corrections for emission off leg 1
  // ME takes values 1..8, denoting standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
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
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsV} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(Bubble1Abb[1]*Bubble1Abb[40]*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[41]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[4]*Bubble1Abb[42]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble1Abb[1]*Bubble1Abb[47]*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[2]*Bubble1Abb[48]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[4]*Bubble1Abb[49]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(Bubble1Abb[0],2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 5) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 6) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} P_i v2
  else if (ME == 7) {
    if (LR == 0) {
      return One*Bubble1Abb[27]/Bubble1Abb[21]

	+(Bubble1Abb[28]*B_0(0.,0.,m_m2,m_mu2))/Bubble1Abb[22]

	+(Bubble1Abb[23]*Bubble1Abb[29]*B_0(m_s1k,0.,m_m2,m_mu2))/Bubble1Abb[21];
    }
    else if (LR == 1) {
      return One*Bubble1Abb[33]/Bubble1Abb[21]

	+(Bubble1Abb[34]*B_0(0.,0.,m_m2,m_mu2))/Bubble1Abb[22]

	+(Bubble1Abb[23]*Bubble1Abb[35]*B_0(m_s1k,0.,m_m2,m_mu2))/Bubble1Abb[21];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*(Bubble1Abb[1]*Bubble1Abb[9])/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[11]*Bubble1Abb[2]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[13]*Bubble1Abb[4]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble1Abb[1]*Bubble1Abb[18])/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[19]*Bubble1Abb[2]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2)

	+(Bubble1Abb[20]*Bubble1Abb[4]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(Bubble1Abb[0],2);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}


DivArrC Z_Decay_RV_Bubbles::RV_Bubble_Insertion_2(const int& ME, const int& LR,
							  const Vec4C& epsV, const Vec4C& epsP)
{
  // Bubble corrections for emission off leg 2
  // ME takes values 1..8, denoting standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
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
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsV} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(Bubble2Abb[1]*Bubble2Abb[40]*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[2]*Bubble2Abb[41]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[4]*Bubble2Abb[42]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble2Abb[1]*Bubble2Abb[47]*(epsP*m_p2))/pow(Bubble2Abb[0],2)
	
	+(Bubble2Abb[2]*Bubble2Abb[48]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[4]*Bubble2Abb[49]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(Bubble2Abb[0],2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 5) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 6) {
    if (LR == 0) {
      return Zero;
    }
    else if (LR == 1) {
      return Zero;
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} P_i v2
  else if (ME == 7) {
    if (LR == 0) {
      return One*Bubble2Abb[27]/Bubble2Abb[21]

	+(Bubble2Abb[28]*B_0(0.,0.,m_m2,m_mu2))/Bubble2Abb[22]

	+(Bubble2Abb[23]*Bubble2Abb[29]*B_0(m_s2k,0.,m_m2,m_mu2))/Bubble2Abb[21];
    }
    else if (LR == 1) {
      return One*Bubble2Abb[33]/Bubble2Abb[21]

	+(Bubble2Abb[34]*B_0(0.,0.,m_m2,m_mu2))/Bubble2Abb[22]

	+(Bubble2Abb[23]*Bubble2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2))/Bubble2Abb[21];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*(Bubble2Abb[1]*Bubble2Abb[9])/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[11]*Bubble2Abb[2]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[13]*Bubble2Abb[4]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2);
    }
    else if (LR == 1) {
      return One*(Bubble2Abb[1]*Bubble2Abb[18])/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[19]*Bubble2Abb[2]*B_0(0.,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2)

	+(Bubble2Abb[20]*Bubble2Abb[4]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(Bubble2Abb[0],2);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}


