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

Higgs_Decay_RV_Vertices::Higgs_Decay_RV_Vertices
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Higgs_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  HVertex1Abb = new Complex[43];
  HVertex2Abb = new Complex[43];
  PVertex1Abb = new Complex[37];
  PVertex2Abb = new Complex[42];
  Init_Coefficients();
}

Higgs_Decay_RV_Vertices::~Higgs_Decay_RV_Vertices()
{
  delete [] HVertex1Abb;
  delete [] PVertex1Abb;
  delete [] HVertex2Abb;
  delete [] PVertex2Abb;
}

void Higgs_Decay_RV_Vertices::Init_Coefficients() 
{
  Init_H_Vertex_1_Coefficients();
  Init_H_Vertex_2_Coefficients();
  Init_P_Vertex_1_Coefficients();
  Init_P_Vertex_2_Coefficients();
  return;
}

// Initialize corrections to Higgs vertex for emission from leg 1
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Vertices::Init_H_Vertex_1_Coefficients() 
{
  HVertex1Abb[0]=m_m2 - m_s1k;

  HVertex1Abb[1]=m_s - m_s1k;

  HVertex1Abb[2]=pow(HVertex1Abb[1],2.) + 3.*m_m4 - 4.*m_m2*m_s1k;

  HVertex1Abb[3]=m_s + m_s1k;

  HVertex1Abb[4]=pow(HVertex1Abb[1],2.) - 2.*HVertex1Abb[3]*m_m2 + m_m4;

  HVertex1Abb[5]=m_m2 + m_s - m_s1k;

  HVertex1Abb[6]=m_m2 - m_s + m_s1k;

  HVertex1Abb[7]=3.*m_s + 8.*m_s1k;

  HVertex1Abb[8]=m_s_2 - m_s*m_s1k + m_s1k_2;

  HVertex1Abb[9]=-4.*HVertex1Abb[8]*m_m2 + HVertex1Abb[7]*m_m4 - 4.*m_m6 + pow(HVertex1Abb[1],2.)*m_s;

  HVertex1Abb[16]=8.*m_m2*m_s;

  HVertex1Abb[23]=3.*m_m2 - m_s + m_s1k;

  HVertex1Abb[28]=8.*m_m3;

  HVertex1Abb[29]=-4.*m_m;

  HVertex1Abb[30]=2.*m_m;

  HVertex1Abb[35]=4.*HVertex1Abb[8]*m_m2 - HVertex1Abb[7]*m_m4 + 4.*m_m6 - pow(HVertex1Abb[1],2.)*m_s;

  HVertex1Abb[42]=16.*m_m2*m_s;
  
  return;
}

// Initialize corrections to photon vertex for emission from leg 1
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Vertices::Init_P_Vertex_1_Coefficients() 
{
  PVertex1Abb[0]=m_m2 - m_s1k;

  PVertex1Abb[1]=2.*m_m2 + m_s1k;

  PVertex1Abb[2]=5.*m_m2 + m_s1k;

  PVertex1Abb[6]=2.*m_m2;

  PVertex1Abb[11]=-2.*m_m;

  PVertex1Abb[12]=2.*m_m;

  PVertex1Abb[14]=-3.*m_m2 + m_s1k;

  PVertex1Abb[15]=-m_m2 + m_s1k;

  PVertex1Abb[16]=m_m2 + m_s1k;

  PVertex1Abb[17]=3.*m_m2 - m_s1k;

  PVertex1Abb[23]=(2.*m_m3)/m_s1k;

  PVertex1Abb[24]=(2.*m_m)/m_s1k;

  PVertex1Abb[30]=m_m4 - 6.*m_m2*m_s1k + m_s1k_2;

  PVertex1Abb[34]=(-2.*m_m2)/m_s1k;

  PVertex1Abb[35]=12.*m_m2;

  PVertex1Abb[36]=2./m_s1k;

  return;
}

// Initialize corrections to Higgs vertex for emission from leg 2
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Vertices::Init_H_Vertex_2_Coefficients() 
{
  HVertex2Abb[0]=m_m2 - m_s2k;

  HVertex2Abb[1]=m_s - m_s2k;

  HVertex2Abb[2]=pow(HVertex2Abb[1],2.) + 3.*m_m4 - 4.*m_m2*m_s2k;

  HVertex2Abb[3]=m_s + m_s2k;

  HVertex2Abb[4]=pow(HVertex2Abb[1],2.) - 2.*HVertex2Abb[3]*m_m2 + m_m4;

  HVertex2Abb[5]=m_m2 + m_s - m_s2k;

  HVertex2Abb[6]=m_m2 - m_s + m_s2k;

  HVertex2Abb[7]=3.*m_s + 8.*m_s2k;

  HVertex2Abb[8]=m_s_2 - m_s*m_s2k + m_s2k_2;

  HVertex2Abb[9]=-4.*HVertex2Abb[8]*m_m2 + HVertex2Abb[7]*m_m4 - 4.*m_m6 + pow(HVertex2Abb[1],2.)*m_s;

  HVertex2Abb[16]=8.*m_m2*m_s;

  HVertex2Abb[23]=3.*m_m2 - m_s + m_s2k;

  HVertex2Abb[28]=8.*m_m3;

  HVertex2Abb[29]=-4.*m_m;

  HVertex2Abb[30]=2.*m_m;

  HVertex2Abb[35]=pow(HVertex2Abb[1],2.) - m_m4;

  HVertex2Abb[36]=4.*HVertex2Abb[8]*m_m2 - HVertex2Abb[7]*m_m4 + 4.*m_m6 - pow(HVertex2Abb[1],2.)*m_s;

  HVertex2Abb[42]=-16.*m_m2*m_s;

  return;
}

// Initialize corrections to photon vertex for emission from leg 2
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Vertices::Init_P_Vertex_2_Coefficients() 
{
  PVertex2Abb[0]=m_m2 - m_s2k;

  PVertex2Abb[1]=2.*m_m2 + m_s2k;

  PVertex2Abb[2]=5.*m_m2 + m_s2k;

  PVertex2Abb[6]=2.*m_m2;

  PVertex2Abb[11]=-2.*m_m;

  PVertex2Abb[12]=2.*m_m;

  PVertex2Abb[14]=-3.*m_m2 + m_s2k;

  PVertex2Abb[15]=-m_m2 + m_s2k;

  PVertex2Abb[16]=m_m2 + m_s2k;

  PVertex2Abb[17]=3.*m_m2 - m_s2k;

  PVertex2Abb[23]=(2.*m_m3)/m_s2k;

  PVertex2Abb[24]=(2.*m_m)/m_s2k;

  PVertex2Abb[30]=m_m4 - 6.*m_m2*m_s2k + m_s2k_2;

  PVertex2Abb[33]=(2.*m_m2)/m_s2k;

  PVertex2Abb[34]=-8.*m_m2;

  PVertex2Abb[35]=-2./m_s2k;

  return;
}


DivArrC Higgs_Decay_RV_Vertices::RV_H_Vertex_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // H Vertex corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,

  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(4.*(epsP*m_p1))/HVertex1Abb[0]

	+(-4.*HVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(HVertex1Abb[42]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(4.*HVertex1Abb[5]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(-4.*HVertex1Abb[35]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4]);
    }
    else if (LR == 1) {
      return One*(4.*(epsP*m_p1))/HVertex1Abb[0]

	+(-4.*HVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(HVertex1Abb[42]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(4.*HVertex1Abb[5]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(-4.*HVertex1Abb[35]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(HVertex1Abb[0]*HVertex1Abb[4]);
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
      return (HVertex1Abb[28]*B_0(m_m2,0.,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[29]*HVertex1Abb[5]*B_0(m_s,m_m2,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[29]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[23]*HVertex1Abb[30]*HVertex1Abb[5]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2))/HVertex1Abb[4];
    }
    else if (LR == 1) {
      return (HVertex1Abb[28]*B_0(m_m2,0.,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[29]*HVertex1Abb[5]*B_0(m_s,m_m2,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[29]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2))/HVertex1Abb[4]

	+(HVertex1Abb[23]*HVertex1Abb[30]*HVertex1Abb[5]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2))/HVertex1Abb[4];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*2./HVertex1Abb[0]

	+(-2.*HVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(HVertex1Abb[16]*B_0(m_s,m_m2,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(2.*HVertex1Abb[5]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(2.*HVertex1Abb[9]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4]);
    }
    else if (LR == 1) {
      return One*2./HVertex1Abb[0]

	+(-2.*HVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(HVertex1Abb[16]*B_0(m_s,m_m2,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(2.*HVertex1Abb[5]*HVertex1Abb[6]*B_0(m_s1k,0.,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4])

	+(2.*HVertex1Abb[9]*C_0(m_m2,m_s,m_s1k,0.,m_m2,m_m2,m_mu2))/(HVertex1Abb[0]*HVertex1Abb[4]);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}

DivArrC Higgs_Decay_RV_Vertices::RV_P_Vertex_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Photon Vertex corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(-2.*PVertex1Abb[16]*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[16]*PVertex1Abb[34]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[35]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[30]*PVertex1Abb[36]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.);
    }
    else if (LR == 1) {
      return One*(-2.*PVertex1Abb[16]*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[16]*PVertex1Abb[34]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[35]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[30]*PVertex1Abb[36]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],2.);
    }
  }
  // ubar1 \slashed{k} P_i v2
  else if (ME == 2) {
    if (LR == 0) {
      return One*(PVertex1Abb[12]*PVertex1Abb[14]*(epsP*m_p1))/pow(PVertex1Abb[0],3.)

	+(PVertex1Abb[16]*PVertex1Abb[23]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[15],3.)

	+(PVertex1Abb[12]*PVertex1Abb[17]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],3.)

	+(PVertex1Abb[24]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/PVertex1Abb[0];
    }
    else if (LR == 1) {
      return One*(PVertex1Abb[12]*PVertex1Abb[14]*(epsP*m_p1))/pow(PVertex1Abb[0],3.)

	+(PVertex1Abb[16]*PVertex1Abb[23]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[15],3.)

	+(PVertex1Abb[12]*PVertex1Abb[17]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[0],3.)

	+(PVertex1Abb[24]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/PVertex1Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} P_i v2
  else if (ME == 3) {
    if (LR == 0) {
      return (PVertex1Abb[11]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex1Abb[0]

	+(PVertex1Abb[12]*B_0(m_s1k,0.,m_m2,m_mu2))/PVertex1Abb[0];
    }
    else if (LR == 1) {
      return (PVertex1Abb[11]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex1Abb[0]

	+(PVertex1Abb[12]*B_0(m_s1k,0.,m_m2,m_mu2))/PVertex1Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*2./PVertex1Abb[0]

	+(2.*PVertex1Abb[1]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex1Abb[0],2.)

	+(-PVertex1Abb[2]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[6]*C_0(m_m2,0.,m_s1k,0.,m_m2,m_m2,m_mu2))/PVertex1Abb[0];
    }
    else if (LR == 1) {
      return One*2./PVertex1Abb[0]

	+(2.*PVertex1Abb[1]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex1Abb[0],2.)

	+(-PVertex1Abb[2]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(PVertex1Abb[0],2.)

	+(PVertex1Abb[6]*C_0(m_m2,0.,m_s1k,0.,m_m2,m_m2,m_mu2))/PVertex1Abb[0];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}

DivArrC Higgs_Decay_RV_Vertices::RV_H_Vertex_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // H Vertex corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(-4.*(epsP*m_p2))/HVertex2Abb[0]

	+(4.*HVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(HVertex2Abb[42]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(4.*HVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(4.*HVertex2Abb[36]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4]);
    }
    else if (LR == 1) {
      return One*(-4.*(epsP*m_p2))/HVertex2Abb[0]

	+(4.*HVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(HVertex2Abb[42]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(4.*HVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(4.*HVertex2Abb[36]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(HVertex2Abb[0]*HVertex2Abb[4]);
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
      return (HVertex2Abb[28]*B_0(m_m2,0.,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[29]*HVertex2Abb[5]*B_0(m_s,m_m2,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[29]*HVertex2Abb[6]*B_0(m_s2k,0.,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[23]*HVertex2Abb[30]*HVertex2Abb[5]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2))/HVertex2Abb[4];
    }
    else if (LR == 1) {
      return (HVertex2Abb[28]*B_0(m_m2,0.,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[29]*HVertex2Abb[5]*B_0(m_s,m_m2,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[29]*HVertex2Abb[6]*B_0(m_s2k,0.,m_m2,m_mu2))/HVertex2Abb[4]

	+(HVertex2Abb[23]*HVertex2Abb[30]*HVertex2Abb[5]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2))/HVertex2Abb[4];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*2./HVertex2Abb[0]

	+(-2.*HVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(HVertex2Abb[16]*B_0(m_s,m_m2,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(2.*HVertex2Abb[5]*HVertex2Abb[6]*B_0(m_s2k,0.,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(2.*HVertex2Abb[9]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4]);
    }
    else if (LR == 1) {
      return One*2./HVertex2Abb[0]

	+(-2.*HVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(HVertex2Abb[16]*B_0(m_s,m_m2,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(2.*HVertex2Abb[5]*HVertex2Abb[6]*B_0(m_s2k,0.,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4])

	+(2.*HVertex2Abb[9]*C_0(m_m2,m_s,m_s2k,0.,m_m2,m_m2,m_mu2))/(HVertex2Abb[0]*HVertex2Abb[4]);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}


DivArrC Higgs_Decay_RV_Vertices::RV_P_Vertex_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Photon Vertex corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(-2.*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[33]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[34]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],2.)

	+(PVertex2Abb[30]*PVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],2.);
    }
    else if (LR == 1) {
      return One*(-2.*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[33]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[34]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],2.)

	+(PVertex2Abb[30]*PVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],2.);
    }
  }
  // ubar1 \slashed{k} P_i v2
  else if (ME == 2) {
    if (LR == 0) {
      return One*(PVertex2Abb[12]*PVertex2Abb[14]*(epsP*m_p2))/pow(PVertex2Abb[0],3.)

	+(PVertex2Abb[16]*PVertex2Abb[23]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[15],3.)

	+(PVertex2Abb[12]*PVertex2Abb[17]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],3.)

	+(PVertex2Abb[24]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[0];
    }
    else if (LR == 1) {
      return One*(PVertex2Abb[12]*PVertex2Abb[14]*(epsP*m_p2))/pow(PVertex2Abb[0],3.)

	+(PVertex2Abb[16]*PVertex2Abb[23]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[15],3.)

	+(PVertex2Abb[12]*PVertex2Abb[17]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[0],3.)

	+(PVertex2Abb[24]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} P_i v2
  else if (ME == 3) {
    if (LR == 0) {
      return (PVertex2Abb[11]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex2Abb[0]

	+(PVertex2Abb[12]*B_0(m_s2k,0.,m_m2,m_mu2))/PVertex2Abb[0];
    }
    else if (LR == 1) {
      return (PVertex2Abb[11]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex2Abb[0]

	+(PVertex2Abb[12]*B_0(m_s2k,0.,m_m2,m_mu2))/PVertex2Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*2./PVertex2Abb[0]

	+(2.*PVertex2Abb[1]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex2Abb[0],2.)

	+(-PVertex2Abb[2]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(PVertex2Abb[0],2.)

	+(PVertex2Abb[6]*C_0(m_m2,0.,m_s2k,0.,m_m2,m_m2,m_mu2))/PVertex2Abb[0];
    }
    else if (LR == 1) {
      return One*2./PVertex2Abb[0]

	+(2.*PVertex2Abb[1]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex2Abb[0],2.)

	+(-PVertex2Abb[2]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(PVertex2Abb[0],2.)

	+(PVertex2Abb[6]*C_0(m_m2,0.,m_s2k,0.,m_m2,m_m2,m_mu2))/PVertex2Abb[0];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 4.";
  }
  return Zero;
}

