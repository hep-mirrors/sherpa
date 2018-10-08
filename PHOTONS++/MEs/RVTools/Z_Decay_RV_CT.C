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

// Caculate counterterms for real-virtual amplitudes in Z-decay
Z_Decay_RV_CT::Z_Decay_RV_CT
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Z_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  p_EW = new EW_One_Loop_Functions_Base(s,mu2,0);
  COne = Complex(1.,0.);
  B = Zero;
  if (m_m2 != 0) {
    B +=  (0.5*(m_s12-2.*m_m2)*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	   +m_m2*C_0(m_m2,m_m2,0.,0.,m_m2,m_m2,m_mu2)
	   +0.25*(B_0(m_s12,m_m2,m_m2,m_mu2)-B_0(0.,m_m2,m_m2,m_mu2)));
  }
  else if (m_m2 == 0) {
    B +=  (0.5*(m_s12-2.*m_m2)*C_0(m_m2,m_m2,m_s12,0.,m_m2,m_m2,m_mu2)
	   +0.5*DivArrC(0.,1.,0.,0.,0.,0.)
	   +0.25*(B_0(m_s12,m_m2,m_m2,m_mu2)-B_0(0.,m_m2,m_m2,m_mu2)));
  }      
}


Z_Decay_RV_CT::~Z_Decay_RV_CT()
{
  delete p_EW;
}

void Z_Decay_RV_CT::Init_Coefficients()
{
  return;
}

DivArrC Z_Decay_RV_CT::Z_Vertex_CT_L()
{  
  // Counterterm for Z-decay vertex - note: no need to include ZA transition in QED, also contains
  // fermion wavefunction counterterm for external fermions, no CT necessary for vertex besides
  // fermion mass. Take out cL as this is taken into account in StandardMEs
  return m_cL*(1./2.*p_EW->dZZZ()
	       +2.*p_EW->dZfermL(m_m2,0.,-1.,-0.5));
}

DivArrC Z_Decay_RV_CT::Z_Vertex_CT_R()
{  
  // Counterterm for Z-decay vertex - note: no need to include ZA transition in QED, also contains
  // fermion wavefunction counterterm for external fermions, no CT necessary for vertex besides
  // fermion mass. Take out cR as this is taken into account in StandardMEs
  return m_cR*(1./2.*p_EW->dZZZ()
	       +2.*p_EW->dZfermR(m_m2,0.,-1.,-0.5));
}

DivArrC Z_Decay_RV_CT::RV_Vertex_CT_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Vertex CT corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME != 4 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return Z_Vertex_CT_L()*2.*(epsP*m_p1)/(m_s1k-m_m2);
    }
    else if (LR == 1) {
      return Z_Vertex_CT_R()*2.*(epsP*m_p1)/(m_s1k-m_m2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return Z_Vertex_CT_L()/(m_s1k-m_m2);
    }
    else if (LR == 1) {
      return Z_Vertex_CT_R()/(m_s1k-m_m2);
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}

DivArrC Z_Decay_RV_CT::RV_Vertex_CT_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Vertex CT corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME != 4 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return -Z_Vertex_CT_L()*2.*(epsP*m_p2)/(m_s2k-m_m2);
    }
    else if (LR == 1) {
      return -Z_Vertex_CT_R()*2.*(epsP*m_p2)/(m_s2k-m_m2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return -Z_Vertex_CT_L()/(m_s2k-m_m2);
    }
    else if (LR == 1) {
      return -Z_Vertex_CT_R()/(m_s2k-m_m2);
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}

DivArrC Z_Decay_RV_CT::RV_Fermion_CT_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Internal fermion CT for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,  
  // ubar1 P_i v2
  if (ME != 4 && ME != 7 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return -p_EW->dZfermL(m_m2,0.,-1.,-0.5)*m_cL*2.*(epsP*m_p1)/(m_s1k-m_m2)
	+2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL*2.*(epsP*m_p1)/(pow(m_s1k-m_m2,2.));
    }
    else if (LR == 1) {
      return -p_EW->dZfermR(m_m2,0.,-1.,-0.5)*m_cR*2.*(epsP*m_p1)/(m_s1k-m_m2)	
	+2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR*2.*(epsP*m_p1)/(pow(m_s1k-m_m2,2.));
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 7) {
    if (LR == 0) {
      return p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL/(m_s1k-m_m2);
    }
    else if (LR == 1) {
      return p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR/(m_s1k-m_m2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return -p_EW->dZfermL(m_m2,0.,-1.,-0.5)*m_cL/(m_s1k-m_m2)
	+2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL/(pow(m_s1k-m_m2,2.));
    }
    else if (LR == 1) {
      return -p_EW->dZfermR(m_m2,0.,-1.,-0.5)*m_cR/(m_s1k-m_m2)	
	+2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR/(pow(m_s1k-m_m2,2.));
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}

DivArrC Z_Decay_RV_CT::RV_Fermion_CT_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Internal fermion CT for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,  
  // ubar1 P_i v2
  if (ME != 4 && ME != 7 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return p_EW->dZfermL(m_m2,0.,-1.,-0.5)*m_cL*2.*(epsP*m_p2)/(m_s2k-m_m2)
	-2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL*2.*(epsP*m_p2)/(pow(m_s2k-m_m2,2.));
    }
    else if (LR == 1) {
      return p_EW->dZfermR(m_m2,0.,-1.,-0.5)*m_cR*2.*(epsP*m_p2)/(m_s2k-m_m2)
	-2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR*2.*(epsP*m_p2)/(pow(m_s2k-m_m2,2.));
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  // cR and cL are exchanged due to having to anticommute GammaV through epsP
  else if (ME == 7) {
    if (LR == 0) {
      return p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR/(m_s2k-m_m2);
    }
    else if (LR == 1) {
      return p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL/(m_s2k-m_m2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return p_EW->dZfermL(m_m2,0.,-1.,-0.5)*m_cL/(m_s2k-m_m2)
	-2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cL/(pow(m_s2k-m_m2,2.));
    }
    else if (LR == 1) {
      return p_EW->dZfermL(m_m2,0.,-1.,-0.5)*m_cR/(m_s2k-m_m2)
	-2.*m_m*p_EW->dm(m_m2,0.,-1.,-0.5)*m_cR/(pow(m_s2k-m_m2,2.));
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}


DivArrC Z_Decay_RV_CT::RV_B_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Contributions due to IR form factor B for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,  
  // ubar1 P_i v2
  if (ME != 4 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsV} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return m_cL*B*2.*(epsP*m_p1)/(m_s1k-m_m2);
    }
    else if (LR == 1) {
      return m_cR*B*2.*(epsP*m_p1)/(m_s1k-m_m2);
    }
  }
  else if (ME == 8) {
    if (LR == 0) {
      return m_cL*B/(m_s1k-m_m2);
    }
    else if (LR == 1) {
      return m_cR*B/(m_s1k-m_m2);
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}


DivArrC Z_Decay_RV_CT::RV_B_2(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Contributions due to IR form factor B for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,  
  // ubar1 P_i v2
  if (ME != 4 && ME != 8) {
    return Zero;
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return -m_cL*B*2.*(epsP*m_p2)/(m_s2k-m_m2);
    }
    else if (LR == 1) {
      return -m_cR*B*2.*(epsP*m_p2)/(m_s2k-m_m2);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return -m_cL*B/(m_s2k-m_m2);
    }
    else if (LR == 1) {
      return -m_cR*B/(m_s2k-m_m2);
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
}


