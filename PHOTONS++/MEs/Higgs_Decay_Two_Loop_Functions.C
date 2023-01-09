#include "PHOTONS++/MEs/Higgs_Decay_Two_Loop_Functions.H"
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

Higgs_Decay_Two_Loop_Functions::Higgs_Decay_Two_Loop_Functions
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  m_m(m), m_s(s), m_mu2(mu2), m_p1(p1), m_p2(p2), m_pP(pP),
  m_cL(cL), m_cR(cR)
{
  // Set parameters in line with p_EW
  p_EW = new EW_One_Loop_Functions_Base(s,mu2,0);
  m_e = p_EW->Get_e();
  m_sw = p_EW->Get_sw();
  m_sw2 = sqr(m_sw);
  m_cw = p_EW->Get_cw();
  m_cw2 = sqr(m_cw);
  MW2 = p_EW->Get_MW2();
  m_s12 = (m_p1+m_p2)*(m_p1+m_p2);
  m_s1k = (m_p1+m_pP)*(m_p1+m_pP);
  m_s2k = (m_p2+m_pP)*(m_p2+m_pP);
  m_p12 = m_p1*m_p1;
  m_p22 = m_p2*m_p2;
  m_pP2 = m_pP*m_pP;
  m_m2 = m_m*m_m;
  Zero = DivArrC(0.,0.,0.,0.,0.,0.);
  One = DivArrC(0.,0.,0.,1.,0.,0.);
  p_bubbles = new Higgs_Decay_RV_Bubbles(m,s,p1,p2,pP,cL,cR,mu2);
  p_vertices = new Higgs_Decay_RV_Vertices(m,s,p1,p2,pP,cL,cR,mu2);
  p_box_1 = new Higgs_Decay_RV_Box_1(m,s,p1,p2,pP,cL,cR,mu2);
  p_box_2 = new Higgs_Decay_RV_Box_2(m,s,p1,p2,pP,cL,cR,mu2);
  p_ct = new Higgs_Decay_RV_CT(m,s,p1,p2,pP,cL,cR,mu2);
  DivArrC p_bubblecoeff_1[4][2];
  DivArrC p_bubblecoeff_2[4][2];
  DivArrC p_H_vertexcoeff_1[4][2];
  DivArrC p_H_vertexcoeff_2[4][2];
  DivArrC p_P_vertexcoeff_1[4][2];
  DivArrC p_P_vertexcoeff_2[4][2];
  DivArrC p_boxcoeff_1[4][2];
  DivArrC p_boxcoeff_2[4][2];
  DivArrC p_fermion_ctcoeff[4][2];
  DivArrC p_vertex_ctcoeff[4][2];
  DivArrC p_bcoeff[4][2];
}

Higgs_Decay_Two_Loop_Functions::~Higgs_Decay_Two_Loop_Functions()
{
  delete p_EW;
  delete p_bubbles;
  delete p_vertices;
  delete p_box_1;
  delete p_box_2;
  delete p_ct;
}


void Higgs_Decay_Two_Loop_Functions::Calculate_RV_Coeffs(const Vec4C& epsP)
{
  // Calculate and store RV coefficients for each standard ME
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      p_bubblecoeff_1[i][j] = p_bubbles->RV_Bubble_Insertion_1(i+1,j,epsP);
      p_bubblecoeff_2[i][j] = p_bubbles->RV_Bubble_Insertion_2(i+1,j,epsP);
      p_H_vertexcoeff_1[i][j] = p_vertices->RV_H_Vertex_1(i+1,j,epsP);
      p_H_vertexcoeff_2[i][j] = p_vertices->RV_H_Vertex_2(i+1,j,epsP);
      p_P_vertexcoeff_1[i][j] = p_vertices->RV_P_Vertex_1(i+1,j,epsP);
      p_P_vertexcoeff_2[i][j] = p_vertices->RV_P_Vertex_2(i+1,j,epsP);
      p_boxcoeff_1[i][j] = p_box_1->RV_Box_1(i+1,j,epsP);
      p_boxcoeff_2[i][j] = p_box_2->RV_Box_2(i+1,j,epsP);
      p_fermion_ctcoeff[i][j] = p_ct->RV_Fermion_CT(i+1,j,epsP);
      p_vertex_ctcoeff[i][j] = p_ct->RV_Vertex_CT(i+1,j,epsP);
      p_bcoeff[i][j] = p_ct->RV_B(i+1,j,epsP);
    }
  }
}

