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

Z_Decay_RV_Vertices::Z_Decay_RV_Vertices
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Z_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  ZVertex1Abb = new Complex[281];
  ZVertex2Abb = new Complex[285];
  PVertex1Abb = new Complex[63];
  PVertex2Abb = new Complex[64];
  Init_Coefficients();
}

Z_Decay_RV_Vertices::~Z_Decay_RV_Vertices()
{
  delete [] ZVertex1Abb;
  delete [] PVertex1Abb;
  delete [] ZVertex2Abb;
  delete [] PVertex2Abb;
}

void Z_Decay_RV_Vertices::Init_Coefficients() 
{
  Init_Z_Vertex_1_Coefficients();
  Init_Z_Vertex_2_Coefficients();
  Init_P_Vertex_1_Coefficients();
  Init_P_Vertex_2_Coefficients();
  return;
}

// Set up coefficients multiplying scalar master integrals in corrections to Z vertex for emission off leg 1
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Vertices::Init_Z_Vertex_1_Coefficients() 
{
  ZVertex1Abb[0]=m_s*m_x - m_s*m_z1k;

  ZVertex1Abb[1]=3. - 5.*m_x - 3.*m_z1k;

  ZVertex1Abb[2]=ZVertex1Abb[1]*m_cL - 4.*m_cR*m_x;

  ZVertex1Abb[3]=m_x - m_z1k;

  ZVertex1Abb[4]=-1. + m_z1k;

  ZVertex1Abb[5]=1. + m_z1k;

  ZVertex1Abb[6]=pow(ZVertex1Abb[4],2.) - 2.*ZVertex1Abb[5]*m_x + m_x_2;

  ZVertex1Abb[7]=-pow(ZVertex1Abb[4],2.) + m_x + m_x_2;

  ZVertex1Abb[8]=1. + m_x - m_z1k;

  ZVertex1Abb[9]=ZVertex1Abb[7]*m_cL + ZVertex1Abb[8]*m_cR*m_x;

  ZVertex1Abb[10]=1. - m_x + m_z1k;

  ZVertex1Abb[11]=5. + 2.*m_z1k;

  ZVertex1Abb[12]=-2. + ZVertex1Abb[11]*m_x - 3.*m_x_2 + m_z1k + m_z1k_2;

  ZVertex1Abb[13]=ZVertex1Abb[12]*m_cL + 2.*ZVertex1Abb[10]*m_cR*m_x;

  ZVertex1Abb[14]=-4. + 3.*m_z1k;

  ZVertex1Abb[15]=pow(ZVertex1Abb[4],2.) + ZVertex1Abb[14]*m_x + 4.*m_x_2;

  ZVertex1Abb[16]=ZVertex1Abb[15]*m_cL + ZVertex1Abb[3]*ZVertex1Abb[8]*m_cR*m_x;

  ZVertex1Abb[24]=2.*m_cL;

  ZVertex1Abb[25]=m_s*m_x;

  ZVertex1Abb[26]=1/m_s;

  ZVertex1Abb[27]=2./m_s;

  ZVertex1Abb[28]=m_s*m_z1k;

  ZVertex1Abb[29]=ZVertex1Abb[1]*m_cR - 4.*m_cL*m_x;

  ZVertex1Abb[30]=ZVertex1Abb[7]*m_cR + ZVertex1Abb[8]*m_cL*m_x;

  ZVertex1Abb[31]=ZVertex1Abb[12]*m_cR + 2.*ZVertex1Abb[10]*m_cL*m_x;

  ZVertex1Abb[32]=ZVertex1Abb[15]*m_cR + ZVertex1Abb[3]*ZVertex1Abb[8]*m_cL*m_x;

  ZVertex1Abb[40]=2.*m_cR;

  ZVertex1Abb[41]=-1. + 3.*m_x + m_z1k;

  ZVertex1Abb[46]=(4.*m_cR*m_m)/m_s;

  ZVertex1Abb[47]=(-2.*m_cR*m_m)/m_s;

  ZVertex1Abb[48]=-2.*m_cR*m_m;

  ZVertex1Abb[53]=(4.*m_cL*m_m)/m_s;

  ZVertex1Abb[54]=(-2.*m_cL*m_m)/m_s;

  ZVertex1Abb[55]=-2.*m_cL*m_m;

  ZVertex1Abb[56]=-1. + m_x + m_z1k;

  ZVertex1Abb[57]=-ZVertex1Abb[41]*m_cR + 2.*m_cL*m_x;

  ZVertex1Abb[58]=m_x + 3.*m_z1k;

  ZVertex1Abb[59]=ZVertex1Abb[41]*m_cL + m_cR - ZVertex1Abb[58]*m_cR;

  ZVertex1Abb[60]=-m_x + m_z1k;

  ZVertex1Abb[61]=2. + 3.*m_z1k;

  ZVertex1Abb[62]=7. - 3.*m_z1k;

  ZVertex1Abb[63]=-4. + ZVertex1Abb[62]*m_z1k;

  ZVertex1Abb[64]=-pow(ZVertex1Abb[4],3.) + ZVertex1Abb[63]*m_x + ZVertex1Abb[61]*m_x_2 + m_x_3;

  ZVertex1Abb[65]=ZVertex1Abb[64]*m_cR + 6.*ZVertex1Abb[60]*m_cL*m_x_2;

  ZVertex1Abb[66]=-5. + 4.*m_x;

  ZVertex1Abb[67]=1. + ZVertex1Abb[66]*m_x + m_z1k - 2.*m_x*m_z1k - 2.*m_z1k_2;

  ZVertex1Abb[68]=-1. + m_x;

  ZVertex1Abb[69]=1. + m_x;

  ZVertex1Abb[70]=-5. + 7.*m_x;

  ZVertex1Abb[71]=7. - 5.*m_x;

  ZVertex1Abb[72]=pow(ZVertex1Abb[68],2.)*ZVertex1Abb[69] + ZVertex1Abb[69]*ZVertex1Abb[70]*m_z1k + ZVertex1Abb[71]*m_z1k_2 - 3.*m_z1k_3;

  ZVertex1Abb[73]=ZVertex1Abb[72]*m_cR - 2.*ZVertex1Abb[67]*m_cL*m_x;

  ZVertex1Abb[74]=9. - 11.*m_z1k;

  ZVertex1Abb[75]=3. + 7.*m_z1k;

  ZVertex1Abb[76]=-pow(ZVertex1Abb[4],3.) + ZVertex1Abb[4]*ZVertex1Abb[75]*m_x + ZVertex1Abb[74]*m_x_2 + 5.*m_x_3;

  ZVertex1Abb[77]=6. + m_z1k;

  ZVertex1Abb[78]=6. - 5.*m_z1k;

  ZVertex1Abb[79]=-9. + ZVertex1Abb[78]*m_z1k;

  ZVertex1Abb[80]=pow(ZVertex1Abb[4],2.)*ZVertex1Abb[61] + ZVertex1Abb[79]*m_x + ZVertex1Abb[77]*m_x_2 + m_x_3;

  ZVertex1Abb[81]=ZVertex1Abb[76]*m_cL - ZVertex1Abb[80]*m_cR;

  ZVertex1Abb[82]=-1. + m_x_2;

  ZVertex1Abb[83]=-2. + 5.*m_x;

  ZVertex1Abb[84]=ZVertex1Abb[83]*m_cL + m_cR - m_cR*m_x;

  ZVertex1Abb[85]=8. + m_x;

  ZVertex1Abb[86]=-5. + ZVertex1Abb[85]*m_x;

  ZVertex1Abb[87]=-6. + 11.*m_x;

  ZVertex1Abb[88]=5. + ZVertex1Abb[87]*m_x;

  ZVertex1Abb[89]=ZVertex1Abb[88]*m_cL + ZVertex1Abb[86]*m_cR;

  ZVertex1Abb[90]=-1. + 5.*m_x;

  ZVertex1Abb[91]=4. + 7.*m_x;

  ZVertex1Abb[92]=ZVertex1Abb[91]*m_cL + ZVertex1Abb[90]*m_cR;

  ZVertex1Abb[93]=m_cL + 3.*m_cR;

  ZVertex1Abb[94]=ZVertex1Abb[82]*ZVertex1Abb[84] - ZVertex1Abb[89]*m_z1k + ZVertex1Abb[92]*m_z1k_2 - ZVertex1Abb[93]*m_z1k_3;

  ZVertex1Abb[95]=-1. + 2.*ZVertex1Abb[69]*m_x + 2.*m_z1k - m_x*m_z1k - m_z1k_2;

  ZVertex1Abb[96]=ZVertex1Abb[56]*ZVertex1Abb[95]*m_cR + 6.*ZVertex1Abb[60]*m_cL*m_x_2;

  ZVertex1Abb[97]=-2. + 4.*m_z1k;

  ZVertex1Abb[98]=pow(ZVertex1Abb[4],2.) + ZVertex1Abb[97]*m_x + m_x_2;

  ZVertex1Abb[99]=7. + 4.*m_z1k;

  ZVertex1Abb[100]=-2. + ZVertex1Abb[99]*m_x - 5.*m_x_2 + m_z1k + m_z1k_2;

  ZVertex1Abb[101]=ZVertex1Abb[3]*ZVertex1Abb[98]*m_cR + ZVertex1Abb[100]*m_cL*m_x;

  ZVertex1Abb[102]=1. + 4.*m_z1k;

  ZVertex1Abb[103]=pow(ZVertex1Abb[4],2.) - 2.*ZVertex1Abb[102]*m_x + m_x_2;

  ZVertex1Abb[104]=-8. + m_z1k;

  ZVertex1Abb[105]=-2. + 5.*m_z1k;

  ZVertex1Abb[106]=m_x + ZVertex1Abb[105]*m_x_2 + m_x_3 - 5.*pow(ZVertex1Abb[4],2.)*m_z1k - ZVertex1Abb[104]*m_x*m_z1k;

  ZVertex1Abb[107]=-ZVertex1Abb[103]*ZVertex1Abb[3]*ZVertex1Abb[56]*m_cL - ZVertex1Abb[106]*m_cR*m_z1k;

  ZVertex1Abb[108]=-2. + 3.*m_x;

  ZVertex1Abb[109]=-4. + 7.*m_x;

  ZVertex1Abb[110]=5. + ZVertex1Abb[109]*m_x;

  ZVertex1Abb[111]=4. + 11.*m_x;

  ZVertex1Abb[112]=ZVertex1Abb[108]*pow(ZVertex1Abb[68],2.) + ZVertex1Abb[110]*m_z1k - ZVertex1Abb[111]*m_z1k_2 + m_z1k_3;

  ZVertex1Abb[113]=ZVertex1Abb[112]*m_cL - 6.*ZVertex1Abb[3]*ZVertex1Abb[56]*m_cR*m_z1k;

  ZVertex1Abb[127]=(-2.*m_cL*m_m*m_x)/(m_s_2*m_z1k);

  ZVertex1Abb[128]=(2.*m_m)/m_s_2;

  ZVertex1Abb[129]=(-4.*m_m)/m_s_2;

  ZVertex1Abb[130]=(-2.*m_m)/m_s_2;

  ZVertex1Abb[131]=(-4.*m_m)/m_s;

  ZVertex1Abb[132]=(-2.*m_m)/(m_s_2*m_z1k);

  ZVertex1Abb[133]=-ZVertex1Abb[41]*m_cL + 2.*m_cR*m_x;

  ZVertex1Abb[134]=-1. + m_x + 3.*m_z1k;

  ZVertex1Abb[135]=-ZVertex1Abb[134]*m_cL + ZVertex1Abb[41]*m_cR;

  ZVertex1Abb[136]=ZVertex1Abb[64]*m_cL + 6.*ZVertex1Abb[60]*m_cR*m_x_2;

  ZVertex1Abb[137]=ZVertex1Abb[72]*m_cL - 2.*ZVertex1Abb[67]*m_cR*m_x;

  ZVertex1Abb[138]=ZVertex1Abb[56]*ZVertex1Abb[95]*m_cL + 6.*ZVertex1Abb[60]*m_cR*m_x_2;

  ZVertex1Abb[139]=ZVertex1Abb[3]*ZVertex1Abb[98]*m_cL + ZVertex1Abb[100]*m_cR*m_x;

  ZVertex1Abb[140]=-9. + 11.*m_z1k;

  ZVertex1Abb[141]=4. - 7.*m_z1k;

  ZVertex1Abb[142]=3. + ZVertex1Abb[141]*m_z1k;

  ZVertex1Abb[143]=pow(ZVertex1Abb[4],3.) + ZVertex1Abb[142]*m_x + ZVertex1Abb[140]*m_x_2 - 5.*m_x_3;

  ZVertex1Abb[144]=ZVertex1Abb[80]*m_cL + ZVertex1Abb[143]*m_cR;

  ZVertex1Abb[145]=2. - 5.*m_x;

  ZVertex1Abb[146]=ZVertex1Abb[68]*m_cL + ZVertex1Abb[145]*m_cR;

  ZVertex1Abb[147]=ZVertex1Abb[86]*m_cL + ZVertex1Abb[88]*m_cR;

  ZVertex1Abb[148]=m_cL - ZVertex1Abb[91]*m_cR - 5.*m_cL*m_x;

  ZVertex1Abb[149]=3.*m_cL + m_cR;

  ZVertex1Abb[150]=ZVertex1Abb[146]*ZVertex1Abb[82] + ZVertex1Abb[147]*m_z1k + ZVertex1Abb[148]*m_z1k_2 + ZVertex1Abb[149]*m_z1k_3;

  ZVertex1Abb[151]=-ZVertex1Abb[103]*ZVertex1Abb[3]*ZVertex1Abb[56]*m_cR - ZVertex1Abb[106]*m_cL*m_z1k;

  ZVertex1Abb[152]=ZVertex1Abb[112]*m_cR - 6.*ZVertex1Abb[3]*ZVertex1Abb[56]*m_cL*m_z1k;

  ZVertex1Abb[166]=(-2.*m_cR*m_m*m_x)/(m_s_2*m_z1k);

  ZVertex1Abb[167]=-3. + 5.*m_x + 3.*m_z1k;

  ZVertex1Abb[168]=ZVertex1Abb[167]*m_cL + 4.*m_cR*m_x;

  ZVertex1Abb[176]=4.*m_cL;

  ZVertex1Abb[177]=-2./m_s;

  ZVertex1Abb[178]=4./m_s;

  ZVertex1Abb[179]=ZVertex1Abb[167]*m_cR + 4.*m_cL*m_x;

  ZVertex1Abb[187]=4.*m_cR;

  ZVertex1Abb[188]=-2. + m_z1k;

  ZVertex1Abb[189]=-pow(ZVertex1Abb[4],2.) + 2.*ZVertex1Abb[188]*m_x - m_x_2;

  ZVertex1Abb[190]=4. - 5.*m_z1k;

  ZVertex1Abb[191]=4. + m_z1k;

  ZVertex1Abb[192]=2. + m_z1k;

  ZVertex1Abb[193]=-11. + ZVertex1Abb[192]*m_z1k;

  ZVertex1Abb[194]=-ZVertex1Abb[191]*pow(ZVertex1Abb[4],2.) - ZVertex1Abb[193]*m_x - ZVertex1Abb[190]*m_x_2 - 3.*m_x_3;

  ZVertex1Abb[195]=-pow(ZVertex1Abb[4],2.) - ZVertex1Abb[188]*m_x + 2.*m_x_2;

  ZVertex1Abb[196]=-3. + m_z1k;

  ZVertex1Abb[197]=3. + ZVertex1Abb[196]*m_z1k;

  ZVertex1Abb[198]=-1. + ZVertex1Abb[188]*ZVertex1Abb[4]*m_x + ZVertex1Abb[192]*m_x_2 - 3.*m_x_3 + ZVertex1Abb[197]*m_z1k;

  ZVertex1Abb[199]=-7. + m_z1k - 3.*m_z1k_2;

  ZVertex1Abb[200]=-5. + m_z1k + 3.*m_z1k_2;

  ZVertex1Abb[201]=-pow(ZVertex1Abb[4],3.) - ZVertex1Abb[200]*ZVertex1Abb[4]*m_x - ZVertex1Abb[199]*m_x_2 - 3.*m_x_4 + 3.*m_x_3*m_z1k;

  ZVertex1Abb[202]=1. - m_x - m_z1k;

  ZVertex1Abb[203]=-4. + 5.*m_z1k;

  ZVertex1Abb[204]=2. + 7.*m_z1k;

  ZVertex1Abb[205]=m_x - ZVertex1Abb[204]*m_x_2 + m_x_3 + pow(ZVertex1Abb[4],2.)*m_z1k + ZVertex1Abb[203]*m_x*m_z1k;

  ZVertex1Abb[206]=1. + 3.*m_x;

  ZVertex1Abb[207]=4. - 3.*m_x;

  ZVertex1Abb[208]=ZVertex1Abb[108]*pow(ZVertex1Abb[68],2.) + ZVertex1Abb[206]*ZVertex1Abb[69]*m_z1k + ZVertex1Abb[207]*m_z1k_2 - 3.*m_z1k_3;

  ZVertex1Abb[221]=(2.*m_cL*m_x)/(m_s*m_z1k);

  ZVertex1Abb[222]=(-4.*m_cL*m_x)/m_s;

  ZVertex1Abb[223]=(-2.*m_cL)/m_s;

  ZVertex1Abb[224]=(-8.*m_cL*m_x)/m_s;

  ZVertex1Abb[225]=(4.*m_cL)/m_s;

  ZVertex1Abb[226]=-8.*m_cL*m_x;

  ZVertex1Abb[227]=(2.*m_cL)/(m_s*m_z1k);

  ZVertex1Abb[228]=(2.*m_cL)/m_s;

  ZVertex1Abb[241]=(2.*m_cR*m_x)/(m_s*m_z1k);

  ZVertex1Abb[242]=(-4.*m_cR*m_x)/m_s;

  ZVertex1Abb[243]=(-2.*m_cR)/m_s;

  ZVertex1Abb[244]=(-8.*m_cR*m_x)/m_s;

  ZVertex1Abb[245]=(4.*m_cR)/m_s;

  ZVertex1Abb[246]=-8.*m_cR*m_x;

  ZVertex1Abb[247]=(2.*m_cR)/(m_s*m_z1k);

  ZVertex1Abb[248]=(2.*m_cR)/m_s;

  ZVertex1Abb[262]=(-4.*m_cL*m_m*m_x)/(m_s_2*m_z1k);

  ZVertex1Abb[263]=(4.*m_m)/m_s_2;

  ZVertex1Abb[264]=(-8.*m_m)/m_s_2;

  ZVertex1Abb[265]=(-8.*m_m)/m_s;

  ZVertex1Abb[266]=(-4.*m_m)/(m_s_2*m_z1k);

  ZVertex1Abb[280]=(-4.*m_cR*m_m*m_x)/(m_s_2*m_z1k);
}

// Set up coefficients multiplying scalar master integrals in corrections to photon vertex for emission off leg 1
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Vertices::Init_P_Vertex_1_Coefficients() 
{
  PVertex1Abb[0]=m_s*m_x - m_s*m_z1k;

  PVertex1Abb[1]=m_x - m_z1k;

  PVertex1Abb[2]=2.*m_x + m_z1k;

  PVertex1Abb[3]=5.*m_x + m_z1k;

  PVertex1Abb[8]=2.*m_cL;

  PVertex1Abb[9]=m_s*m_x;

  PVertex1Abb[10]=(2.*m_cL)/m_s;

  PVertex1Abb[11]=m_s*m_z1k;

  PVertex1Abb[12]=(-m_cL)/m_s;

  PVertex1Abb[13]=2.*m_cL*m_x;

  PVertex1Abb[18]=2.*m_cR;

  PVertex1Abb[19]=(2.*m_cR)/m_s;

  PVertex1Abb[20]=(-m_cR)/m_s;

  PVertex1Abb[21]=2.*m_cR*m_x;

  PVertex1Abb[23]=-2.*m_cL*m_m;

  PVertex1Abb[24]=2.*m_cL*m_m;

  PVertex1Abb[26]=-2.*m_cR*m_m;

  PVertex1Abb[27]=2.*m_cR*m_m;

  PVertex1Abb[28]=-3.*m_x + m_z1k;

  PVertex1Abb[29]=-m_x + m_z1k;

  PVertex1Abb[30]=m_x + m_z1k;

  PVertex1Abb[31]=3.*m_x - m_z1k;

  PVertex1Abb[37]=(2.*m_cL*m_m)/m_s_2;

  PVertex1Abb[38]=(2.*m_cL*m_m*m_x)/(m_s_2*m_z1k);

  PVertex1Abb[39]=(2.*m_cL*m_m)/(m_s_2*m_z1k);

  PVertex1Abb[45]=(2.*m_cR*m_m)/m_s_2;

  PVertex1Abb[46]=(2.*m_cR*m_m*m_x)/(m_s_2*m_z1k);

  PVertex1Abb[47]=(2.*m_cR*m_m)/(m_s_2*m_z1k);

  PVertex1Abb[48]=m_x_2 - 6.*m_x*m_z1k + m_z1k_2;

  PVertex1Abb[52]=(-2.*m_cL)/m_s;

  PVertex1Abb[53]=(-2.*m_cL*m_x)/(m_s*m_z1k);

  PVertex1Abb[54]=(12.*m_cL*m_x)/m_s;

  PVertex1Abb[55]=(2.*m_cL)/(m_s*m_z1k);

  PVertex1Abb[59]=(-2.*m_cR)/m_s;

  PVertex1Abb[60]=(-2.*m_cR*m_x)/(m_s*m_z1k);

  PVertex1Abb[61]=(12.*m_cR*m_x)/m_s;

  PVertex1Abb[62]=(2.*m_cR)/(m_s*m_z1k);
}

// Set up coefficients multiplying scalar master integrals in corrections to Z vertex for emission off leg 2
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Vertices::Init_Z_Vertex_2_Coefficients() 
{
  ZVertex2Abb[0]=m_s*m_x - m_s*m_z2k;

  ZVertex2Abb[1]=m_x - m_z2k;

  ZVertex2Abb[2]=-1. + m_z2k;

  ZVertex2Abb[3]=1. + m_z2k;

  ZVertex2Abb[4]=pow(ZVertex2Abb[2],2.) - 2.*ZVertex2Abb[3]*m_x + m_x_2;

  ZVertex2Abb[5]=-3. + 5.*m_x + 3.*m_z2k;

  ZVertex2Abb[6]=ZVertex2Abb[5]*m_cL + 4.*m_cR*m_x;

  ZVertex2Abb[7]=-pow(ZVertex2Abb[2],2.) + m_x + m_x_2;

  ZVertex2Abb[8]=1. + m_x - m_z2k;

  ZVertex2Abb[9]=ZVertex2Abb[7]*m_cL + ZVertex2Abb[8]*m_cR*m_x;

  ZVertex2Abb[10]=-1. + m_x;

  ZVertex2Abb[11]=-2. + 3.*m_x;

  ZVertex2Abb[12]=ZVertex2Abb[11]*m_cL + 2.*m_cR*m_x;

  ZVertex2Abb[13]=m_cL + m_cR;

  ZVertex2Abb[14]=m_cL + 2.*ZVertex2Abb[13]*m_x;

  ZVertex2Abb[15]=ZVertex2Abb[10]*ZVertex2Abb[12] - ZVertex2Abb[14]*m_z2k - m_cL*m_z2k_2;

  ZVertex2Abb[16]=-4. + 3.*m_z2k;

  ZVertex2Abb[17]=pow(ZVertex2Abb[2],2.) + ZVertex2Abb[16]*m_x + 4.*m_x_2;

  ZVertex2Abb[18]=ZVertex2Abb[17]*m_cL + ZVertex2Abb[1]*ZVertex2Abb[8]*m_cR*m_x;

  ZVertex2Abb[26]=-2.*m_cL;

  ZVertex2Abb[27]=m_s*m_x;

  ZVertex2Abb[28]=1./m_s;

  ZVertex2Abb[29]=-2./m_s;

  ZVertex2Abb[30]=m_s*m_z2k;

  ZVertex2Abb[31]=ZVertex2Abb[5]*m_cR + 4.*m_cL*m_x;

  ZVertex2Abb[32]=ZVertex2Abb[7]*m_cR + ZVertex2Abb[8]*m_cL*m_x;

  ZVertex2Abb[33]=ZVertex2Abb[11]*m_cR + 2.*m_cL*m_x;

  ZVertex2Abb[34]=m_cR + 2.*ZVertex2Abb[13]*m_x;

  ZVertex2Abb[35]=ZVertex2Abb[10]*ZVertex2Abb[33] - ZVertex2Abb[34]*m_z2k - m_cR*m_z2k_2;

  ZVertex2Abb[36]=ZVertex2Abb[17]*m_cR + ZVertex2Abb[1]*ZVertex2Abb[8]*m_cL*m_x;

  ZVertex2Abb[44]=-2.*m_cR;

  ZVertex2Abb[45]=1. - m_x + m_z2k;

  ZVertex2Abb[46]=-1. + 3.*m_x + m_z2k;

  ZVertex2Abb[51]=(4.*m_cL*m_m)/m_s;

  ZVertex2Abb[52]=(-2.*m_cL*m_m)/m_s;

  ZVertex2Abb[53]=-2.*m_cL*m_m;

  ZVertex2Abb[58]=(4.*m_cR*m_m)/m_s;

  ZVertex2Abb[59]=(-2.*m_cR*m_m)/m_s;

  ZVertex2Abb[60]=-2.*m_cR*m_m;

  ZVertex2Abb[61]=-1. + m_x + m_z2k;

  ZVertex2Abb[62]=ZVertex2Abb[61]*m_cL - 2.*m_cR*m_x;

  ZVertex2Abb[63]=m_x + m_z2k;

  ZVertex2Abb[64]=-1. + m_x + 3.*m_z2k;

  ZVertex2Abb[65]=ZVertex2Abb[64]*m_cL + m_cR - ZVertex2Abb[63]*m_cR;

  ZVertex2Abb[66]=ZVertex2Abb[61]*m_cR - 2.*m_cL*m_z2k;

  ZVertex2Abb[67]=-m_x + m_z2k;

  ZVertex2Abb[68]=1. + m_x;

  ZVertex2Abb[69]=-1. + 2.*ZVertex2Abb[68]*m_x + 2.*m_z2k - m_x*m_z2k - m_z2k_2;

  ZVertex2Abb[70]=ZVertex2Abb[61]*ZVertex2Abb[69]*m_cL + 6.*ZVertex2Abb[67]*m_cR*m_x_2;

  ZVertex2Abb[71]=-2. + 7.*m_x - 5.*m_x_2 + m_z2k + 4.*m_x*m_z2k + m_z2k_2;

  ZVertex2Abb[72]=-5. + 7.*m_x;

  ZVertex2Abb[73]=7. - 5.*m_x;

  ZVertex2Abb[74]=pow(ZVertex2Abb[10],2.)*ZVertex2Abb[68] + ZVertex2Abb[68]*ZVertex2Abb[72]*m_z2k + ZVertex2Abb[73]*m_z2k_2 - 3.*m_z2k_3;

  ZVertex2Abb[75]=ZVertex2Abb[74]*m_cL + 2.*ZVertex2Abb[71]*m_cR*m_x;

  ZVertex2Abb[76]=-2. + 4.*m_z2k;

  ZVertex2Abb[77]=pow(ZVertex2Abb[2],2.) + ZVertex2Abb[76]*m_x + m_x_2;

  ZVertex2Abb[78]=7. + 4.*m_z2k;

  ZVertex2Abb[79]=-2. + ZVertex2Abb[78]*m_x - 5.*m_x_2 + m_z2k + m_z2k_2;

  ZVertex2Abb[80]=ZVertex2Abb[1]*ZVertex2Abb[77]*m_cL + ZVertex2Abb[79]*m_cR*m_x;

  ZVertex2Abb[81]=-9. + 11.*m_z2k;

  ZVertex2Abb[82]=4. - 7.*m_z2k;

  ZVertex2Abb[83]=3. + ZVertex2Abb[82]*m_z2k;

  ZVertex2Abb[84]=pow(ZVertex2Abb[2],3.) + ZVertex2Abb[83]*m_x + ZVertex2Abb[81]*m_x_2 - 5.*m_x_3;

  ZVertex2Abb[85]=6. + m_z2k;

  ZVertex2Abb[86]=2. + 3.*m_z2k;

  ZVertex2Abb[87]=6. - 5.*m_z2k;

  ZVertex2Abb[88]=-9. + ZVertex2Abb[87]*m_z2k;

  ZVertex2Abb[89]=pow(ZVertex2Abb[2],2.)*ZVertex2Abb[86] + ZVertex2Abb[88]*m_x + ZVertex2Abb[85]*m_x_2 + m_x_3;

  ZVertex2Abb[90]=ZVertex2Abb[89]*m_cL + ZVertex2Abb[84]*m_cR;

  ZVertex2Abb[91]=2. - 5.*m_x;

  ZVertex2Abb[92]=ZVertex2Abb[10]*m_cL + ZVertex2Abb[91]*m_cR;

  ZVertex2Abb[93]=-1. + m_x_2;

  ZVertex2Abb[94]=8. + m_x;

  ZVertex2Abb[95]=-5. + ZVertex2Abb[94]*m_x;

  ZVertex2Abb[96]=-6. + 11.*m_x;

  ZVertex2Abb[97]=5. + ZVertex2Abb[96]*m_x;

  ZVertex2Abb[98]=ZVertex2Abb[95]*m_cL + ZVertex2Abb[97]*m_cR;

  ZVertex2Abb[99]=4. + 7.*m_x;

  ZVertex2Abb[100]=m_cL - ZVertex2Abb[99]*m_cR - 5.*m_cL*m_x;

  ZVertex2Abb[101]=3.*m_cL + m_cR;

  ZVertex2Abb[102]=ZVertex2Abb[92]*ZVertex2Abb[93] + ZVertex2Abb[98]*m_z2k + ZVertex2Abb[100]*m_z2k_2 + ZVertex2Abb[101]*m_z2k_3;

  ZVertex2Abb[103]=1. + 4.*m_z2k;

  ZVertex2Abb[104]=pow(ZVertex2Abb[2],2.) - 2.*ZVertex2Abb[103]*m_x + m_x_2;

  ZVertex2Abb[105]=-8. + m_z2k;

  ZVertex2Abb[106]=-2. + 5.*m_z2k;

  ZVertex2Abb[107]=m_x + ZVertex2Abb[106]*m_x_2 + m_x_3 - 5.*pow(ZVertex2Abb[2],2.)*m_z2k - ZVertex2Abb[105]*m_x*m_z2k;

  ZVertex2Abb[108]=-ZVertex2Abb[1]*ZVertex2Abb[104]*ZVertex2Abb[61]*m_cR - ZVertex2Abb[107]*m_cL*m_z2k;

  ZVertex2Abb[109]=-4. + 7.*m_x;

  ZVertex2Abb[110]=5. + ZVertex2Abb[109]*m_x;

  ZVertex2Abb[111]=4. + 11.*m_x;

  ZVertex2Abb[112]=pow(ZVertex2Abb[10],2.)*ZVertex2Abb[11] + ZVertex2Abb[110]*m_z2k - ZVertex2Abb[111]*m_z2k_2 + m_z2k_3;

  ZVertex2Abb[113]=ZVertex2Abb[112]*m_cR - 6.*ZVertex2Abb[1]*ZVertex2Abb[61]*m_cL*m_z2k;

  ZVertex2Abb[126]=(-2.*m_m)/m_s_2;

  ZVertex2Abb[127]=(-2.*m_m*m_x)/(m_s_2*m_z2k);

  ZVertex2Abb[128]=(-4.*m_cR*m_m*m_x)/m_s_2;

  ZVertex2Abb[129]=(-4.*m_m)/m_s_2;

  ZVertex2Abb[130]=(-4.*m_m)/m_s;

  ZVertex2Abb[131]=(2.*m_m)/m_s_2;

  ZVertex2Abb[132]=(-2.*m_m)/(m_s_2*m_z2k);

  ZVertex2Abb[133]=-ZVertex2Abb[61]*m_cR + 2.*m_cL*m_x;

  ZVertex2Abb[134]=m_x + 3.*m_z2k;

  ZVertex2Abb[135]=ZVertex2Abb[61]*m_cL + m_cR - ZVertex2Abb[134]*m_cR;

  ZVertex2Abb[136]=ZVertex2Abb[61]*m_cL - 2.*m_cR*m_z2k;

  ZVertex2Abb[137]=9. - 11.*m_z2k;

  ZVertex2Abb[138]=3. + 7.*m_z2k;

  ZVertex2Abb[139]=-pow(ZVertex2Abb[2],3.) + ZVertex2Abb[138]*ZVertex2Abb[2]*m_x + ZVertex2Abb[137]*m_x_2 + 5.*m_x_3;

  ZVertex2Abb[140]=ZVertex2Abb[139]*m_cL - ZVertex2Abb[89]*m_cR;

  ZVertex2Abb[141]=-2. + 5.*m_x;

  ZVertex2Abb[142]=ZVertex2Abb[141]*m_cL + m_cR - m_cR*m_x;

  ZVertex2Abb[143]=ZVertex2Abb[97]*m_cL + ZVertex2Abb[95]*m_cR;

  ZVertex2Abb[144]=-1. + 5.*m_x;

  ZVertex2Abb[145]=ZVertex2Abb[99]*m_cL + ZVertex2Abb[144]*m_cR;

  ZVertex2Abb[146]=m_cL + 3.*m_cR;

  ZVertex2Abb[147]=ZVertex2Abb[142]*ZVertex2Abb[93] - ZVertex2Abb[143]*m_z2k + ZVertex2Abb[145]*m_z2k_2 - ZVertex2Abb[146]*m_z2k_3;

  ZVertex2Abb[148]=ZVertex2Abb[61]*ZVertex2Abb[69]*m_cR + 6.*ZVertex2Abb[67]*m_cL*m_x_2;

  ZVertex2Abb[149]=ZVertex2Abb[74]*m_cR + 2.*ZVertex2Abb[71]*m_cL*m_x;

  ZVertex2Abb[150]=ZVertex2Abb[1]*ZVertex2Abb[77]*m_cR + ZVertex2Abb[79]*m_cL*m_x;

  ZVertex2Abb[151]=-ZVertex2Abb[1]*ZVertex2Abb[104]*ZVertex2Abb[61]*m_cL - ZVertex2Abb[107]*m_cR*m_z2k;

  ZVertex2Abb[152]=ZVertex2Abb[112]*m_cL - 6.*ZVertex2Abb[1]*ZVertex2Abb[61]*m_cR*m_z2k;

  ZVertex2Abb[165]=(-4.*m_cL*m_m*m_x)/m_s_2;

  ZVertex2Abb[173]=-4.*m_cL;

  ZVertex2Abb[174]=2./m_s;

  ZVertex2Abb[175]=-4./m_s;

  ZVertex2Abb[183]=-4.*m_cR;

  ZVertex2Abb[184]=-2. + m_z2k;

  ZVertex2Abb[185]=pow(ZVertex2Abb[2],2.) + ZVertex2Abb[184]*m_x - 2.*m_x_2;

  ZVertex2Abb[186]=4. + 3.*m_z2k;

  ZVertex2Abb[187]=pow(ZVertex2Abb[2],3.) + m_x + ZVertex2Abb[186]*m_x_2 - 4.*m_x_3 - m_x*m_z2k;

  ZVertex2Abb[188]=-pow(ZVertex2Abb[2],2.) + 2.*ZVertex2Abb[184]*m_x - m_x_2;

  ZVertex2Abb[189]=4. - 5.*m_z2k;

  ZVertex2Abb[190]=4. + m_z2k;

  ZVertex2Abb[191]=2. + m_z2k;

  ZVertex2Abb[192]=-11. + ZVertex2Abb[191]*m_z2k;

  ZVertex2Abb[193]=-ZVertex2Abb[190]*pow(ZVertex2Abb[2],2.) - ZVertex2Abb[192]*m_x - ZVertex2Abb[189]*m_x_2 - 3.*m_x_3;

  ZVertex2Abb[194]=-pow(ZVertex2Abb[2],2.) - ZVertex2Abb[184]*m_x + 2.*m_x_2;

  ZVertex2Abb[195]=-7. + m_z2k - 3.*m_z2k_2;

  ZVertex2Abb[196]=-5. + m_z2k + 3.*m_z2k_2;

  ZVertex2Abb[197]=pow(ZVertex2Abb[2],3.) + ZVertex2Abb[196]*ZVertex2Abb[2]*m_x + ZVertex2Abb[195]*m_x_2 + 3.*m_x_4 - 3.*m_x_3*m_z2k;

  ZVertex2Abb[198]=1. - m_x - m_z2k;

  ZVertex2Abb[199]=-4. + 5.*m_z2k;

  ZVertex2Abb[200]=2. + 7.*m_z2k;

  ZVertex2Abb[201]=m_x - ZVertex2Abb[200]*m_x_2 + m_x_3 + pow(ZVertex2Abb[2],2.)*m_z2k + ZVertex2Abb[199]*m_x*m_z2k;

  ZVertex2Abb[202]=1. + 3.*m_x;

  ZVertex2Abb[203]=4. - 3.*m_x;

  ZVertex2Abb[204]=pow(ZVertex2Abb[10],2.)*ZVertex2Abb[11] + ZVertex2Abb[202]*ZVertex2Abb[68]*m_z2k + ZVertex2Abb[203]*m_z2k_2 - 3.*m_z2k_3;

  ZVertex2Abb[205]=pow(ZVertex2Abb[10],2.)*m_x + m_z2k - ZVertex2Abb[99]*m_x*m_z2k + ZVertex2Abb[141]*m_z2k_2 + m_z2k_3;

  ZVertex2Abb[220]=(4.*m_cL*m_x)/m_s;

  ZVertex2Abb[221]=(2.*m_cL)/m_s;

  ZVertex2Abb[222]=(-2.*m_cL*m_x)/(m_s*m_z2k);

  ZVertex2Abb[223]=(-4.*m_cL*m_x)/m_s;

  ZVertex2Abb[224]=(-8.*m_cL*m_x)/m_s;

  ZVertex2Abb[225]=(-4.*m_cL)/m_s;

  ZVertex2Abb[226]=8.*m_cL*m_x;

  ZVertex2Abb[227]=4.*m_cL;

  ZVertex2Abb[228]=(-2.*m_cL)/(m_s*m_z2k);

  ZVertex2Abb[229]=(-2.*m_cL)/m_s;

  ZVertex2Abb[244]=(4.*m_cR*m_x)/m_s;

  ZVertex2Abb[245]=(2.*m_cR)/m_s;

  ZVertex2Abb[246]=(-2.*m_cR*m_x)/(m_s*m_z2k);

  ZVertex2Abb[247]=(-4.*m_cR*m_x)/m_s;

  ZVertex2Abb[248]=(-8.*m_cR*m_x)/m_s;

  ZVertex2Abb[249]=(-4.*m_cR)/m_s;

  ZVertex2Abb[250]=8.*m_cR*m_x;

  ZVertex2Abb[251]=4.*m_cR;

  ZVertex2Abb[252]=(-2.*m_cR)/(m_s*m_z2k);

  ZVertex2Abb[253]=(-2.*m_cR)/m_s;

  ZVertex2Abb[266]=(-4.*m_m*m_x)/(m_s_2*m_z2k);

  ZVertex2Abb[267]=(-8.*m_cR*m_m*m_x)/m_s_2;

  ZVertex2Abb[268]=(-8.*m_m)/m_s_2;

  ZVertex2Abb[269]=(-8.*m_m)/m_s;

  ZVertex2Abb[270]=(4.*m_m)/m_s_2;

  ZVertex2Abb[271]=(-4.*m_m)/(m_s_2*m_z2k);

  ZVertex2Abb[284]=(-8.*m_cL*m_m*m_x)/m_s_2;
}

// Set up coefficients multiplying scalar master integrals in corrections to photon vertex for emission off leg 2
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Vertices::Init_P_Vertex_2_Coefficients() 
{
  PVertex2Abb[0]=m_s*m_x - m_s*m_z2k;

  PVertex2Abb[1]=m_x - m_z2k;

  PVertex2Abb[2]=2.*m_x + m_z2k;

  PVertex2Abb[3]=5.*m_x + m_z2k;

  PVertex2Abb[8]=-2.*m_cL;

  PVertex2Abb[9]=m_s*m_x;

  PVertex2Abb[10]=(-2.*m_cL)/m_s;

  PVertex2Abb[11]=m_s*m_z2k;

  PVertex2Abb[12]=m_cL/m_s;

  PVertex2Abb[13]=-2.*m_cL*m_x;

  PVertex2Abb[18]=-2.*m_cR;

  PVertex2Abb[19]=(-2.*m_cR)/m_s;

  PVertex2Abb[20]=m_cR/m_s;

  PVertex2Abb[21]=-2.*m_cR*m_x;

  PVertex2Abb[23]=-2.*m_cR*m_m;

  PVertex2Abb[24]=2.*m_cR*m_m;

  PVertex2Abb[26]=-2.*m_cL*m_m;

  PVertex2Abb[27]=2.*m_cL*m_m;

  PVertex2Abb[28]=-3.*m_x + m_z2k;

  PVertex2Abb[29]=-m_x + m_z2k;

  PVertex2Abb[30]=m_x + m_z2k;

  PVertex2Abb[31]=3.*m_x - m_z2k;

  PVertex2Abb[37]=(2.*m_cR*m_m)/m_s_2;

  PVertex2Abb[38]=(2.*m_cR*m_m*m_x)/(m_s_2*m_z2k);

  PVertex2Abb[39]=(2.*m_cR*m_m)/(m_s_2*m_z2k);

  PVertex2Abb[45]=(2.*m_cL*m_m)/m_s_2;

  PVertex2Abb[46]=(2.*m_cL*m_m*m_x)/(m_s_2*m_z2k);

  PVertex2Abb[47]=(2.*m_cL*m_m)/(m_s_2*m_z2k);

  PVertex2Abb[48]=m_s*m_x*m_z2k - m_s*m_z2k_2;

  PVertex2Abb[49]=m_x_2 - 6.*m_x*m_z2k + m_z2k_2;

  PVertex2Abb[54]=2.*m_cL*m_x;

  PVertex2Abb[55]=(-8.*m_cL*m_x)/m_s;

  PVertex2Abb[56]=(-2.*m_cL)/(m_s*m_z2k);

  PVertex2Abb[61]=2.*m_cR*m_x;

  PVertex2Abb[62]=(-8.*m_cR*m_x)/m_s;

  PVertex2Abb[63]=(-2.*m_cR)/(m_s*m_z2k);
}


DivArrC Z_Decay_RV_Vertices::RV_Z_Vertex_1(const int& ME, const int& LR,
					   const Vec4C& epsV, const Vec4C& epsP)
{
  // Z Vertex corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,

  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*(epsP*m_p1)*((ZVertex1Abb[263]*ZVertex1Abb[57]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[263]*ZVertex1Abb[57]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[263]*ZVertex1Abb[59]*(epsV*m_p2))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[262]*ZVertex1Abb[56]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[262]*ZVertex1Abb[56]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[264]*ZVertex1Abb[65]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[264]*ZVertex1Abb[65]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[73]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[129]*ZVertex1Abb[81]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[81]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[94]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[107]*ZVertex1Abb[266]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[107]*ZVertex1Abb[266]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[113]*ZVertex1Abb[129]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[265]*ZVertex1Abb[96]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[265]*ZVertex1Abb[96]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[101]*ZVertex1Abb[265]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
    }
    else if (LR == 1) {
      return One*(epsP*m_p1)*((ZVertex1Abb[133]*ZVertex1Abb[263]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[133]*ZVertex1Abb[263]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[135]*ZVertex1Abb[263]*(epsV*m_p2))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[280]*ZVertex1Abb[56]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[280]*ZVertex1Abb[56]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[136]*ZVertex1Abb[264]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[136]*ZVertex1Abb[264]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[137]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[144]*ZVertex1Abb[263]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[144]*ZVertex1Abb[263]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[150]*ZVertex1Abb[263]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[151]*ZVertex1Abb[266]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[151]*ZVertex1Abb[266]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[152]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p1)*((ZVertex1Abb[138]*ZVertex1Abb[265]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[138]*ZVertex1Abb[265]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[139]*ZVertex1Abb[265]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
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
      return One
	*((ZVertex1Abb[222]*(epsV*m_pP))/ZVertex1Abb[6] + (ZVertex1Abb[222]*(epsV*m_p1))/ZVertex1Abb[6] + (ZVertex1Abb[223]*ZVertex1Abb[41]*(epsV*m_p2))/ZVertex1Abb[6])

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex1Abb[221]*ZVertex1Abb[56]*(epsV*m_pP))/ZVertex1Abb[6] + (ZVertex1Abb[221]*ZVertex1Abb[56]*(epsV*m_p1))/ZVertex1Abb[6])

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex1Abb[189]*ZVertex1Abb[223]*ZVertex1Abb[41]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[189]*ZVertex1Abb[223]*ZVertex1Abb[41]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[194]*ZVertex1Abb[223]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex1Abb[195]*ZVertex1Abb[224]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[195]*ZVertex1Abb[224]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[198]*ZVertex1Abb[225]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((ZVertex1Abb[202]*ZVertex1Abb[205]*ZVertex1Abb[227]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[202]*ZVertex1Abb[205]*ZVertex1Abb[227]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[208]*ZVertex1Abb[228]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((ZVertex1Abb[195]*ZVertex1Abb[226]*ZVertex1Abb[3]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[195]*ZVertex1Abb[226]*ZVertex1Abb[3]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[176]*ZVertex1Abb[201]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
    }
    else if (LR == 1) {
      return One
	*((ZVertex1Abb[242]*(epsV*m_pP))/ZVertex1Abb[6] + (ZVertex1Abb[242]*(epsV*m_p1))/ZVertex1Abb[6] + (ZVertex1Abb[243]*ZVertex1Abb[41]*(epsV*m_p2))/ZVertex1Abb[6])

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex1Abb[241]*ZVertex1Abb[56]*(epsV*m_pP))/ZVertex1Abb[6] + (ZVertex1Abb[241]*ZVertex1Abb[56]*(epsV*m_p1))/ZVertex1Abb[6])

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex1Abb[189]*ZVertex1Abb[243]*ZVertex1Abb[41]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[189]*ZVertex1Abb[243]*ZVertex1Abb[41]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[194]*ZVertex1Abb[243]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex1Abb[195]*ZVertex1Abb[244]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[195]*ZVertex1Abb[244]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[198]*ZVertex1Abb[245]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((ZVertex1Abb[202]*ZVertex1Abb[205]*ZVertex1Abb[247]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[202]*ZVertex1Abb[205]*ZVertex1Abb[247]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[208]*ZVertex1Abb[248]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((ZVertex1Abb[195]*ZVertex1Abb[246]*ZVertex1Abb[3]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[195]*ZVertex1Abb[246]*ZVertex1Abb[3]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[187]*ZVertex1Abb[201]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
    }
  }
  // ubar1 \slashed{m_epsV} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(ZVertex1Abb[176]*(epsP*m_p1))/ZVertex1Abb[0]

	+(ZVertex1Abb[168]*ZVertex1Abb[177]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[178]*ZVertex1Abb[9]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[13]*ZVertex1Abb[27]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(4.*ZVertex1Abb[16]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]);
    }
    else if (LR == 1) {
      return One*(ZVertex1Abb[187]*(epsP*m_p1))/ZVertex1Abb[0]

	+(ZVertex1Abb[177]*ZVertex1Abb[179]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[178]*ZVertex1Abb[30]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[27]*ZVertex1Abb[31]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(4.*ZVertex1Abb[32]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 5) {
    if (LR == 0) {
      return One
	*((ZVertex1Abb[128]*ZVertex1Abb[57]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[128]*ZVertex1Abb[57]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[128]*ZVertex1Abb[59]*(epsV*m_p2))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex1Abb[127]*ZVertex1Abb[56]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[127]*ZVertex1Abb[56]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex1Abb[129]*ZVertex1Abb[65]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[65]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[130]*ZVertex1Abb[73]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex1Abb[130]*ZVertex1Abb[81]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[130]*ZVertex1Abb[81]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[130]*ZVertex1Abb[94]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((ZVertex1Abb[107]*ZVertex1Abb[132]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[107]*ZVertex1Abb[132]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[113]*ZVertex1Abb[130]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((ZVertex1Abb[131]*ZVertex1Abb[96]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[131]*ZVertex1Abb[96]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[101]*ZVertex1Abb[131]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
    }
    else if (LR == 1) {
      return One
	*((ZVertex1Abb[128]*ZVertex1Abb[133]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[128]*ZVertex1Abb[133]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[128]*ZVertex1Abb[135]*(epsV*m_p2))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex1Abb[166]*ZVertex1Abb[56]*(epsV*m_pP))/(ZVertex1Abb[3]*ZVertex1Abb[6]) + (ZVertex1Abb[166]*ZVertex1Abb[56]*(epsV*m_p1))/(ZVertex1Abb[3]*ZVertex1Abb[6]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex1Abb[129]*ZVertex1Abb[136]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[129]*ZVertex1Abb[136]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[130]*ZVertex1Abb[137]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex1Abb[128]*ZVertex1Abb[144]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[128]*ZVertex1Abb[144]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[128]*ZVertex1Abb[150]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((ZVertex1Abb[132]*ZVertex1Abb[151]*(epsV*m_pP))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[132]*ZVertex1Abb[151]*(epsV*m_p1))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)) + (ZVertex1Abb[130]*ZVertex1Abb[152]*(epsV*m_p2))/(ZVertex1Abb[3]*pow(ZVertex1Abb[6],2.)))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((ZVertex1Abb[131]*ZVertex1Abb[138]*(epsV*m_pP))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[131]*ZVertex1Abb[138]*(epsV*m_p1))/pow(ZVertex1Abb[6],2.) + (ZVertex1Abb[131]*ZVertex1Abb[139]*(epsV*m_p2))/pow(ZVertex1Abb[6],2.));
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
      return (ZVertex1Abb[46]*B_0(m_s,m_m2,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[47]*ZVertex1Abb[8]*B_0(m_m2,0.,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[10]*ZVertex1Abb[47]*B_0(m_s1k,0.,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[41]*ZVertex1Abb[48]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/ZVertex1Abb[6];
    }
    else if (LR == 1) {
      return (ZVertex1Abb[53]*B_0(m_s,m_m2,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[54]*ZVertex1Abb[8]*B_0(m_m2,0.,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[10]*ZVertex1Abb[54]*B_0(m_s1k,0.,m_m2,m_mu2))/ZVertex1Abb[6]

	+(ZVertex1Abb[41]*ZVertex1Abb[55]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/ZVertex1Abb[6];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*ZVertex1Abb[24]/ZVertex1Abb[0]

	+(ZVertex1Abb[2]*ZVertex1Abb[26]*B_0(m_s,m_m2,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[27]*ZVertex1Abb[9]*B_0(m_m2,0.,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[13]*ZVertex1Abb[26]*B_0(m_s1k,0.,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(2.*ZVertex1Abb[16]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6]);
    }
    else if (LR == 1) {
      return One*ZVertex1Abb[40]/ZVertex1Abb[0]

	+(ZVertex1Abb[26]*ZVertex1Abb[29]*B_0(m_s,m_m2,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[27]*ZVertex1Abb[30]*B_0(m_m2,0.,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(ZVertex1Abb[26]*ZVertex1Abb[31]*B_0(m_s1k,0.,m_m2,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6])

	+(2.*ZVertex1Abb[32]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/(ZVertex1Abb[3]*ZVertex1Abb[6]);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}

DivArrC Z_Decay_RV_Vertices::RV_P_Vertex_1(const int& ME, const int& LR,
					   const Vec4C& epsV, const Vec4C& epsP)
{
  // Photon Vertex corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
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
      return One*(PVertex1Abb[30]*PVertex1Abb[52]*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[30]*PVertex1Abb[53]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[54]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[48]*PVertex1Abb[55]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.);
    }
    else if (LR == 1) {
      return One*(PVertex1Abb[30]*PVertex1Abb[59]*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[30]*PVertex1Abb[60]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[61]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[48]*PVertex1Abb[62]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],2.);
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
      return One*(PVertex1Abb[28]*PVertex1Abb[37]*(epsP*m_p1))/pow(PVertex1Abb[1],3.)

	+(PVertex1Abb[30]*PVertex1Abb[38]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[29],3.)

	+(PVertex1Abb[31]*PVertex1Abb[37]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],3.)

	+(PVertex1Abb[39]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/PVertex1Abb[1];
    }
    else if (LR == 1) {
      return One*(PVertex1Abb[28]*PVertex1Abb[45]*(epsP*m_p1))/pow(PVertex1Abb[1],3.)

	+(PVertex1Abb[30]*PVertex1Abb[46]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[29],3.)

	+(PVertex1Abb[31]*PVertex1Abb[45]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p1))/pow(PVertex1Abb[1],3.)

	+(PVertex1Abb[47]*B_0(m_s1k,0.,m_m2,m_mu2)*(epsP*m_p1))/PVertex1Abb[1];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} P_i v2
  else if (ME == 7) {
    if (LR == 0) {
      return (PVertex1Abb[23]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex1Abb[0]

	+(PVertex1Abb[24]*B_0(m_s1k,0.,m_m2,m_mu2))/PVertex1Abb[0];
    }
    else if (LR == 1) {
      return (PVertex1Abb[26]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex1Abb[0]

	+(PVertex1Abb[27]*B_0(m_s1k,0.,m_m2,m_mu2))/PVertex1Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*PVertex1Abb[8]/PVertex1Abb[0]

	+(PVertex1Abb[10]*PVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[12]*PVertex1Abb[3]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[13]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/PVertex1Abb[1];
    }
    else if (LR == 1) {
      return One*PVertex1Abb[18]/PVertex1Abb[0]

	+(PVertex1Abb[19]*PVertex1Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[20]*PVertex1Abb[3]*B_0(m_s1k,0.,m_m2,m_mu2))/pow(PVertex1Abb[1],2.)

	+(PVertex1Abb[21]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/PVertex1Abb[1];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}

DivArrC Z_Decay_RV_Vertices::RV_Z_Vertex_2(const int& ME, const int& LR,
					   const Vec4C& epsV, const Vec4C& epsP)
{
  // Z Vertex corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,

  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One*((ZVertex2Abb[129]*ZVertex2Abb[62]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + ((ZVertex2Abb[129]*ZVertex2Abb[65]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[129]*ZVertex2Abb[62]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))*(epsP*m_p2))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[266]*ZVertex2Abb[66]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + ((ZVertex2Abb[267]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[266]*ZVertex2Abb[66]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))*(epsP*m_p2))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[268]*ZVertex2Abb[70]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[129]*ZVertex2Abb[75]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[268]*ZVertex2Abb[70]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[269]*ZVertex2Abb[70]*(epsV*m_pP)*(epsP*m_p2))/pow(ZVertex2Abb[4],2.) + ((ZVertex2Abb[269]*ZVertex2Abb[80]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[269]*ZVertex2Abb[70]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))*(epsP*m_p2))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[270]*ZVertex2Abb[90]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[102]*ZVertex2Abb[270]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[270]*ZVertex2Abb[90]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[108]*ZVertex2Abb[271]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[113]*ZVertex2Abb[129]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[108]*ZVertex2Abb[271]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2));
    }
    else if (LR == 1) {
      return One*((ZVertex2Abb[133]*ZVertex2Abb[270]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + ((ZVertex2Abb[135]*ZVertex2Abb[270]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[133]*ZVertex2Abb[270]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))*(epsP*m_p2))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[136]*ZVertex2Abb[266]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + ((ZVertex2Abb[284]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[136]*ZVertex2Abb[266]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))*(epsP*m_p2))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[129]*ZVertex2Abb[140]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[129]*ZVertex2Abb[147]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[129]*ZVertex2Abb[140]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[148]*ZVertex2Abb[268]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[129]*ZVertex2Abb[149]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[148]*ZVertex2Abb[268]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[148]*ZVertex2Abb[269]*(epsV*m_pP)*(epsP*m_p2))/pow(ZVertex2Abb[4],2.) + ((ZVertex2Abb[150]*ZVertex2Abb[269]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[148]*ZVertex2Abb[269]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))*(epsP*m_p2))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[151]*ZVertex2Abb[271]*(epsV*m_pP)*(epsP*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + ((ZVertex2Abb[129]*ZVertex2Abb[152]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[151]*ZVertex2Abb[271]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))*(epsP*m_p2));
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
      return One*((ZVertex2Abb[220]*(epsV*m_pP))/ZVertex2Abb[4] + (ZVertex2Abb[221]*ZVertex2Abb[61]*(epsV*m_p1))/ZVertex2Abb[4] + (ZVertex2Abb[220]*(epsV*m_p2))/ZVertex2Abb[4])

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[222]*ZVertex2Abb[61]*(epsV*m_pP))/ZVertex2Abb[4] + (ZVertex2Abb[223]*(epsV*m_p1))/ZVertex2Abb[4] + (ZVertex2Abb[222]*ZVertex2Abb[61]*(epsV*m_p2))/ZVertex2Abb[4])

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[185]*ZVertex2Abb[224]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[187]*ZVertex2Abb[225]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[185]*ZVertex2Abb[224]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[188]*ZVertex2Abb[221]*ZVertex2Abb[46]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[193]*ZVertex2Abb[221]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[188]*ZVertex2Abb[221]*ZVertex2Abb[46]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[1]*ZVertex2Abb[194]*ZVertex2Abb[226]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[197]*ZVertex2Abb[227]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[1]*ZVertex2Abb[194]*ZVertex2Abb[226]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[198]*ZVertex2Abb[201]*ZVertex2Abb[228]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[204]*ZVertex2Abb[229]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[198]*ZVertex2Abb[205]*ZVertex2Abb[228]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.));
    }
    else if (LR == 1) {
      return One*((ZVertex2Abb[244]*(epsV*m_pP))/ZVertex2Abb[4] + (ZVertex2Abb[245]*ZVertex2Abb[61]*(epsV*m_p1))/ZVertex2Abb[4] + (ZVertex2Abb[244]*(epsV*m_p2))/ZVertex2Abb[4])

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[246]*ZVertex2Abb[61]*(epsV*m_pP))/ZVertex2Abb[4] + (ZVertex2Abb[247]*(epsV*m_p1))/ZVertex2Abb[4] + (ZVertex2Abb[246]*ZVertex2Abb[61]*(epsV*m_p2))/ZVertex2Abb[4])

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[185]*ZVertex2Abb[248]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[187]*ZVertex2Abb[249]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[185]*ZVertex2Abb[248]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[188]*ZVertex2Abb[245]*ZVertex2Abb[46]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[193]*ZVertex2Abb[245]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[188]*ZVertex2Abb[245]*ZVertex2Abb[46]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[1]*ZVertex2Abb[194]*ZVertex2Abb[250]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[197]*ZVertex2Abb[251]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[1]*ZVertex2Abb[194]*ZVertex2Abb[250]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[198]*ZVertex2Abb[201]*ZVertex2Abb[252]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[204]*ZVertex2Abb[253]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[198]*ZVertex2Abb[205]*ZVertex2Abb[252]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.));
    }
  }
  // ubar1 \slashed{m_epsV} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return One*(ZVertex2Abb[173]*(epsP*m_p2))/ZVertex2Abb[0]
	
	+(ZVertex2Abb[174]*ZVertex2Abb[6]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[175]*ZVertex2Abb[9]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[15]*ZVertex2Abb[174]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(-4.*ZVertex2Abb[18]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]);
    }
    else if (LR == 1) {
      return One*(ZVertex2Abb[183]*(epsP*m_p2))/ZVertex2Abb[0]

	+(ZVertex2Abb[174]*ZVertex2Abb[31]*B_0(m_s,m_m2,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[175]*ZVertex2Abb[32]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[174]*ZVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(-4.*ZVertex2Abb[36]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*(epsP*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]);
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 5) {
    if (LR == 0) {
      return One
	*((ZVertex2Abb[126]*ZVertex2Abb[62]*(epsV*m_pP))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[126]*ZVertex2Abb[65]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[126]*ZVertex2Abb[62]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[127]*ZVertex2Abb[66]*(epsV*m_pP))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[128]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[127]*ZVertex2Abb[66]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[129]*ZVertex2Abb[70]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[126]*ZVertex2Abb[75]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[129]*ZVertex2Abb[70]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[130]*ZVertex2Abb[70]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[130]*ZVertex2Abb[80]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[130]*ZVertex2Abb[70]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[131]*ZVertex2Abb[90]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[102]*ZVertex2Abb[131]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[131]*ZVertex2Abb[90]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[108]*ZVertex2Abb[132]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[113]*ZVertex2Abb[126]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[108]*ZVertex2Abb[132]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)));
    }
    else if (LR == 1) {
      return One
	*((ZVertex2Abb[131]*ZVertex2Abb[133]*(epsV*m_pP))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[131]*ZVertex2Abb[135]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[131]*ZVertex2Abb[133]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))

	+B_0(0.,0.,m_m2,m_mu2)*((ZVertex2Abb[127]*ZVertex2Abb[136]*(epsV*m_pP))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[165]*(epsV*m_p1))/(ZVertex2Abb[1]*ZVertex2Abb[4]) + (ZVertex2Abb[127]*ZVertex2Abb[136]*(epsV*m_p2))/(ZVertex2Abb[1]*ZVertex2Abb[4]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((ZVertex2Abb[126]*ZVertex2Abb[140]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[126]*ZVertex2Abb[147]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[126]*ZVertex2Abb[140]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))

	+B_0(m_m2,0.,m_m2,m_mu2)*((ZVertex2Abb[129]*ZVertex2Abb[148]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[126]*ZVertex2Abb[149]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[129]*ZVertex2Abb[148]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((ZVertex2Abb[130]*ZVertex2Abb[148]*(epsV*m_pP))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[130]*ZVertex2Abb[150]*(epsV*m_p1))/pow(ZVertex2Abb[4],2.) + (ZVertex2Abb[130]*ZVertex2Abb[148]*(epsV*m_p2))/pow(ZVertex2Abb[4],2.))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((ZVertex2Abb[132]*ZVertex2Abb[151]*(epsV*m_pP))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[126]*ZVertex2Abb[152]*(epsV*m_p1))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)) + (ZVertex2Abb[132]*ZVertex2Abb[151]*(epsV*m_p2))/(ZVertex2Abb[1]*pow(ZVertex2Abb[4],2.)));
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
      return (ZVertex2Abb[51]*B_0(m_s,m_m2,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[52]*ZVertex2Abb[8]*B_0(m_m2,0.,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[45]*ZVertex2Abb[52]*B_0(m_s2k,0.,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[46]*ZVertex2Abb[53]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/ZVertex2Abb[4];
    }
    else if (LR == 1) {
      return (ZVertex2Abb[58]*B_0(m_s,m_m2,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[59]*ZVertex2Abb[8]*B_0(m_m2,0.,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[45]*ZVertex2Abb[59]*B_0(m_s2k,0.,m_m2,m_mu2))/ZVertex2Abb[4]

	+(ZVertex2Abb[46]*ZVertex2Abb[60]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/ZVertex2Abb[4];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*ZVertex2Abb[26]/ZVertex2Abb[0]

	+(ZVertex2Abb[28]*ZVertex2Abb[6]*B_0(m_s,m_m2,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[29]*ZVertex2Abb[9]*B_0(m_m2,0.,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[15]*ZVertex2Abb[28]*B_0(m_s2k,0.,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(-2.*ZVertex2Abb[18]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4]);
    }
    else if (LR == 1) {
      return One*ZVertex2Abb[44]/ZVertex2Abb[0]

	+(ZVertex2Abb[28]*ZVertex2Abb[31]*B_0(m_s,m_m2,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[29]*ZVertex2Abb[32]*B_0(m_m2,0.,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(ZVertex2Abb[28]*ZVertex2Abb[35]*B_0(m_s2k,0.,m_m2,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4])

	+(-2.*ZVertex2Abb[36]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/(ZVertex2Abb[1]*ZVertex2Abb[4]);
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}


DivArrC Z_Decay_RV_Vertices::RV_P_Vertex_2(const int& ME, const int& LR,
					   const Vec4C& epsV, const Vec4C& epsP)
{
  // Photon Vertex corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
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
      return One*(PVertex2Abb[8]*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[54]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[48]

	+(PVertex2Abb[55]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],2.)

	+(PVertex2Abb[49]*PVertex2Abb[56]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],2.);
    }
    else if (LR == 1) {
      return One*(PVertex2Abb[18]*(epsP*m_p2))/PVertex2Abb[0]

	+(PVertex2Abb[61]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[48]

	+(PVertex2Abb[62]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],2.)

	+(PVertex2Abb[49]*PVertex2Abb[63]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],2.);
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
      return One*(PVertex2Abb[28]*PVertex2Abb[37]*(epsP*m_p2))/pow(PVertex2Abb[1],3.)

	+(PVertex2Abb[30]*PVertex2Abb[38]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[29],3.)

	+(PVertex2Abb[31]*PVertex2Abb[37]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],3.)

	+(PVertex2Abb[39]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[1];
    }
    else if (LR == 1) {
      return One*(PVertex2Abb[28]*PVertex2Abb[45]*(epsP*m_p2))/pow(PVertex2Abb[1],3.)

	+(PVertex2Abb[30]*PVertex2Abb[46]*B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[29],3.)

	+(PVertex2Abb[31]*PVertex2Abb[45]*B_0(m_m2,0.,m_m2,m_mu2)*(epsP*m_p2))/pow(PVertex2Abb[1],3.)

	+(PVertex2Abb[47]*B_0(m_s2k,0.,m_m2,m_mu2)*(epsP*m_p2))/PVertex2Abb[1];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} P_i v2
  else if (ME == 7) {
    if (LR == 0) {
      return (PVertex2Abb[23]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex2Abb[0]
	
	+(PVertex2Abb[24]*B_0(m_s2k,0.,m_m2,m_mu2))/PVertex2Abb[0];
    }
    else if (LR == 1) {
      return (PVertex2Abb[26]*B_0(m_m2,0.,m_m2,m_mu2))/PVertex2Abb[0]

	+(PVertex2Abb[27]*B_0(m_s2k,0.,m_m2,m_mu2))/PVertex2Abb[0];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{m_epsV} \slashed{k} P_i v2
  else if (ME == 8) {
    if (LR == 0) {
      return One*PVertex2Abb[8]/PVertex2Abb[0]

	+(PVertex2Abb[10]*PVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex2Abb[1],2.)

	+(PVertex2Abb[12]*PVertex2Abb[3]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(PVertex2Abb[1],2.)

	+(PVertex2Abb[13]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/PVertex2Abb[1];
    }
    else if (LR == 1) {
      return One*PVertex2Abb[18]/PVertex2Abb[0]

	+(PVertex2Abb[19]*PVertex2Abb[2]*B_0(m_m2,0.,m_m2,m_mu2))/pow(PVertex2Abb[1],2.)

	+(PVertex2Abb[20]*PVertex2Abb[3]*B_0(m_s2k,0.,m_m2,m_mu2))/pow(PVertex2Abb[1],2.)
	
	+(PVertex2Abb[21]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/PVertex2Abb[1];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	      << "Values range from 1 to 8.";
  }
  return Zero;
}

