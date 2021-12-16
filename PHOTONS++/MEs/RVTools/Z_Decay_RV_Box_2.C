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

Z_Decay_RV_Box_2::Z_Decay_RV_Box_2
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Z_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Box2Abb = new Complex[6327];
  Init_Coefficients();
}

Z_Decay_RV_Box_2::~Z_Decay_RV_Box_2()
{
  delete [] Box2Abb;
}

void Z_Decay_RV_Box_2::Init_Coefficients() 
{
  Init_Box_2_Coefficients();
  return;
}

// Set up coefficients multiplying master integrals for emission from leg 2
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Box_2::Init_Box_2_Coefficients() 
{
  Box2Abb[0]=1. - m_z12;

  Box2Abb[1]=-1. + m_x;

  Box2Abb[2]=-1. - 2.*m_x + m_z12;

  Box2Abb[3]=m_x + Box2Abb[1]*m_x*m_z12 + Box2Abb[2]*m_z12*m_z2k + m_z12*m_z2k_2;

  Box2Abb[4]=-1. + m_z12;

  Box2Abb[5]=-1. + m_z12 + m_z2k;

  Box2Abb[6]=m_x + Box2Abb[5]*m_z12 - 2.*m_x*m_z12;

  Box2Abb[7]=Box2Abb[6]*m_cL - Box2Abb[4]*m_cR*m_x;

  Box2Abb[8]=m_x - m_z2k;

  Box2Abb[9]=1. + m_z12;

  Box2Abb[10]=1. + m_x - m_z12 + m_x*m_z12 - Box2Abb[9]*m_z2k;

  Box2Abb[11]=1. + 2.*m_z12;

  Box2Abb[12]=1. + m_z2k;

  Box2Abb[13]=3. + m_z12 + 3.*m_z2k;

  Box2Abb[14]=-3. + Box2Abb[13]*m_z12 + m_z2k;

  Box2Abb[15]=-Box2Abb[14]*m_x + Box2Abb[11]*m_x_2 + Box2Abb[12]*Box2Abb[5]*m_z12;

  Box2Abb[16]=Box2Abb[15]*m_cL + Box2Abb[10]*m_cR*m_x;

  Box2Abb[17]=1. + 2.*m_x - m_z12 - 2.*m_z2k;

  Box2Abb[18]=3. - 4.*m_z12 - 2.*m_z2k;

  Box2Abb[19]=Box2Abb[18]*m_x + 2.*m_x_2 + Box2Abb[5]*m_z12;

  Box2Abb[20]=Box2Abb[19]*m_cL + Box2Abb[17]*m_cR*m_x;

  Box2Abb[21]=2. - 3.*Box2Abb[12]*m_z12 + m_z12_2 + m_z2k - 2.*m_z2k_2;

  Box2Abb[22]=7. - 2.*m_z12 + 4.*m_z2k;

  Box2Abb[23]=-4. + Box2Abb[22]*m_z12;

  Box2Abb[24]=Box2Abb[23]*m_x_2 + Box2Abb[21]*m_x*m_z12 - 2.*m_x_3*m_z12 + Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[25]=Box2Abb[24]*m_cL + Box2Abb[17]*Box2Abb[8]*m_cR*m_x*m_z12;

  Box2Abb[33]=m_s*m_z12;

  Box2Abb[34]=m_s*m_x;

  Box2Abb[35]=m_s*m_z2k;

  Box2Abb[36]=m_z12;

  Box2Abb[37]=m_s;

  Box2Abb[38]=m_x - m_x*m_z12;

  Box2Abb[39]=Box2Abb[38]*m_cL + Box2Abb[6]*m_cR;

  Box2Abb[40]=Box2Abb[15]*m_cR + Box2Abb[10]*m_cL*m_x;

  Box2Abb[41]=Box2Abb[19]*m_cR + Box2Abb[17]*m_cL*m_x;

  Box2Abb[42]=Box2Abb[24]*m_cR + Box2Abb[17]*Box2Abb[8]*m_cL*m_x*m_z12;

  Box2Abb[50]=m_cL + m_cR;

  Box2Abb[51]=4.*m_x - m_z12;

  Box2Abb[52]=1. - 2.*m_z12 - 4.*m_z2k;

  Box2Abb[53]=-3. + 4.*m_z12;

  Box2Abb[54]=pow(Box2Abb[4],2.) + Box2Abb[52]*m_x + 2.*m_x_2 + Box2Abb[53]*m_z2k + 2.*m_z2k_2;

  Box2Abb[55]=-2. + m_z12 - 4.*m_z2k;

  Box2Abb[56]=m_z12 + 2.*m_z2k;

  Box2Abb[57]=Box2Abb[55]*m_x + 2.*m_x_2 + Box2Abb[56]*m_z2k;

  Box2Abb[58]=-Box2Abb[57]*m_cL + Box2Abb[54]*m_cR;

  Box2Abb[59]=1. + 6.*m_z2k;

  Box2Abb[60]=-2. + m_z12_2 + 4.*m_z2k - 2.*m_z12*m_z2k - 6.*m_z2k_2;

  Box2Abb[61]=-1. + m_z2k;

  Box2Abb[62]=-3. + 2.*m_z12;

  Box2Abb[63]=pow(Box2Abb[4],2.) + Box2Abb[62]*m_z2k + 2.*m_z2k_2;

  Box2Abb[64]=Box2Abb[61]*Box2Abb[63] + Box2Abb[60]*m_x + Box2Abb[59]*m_x_2 - 2.*m_x_3;

  Box2Abb[65]=-1. + 2.*m_z12 + m_z2k;

  Box2Abb[66]=4. + 3.*m_z12 + 6.*m_z2k;

  Box2Abb[67]=-2. + 3.*m_z12 + 2.*m_z2k + 8.*m_z12*m_z2k + 6.*m_z2k_2;

  Box2Abb[68]=Box2Abb[67]*m_x - Box2Abb[66]*m_x_2 + 2.*m_x_3 - Box2Abb[56]*Box2Abb[65]*m_z2k;

  Box2Abb[69]=Box2Abb[68]*m_cL + Box2Abb[64]*m_cR;

  Box2Abb[70]=2. + m_z12;

  Box2Abb[71]=2.*Box2Abb[70]*m_x - Box2Abb[56]*m_z12;

  Box2Abb[72]=-2. + m_z12;

  Box2Abb[73]=7. - 4.*m_z12 - 4.*m_z2k;

  Box2Abb[74]=-2. + Box2Abb[73]*m_z12 + 4.*m_z2k;

  Box2Abb[75]=Box2Abb[74]*m_x + 2.*Box2Abb[72]*m_x_2 + Box2Abb[63]*m_z12;

  Box2Abb[76]=Box2Abb[71]*Box2Abb[8]*m_cL + Box2Abb[75]*m_cR;

  Box2Abb[85]=(-2.*m_m)/m_s;

  Box2Abb[86]=(2.*m_m)/m_s;

  Box2Abb[87]=m_m;

  Box2Abb[88]=m_m*m_s;

  Box2Abb[89]=Box2Abb[54]*m_cL - Box2Abb[57]*m_cR;

  Box2Abb[90]=Box2Abb[64]*m_cL + Box2Abb[68]*m_cR;

  Box2Abb[91]=Box2Abb[75]*m_cL + Box2Abb[71]*Box2Abb[8]*m_cR;

  Box2Abb[130]=m_x + Box2Abb[1]*m_z12;

  Box2Abb[131]=1. - 2.*Box2Abb[9]*m_x + Box2Abb[4]*m_z12;

  Box2Abb[132]=Box2Abb[131]*m_cL + m_cR + 2.*Box2Abb[130]*m_cR;

  Box2Abb[133]=-2.*Box2Abb[9]*m_x + m_z12_2;

  Box2Abb[134]=Box2Abb[133]*m_cR + 2.*Box2Abb[9]*m_cL*m_x - m_cL*m_z12;

  Box2Abb[135]=Box2Abb[132]*m_x + Box2Abb[134]*m_z2k;

  Box2Abb[136]=m_x + 2.*Box2Abb[9]*m_x_2 - 2.*Box2Abb[65]*m_x*m_z12 + Box2Abb[5]*m_z12_2 - 2.*m_x*m_z2k;

  Box2Abb[137]=m_z12 + m_z2k;

  Box2Abb[138]=3. + m_z12 + 2.*m_z2k;

  Box2Abb[139]=-3. + Box2Abb[138]*m_z12 + 2.*m_z2k;

  Box2Abb[140]=Box2Abb[139]*m_x - 2.*Box2Abb[9]*m_x_2 + m_z12 - Box2Abb[137]*m_z12;

  Box2Abb[141]=-Box2Abb[136]*m_cL - Box2Abb[140]*m_cR;

  Box2Abb[142]=pow(Box2Abb[61],2.) - 2.*Box2Abb[12]*m_x + m_x_2;

  Box2Abb[143]=-2. + m_x + m_z12;

  Box2Abb[144]=-2.*Box2Abb[4]*m_x - 3.*m_x_2 + Box2Abb[62]*m_z12;

  Box2Abb[145]=1. + 3.*m_x + m_z12;

  Box2Abb[146]=Box2Abb[1]*Box2Abb[143]*m_x + Box2Abb[144]*m_z2k + Box2Abb[145]*m_z2k_2 - m_z2k_3;

  Box2Abb[147]=Box2Abb[146]*m_cL - Box2Abb[142]*Box2Abb[8]*m_cR;

  Box2Abb[148]=1. + m_x - m_z12 - m_z2k;

  Box2Abb[149]=m_z12 - m_z2k;

  Box2Abb[150]=2. + 3.*m_z2k;

  Box2Abb[151]=-1. + m_z12 + m_z12_2 + 3.*m_z2k_2;

  Box2Abb[152]=Box2Abb[149]*Box2Abb[5]*Box2Abb[61] + Box2Abb[151]*m_x - Box2Abb[150]*m_x_2 + m_x_3;

  Box2Abb[153]=Box2Abb[152]*m_cL - Box2Abb[142]*Box2Abb[148]*m_cR;

  Box2Abb[154]=m_z12_2 + 4.*m_z2k + 3.*m_z12*m_z2k;

  Box2Abb[155]=-2. + 3.*m_z12;

  Box2Abb[156]=pow(Box2Abb[4],2.) + Box2Abb[155]*m_z2k + m_z2k_2;

  Box2Abb[157]=-1. + m_z12 + m_z12_2;

  Box2Abb[158]=2. + 3.*m_z12;

  Box2Abb[159]=pow(Box2Abb[4],2.) + 2.*Box2Abb[157]*m_z2k + Box2Abb[158]*m_z2k_2;

  Box2Abb[160]=Box2Abb[159]*m_x - Box2Abb[154]*m_x_2 + Box2Abb[70]*m_x_3 - Box2Abb[156]*m_z12*m_z2k;

  Box2Abb[161]=2.*m_z12 + 2.*m_z2k + 3.*m_z12*m_z2k;

  Box2Abb[162]=-2. + m_z12 + 3.*m_z2k;

  Box2Abb[163]=2. + Box2Abb[162]*m_z12 + 4.*m_z2k;

  Box2Abb[164]=Box2Abb[163]*m_x_2 - Box2Abb[70]*m_x_3 - Box2Abb[161]*m_x*m_z2k + Box2Abb[137]*m_z12*m_z2k_2;

  Box2Abb[165]=Box2Abb[164]*m_cL + Box2Abb[160]*m_cR;

  Box2Abb[166]=-4. + 3.*m_z12;

  Box2Abb[167]=6. - 4.*m_z12 - 3.*m_z2k;

  Box2Abb[168]=2.*Box2Abb[61] + Box2Abb[167]*m_z12;

  Box2Abb[169]=-Box2Abb[72]*m_x_3 + Box2Abb[168]*m_x*m_z2k + Box2Abb[166]*m_x_2*m_z2k + pow(Box2Abb[5],2.)*m_z12*m_z2k;

  Box2Abb[170]=-5. + 3.*m_z12 - 3.*m_z2k;

  Box2Abb[171]=2. + Box2Abb[170]*m_z12 + 4.*m_z2k;

  Box2Abb[172]=Box2Abb[171]*m_x_2 + Box2Abb[72]*m_x_3 - pow(Box2Abb[4],2.)*m_x*m_z12 + Box2Abb[155]*m_x*m_z2k_2 - Box2Abb[5]*m_z12*m_z2k_2;

  Box2Abb[173]=Box2Abb[172]*m_cL + Box2Abb[169]*m_cR;

  Box2Abb[174]=-6. + 3.*m_z12 - 4.*m_z2k;

  Box2Abb[175]=-2. + Box2Abb[174]*m_z12;

  Box2Abb[176]=2. + 4.*m_z2k;

  Box2Abb[177]=Box2Abb[176]*m_z12 - m_z12_2 + 2.*Box2Abb[12]*m_z2k;

  Box2Abb[178]=Box2Abb[175]*m_x_2 + Box2Abb[177]*m_x*m_z12 + 2.*m_x_3*m_z12 - Box2Abb[137]*m_z12_2*m_z2k;

  Box2Abb[179]=-2. + m_z2k;

  Box2Abb[180]=pow(Box2Abb[61],2.) + Box2Abb[179]*m_z12 + m_z12_2;

  Box2Abb[181]=4. - 3.*m_z12 + 6.*m_z2k;

  Box2Abb[182]=2. + Box2Abb[181]*m_z12;

  Box2Abb[183]=1. + 4.*m_z2k;

  Box2Abb[184]=6. + m_z2k;

  Box2Abb[185]=2. + Box2Abb[184]*m_z2k;

  Box2Abb[186]=-1. + Box2Abb[185]*m_z12 - Box2Abb[183]*m_z12_2 + 2.*m_z2k_3;

  Box2Abb[187]=2. + m_z2k;

  Box2Abb[188]=-6. + Box2Abb[187]*m_z12 + m_z12_2 - 2.*Box2Abb[150]*m_z2k;

  Box2Abb[189]=2. + Box2Abb[188]*m_z12 - 2.*m_z2k;

  Box2Abb[190]=Box2Abb[189]*m_x_2 + Box2Abb[182]*m_x_3 + Box2Abb[186]*m_x*m_z12 - 2.*m_x_4*m_z12 + Box2Abb[180]*m_z12_2*m_z2k;

  Box2Abb[191]=Box2Abb[178]*Box2Abb[8]*m_cL + Box2Abb[190]*m_cR;

  Box2Abb[192]=-5. + 3.*m_z12;

  Box2Abb[193]=4. + 2.*Box2Abb[192]*m_z12 - 6.*m_z2k + 9.*m_z12*m_z2k + 2.*m_z2k_2;

  Box2Abb[194]=-2. + m_z12 + 6.*m_z2k;

  Box2Abb[195]=2. + Box2Abb[194]*m_z12;

  Box2Abb[196]=-8. + 9.*m_z12 + 6.*m_z2k;

  Box2Abb[197]=2. + Box2Abb[196]*m_z12;

  Box2Abb[198]=Box2Abb[195]*m_x_3 - 2.*m_x_4*m_z12 - Box2Abb[197]*m_x_2*m_z2k + Box2Abb[193]*m_x*m_z12*m_z2k - pow(Box2Abb[5],2.)*m_z12_2*m_z2k;

  Box2Abb[199]=m_z12 - 2.*m_z2k;

  Box2Abb[200]=-2. + 3.*Box2Abb[199]*m_z12;

  Box2Abb[201]=-4. + m_z12;

  Box2Abb[202]=pow(Box2Abb[4],2.)*m_z12 + 4.*Box2Abb[4]*m_z2k - Box2Abb[201]*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[203]=4. + m_z12;

  Box2Abb[204]=6. - 9.*m_z12 + 5.*m_z12_2 + Box2Abb[203]*m_z2k - 6.*m_z2k_2;

  Box2Abb[205]=2.*Box2Abb[12] - Box2Abb[204]*m_z12;

  Box2Abb[206]=Box2Abb[205]*m_x_2 + Box2Abb[200]*m_x_3 + Box2Abb[202]*m_x*m_z12 + 2.*m_x_4*m_z12 - Box2Abb[5]*m_z12_2*m_z2k_2;

  Box2Abb[207]=Box2Abb[206]*m_cL + Box2Abb[198]*m_cR;

  Box2Abb[208]=-Box2Abb[138]*m_x + m_x_2 + m_z12 + Box2Abb[4]*m_z2k + m_z2k_2;

  Box2Abb[209]=-Box2Abb[208]*m_cL + Box2Abb[142]*m_cR;

  Box2Abb[210]=1. + m_x;

  Box2Abb[211]=pow(Box2Abb[1],2.) - 2.*Box2Abb[210]*m_z2k + m_z2k_2;

  Box2Abb[212]=1. + Box2Abb[1]*m_z12;

  Box2Abb[213]=-2. + Box2Abb[143]*m_x + m_z12;

  Box2Abb[214]=1. + 3.*m_x_2 - Box2Abb[210]*m_z12;

  Box2Abb[215]=-2. - 3.*m_x + m_z12;

  Box2Abb[216]=Box2Abb[1]*Box2Abb[212]*m_x - Box2Abb[213]*m_x*m_z2k + Box2Abb[214]*m_z2k_2 + Box2Abb[215]*m_z2k_3 + m_z2k_4;

  Box2Abb[217]=Box2Abb[216]*m_cL + Box2Abb[148]*Box2Abb[211]*m_cR*m_z2k;

  Box2Abb[218]=-1. + 2.*m_x + 2.*m_z2k;

  Box2Abb[219]=7. + m_z12 + 2.*m_z2k;

  Box2Abb[220]=-3. + 5.*m_z12;

  Box2Abb[221]=1. - 2.*m_z12 + Box2Abb[4]*m_z12*m_z2k + Box2Abb[220]*m_z2k_2 + 2.*m_z2k_3;

  Box2Abb[222]=3. + m_z2k;

  Box2Abb[223]=4. - 3.*Box2Abb[61]*m_z12 - 2.*Box2Abb[222]*m_z2k;

  Box2Abb[224]=Box2Abb[221]*m_x + Box2Abb[223]*m_x_2 - Box2Abb[219]*m_x_3 + 2.*m_x_4 - Box2Abb[5]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[225]=Box2Abb[224]*m_cL - Box2Abb[142]*Box2Abb[218]*m_cR*m_x;

  Box2Abb[226]=-Box2Abb[4]*m_x_2 + 2.*m_x_3 + Box2Abb[18]*m_x*m_z2k + Box2Abb[5]*m_z12*m_z2k;

  Box2Abb[227]=5. + 8.*m_z12 + 4.*m_z2k;

  Box2Abb[228]=6. + m_z12 + 8.*m_z2k;

  Box2Abb[229]=-8. + 2.*Box2Abb[228]*m_z12 - m_z2k;

  Box2Abb[230]=-1. + 2.*m_z2k;

  Box2Abb[231]=3. - 2.*m_z2k;

  Box2Abb[232]=-1. + Box2Abb[231]*m_z2k;

  Box2Abb[233]=m_z12 + Box2Abb[232]*m_z2k + Box2Abb[230]*m_z12*m_z2k;

  Box2Abb[234]=3. + 4.*m_z2k;

  Box2Abb[235]=-1. + 4.*m_z2k - 8.*m_z2k_2;

  Box2Abb[236]=-4. + m_z2k + 4.*m_z2k_2;

  Box2Abb[237]=3. + Box2Abb[235]*m_z12 - Box2Abb[234]*m_z12_2 + Box2Abb[236]*m_z2k;

  Box2Abb[238]=Box2Abb[233]*Box2Abb[5] + Box2Abb[237]*m_x + Box2Abb[229]*m_x_2 - Box2Abb[227]*m_x_3 + 2.*m_x_4;

  Box2Abb[239]=-Box2Abb[142]*Box2Abb[226]*m_cR + Box2Abb[238]*m_cL*m_x;

  Box2Abb[240]=-1. + m_z12 - 4.*m_z2k;

  Box2Abb[241]=-2. + Box2Abb[240]*m_z12 - 6.*m_z2k;

  Box2Abb[242]=3. + m_z12;

  Box2Abb[243]=1. - 2.*Box2Abb[4]*m_z12 - 2.*m_z2k + Box2Abb[242]*m_z12*m_z2k + 6.*Box2Abb[9]*m_z2k_2;

  Box2Abb[244]=-3. + m_z12 - m_z12_2 + m_z12_3;

  Box2Abb[245]=4. + m_z12 - 3.*m_z12_2;

  Box2Abb[246]=pow(Box2Abb[4],2.) + Box2Abb[244]*m_z2k + Box2Abb[245]*m_z2k_2 - 2.*Box2Abb[11]*m_z2k_3;

  Box2Abb[247]=Box2Abb[246]*m_x + Box2Abb[243]*m_x_2 + Box2Abb[241]*m_x_3 + Box2Abb[70]*m_x_4 + Box2Abb[180]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[248]=5. + Box2Abb[70]*m_z12;

  Box2Abb[249]=6. + 5.*m_z12;

  Box2Abb[250]=12. + m_z12 - 3.*m_z12_2 + 2.*Box2Abb[248]*m_z2k + 2.*Box2Abb[249]*m_z2k_2;

  Box2Abb[251]=-2. + m_z12 + 5.*m_z12_2;

  Box2Abb[252]=4. + 5.*m_z12;

  Box2Abb[253]=-3.*Box2Abb[72]*Box2Abb[9] + Box2Abb[155]*Box2Abb[242]*m_z2k + 2.*Box2Abb[251]*m_z2k_2 + 2.*Box2Abb[252]*m_z2k_3;

  Box2Abb[254]=-1. + m_z12 - 5.*m_z2k;

  Box2Abb[255]=-8.*Box2Abb[12] + Box2Abb[254]*m_z12;

  Box2Abb[256]=1. + 3.*m_z2k;

  Box2Abb[257]=Box2Abb[256]*Box2Abb[61]*m_z12 + pow(Box2Abb[61],2.)*m_z2k + 2.*m_z12_2*m_z2k;

  Box2Abb[258]=-4. + 5.*m_z2k;

  Box2Abb[259]=-2. + Box2Abb[258]*m_z2k;

  Box2Abb[260]=-7. + 10.*m_z2k;

  Box2Abb[261]=2. + Box2Abb[260]*m_z2k;

  Box2Abb[262]=-1. + Box2Abb[261]*m_z2k;

  Box2Abb[263]=Box2Abb[12]*Box2Abb[259]*Box2Abb[61]*m_z12 + Box2Abb[262]*m_z12_2 + 2.*pow(Box2Abb[61],3.)*m_z2k + 2.*m_z12_3*m_z2k_2;

  Box2Abb[264]=Box2Abb[263]*m_x - Box2Abb[253]*m_x_2 + Box2Abb[250]*m_x_3 + Box2Abb[255]*m_x_4 + Box2Abb[70]*m_x_5 - Box2Abb[257]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[265]=-Box2Abb[264]*Box2Abb[8]*m_cL + Box2Abb[142]*Box2Abb[247]*m_cR;

  Box2Abb[266]=2. - 5.*m_z12 + 6.*Box2Abb[4]*m_z2k;

  Box2Abb[267]=-2. + m_z12 - 6.*m_z2k + 4.*m_z12*m_z2k;

  Box2Abb[268]=9. - 4.*m_z2k;

  Box2Abb[269]=-5. + Box2Abb[268]*m_z2k;

  Box2Abb[270]=2.*pow(Box2Abb[61],2.) + Box2Abb[269]*m_z12 - 2.*Box2Abb[61]*m_z12_2 + m_z12_3;

  Box2Abb[271]=-Box2Abb[267]*m_x_3 + Box2Abb[72]*m_x_4 + Box2Abb[270]*m_x*m_z2k + Box2Abb[266]*m_x_2*m_z2k + pow(Box2Abb[5],2.)*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[272]=-4. + 7.*m_z12;

  Box2Abb[273]=27. + 5.*m_z12;

  Box2Abb[274]=-1. + 6.*m_z12;

  Box2Abb[275]=5. + 2.*Box2Abb[274]*m_z12;

  Box2Abb[276]=3.*Box2Abb[4]*m_z12_2 + Box2Abb[272]*m_z12*m_z2k + Box2Abb[273]*Box2Abb[4]*m_z12*m_z2k_2 + 2.*Box2Abb[275]*m_z2k_3 + 5.*Box2Abb[155]*m_z2k_4;

  Box2Abb[277]=8. + 3.*m_z12 + 6.*m_z2k;

  Box2Abb[278]=8. - Box2Abb[277]*m_z12 + 10.*m_z2k;

  Box2Abb[279]=7. + 13.*m_z2k;

  Box2Abb[280]=8. + 5.*m_z2k;

  Box2Abb[281]=5. + 3.*Box2Abb[280]*m_z2k;

  Box2Abb[282]=9. + 10.*m_z2k;

  Box2Abb[283]=6. + Box2Abb[282]*m_z2k;

  Box2Abb[284]=-2.*Box2Abb[283] + Box2Abb[281]*m_z12 + Box2Abb[279]*m_z12_2 + m_z12_3;

  Box2Abb[285]=3. + 10.*m_z2k;

  Box2Abb[286]=3. + Box2Abb[285]*m_z2k;

  Box2Abb[287]=17. - 20.*Box2Abb[12]*m_z2k;

  Box2Abb[288]=3. + Box2Abb[287]*m_z2k;

  Box2Abb[289]=2. + Box2Abb[288]*m_z12 - Box2Abb[234]*Box2Abb[59]*m_z12_2 - 3.*Box2Abb[12]*m_z12_3 + 2.*Box2Abb[286]*m_z2k;

  Box2Abb[290]=1. - 5.*m_z2k + 6.*m_z2k_3;

  Box2Abb[291]=9. + 5.*m_z2k;

  Box2Abb[292]=-3. + Box2Abb[291]*m_z2k;

  Box2Abb[293]=1. + Box2Abb[292]*m_z2k;

  Box2Abb[294]=11. + 13.*m_z2k;

  Box2Abb[295]=-8. + Box2Abb[294]*m_z2k;

  Box2Abb[296]=2. + Box2Abb[295]*m_z2k;

  Box2Abb[297]=Box2Abb[290]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[296]*Box2Abb[61]*m_z12_2 + Box2Abb[293]*m_z12_3 - 2.*pow(Box2Abb[61],3.)*m_z2k_2;

  Box2Abb[298]=-Box2Abb[297]*m_x + Box2Abb[276]*m_x_2 + Box2Abb[289]*m_x_3 + Box2Abb[284]*m_x_4 + Box2Abb[278]*m_x_5 + Box2Abb[72]*m_x_6 + Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[65]*m_z12*m_z2k_2;

  Box2Abb[299]=Box2Abb[298]*m_cL - Box2Abb[142]*Box2Abb[271]*m_cR;

  Box2Abb[319]=(-4.*m_m)/m_s_2;

  Box2Abb[320]=2.*m_m;

  Box2Abb[321]=(4.*m_m*m_z2k)/m_s_2;

  Box2Abb[322]=(4.*m_m)/m_s_2;

  Box2Abb[323]=2.*Box2Abb[9]*m_x - m_z12;

  Box2Abb[324]=m_x + 2.*Box2Abb[9]*m_x_2 - 2.*m_x*m_z12 - 2.*Box2Abb[9]*m_x*m_z2k + m_z12_2*m_z2k;

  Box2Abb[325]=Box2Abb[324]*m_cL + Box2Abb[131]*m_cR*m_x + Box2Abb[323]*m_cR*m_z2k;

  Box2Abb[326]=-Box2Abb[139]*m_x + 2.*Box2Abb[9]*m_x_2 + Box2Abb[5]*m_z12;

  Box2Abb[327]=Box2Abb[326]*m_cL - Box2Abb[136]*m_cR;

  Box2Abb[328]=3. - 2.*m_z12;

  Box2Abb[329]=2.*Box2Abb[4]*m_x + 3.*m_x_2 + Box2Abb[328]*m_z12;

  Box2Abb[330]=-Box2Abb[1]*Box2Abb[143]*m_x + Box2Abb[329]*m_z2k - Box2Abb[145]*m_z2k_2 + m_z2k_3;

  Box2Abb[331]=Box2Abb[142]*Box2Abb[8]*m_cL + Box2Abb[330]*m_cR;

  Box2Abb[332]=-Box2Abb[149]*Box2Abb[5]*Box2Abb[61] - Box2Abb[151]*m_x + Box2Abb[150]*m_x_2 - m_x_3;

  Box2Abb[333]=Box2Abb[142]*Box2Abb[148]*m_cL + Box2Abb[332]*m_cR;

  Box2Abb[334]=Box2Abb[160]*m_cL + Box2Abb[164]*m_cR;

  Box2Abb[335]=Box2Abb[169]*m_cL + Box2Abb[172]*m_cR;

  Box2Abb[336]=Box2Abb[190]*m_cL + Box2Abb[178]*Box2Abb[8]*m_cR;

  Box2Abb[337]=Box2Abb[198]*m_cL + Box2Abb[206]*m_cR;

  Box2Abb[338]=Box2Abb[142]*m_cL - Box2Abb[208]*m_cR;

  Box2Abb[339]=Box2Abb[216]*m_cR + Box2Abb[148]*Box2Abb[211]*m_cL*m_z2k;

  Box2Abb[340]=-m_x + m_z2k;

  Box2Abb[341]=3. - 5.*m_z12;

  Box2Abb[342]=-1. + 2.*m_z12 - Box2Abb[4]*m_z12*m_z2k + Box2Abb[341]*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[343]=-4. + 3.*Box2Abb[61]*m_z12 + 2.*Box2Abb[222]*m_z2k;

  Box2Abb[344]=Box2Abb[342]*m_x + Box2Abb[343]*m_x_2 + Box2Abb[219]*m_x_3 - 2.*m_x_4 + Box2Abb[5]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[345]=Box2Abb[344]*m_cR + Box2Abb[142]*Box2Abb[218]*m_cL*m_x;

  Box2Abb[346]=-Box2Abb[142]*Box2Abb[226]*m_cL + Box2Abb[238]*m_cR*m_x;

  Box2Abb[347]=Box2Abb[142]*Box2Abb[247]*m_cL - Box2Abb[264]*Box2Abb[8]*m_cR;

  Box2Abb[348]=-Box2Abb[142]*Box2Abb[271]*m_cL + Box2Abb[298]*m_cR;

  Box2Abb[369]=Box2Abb[340]*m_cL + Box2Abb[148]*m_cR;

  Box2Abb[370]=m_x - 2.*m_x_2 + 2.*m_x*m_z2k - m_z12*m_z2k;

  Box2Abb[371]=Box2Abb[370]*m_cL + Box2Abb[17]*m_cR*m_x;

  Box2Abb[372]=-1. - 2.*m_x + m_z12 + 2.*m_z2k;

  Box2Abb[373]=Box2Abb[19]*m_cR + Box2Abb[372]*m_cL*m_x;

  Box2Abb[374]=pow(Box2Abb[4],2.) + 2.*Box2Abb[4]*m_z2k + 2.*m_z2k_2;

  Box2Abb[375]=-3. + 2.*m_z12 + 2.*m_z2k;

  Box2Abb[376]=1. + Box2Abb[375]*m_z12;

  Box2Abb[377]=-2.*Box2Abb[376]*m_x + Box2Abb[374]*m_z12 + 2.*m_x_2*m_z12;

  Box2Abb[378]=-2. + m_z12 + 4.*m_z2k;

  Box2Abb[379]=2. + Box2Abb[378]*m_z12;

  Box2Abb[380]=Box2Abb[379]*m_x_2 - 2.*m_x_3*m_z12 - 2.*Box2Abb[12]*m_x*m_z12*m_z2k + m_z12_2*m_z2k_2;

  Box2Abb[381]=Box2Abb[380]*m_cL + Box2Abb[377]*m_cR*m_x;

  Box2Abb[382]=Box2Abb[72]*m_x + m_z12*m_z2k;

  Box2Abb[383]=-2. - Box2Abb[179]*m_z12 + 4.*m_z2k;

  Box2Abb[384]=-1. + m_z12_2;

  Box2Abb[385]=pow(Box2Abb[4],2.) + 2.*Box2Abb[384]*m_z2k + Box2Abb[70]*m_z2k_2;

  Box2Abb[386]=-Box2Abb[385]*m_x + Box2Abb[383]*m_x_2 + Box2Abb[72]*m_x_3 + pow(Box2Abb[5],2.)*m_z12*m_z2k;

  Box2Abb[387]=-Box2Abb[382]*pow(Box2Abb[8],2.)*m_cL + Box2Abb[386]*m_cR;

  Box2Abb[388]=Box2Abb[155]*m_x + m_z12 - Box2Abb[137]*m_z12;

  Box2Abb[389]=-2. + m_z12 + m_z2k;

  Box2Abb[390]=pow(Box2Abb[4],2.) + 3.*m_z12*m_z2k;

  Box2Abb[391]=1. + m_z12 + 3.*m_z12*m_z2k;

  Box2Abb[392]=m_x - Box2Abb[391]*m_x_2 + m_x_3*m_z12 + Box2Abb[390]*m_x*m_z2k - Box2Abb[389]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[393]=1. + m_z12 + Box2Abb[72]*m_z2k + 3.*m_z2k_2;

  Box2Abb[394]=-1. + Box2Abb[393]*m_z12 + m_z2k;

  Box2Abb[395]=-Box2Abb[394]*m_x + Box2Abb[391]*m_x_2 + Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12 - m_x_3*m_z12;

  Box2Abb[396]=Box2Abb[392]*m_cL + Box2Abb[395]*m_cR;

  Box2Abb[397]=-1. + m_z12 + 3.*m_z2k;

  Box2Abb[398]=Box2Abb[256]*m_x_2 - m_x_3 + Box2Abb[5]*Box2Abb[61]*m_z2k - Box2Abb[397]*m_x*m_z2k;

  Box2Abb[399]=1. + m_z12 + 3.*m_z2k;

  Box2Abb[400]=m_z12 + Box2Abb[4]*m_z2k + m_z2k_2;

  Box2Abb[401]=m_x - Box2Abb[150]*m_x_2 + m_x_3 - Box2Abb[400]*m_z2k + Box2Abb[399]*m_x*m_z2k;

  Box2Abb[402]=Box2Abb[401]*m_cL + Box2Abb[398]*m_cR;

  Box2Abb[403]=-Box2Abb[5]*pow(Box2Abb[61],2.) + m_x_3 + Box2Abb[12]*m_x*m_z12 + 3.*Box2Abb[61]*m_x*m_z2k - 3.*m_x_2*m_z2k;

  Box2Abb[404]=Box2Abb[398]*m_cL + Box2Abb[403]*m_cR;

  Box2Abb[405]=-4.*m_x + m_z12;

  Box2Abb[406]=Box2Abb[382]*Box2Abb[8]*m_cL + Box2Abb[372]*Box2Abb[4]*m_cR*m_x;

  Box2Abb[407]=-Box2Abb[379]*m_x_2 + 2.*m_x_3*m_z12 + 2.*Box2Abb[12]*m_x*m_z12*m_z2k - m_z12_2*m_z2k_2;

  Box2Abb[408]=pow(Box2Abb[4],2.) + 4.*Box2Abb[4]*m_z12*m_z2k + 3.*m_z12*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[409]=-4. + m_z12 + 6.*m_z2k;

  Box2Abb[410]=2. + Box2Abb[409]*m_z12;

  Box2Abb[411]=4. + m_z2k;

  Box2Abb[412]=-6. + Box2Abb[411]*m_z12 + 4.*m_z2k - 6.*m_z2k_2;

  Box2Abb[413]=2. + Box2Abb[412]*m_z12 - 2.*m_z2k;

  Box2Abb[414]=-Box2Abb[413]*m_x_2 - Box2Abb[410]*m_x_3 + Box2Abb[408]*m_x*m_z12 + 2.*m_x_4*m_z12 - pow(Box2Abb[5],2.)*m_z12_2*m_z2k;

  Box2Abb[415]=-Box2Abb[407]*Box2Abb[8]*m_cL + Box2Abb[414]*m_cR;

  Box2Abb[416]=-2. + 3.*m_z12 + m_z2k;

  Box2Abb[417]=-10. + 9.*m_z12 + 4.*m_z2k;

  Box2Abb[418]=2. + Box2Abb[417]*m_z12;

  Box2Abb[419]=-Box2Abb[418]*m_x_2 + 2.*Box2Abb[416]*Box2Abb[5]*m_x*m_z12 + 2.*m_x_3*m_z12 - pow(Box2Abb[5],2.)*m_z12_2;

  Box2Abb[420]=-6. + m_z2k;

  Box2Abb[421]=4. + Box2Abb[420]*m_z12 + 2.*m_z12_2 + 2.*Box2Abb[179]*m_z2k;

  Box2Abb[422]=-8. + 5.*m_z12 + 6.*m_z2k;

  Box2Abb[423]=2. + Box2Abb[422]*m_z12;

  Box2Abb[424]=-1. + 7.*m_z2k;

  Box2Abb[425]=-2. + Box2Abb[424]*m_z12 + m_z12_2 + 6.*Box2Abb[179]*m_z2k;

  Box2Abb[426]=2.*Box2Abb[12] + Box2Abb[425]*m_z12;

  Box2Abb[427]=-Box2Abb[426]*m_x_2 + Box2Abb[423]*m_x_3 - 2.*m_x_4*m_z12 + Box2Abb[421]*m_x*m_z12*m_z2k + Box2Abb[5]*m_z12_2*m_z2k_2;

  Box2Abb[428]=Box2Abb[427]*m_cL + Box2Abb[419]*Box2Abb[8]*m_cR;

  Box2Abb[429]=Box2Abb[70]*m_x + Box2Abb[61]*m_z12;

  Box2Abb[430]=3. + 2.*m_x - 2.*m_z12;

  Box2Abb[431]=-1. + m_x + pow(Box2Abb[1],2.)*m_z12 + Box2Abb[430]*m_z2k - Box2Abb[70]*m_z2k_2;

  Box2Abb[432]=-Box2Abb[431]*m_cR*m_x + Box2Abb[429]*Box2Abb[8]*m_cL*m_z2k;

  Box2Abb[433]=Box2Abb[5]*Box2Abb[61] + Box2Abb[199]*m_x + m_x_2;

  Box2Abb[434]=2.*m_z12 + m_z2k;

  Box2Abb[435]=m_x - m_x_2 + m_z2k - Box2Abb[434]*m_z2k + 2.*m_x*m_z2k;

  Box2Abb[436]=Box2Abb[435]*m_cL + Box2Abb[433]*m_cR;

  Box2Abb[437]=1. - 2.*m_z12 + m_z2k - 2.*m_z2k_2;

  Box2Abb[438]=Box2Abb[5]*pow(Box2Abb[61],2.) + Box2Abb[437]*m_x + Box2Abb[137]*m_x_2;

  Box2Abb[439]=Box2Abb[438]*m_cR - Box2Abb[433]*m_cL*m_z2k;

  Box2Abb[440]=-1. + m_x + m_z2k + 2.*m_z12*m_z2k;

  Box2Abb[441]=2. - 3.*m_z12 - 2.*m_z12*m_z2k;

  Box2Abb[442]=-Box2Abb[5]*Box2Abb[61] + Box2Abb[441]*m_x + Box2Abb[11]*m_x_2;

  Box2Abb[443]=-Box2Abb[440]*Box2Abb[8]*m_cL + Box2Abb[442]*m_cR;

  Box2Abb[444]=1. + 2.*m_z2k;

  Box2Abb[445]=3. + 2.*m_z12 - 2.*m_z2k;

  Box2Abb[446]=-1. + Box2Abb[445]*m_z2k;

  Box2Abb[447]=m_z12 + 2.*Box2Abb[12]*m_z2k - 8.*m_z12*m_z2k;

  Box2Abb[448]=Box2Abb[446]*Box2Abb[5] + Box2Abb[447]*m_x + Box2Abb[444]*m_x_2 - 2.*m_x_3;

  Box2Abb[449]=-1. + m_z12 + Box2Abb[70]*m_z12*m_z2k + Box2Abb[341]*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[450]=2. + 5.*m_z2k;

  Box2Abb[451]=Box2Abb[450]*m_z12 - 2.*Box2Abb[12]*m_z2k;

  Box2Abb[452]=-Box2Abb[449]*m_x + Box2Abb[451]*m_x_2 - Box2Abb[138]*m_x_3 + 2.*m_x_4 - Box2Abb[5]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[453]=Box2Abb[452]*m_cL + Box2Abb[448]*m_cR*m_x;

  Box2Abb[454]=1. + m_x + 2.*m_x_2 - m_z12;

  Box2Abb[455]=m_x + m_x_2 - 4.*m_x*m_z12 + m_z12_2;

  Box2Abb[456]=-4. + 2.*Box2Abb[455] + m_z12;

  Box2Abb[457]=5. + 2.*m_x;

  Box2Abb[458]=Box2Abb[1]*Box2Abb[454] - Box2Abb[456]*m_z2k - Box2Abb[457]*m_z2k_2 + 2.*m_z2k_3;

  Box2Abb[459]=-1. - 7.*m_z12 + 2.*m_z2k;

  Box2Abb[460]=-3. + 2.*m_z2k;

  Box2Abb[461]=-1. + 5.*m_z2k;

  Box2Abb[462]=Box2Abb[460]*pow(Box2Abb[61],2.) + Box2Abb[461]*Box2Abb[61]*m_z12 + Box2Abb[222]*m_z12_2;

  Box2Abb[463]=11. + 3.*m_z2k;

  Box2Abb[464]=-4. + m_z2k + m_z2k_2;

  Box2Abb[465]=2.*Box2Abb[464] + Box2Abb[463]*m_z12 + 2.*m_z12_2;

  Box2Abb[466]=-Box2Abb[462]*m_x + Box2Abb[465]*m_x_2 + Box2Abb[459]*m_x_3 - 2.*m_x_4 + Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12;

  Box2Abb[467]=Box2Abb[466]*m_cR + Box2Abb[458]*m_cL*m_x;

  Box2Abb[468]=-4. + m_z12 - 8.*m_z2k;

  Box2Abb[469]=-2. + Box2Abb[468]*m_z12;

  Box2Abb[470]=-1. + 3.*m_z2k;

  Box2Abb[471]=1. + m_z2k_2;

  Box2Abb[472]=2.*Box2Abb[471]*pow(Box2Abb[61],2.) + Box2Abb[230]*Box2Abb[470]*Box2Abb[61]*m_z12 + 4.*m_z12_2*m_z2k_2;

  Box2Abb[473]=-3. + 6.*m_z2k;

  Box2Abb[474]=1. + m_z2k + 3.*m_z2k_2;

  Box2Abb[475]=6. + 4.*Box2Abb[474]*m_z12 + Box2Abb[473]*m_z12_2 + 4.*m_z2k;

  Box2Abb[476]=-3. + 14.*m_z2k;

  Box2Abb[477]=4. + 4.*Box2Abb[444]*Box2Abb[61]*m_z2k;

  Box2Abb[478]=2.*Box2Abb[179]*Box2Abb[61] + Box2Abb[477]*m_z12 + Box2Abb[12]*Box2Abb[476]*m_z12_2;

  Box2Abb[479]=-Box2Abb[478]*m_x_2 + Box2Abb[475]*m_x_3 + Box2Abb[469]*m_x_4 + Box2Abb[472]*m_x*m_z12 + 2.*m_x_5*m_z12 + pow(Box2Abb[61],3.)*m_z12_2*m_z2k;

  Box2Abb[480]=-1. + 2.*m_z12_2 - 5.*m_z12*m_z2k;

  Box2Abb[481]=-11. + 3.*Box2Abb[70]*m_z12;

  Box2Abb[482]=11. + 3.*m_z12;

  Box2Abb[483]=1. + Box2Abb[201]*m_z12;

  Box2Abb[484]=2. + Box2Abb[481]*m_z12 + Box2Abb[4]*Box2Abb[482]*m_z12*m_z2k - 6.*Box2Abb[483]*m_z2k_2 - 20.*m_z12*m_z2k_3;

  Box2Abb[485]=-7. + 3.*m_z12;

  Box2Abb[486]=pow(Box2Abb[4],2.) - 2.*Box2Abb[4]*Box2Abb[72]*m_z2k + Box2Abb[4]*Box2Abb[485]*m_z2k_2 + 6.*Box2Abb[4]*m_z2k_3 + 2.*m_z2k_4;

  Box2Abb[487]=2. - 5.*m_z2k;

  Box2Abb[488]=5. + 3.*m_z2k;

  Box2Abb[489]=-5. + 2.*Box2Abb[488]*m_z12 + m_z12_2 + 4.*Box2Abb[487]*m_z2k;

  Box2Abb[490]=-6.*Box2Abb[12] + Box2Abb[489]*m_z12;

  Box2Abb[491]=-3. + m_z2k_2;

  Box2Abb[492]=1. + m_z2k + 7.*Box2Abb[61]*m_z2k_2;

  Box2Abb[493]=-7. + 5.*m_z2k;

  Box2Abb[494]=5. + 2.*Box2Abb[493]*m_z2k;

  Box2Abb[495]=-3. + Box2Abb[494]*m_z2k;

  Box2Abb[496]=2.*pow(Box2Abb[61],3.) + Box2Abb[495]*Box2Abb[61]*m_z12 + 2.*Box2Abb[492]*m_z12_2 + Box2Abb[491]*m_z12_3;

  Box2Abb[497]=Box2Abb[496]*m_x + Box2Abb[484]*m_x_2 - Box2Abb[490]*m_x_3 + 2.*Box2Abb[480]*m_x_4 - Box2Abb[486]*Box2Abb[61]*m_z12 + 2.*m_x_5*m_z12;

  Box2Abb[498]=Box2Abb[479]*Box2Abb[8]*m_cL - Box2Abb[497]*m_cR*m_x;

  Box2Abb[499]=-4. + m_z12 - 10.*m_z2k + 4.*m_z12*m_z2k;

  Box2Abb[500]=7. + m_z12;

  Box2Abb[501]=-11. + Box2Abb[500]*m_z12;

  Box2Abb[502]=-5. + 8.*m_z12;

  Box2Abb[503]=6. + Box2Abb[502]*m_z12;

  Box2Abb[504]=1. + 3.*Box2Abb[4]*m_z12 - 3.*m_z2k + Box2Abb[501]*m_z12*m_z2k + 2.*Box2Abb[503]*m_z2k_2 - 20.*m_z2k_3;

  Box2Abb[505]=4. - 3.*m_z12;

  Box2Abb[506]=15. + m_z12 + m_z12_2 - 5.*m_z12_3;

  Box2Abb[507]=-10. + m_z12 - 6.*m_z12_2;

  Box2Abb[508]=-1. + Box2Abb[505]*m_z12 - 4.*m_z2k + Box2Abb[505]*m_z12_2*m_z2k + Box2Abb[506]*m_z2k_2 + 2.*Box2Abb[507]*m_z2k_3 + 5.*Box2Abb[70]*m_z2k_4;

  Box2Abb[509]=3. - 5.*m_z2k;

  Box2Abb[510]=1. + 5.*m_z2k;

  Box2Abb[511]=1. + m_z12 - Box2Abb[59]*m_z12_2 + 4.*Box2Abb[510]*m_z2k + Box2Abb[509]*m_z12*m_z2k;

  Box2Abb[512]=1. + 2.*Box2Abb[61]*m_z2k;

  Box2Abb[513]=-1. + Box2Abb[184]*m_z2k;

  Box2Abb[514]=-2. + 3.*m_z2k;

  Box2Abb[515]=3. + Box2Abb[514]*m_z2k;

  Box2Abb[516]=-2. + Box2Abb[234]*m_z2k;

  Box2Abb[517]=-Box2Abb[512]*pow(Box2Abb[61],3.) - Box2Abb[516]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[513]*Box2Abb[61]*m_z12_2 + Box2Abb[515]*m_z12_3*m_z2k;

  Box2Abb[518]=Box2Abb[517]*m_x + Box2Abb[508]*m_x_2 + Box2Abb[504]*m_x_3 + Box2Abb[511]*m_x_4 + Box2Abb[499]*m_x_5 - Box2Abb[72]*m_x_6 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[519]=pow(Box2Abb[61],2.) + 4.*Box2Abb[61]*m_z12 + 2.*m_z12_2;

  Box2Abb[520]=-6. + m_z12;

  Box2Abb[521]=-6.*Box2Abb[12] + 3.*m_z12 + 4.*m_z12_2*m_z2k + 2.*Box2Abb[520]*m_z2k_2;

  Box2Abb[522]=6. - 3.*Box2Abb[12]*m_z12 + 8.*m_z2k;

  Box2Abb[523]=-1. + m_z2k + m_z2k_2;

  Box2Abb[524]=-8. + 11.*m_z2k - 3.*m_z2k_3;

  Box2Abb[525]=-2.*pow(Box2Abb[61],3.) + Box2Abb[524]*m_z12 - 4.*Box2Abb[523]*m_z12_2 + 2.*m_z12_3*m_z2k;

  Box2Abb[526]=11. + 2.*Box2Abb[222]*m_z2k;

  Box2Abb[527]=-1. + Box2Abb[526]*m_z2k;

  Box2Abb[528]=2. + Box2Abb[527]*m_z12 - 4.*Box2Abb[187]*m_z12_2*m_z2k - 6.*m_z2k_2 + 8.*m_z2k_3;

  Box2Abb[529]=Box2Abb[528]*m_x_2 + Box2Abb[521]*m_x_3 + Box2Abb[522]*m_x_4 + Box2Abb[72]*m_x_5 + Box2Abb[525]*m_x*m_z2k + Box2Abb[519]*Box2Abb[61]*m_z12*m_z2k_2;

  Box2Abb[530]=Box2Abb[529]*Box2Abb[8]*m_cL + Box2Abb[518]*m_cR;

  Box2Abb[531]=-3. + 4.*m_z12 - 13.*m_z2k;

  Box2Abb[532]=4. + Box2Abb[531]*m_z12 + 8.*m_z2k;

  Box2Abb[533]=10. + m_z12;

  Box2Abb[534]=-8. + Box2Abb[533]*m_z12;

  Box2Abb[535]=6. - 11.*m_z12;

  Box2Abb[536]=Box2Abb[534]*m_z12 + 2.*Box2Abb[203]*m_z12*m_z2k + 2.*Box2Abb[535]*m_z2k_2;

  Box2Abb[537]=-4. + Box2Abb[70]*m_z12;

  Box2Abb[538]=-3. + m_z12;

  Box2Abb[539]=6. + 5.*Box2Abb[538]*m_z12;

  Box2Abb[540]=4. - 9.*m_z12;

  Box2Abb[541]=4. + 3.*Box2Abb[537]*m_z12 + 3.*m_z12_3*m_z2k - 2.*Box2Abb[539]*m_z2k_2 + 2.*Box2Abb[540]*m_z2k_3;

  Box2Abb[542]=-2.*pow(Box2Abb[61],2.) + Box2Abb[424]*Box2Abb[61]*m_z12 + 3.*Box2Abb[12]*m_z12_2;

  Box2Abb[543]=Box2Abb[5]*Box2Abb[542]*Box2Abb[61]*m_x + Box2Abb[541]*m_x_2 - Box2Abb[536]*m_x_3 + Box2Abb[532]*m_x_4 + Box2Abb[155]*m_x_5 - pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12;

  Box2Abb[544]=-8. + m_z12 - 16.*m_z2k;

  Box2Abb[545]=6. + Box2Abb[544]*m_z12 + 10.*m_z2k;

  Box2Abb[546]=7. + 3.*m_z12;

  Box2Abb[547]=3. + Box2Abb[546]*m_z12;

  Box2Abb[548]=-14. + 19.*m_z12;

  Box2Abb[549]=6. + Box2Abb[548]*m_z12;

  Box2Abb[550]=-2. + 5.*m_z12;

  Box2Abb[551]=-2. - Box2Abb[538]*m_z12 + 6.*m_z2k + 3.*Box2Abb[166]*m_z12*m_z2k + 2.*Box2Abb[4]*Box2Abb[547]*m_z2k_2 + 2.*Box2Abb[549]*m_z2k_3 + 5.*Box2Abb[550]*m_z2k_4;

  Box2Abb[552]=-1. - 2.*m_z2k + 4.*m_z2k_2;

  Box2Abb[553]=-2. + m_z2k + 19.*m_z2k_2;

  Box2Abb[554]=2.*Box2Abb[552]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[553]*Box2Abb[61]*m_z12_2 - 2.*pow(Box2Abb[61],3.)*m_z2k + 4.*Box2Abb[444]*m_z12_3*m_z2k;

  Box2Abb[555]=-8. + 6.*m_z2k + 4.*m_z2k_2 - 40.*m_z2k_3;

  Box2Abb[556]=7. + 15.*m_z2k;

  Box2Abb[557]=3. - 2.*Box2Abb[556]*m_z2k;

  Box2Abb[558]=4. + Box2Abb[555]*m_z12 + Box2Abb[557]*m_z12_2 + 6.*m_z2k + 20.*m_z2k_3;

  Box2Abb[559]=3. + 5.*m_z2k;

  Box2Abb[560]=-3. + 7.*m_z2k;

  Box2Abb[561]=16. + 35.*m_z2k;

  Box2Abb[562]=10. + Box2Abb[561]*m_z2k;

  Box2Abb[563]=-6. + Box2Abb[562]*m_z12 + Box2Abb[560]*m_z12_2 - 4.*Box2Abb[559]*m_z2k;

  Box2Abb[564]=Box2Abb[551]*m_x_2 + Box2Abb[558]*m_x_3 + Box2Abb[563]*m_x_4 + Box2Abb[545]*m_x_5 + Box2Abb[155]*m_x_6 - Box2Abb[554]*m_x*m_z2k + Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[65]*m_z12*m_z2k_2;

  Box2Abb[565]=Box2Abb[564]*m_cL - Box2Abb[543]*Box2Abb[8]*m_cR;

  Box2Abb[595]=(4.*m_m*m_z12)/m_s_2;

  Box2Abb[596]=(2.*m_m*m_z12)/m_s;

  Box2Abb[597]=2.*m_m*m_z12;

  Box2Abb[598]=-2.*m_m;

  Box2Abb[599]=(-4.*m_m*m_z2k)/m_s_2;

  Box2Abb[600]=(-4.*m_m*m_x)/m_s_2;

  Box2Abb[601]=Box2Abb[148]*m_cL + Box2Abb[340]*m_cR;

  Box2Abb[602]=Box2Abb[370]*m_cR + Box2Abb[17]*m_cL*m_x;

  Box2Abb[603]=Box2Abb[19]*m_cL + Box2Abb[372]*m_cR*m_x;

  Box2Abb[604]=2. - 2.*Box2Abb[210]*m_z12 + m_z12_2;

  Box2Abb[605]=3. + m_x;

  Box2Abb[606]=1. + 2.*m_x;

  Box2Abb[607]=-2.*m_x + m_z12 + 2.*Box2Abb[605]*m_x*m_z12 - 2.*Box2Abb[606]*m_z12_2 + m_z12_3;

  Box2Abb[608]=Box2Abb[607]*m_cL + Box2Abb[604]*m_cR*m_x;

  Box2Abb[609]=m_cL + m_cR + 2.*m_cL*m_x - 2.*m_cR*m_x - m_cL*m_z12;

  Box2Abb[610]=-2.*m_x + m_z12;

  Box2Abb[611]=Box2Abb[610]*m_cR + 2.*m_cL*m_x;

  Box2Abb[612]=-Box2Abb[608]*m_x + 2.*Box2Abb[609]*m_x*m_z12*m_z2k - Box2Abb[611]*m_z12*m_z2k_2;

  Box2Abb[613]=2. + Box2Abb[179]*m_z12 - 4.*m_z2k;

  Box2Abb[614]=Box2Abb[385]*m_x + Box2Abb[613]*m_x_2 - Box2Abb[72]*m_x_3 - pow(Box2Abb[5],2.)*m_z12*m_z2k;

  Box2Abb[615]=Box2Abb[614]*m_cL + Box2Abb[382]*pow(Box2Abb[8],2.)*m_cR;

  Box2Abb[616]=Box2Abb[395]*m_cL + Box2Abb[392]*m_cR;

  Box2Abb[617]=Box2Abb[398]*m_cL + Box2Abb[401]*m_cR;

  Box2Abb[618]=Box2Abb[403]*m_cL + Box2Abb[398]*m_cR;

  Box2Abb[619]=-Box2Abb[382]*Box2Abb[8]*m_cR + Box2Abb[17]*Box2Abb[4]*m_cL*m_x;

  Box2Abb[620]=Box2Abb[414]*m_cL - Box2Abb[407]*Box2Abb[8]*m_cR;

  Box2Abb[621]=Box2Abb[419]*Box2Abb[8]*m_cL + Box2Abb[427]*m_cR;

  Box2Abb[622]=-Box2Abb[431]*m_cL*m_x + Box2Abb[429]*Box2Abb[8]*m_cR*m_z2k;

  Box2Abb[623]=Box2Abb[433]*m_cL + Box2Abb[435]*m_cR;

  Box2Abb[624]=Box2Abb[438]*m_cL - Box2Abb[433]*m_cR*m_z2k;

  Box2Abb[625]=Box2Abb[442]*m_cL - Box2Abb[440]*Box2Abb[8]*m_cR;

  Box2Abb[626]=-Box2Abb[450]*m_z12 + 2.*Box2Abb[12]*m_z2k;

  Box2Abb[627]=Box2Abb[449]*m_x + Box2Abb[626]*m_x_2 + Box2Abb[138]*m_x_3 - 2.*m_x_4 + Box2Abb[5]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[628]=Box2Abb[627]*m_cR + Box2Abb[458]*m_cL*m_x;

  Box2Abb[629]=1. + 7.*m_z12 - 2.*m_z2k;

  Box2Abb[630]=Box2Abb[462]*m_x - Box2Abb[465]*m_x_2 + Box2Abb[629]*m_x_3 + 2.*m_x_4 - Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12;

  Box2Abb[631]=Box2Abb[630]*m_cL + Box2Abb[448]*m_cR*m_x;

  Box2Abb[632]=Box2Abb[479]*Box2Abb[8]*m_cR - Box2Abb[497]*m_cL*m_x;

  Box2Abb[633]=Box2Abb[518]*m_cL + Box2Abb[529]*Box2Abb[8]*m_cR;

  Box2Abb[634]=-Box2Abb[543]*Box2Abb[8]*m_cL + Box2Abb[564]*m_cR;

  Box2Abb[664]=(-2.*m_m*m_z12)/m_s;

  Box2Abb[665]=-2.*m_m*m_z12;

  Box2Abb[666]=pow(Box2Abb[4],2.) + 2.*m_z12*m_z2k;

  Box2Abb[667]=-Box2Abb[666]*m_x + m_x_2*m_z12 + m_z12*m_z2k_2;

  Box2Abb[668]=-2. + 2.*m_z12 + m_z2k;

  Box2Abb[669]=6. + m_z12;

  Box2Abb[670]=-6. + Box2Abb[669]*m_z12 + 8.*m_z2k;

  Box2Abb[671]=Box2Abb[12]*m_z12 - Box2Abb[256]*m_z12_2 + m_z2k - 2.*m_z2k_2;

  Box2Abb[672]=2.*Box2Abb[671]*m_x + Box2Abb[670]*m_x_2 - 4.*m_x_3 + Box2Abb[668]*m_z12_2*m_z2k;

  Box2Abb[673]=Box2Abb[672]*m_cL + 2.*Box2Abb[667]*m_cR*m_x;

  Box2Abb[674]=-1. - m_x + m_z12 + m_z2k;

  Box2Abb[675]=pow(Box2Abb[4],2.) + Box2Abb[505]*m_x;

  Box2Abb[676]=Box2Abb[675]*m_x + 2.*Box2Abb[72]*m_x*m_z2k + m_z12*m_z2k_2;

  Box2Abb[677]=-2. + m_z12 + 8.*m_z2k;

  Box2Abb[678]=2. + Box2Abb[677]*m_z12;

  Box2Abb[679]=8. + 7.*m_z2k;

  Box2Abb[680]=-1. + Box2Abb[61]*m_z2k;

  Box2Abb[681]=6.*Box2Abb[680]*m_z12 + Box2Abb[679]*m_z12_2 - 2.*m_z12_3 + 2.*Box2Abb[512]*m_z2k;

  Box2Abb[682]=3. + Box2Abb[679]*m_z2k;

  Box2Abb[683]=5. + 6.*m_z2k - 8.*m_z2k_2;

  Box2Abb[684]=2. + Box2Abb[683]*m_z2k;

  Box2Abb[685]=2.*Box2Abb[684]*m_z12 - 2.*Box2Abb[682]*m_z12_2 + Box2Abb[150]*m_z12_3 - 4.*m_z2k_2;

  Box2Abb[686]=2. - 10.*m_z2k;

  Box2Abb[687]=-1. + m_z2k - 2.*m_z2k_2;

  Box2Abb[688]=12.*Box2Abb[687] + Box2Abb[686]*m_z12 + m_z12_2;

  Box2Abb[689]=8. + Box2Abb[688]*m_z12 - 8.*m_z2k;

  Box2Abb[690]=Box2Abb[685]*m_x_2 - Box2Abb[689]*m_x_3 - 2.*Box2Abb[678]*m_x_4 + 4.*m_x_5*m_z12 + Box2Abb[681]*m_x*m_z12*m_z2k - Box2Abb[668]*m_z12_3*m_z2k_2;

  Box2Abb[691]=Box2Abb[690]*m_cL + 2.*Box2Abb[676]*Box2Abb[8]*m_cR*m_x*m_z12;

  Box2Abb[692]=3. - 3.*m_z12 - 2.*m_z2k;

  Box2Abb[693]=-1. + 2.*m_z12;

  Box2Abb[694]=3.*Box2Abb[4]*m_z12 - 2.*m_z2k + 6.*Box2Abb[693]*m_z12*m_z2k + 4.*Box2Abb[9]*m_z2k_2;

  Box2Abb[695]=5. + m_z12 + 4.*m_z2k;

  Box2Abb[696]=-5. + Box2Abb[695]*m_z12 + 4.*m_z2k;

  Box2Abb[697]=Box2Abb[694]*m_x - 2.*Box2Abb[696]*m_x_2 + 4.*Box2Abb[9]*m_x_3 + Box2Abb[692]*m_z12_2*m_z2k;

  Box2Abb[698]=Box2Abb[697]*m_cL + 2.*Box2Abb[382]*Box2Abb[4]*m_cR*m_x;

  Box2Abb[699]=2. - 3.*m_z12;

  Box2Abb[700]=Box2Abb[699]*m_x + Box2Abb[5]*m_z12;

  Box2Abb[701]=7. - 9.*m_z12 + m_z12_2 - 4.*Box2Abb[9]*m_z2k;

  Box2Abb[702]=7. + Box2Abb[520]*m_z12;

  Box2Abb[703]=1. + 4.*m_z12;

  Box2Abb[704]=2. - Box2Abb[702]*m_z12 - 6.*m_z2k + 2.*Box2Abb[703]*m_z12*m_z2k + 4.*Box2Abb[9]*m_z2k_2;

  Box2Abb[705]=Box2Abb[704]*m_x + 2.*Box2Abb[701]*m_x_2 + 4.*Box2Abb[9]*m_x_3 - 2.*Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[706]=Box2Abb[705]*m_cL - 2.*Box2Abb[4]*Box2Abb[700]*m_cR*m_x;

  Box2Abb[707]=-Box2Abb[444]*m_x + m_x_2 + Box2Abb[65]*m_z2k;

  Box2Abb[708]=-3.*Box2Abb[4]*m_z12 + 2.*Box2Abb[61]*m_z2k;

  Box2Abb[709]=-5. + 7.*m_z12 - 4.*m_z2k + 8.*m_z12*m_z2k + 12.*m_z2k_2;

  Box2Abb[710]=-6. + 5.*m_z12;

  Box2Abb[711]=-3. + 3.*m_z12 + 5.*m_z2k + Box2Abb[710]*m_z12*m_z2k + 4.*Box2Abb[72]*m_z2k_2 + 8.*m_z2k_3;

  Box2Abb[712]=-Box2Abb[711]*m_x + Box2Abb[709]*m_x_2 - 4.*Box2Abb[56]*m_x_3 + 2.*m_x_4 + Box2Abb[61]*Box2Abb[708]*m_z2k;

  Box2Abb[713]=Box2Abb[712]*m_cL - 2.*Box2Abb[4]*Box2Abb[707]*m_cR*m_x;

  Box2Abb[714]=2. - 6.*m_z12 - 8.*m_z2k;

  Box2Abb[715]=12. - 7.*m_z12;

  Box2Abb[716]=7. - 5.*m_z12;

  Box2Abb[717]=3. + Box2Abb[201]*m_z12 - 9.*m_z2k + Box2Abb[715]*m_z12*m_z2k + 2.*Box2Abb[716]*m_z2k_2 - 8.*m_z2k_3;

  Box2Abb[718]=5. - 6.*m_z2k;

  Box2Abb[719]=4. + 7.*m_z2k;

  Box2Abb[720]=5. - 2.*Box2Abb[719]*m_z12 + m_z12_2 + 2.*Box2Abb[718]*m_z2k;

  Box2Abb[721]=Box2Abb[717]*m_x - Box2Abb[720]*m_x_2 + Box2Abb[714]*m_x_3 + 2.*m_x_4 + 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[722]=Box2Abb[721]*m_cL - 2.*Box2Abb[4]*Box2Abb[433]*m_cR*m_x;

  Box2Abb[723]=1. + m_z2k - m_z2k_2;

  Box2Abb[724]=6.*Box2Abb[723]*m_z12 - Box2Abb[679]*m_z12_2 + 2.*m_z12_3 - 2.*Box2Abb[512]*m_z2k;

  Box2Abb[725]=-5. + 4.*m_z2k;

  Box2Abb[726]=-2. + Box2Abb[444]*Box2Abb[725]*m_z2k;

  Box2Abb[727]=2.*Box2Abb[726]*m_z12 + 2.*Box2Abb[682]*m_z12_2 - Box2Abb[150]*m_z12_3 + 4.*m_z2k_2;

  Box2Abb[728]=Box2Abb[727]*m_x_2 + Box2Abb[689]*m_x_3 + 2.*Box2Abb[678]*m_x_4 - 4.*m_x_5*m_z12 + Box2Abb[724]*m_x*m_z12*m_z2k + Box2Abb[668]*m_z12_3*m_z2k_2;

  Box2Abb[729]=Box2Abb[728]*m_cL - 2.*Box2Abb[676]*Box2Abb[8]*m_cR*m_x*m_z12;

  Box2Abb[730]=-pow(Box2Abb[4],2.)*m_x + Box2Abb[166]*m_x_2 - 2.*Box2Abb[72]*m_x*m_z2k - m_z12*m_z2k_2;

  Box2Abb[731]=4. + m_z12 - 8.*m_z2k;

  Box2Abb[732]=-2. + Box2Abb[731]*m_z12;

  Box2Abb[733]=-2. + 2.*m_z12 + 3.*m_z2k;

  Box2Abb[734]=8. + m_z12;

  Box2Abb[735]=12. + m_z12;

  Box2Abb[736]=-36. + 3.*Box2Abb[734]*m_z12 + 2.*Box2Abb[735]*m_z2k - 24.*m_z2k_2;

  Box2Abb[737]=12. + Box2Abb[736]*m_z12 - 8.*m_z2k;

  Box2Abb[738]=4. - 5.*m_z2k;

  Box2Abb[739]=4. + 3.*Box2Abb[738]*m_z2k;

  Box2Abb[740]=-1. + m_z2k + m_z2k_3;

  Box2Abb[741]=2.*Box2Abb[740]*m_z12 + Box2Abb[739]*m_z12_2 - 2.*Box2Abb[59]*m_z12_3 + 2.*Box2Abb[512]*Box2Abb[61]*m_z2k;

  Box2Abb[742]=11. + 23.*m_z2k;

  Box2Abb[743]=10. + m_z2k_2;

  Box2Abb[744]=-9. + 4.*Box2Abb[231]*m_z2k;

  Box2Abb[745]=4. + Box2Abb[744]*m_z2k;

  Box2Abb[746]=2.*Box2Abb[745]*m_z12 - 2.*Box2Abb[743]*m_z12_2 + Box2Abb[742]*m_z12_3 + m_z12_4 - 4.*Box2Abb[61]*m_z2k;

  Box2Abb[747]=Box2Abb[746]*m_x_2 - Box2Abb[737]*m_x_3 + 2.*Box2Abb[732]*m_x_4 + Box2Abb[741]*m_x*m_z12 + 4.*m_x_5*m_z12 + Box2Abb[5]*Box2Abb[733]*m_z12_3*m_z2k;

  Box2Abb[748]=-Box2Abb[747]*m_cL + 2.*Box2Abb[148]*Box2Abb[730]*m_cR*m_x*m_z12;

  Box2Abb[749]=-2. + 2.*Box2Abb[215]*m_x + 3.*m_z12;

  Box2Abb[750]=4. + 6.*m_x - 3.*m_z12;

  Box2Abb[751]=2.*m_x_3 + Box2Abb[1]*m_x*m_z12 + Box2Abb[749]*m_z2k + Box2Abb[750]*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[752]=Box2Abb[751]*m_cL + 2.*Box2Abb[10]*m_cR*m_x;

  Box2Abb[753]=1. - Box2Abb[187]*m_z12 + m_z2k;

  Box2Abb[754]=-Box2Abb[5]*Box2Abb[61] + Box2Abb[753]*m_x + m_x_2*m_z12;

  Box2Abb[755]=3.*m_z12 + 2.*m_z2k;

  Box2Abb[756]=8. + 9.*m_z2k;

  Box2Abb[757]=3. - Box2Abb[756]*m_z12 + 2.*m_z2k - 8.*m_z2k_2;

  Box2Abb[758]=-7. + 4.*m_z2k;

  Box2Abb[759]=7. + 2.*Box2Abb[758]*m_z2k;

  Box2Abb[760]=1. + 11.*m_z2k;

  Box2Abb[761]=2. + Box2Abb[61]*Box2Abb[760]*m_z2k;

  Box2Abb[762]=-2. + Box2Abb[761]*m_z12 + m_z2k + 4.*Box2Abb[12]*m_z12_2*m_z2k + Box2Abb[759]*m_z2k_2;

  Box2Abb[763]=8. + 13.*m_z2k;

  Box2Abb[764]=7. + Box2Abb[763]*m_z2k;

  Box2Abb[765]=-5. + Box2Abb[764]*m_z12 + 2.*Box2Abb[59]*Box2Abb[61]*m_z2k + 2.*m_z12_2*m_z2k;

  Box2Abb[766]=-Box2Abb[762]*m_x + Box2Abb[765]*m_x_2 + Box2Abb[757]*m_x_3 + Box2Abb[755]*m_x_4 + 2.*Box2Abb[137]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[767]=Box2Abb[766]*m_cL + 2.*Box2Abb[754]*Box2Abb[8]*m_cR*m_x;

  Box2Abb[768]=-1. + 2.*m_x;

  Box2Abb[769]=2. + 2.*m_x - m_z12;

  Box2Abb[770]=-2. + m_x + 4.*m_x_2;

  Box2Abb[771]=-1. + 9.*m_x;

  Box2Abb[772]=1. - 2.*m_x;

  Box2Abb[773]=2. - 2.*Box2Abb[770]*m_x - 5.*m_z12 + 2.*Box2Abb[771]*m_x*m_z12 + 2.*Box2Abb[772]*m_z12_2;

  Box2Abb[774]=7. + 4.*m_z12;

  Box2Abb[775]=3. + 7.*m_z12;

  Box2Abb[776]=-8. - 2.*Box2Abb[775]*m_x + Box2Abb[774]*m_z12;

  Box2Abb[777]=5. + 4.*m_x - m_z12;

  Box2Abb[778]=Box2Abb[1]*Box2Abb[768]*Box2Abb[769]*m_x + Box2Abb[773]*m_z2k + Box2Abb[776]*m_z2k_2 + 2.*Box2Abb[777]*m_z2k_3 - 4.*m_z2k_4;

  Box2Abb[779]=2. + m_z2k - 5.*m_z2k_2;

  Box2Abb[780]=-2.*pow(Box2Abb[61],2.) + Box2Abb[779]*m_z12 - m_z12_2*m_z2k;

  Box2Abb[781]=Box2Abb[780]*m_x + Box2Abb[70]*m_x_3 + 3.*Box2Abb[61]*m_x_2*m_z12 + Box2Abb[5]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[782]=Box2Abb[778]*m_cL + 2.*Box2Abb[781]*m_cR;

  Box2Abb[783]=2. + 7.*m_z12;

  Box2Abb[784]=-2.*pow(Box2Abb[61],2.) + 3.*Box2Abb[12]*m_z12_2 + 5.*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[785]=10. + 2.*m_z12 + 11.*m_z2k;

  Box2Abb[786]=-8. + Box2Abb[785]*m_z12;

  Box2Abb[787]=Box2Abb[784]*m_x - Box2Abb[786]*m_x_2 + Box2Abb[783]*m_x_3 - Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12;

  Box2Abb[788]=1. + 7.*m_z12 - 6.*m_z2k;

  Box2Abb[789]=20. + 9.*m_z12;

  Box2Abb[790]=43. + 18.*m_z12;

  Box2Abb[791]=-1. + 10.*m_z12;

  Box2Abb[792]=-18.*Box2Abb[12] + Box2Abb[789]*m_z12 + Box2Abb[790]*m_z12*m_z2k + 4.*Box2Abb[791]*m_z2k_2 + 8.*m_z2k_3;

  Box2Abb[793]=33. + 36.*m_z2k;

  Box2Abb[794]=-14. + Box2Abb[793]*m_z12 + 4.*m_z12_2 - 8.*Box2Abb[61]*m_z2k;

  Box2Abb[795]=3. + m_z2k + 6.*Box2Abb[61]*m_z2k_2;

  Box2Abb[796]=7. + 26.*Box2Abb[12]*m_z2k;

  Box2Abb[797]=-3. + 4.*m_z2k;

  Box2Abb[798]=-2. + Box2Abb[797]*m_z2k;

  Box2Abb[799]=2.*Box2Abb[61]*Box2Abb[795] + m_z12 + Box2Abb[796]*m_z12_2 + 7.*Box2Abb[798]*m_z12*m_z2k + 2.*m_z12_3*m_z2k;

  Box2Abb[800]=-2. + 7.*Box2Abb[230]*m_z2k;

  Box2Abb[801]=-11. + 10.*m_z2k;

  Box2Abb[802]=-2. + Box2Abb[801]*m_z2k;

  Box2Abb[803]=pow(Box2Abb[61],2.)*Box2Abb[802]*m_z12 + Box2Abb[61]*Box2Abb[800]*m_z12_2 + 2.*Box2Abb[230]*pow(Box2Abb[61],3.)*m_z2k + 4.*Box2Abb[12]*m_z12_3*m_z2k;

  Box2Abb[804]=Box2Abb[803]*m_x - Box2Abb[799]*m_x_2 + Box2Abb[792]*m_x_3 - Box2Abb[794]*m_x_4 + 2.*Box2Abb[788]*m_x_5 + 4.*m_x_6 - 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k;

  Box2Abb[805]=Box2Abb[804]*m_cL + 2.*Box2Abb[787]*Box2Abb[8]*m_cR*m_x;

  Box2Abb[806]=-6. - 6.*m_z12 + m_z12_2 - 24.*m_z2k;

  Box2Abb[807]=-31. + 6.*m_z12;

  Box2Abb[808]=36. + Box2Abb[807]*m_z12;

  Box2Abb[809]=16. + Box2Abb[166]*m_z12;

  Box2Abb[810]=-8. + Box2Abb[693]*Box2Abb[70]*m_z12;

  Box2Abb[811]=2. + 31.*m_z12;

  Box2Abb[812]=50. + Box2Abb[811]*m_z12;

  Box2Abb[813]=4. + 2.*Box2Abb[538]*m_z12 - 14.*m_z2k + Box2Abb[808]*m_z12*m_z2k - 2.*Box2Abb[4]*Box2Abb[809]*m_z2k_2 + 6.*Box2Abb[810]*m_z2k_3 + Box2Abb[812]*m_z2k_4 - 24.*m_z2k_5;

  Box2Abb[814]=-20. + 7.*m_z12;

  Box2Abb[815]=-5. + m_z12;

  Box2Abb[816]=2. + Box2Abb[815]*m_z12;

  Box2Abb[817]=-5. + 4.*m_z12;

  Box2Abb[818]=-18. + Box2Abb[70]*Box2Abb[817]*m_z12;

  Box2Abb[819]=-3. + 11.*m_z12;

  Box2Abb[820]=15. + Box2Abb[819]*m_z12;

  Box2Abb[821]=6.*Box2Abb[187] + Box2Abb[814]*m_z12 + 6.*Box2Abb[816]*m_z12*m_z2k + 2.*Box2Abb[818]*m_z2k_2 + 4.*Box2Abb[820]*m_z2k_3 - 60.*m_z2k_4;

  Box2Abb[822]=20. + 22.*m_z2k;

  Box2Abb[823]=-6. + Box2Abb[822]*m_z12 - 5.*Box2Abb[444]*m_z12_2 + 10.*Box2Abb[59]*m_z2k;

  Box2Abb[824]=1. - 4.*m_z2k;

  Box2Abb[825]=9. + 7.*m_z2k;

  Box2Abb[826]=7. + Box2Abb[825]*m_z2k;

  Box2Abb[827]=3. + 31.*m_z2k;

  Box2Abb[828]=9. + Box2Abb[827]*m_z2k;

  Box2Abb[829]=16. - 4.*Box2Abb[826]*m_z12 + Box2Abb[828]*m_z12_2 + 2.*m_z12_3*m_z2k + 20.*Box2Abb[824]*m_z2k_2;

  Box2Abb[830]=3. + Box2Abb[179]*m_z2k;

  Box2Abb[831]=-1. + 4.*m_z2k;

  Box2Abb[832]=-1. + Box2Abb[187]*m_z2k;

  Box2Abb[833]=3. + 2.*m_z2k;

  Box2Abb[834]=-19. + 5.*Box2Abb[833]*m_z2k;

  Box2Abb[835]=8. + Box2Abb[834]*m_z2k;

  Box2Abb[836]=2.*pow(Box2Abb[61],2.)*Box2Abb[830]*m_z12 + Box2Abb[61]*Box2Abb[835]*m_z12_2 + 2.*Box2Abb[831]*Box2Abb[832]*m_z12_3 - 2.*Box2Abb[230]*pow(Box2Abb[61],3.)*m_z2k;

  Box2Abb[837]=Box2Abb[813]*m_x_2 - Box2Abb[821]*m_x_3 + Box2Abb[829]*m_x_4 + Box2Abb[823]*m_x_5 + Box2Abb[806]*m_x_6 + 4.*m_x_7 - Box2Abb[836]*m_x*m_z2k + pow(Box2Abb[61],3.)*Box2Abb[668]*m_z12_2*m_z2k_2;

  Box2Abb[838]=-1. + Box2Abb[254]*m_z12;

  Box2Abb[839]=1. + Box2Abb[155]*m_z12;

  Box2Abb[840]=-1. + 3.*Box2Abb[4]*m_z12 + m_z12*m_z2k - 3.*Box2Abb[839]*m_z2k_2 - 10.*m_z12*m_z2k_3;

  Box2Abb[841]=-3. + m_z2k;

  Box2Abb[842]=3.*Box2Abb[12] + m_z12 + Box2Abb[841]*m_z12_2 + 10.*m_z12*m_z2k_2;

  Box2Abb[843]=2. + m_z2k_2 - 8.*m_z2k_3 + 5.*m_z2k_4;

  Box2Abb[844]=-5. + 11.*m_z2k;

  Box2Abb[845]=-1. + Box2Abb[844]*m_z2k;

  Box2Abb[846]=-1. + Box2Abb[845]*m_z2k;

  Box2Abb[847]=pow(Box2Abb[61],3.) + Box2Abb[843]*m_z12 + Box2Abb[846]*m_z12_2 + 2.*m_z12_3*m_z2k_2;

  Box2Abb[848]=Box2Abb[847]*m_x + Box2Abb[840]*m_x_2 + Box2Abb[842]*m_x_3 + Box2Abb[838]*m_x_4 + m_x_5*m_z12 - Box2Abb[519]*Box2Abb[61]*m_z12*m_z2k_2;

  Box2Abb[849]=Box2Abb[837]*m_cL - 2.*Box2Abb[8]*Box2Abb[848]*m_cR*m_x;

  Box2Abb[850]=1. - 2.*m_z12_2 + 6.*m_z12*m_z2k;

  Box2Abb[851]=4. - 15.*m_z2k;

  Box2Abb[852]=-5. + Box2Abb[559]*m_z12 + m_z12_2 + Box2Abb[851]*m_z2k;

  Box2Abb[853]=-2. + Box2Abb[852]*m_z12 - 4.*m_z2k;

  Box2Abb[854]=-3. + 5.*m_z2k;

  Box2Abb[855]=1. + 2.*m_z2k_2;

  Box2Abb[856]=3. - 4.*m_z2k + 6.*m_z2k_2;

  Box2Abb[857]=-3. + Box2Abb[856]*m_z2k;

  Box2Abb[858]=pow(Box2Abb[61],4.) + pow(Box2Abb[61],2.)*Box2Abb[857]*m_z12 + Box2Abb[61]*Box2Abb[854]*Box2Abb[855]*m_z12_2 + Box2Abb[470]*Box2Abb[471]*m_z12_3;

  Box2Abb[859]=3. - 4.*m_z2k;

  Box2Abb[860]=15. + 4.*Box2Abb[258]*m_z2k;

  Box2Abb[861]=5. + Box2Abb[860]*m_z2k;

  Box2Abb[862]=-4. + Box2Abb[861]*m_z12 + Box2Abb[859]*m_z12_2 - 3.*Box2Abb[12]*m_z12_3 + 6.*m_z2k_2;

  Box2Abb[863]=3. + m_z2k_2;

  Box2Abb[864]=5. + 2.*Box2Abb[258]*m_z2k;

  Box2Abb[865]=7. + Box2Abb[864]*m_z2k;

  Box2Abb[866]=-8. + 5.*m_z2k;

  Box2Abb[867]=7. + Box2Abb[866]*m_z2k;

  Box2Abb[868]=4. - 3.*Box2Abb[867]*m_z2k;

  Box2Abb[869]=4. + Box2Abb[868]*m_z2k;

  Box2Abb[870]=Box2Abb[869]*m_z12 - Box2Abb[865]*m_z12_2 + Box2Abb[863]*m_z12_3 - 2.*Box2Abb[230]*Box2Abb[61]*m_z2k;

  Box2Abb[871]=Box2Abb[858]*m_x + Box2Abb[870]*m_x_2 + Box2Abb[862]*m_x_3 + Box2Abb[853]*m_x_4 + Box2Abb[850]*m_x_5 - m_x_6*m_z12 - Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[65]*m_z12*m_z2k_2;

  Box2Abb[872]=-10. + m_z12;

  Box2Abb[873]=-2. + Box2Abb[872]*m_z12 - 24.*m_z2k;

  Box2Abb[874]=22. + 9.*m_z12;

  Box2Abb[875]=40. - Box2Abb[874]*m_z12;

  Box2Abb[876]=17. + m_z12;

  Box2Abb[877]=-14. + Box2Abb[876]*m_z12;

  Box2Abb[878]=7. + Box2Abb[877]*m_z12;

  Box2Abb[879]=2. - 5.*m_z12;

  Box2Abb[880]=-8. + 3.*Box2Abb[879]*m_z12;

  Box2Abb[881]=30. + Box2Abb[880]*m_z12;

  Box2Abb[882]=-13. + 11.*m_z12;

  Box2Abb[883]=25. + Box2Abb[882]*m_z12;

  Box2Abb[884]=-12.*Box2Abb[12] + Box2Abb[875]*m_z12 - 2.*Box2Abb[878]*m_z12*m_z2k + 2.*Box2Abb[881]*m_z2k_2 - 4.*Box2Abb[883]*m_z2k_3 + 60.*m_z2k_4;

  Box2Abb[885]=-2. + 2.*Box2Abb[12]*m_z12 + m_z2k + m_z2k_2;

  Box2Abb[886]=1. - 6.*m_z2k;

  Box2Abb[887]=7. + 5.*m_z2k;

  Box2Abb[888]=10. + 21.*m_z2k;

  Box2Abb[889]=12. - 2.*Box2Abb[888]*m_z12 + 2.*Box2Abb[887]*m_z12_2 + m_z12_3 + 10.*Box2Abb[886]*m_z2k;

  Box2Abb[890]=1. + m_z2k + m_z2k_2;

  Box2Abb[891]=2. + m_z2k - 17.*m_z2k_2 + 18.*m_z2k_3 + m_z2k_4 - 5.*m_z2k_5;

  Box2Abb[892]=3. + Box2Abb[150]*m_z2k;

  Box2Abb[893]=2. + 17.*m_z2k;

  Box2Abb[894]=-2. + Box2Abb[893]*m_z2k;

  Box2Abb[895]=2.*pow(Box2Abb[61],3.)*Box2Abb[890]*m_z12 + 2.*Box2Abb[891]*m_z12_2 - Box2Abb[12]*Box2Abb[61]*Box2Abb[894]*m_z12_3 + 2.*Box2Abb[230]*pow(Box2Abb[61],4.)*m_z2k - 2.*Box2Abb[892]*m_z12_4*m_z2k;

  Box2Abb[896]=5. + 11.*m_z2k;

  Box2Abb[897]=1. + 5.*Box2Abb[859]*m_z2k;

  Box2Abb[898]=6. + 17.*m_z2k;

  Box2Abb[899]=17. + 2.*Box2Abb[898]*m_z2k;

  Box2Abb[900]=14. + 31.*m_z2k;

  Box2Abb[901]=32. + Box2Abb[900]*m_z2k;

  Box2Abb[902]=10. - 2.*Box2Abb[899]*m_z12 + Box2Abb[901]*m_z12_2 + Box2Abb[896]*m_z12_3 + 4.*Box2Abb[897]*m_z2k;

  Box2Abb[903]=2. + Box2Abb[831]*m_z2k;

  Box2Abb[904]=22. - 3.*m_z2k + 9.*m_z2k_2;

  Box2Abb[905]=-7. + Box2Abb[904]*m_z2k;

  Box2Abb[906]=25. + 4.*m_z2k + 34.*m_z2k_2;

  Box2Abb[907]=7. + Box2Abb[906]*m_z2k;

  Box2Abb[908]=-12. + 31.*m_z2k;

  Box2Abb[909]=2. + Box2Abb[908]*m_z2k;

  Box2Abb[910]=-48. + Box2Abb[909]*m_z2k;

  Box2Abb[911]=-1. + Box2Abb[910]*m_z2k;

  Box2Abb[912]=-2.*Box2Abb[514]*pow(Box2Abb[61],2.)*Box2Abb[903] - 2.*Box2Abb[61]*Box2Abb[905]*m_z12 + Box2Abb[911]*m_z12_2 + Box2Abb[907]*m_z12_3 + 6.*Box2Abb[12]*m_z12_4*m_z2k;

  Box2Abb[913]=Box2Abb[895]*m_x + Box2Abb[912]*m_x_2 + Box2Abb[884]*m_x_3 + Box2Abb[902]*m_x_4 - Box2Abb[889]*m_x_5 + Box2Abb[873]*m_x_6 + 4.*m_x_7 + Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[885]*m_z12_2*m_z2k;

  Box2Abb[914]=Box2Abb[913]*m_cL + 2.*Box2Abb[871]*m_cR*m_x;

  Box2Abb[938]=2./m_s;

  Box2Abb[939]=-2./m_s;

  Box2Abb[940]=(-2.*m_z2k)/m_s;

  Box2Abb[941]=(2.*m_x)/m_s;

  Box2Abb[942]=Box2Abb[672]*m_cR + 2.*Box2Abb[667]*m_cL*m_x;

  Box2Abb[943]=Box2Abb[690]*m_cR + 2.*Box2Abb[676]*Box2Abb[8]*m_cL*m_x*m_z12;

  Box2Abb[944]=Box2Abb[697]*m_cR + 2.*Box2Abb[382]*Box2Abb[4]*m_cL*m_x;

  Box2Abb[945]=Box2Abb[705]*m_cR + 2.*Box2Abb[388]*Box2Abb[4]*m_cL*m_x;

  Box2Abb[946]=Box2Abb[712]*m_cR - 2.*Box2Abb[4]*Box2Abb[707]*m_cL*m_x;

  Box2Abb[947]=Box2Abb[721]*m_cR - 2.*Box2Abb[4]*Box2Abb[433]*m_cL*m_x;

  Box2Abb[948]=Box2Abb[728]*m_cR - 2.*Box2Abb[676]*Box2Abb[8]*m_cL*m_x*m_z12;

  Box2Abb[949]=-Box2Abb[747]*m_cR + 2.*Box2Abb[148]*Box2Abb[730]*m_cL*m_x*m_z12;

  Box2Abb[950]=Box2Abb[751]*m_cR + 2.*Box2Abb[10]*m_cL*m_x;

  Box2Abb[951]=Box2Abb[766]*m_cR + 2.*Box2Abb[754]*Box2Abb[8]*m_cL*m_x;

  Box2Abb[952]=2.*Box2Abb[781]*m_cL + Box2Abb[778]*m_cR;

  Box2Abb[953]=43. + 40.*m_z2k;

  Box2Abb[954]=-9. - 2.*m_z2k + 4.*m_z2k_2;

  Box2Abb[955]=-18. + 20.*m_z12 + 9.*Box2Abb[444]*m_z12_2 + 2.*Box2Abb[954]*m_z2k + Box2Abb[953]*m_z12*m_z2k;

  Box2Abb[956]=-Box2Abb[803]*m_x + Box2Abb[799]*m_x_2 - Box2Abb[955]*m_x_3 + Box2Abb[794]*m_x_4 - 2.*Box2Abb[788]*m_x_5 - 4.*m_x_6 + 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k;

  Box2Abb[957]=Box2Abb[956]*m_cR - 2.*Box2Abb[787]*Box2Abb[8]*m_cL*m_x;

  Box2Abb[958]=Box2Abb[837]*m_cR - 2.*Box2Abb[8]*Box2Abb[848]*m_cL*m_x;

  Box2Abb[959]=2. - Box2Abb[872]*m_z12 + 24.*m_z2k;

  Box2Abb[960]=-40. + Box2Abb[874]*m_z12;

  Box2Abb[961]=8. + 3.*Box2Abb[550]*m_z12;

  Box2Abb[962]=-30. + Box2Abb[961]*m_z12;

  Box2Abb[963]=12.*Box2Abb[12] + Box2Abb[960]*m_z12 + 2.*Box2Abb[878]*m_z12*m_z2k + 2.*Box2Abb[962]*m_z2k_2 + 4.*Box2Abb[883]*m_z2k_3 - 60.*m_z2k_4;

  Box2Abb[964]=-5. + Box2Abb[291]*m_z2k;

  Box2Abb[965]=-2. + Box2Abb[964]*m_z2k;

  Box2Abb[966]=-2.*pow(Box2Abb[61],3.)*Box2Abb[890]*m_z12 + 2.*pow(Box2Abb[61],2.)*Box2Abb[965]*m_z12_2 + Box2Abb[12]*Box2Abb[61]*Box2Abb[894]*m_z12_3 - 2.*Box2Abb[230]*pow(Box2Abb[61],4.)*m_z2k + 2.*Box2Abb[892]*m_z12_4*m_z2k;

  Box2Abb[967]=12. - 31.*m_z2k;

  Box2Abb[968]=-2. + Box2Abb[967]*m_z2k;

  Box2Abb[969]=48. + Box2Abb[968]*m_z2k;

  Box2Abb[970]=1. + Box2Abb[969]*m_z2k;

  Box2Abb[971]=2.*Box2Abb[514]*pow(Box2Abb[61],2.)*Box2Abb[903] + 2.*Box2Abb[61]*Box2Abb[905]*m_z12 + Box2Abb[970]*m_z12_2 - Box2Abb[907]*m_z12_3 - 6.*Box2Abb[12]*m_z12_4*m_z2k;

  Box2Abb[972]=Box2Abb[966]*m_x + Box2Abb[971]*m_x_2 + Box2Abb[963]*m_x_3 - Box2Abb[902]*m_x_4 + Box2Abb[889]*m_x_5 + Box2Abb[959]*m_x_6 - 4.*m_x_7 - Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[885]*m_z12_2*m_z2k;

  Box2Abb[973]=m_z12 - 3.*m_z2k;

  Box2Abb[974]=-1. + 2.*Box2Abb[973]*m_z12;

  Box2Abb[975]=-5. + 3.*Box2Abb[4]*m_z12;

  Box2Abb[976]=-3. + 8.*m_z12;

  Box2Abb[977]=4. + Box2Abb[975]*m_z12 + Box2Abb[192]*Box2Abb[242]*m_z12*m_z2k + 2.*Box2Abb[976]*m_z2k_2 - 20.*m_z12*m_z2k_3;

  Box2Abb[978]=-pow(Box2Abb[61],4.) - pow(Box2Abb[61],2.)*Box2Abb[857]*m_z12 - Box2Abb[61]*Box2Abb[854]*Box2Abb[855]*m_z12_2 - Box2Abb[470]*Box2Abb[471]*m_z12_3;

  Box2Abb[979]=-4. + 3.*Box2Abb[867]*m_z2k;

  Box2Abb[980]=-4. + Box2Abb[979]*m_z2k;

  Box2Abb[981]=Box2Abb[980]*m_z12 + Box2Abb[865]*m_z12_2 - Box2Abb[863]*m_z12_3 + 2.*Box2Abb[230]*Box2Abb[61]*m_z2k;

  Box2Abb[982]=Box2Abb[978]*m_x + Box2Abb[981]*m_x_2 + Box2Abb[977]*m_x_3 - Box2Abb[853]*m_x_4 + Box2Abb[974]*m_x_5 + m_x_6*m_z12 + Box2Abb[5]*pow(Box2Abb[61],2.)*Box2Abb[65]*m_z12*m_z2k_2;

  Box2Abb[983]=-Box2Abb[972]*m_cR - 2.*Box2Abb[982]*m_cL*m_x;

  Box2Abb[1007]=2. - 3.*m_z12 - 2.*m_z2k;

  Box2Abb[1008]=m_z12 + 4.*m_z2k;

  Box2Abb[1009]=Box2Abb[1008]*m_x - 2.*m_x_2 + Box2Abb[1007]*m_z2k;

  Box2Abb[1010]=1. + m_z12 + 4.*m_z2k;

  Box2Abb[1011]=m_z12 - 2.*m_z2k + 6.*m_z12*m_z2k + 4.*m_z2k_2;

  Box2Abb[1012]=Box2Abb[1011]*m_x - 2.*Box2Abb[1010]*m_x_2 + 4.*m_x_3 - m_z12_2*m_z2k;

  Box2Abb[1013]=-1. + m_z12 + 4.*m_z2k;

  Box2Abb[1014]=2. - 2.*Box2Abb[1013]*m_x + 4.*m_x_2 - 3.*m_z12 + m_z12_2 + 6.*Box2Abb[4]*m_z2k + 4.*m_z2k_2;

  Box2Abb[1015]=2. + Box2Abb[201]*m_z12 - 4.*m_z2k;

  Box2Abb[1016]=Box2Abb[1015]*m_x + 4.*m_x_2 + m_z12_2*m_z2k;

  Box2Abb[1017]=-1. + m_z12 + m_z12_2 - 3.*m_z12*m_z2k;

  Box2Abb[1018]=-1. + Box2Abb[220]*m_z12;

  Box2Abb[1019]=Box2Abb[4]*m_z12 + Box2Abb[1018]*m_z2k + 2.*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[1020]=-10. + 10.*m_z12 + m_z12_2 + 4.*Box2Abb[70]*m_z2k - 12.*m_z2k_2;

  Box2Abb[1021]=Box2Abb[1020]*m_z12 - 4.*m_z2k;

  Box2Abb[1022]=-Box2Abb[1021]*m_x_2 + 4.*Box2Abb[1017]*m_x_3 + 2.*Box2Abb[1019]*m_x*m_z12 + 4.*m_x_4*m_z12 - Box2Abb[668]*m_z12_3*m_z2k;

  Box2Abb[1023]=2.*m_x - m_z12;

  Box2Abb[1024]=-2. + Box2Abb[769]*m_z12;

  Box2Abb[1025]=4. - 8.*m_z12;

  Box2Abb[1026]=Box2Abb[1025]*m_x + 8.*m_x_2 + m_z12_2;

  Box2Abb[1027]=Box2Abb[1023]*Box2Abb[1024]*m_x - Box2Abb[1026]*m_z12*m_z2k + 4.*m_x*m_z12*m_z2k_2;

  Box2Abb[1028]=6. + Box2Abb[815]*m_z12 - 6.*m_z2k;

  Box2Abb[1029]=-2. + Box2Abb[1028]*m_z12;

  Box2Abb[1030]=-1. + m_z12 + Box2Abb[4]*Box2Abb[72]*m_z2k - 2.*m_z2k_2;

  Box2Abb[1031]=Box2Abb[4]*pow(Box2Abb[72],2.) + 4.*pow(Box2Abb[72],2.)*m_z2k - 12.*m_z2k_2;

  Box2Abb[1032]=Box2Abb[1031]*m_z12 - 4.*m_z2k;

  Box2Abb[1033]=Box2Abb[1032]*m_x_2 - 2.*Box2Abb[1029]*m_x_3 - 4.*m_x_4*m_z12 - 2.*Box2Abb[1030]*m_x*m_z12*m_z2k + Box2Abb[4]*m_z12_3*m_z2k_2;

  Box2Abb[1034]=2.*m_x + m_z12;

  Box2Abb[1035]=2. + Box2Abb[1034]*m_x - 2.*m_z12;

  Box2Abb[1036]=2. + m_z12_2;

  Box2Abb[1037]=Box2Abb[1036]*m_x - 6.*Box2Abb[72]*m_x_2 + 2.*Box2Abb[4]*Box2Abb[538]*m_z12;

  Box2Abb[1038]=-2. + Box2Abb[699]*m_z12;

  Box2Abb[1039]=Box2Abb[1038]*m_x + 6.*Box2Abb[72]*m_x_2 + 2.*Box2Abb[4]*m_z12_2;

  Box2Abb[1040]=-2.*Box2Abb[72]*m_x + m_z12_2;

  Box2Abb[1041]=Box2Abb[1035]*Box2Abb[72]*m_x_2 + Box2Abb[1037]*m_x*m_z2k + Box2Abb[1039]*m_z2k_2 + Box2Abb[1040]*m_z2k_3;

  Box2Abb[1042]=-3. + m_z12 + m_z2k;

  Box2Abb[1043]=1. + Box2Abb[1042]*m_z12 - 2.*m_z2k;

  Box2Abb[1044]=-Box2Abb[1008]*Box2Abb[72]*m_x_2 + 2.*Box2Abb[72]*m_x_3 + 2.*Box2Abb[1043]*m_x*m_z2k + m_z12_2*m_z2k_2;

  Box2Abb[1045]=2. + m_z12 + 2.*m_z12_2 + 8.*m_z12*m_z2k;

  Box2Abb[1046]=5. - 2.*m_z12;

  Box2Abb[1047]=-8. + Box2Abb[1046]*m_z12;

  Box2Abb[1048]=13. - 10.*m_z12;

  Box2Abb[1049]=-2. + Box2Abb[1048]*m_z12;

  Box2Abb[1050]=1. + Box2Abb[538]*m_z12 + m_z2k + Box2Abb[1047]*m_z12*m_z2k + Box2Abb[1049]*m_z2k_2 - 8.*m_z12*m_z2k_3;

  Box2Abb[1051]=Box2Abb[538]*m_z12 + 2.*m_z12*m_z2k + m_z2k_2;

  Box2Abb[1052]=3. + 2.*Box2Abb[1051] - 5.*m_z2k;

  Box2Abb[1053]=-5. + 12.*m_z2k;

  Box2Abb[1054]=m_z12 + Box2Abb[1053]*m_z2k + 8.*m_z12*m_z2k;

  Box2Abb[1055]=3. + Box2Abb[1054]*m_z12 + 4.*m_z2k;

  Box2Abb[1056]=Box2Abb[1050]*m_x + Box2Abb[1055]*m_x_2 - Box2Abb[1045]*m_x_3 + 2.*m_x_4*m_z12 + Box2Abb[1052]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1057]=2. + m_z12 + 4.*m_z2k;

  Box2Abb[1058]=5. + 2.*m_z12;

  Box2Abb[1059]=1. - 2.*m_z12 + m_z2k - Box2Abb[1058]*m_z12*m_z2k + 2.*Box2Abb[879]*m_z2k_2 - 8.*m_z2k_3;

  Box2Abb[1060]=Box2Abb[183]*Box2Abb[61]*m_z12 + 2.*Box2Abb[12]*m_z12_2 + 2.*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[1061]=4. + 8.*m_z2k;

  Box2Abb[1062]=1. + Box2Abb[1061]*m_z12 + 4.*Box2Abb[256]*m_z2k;

  Box2Abb[1063]=Box2Abb[1059]*m_x + Box2Abb[1062]*m_x_2 - 2.*Box2Abb[1057]*m_x_3 + 2.*m_x_4 + Box2Abb[1060]*m_z2k;

  Box2Abb[1064]=1. - m_z12 + m_z2k + Box2Abb[62]*m_z12*m_z2k + 10.*Box2Abb[4]*m_z2k_2 + 8.*m_z2k_3;

  Box2Abb[1065]=-1. + 6.*m_z2k;

  Box2Abb[1066]=-1. + m_z12 + 2.*Box2Abb[1065]*m_z2k + 8.*m_z12*m_z2k;

  Box2Abb[1067]=-Box2Abb[1064]*m_x + Box2Abb[1066]*m_x_2 - 2.*Box2Abb[1010]*m_x_3 + 2.*m_x_4 + 2.*pow(Box2Abb[5],2.)*Box2Abb[61]*m_z2k;

  Box2Abb[1068]=2. + m_z12_2 + 8.*m_z12*m_z2k;

  Box2Abb[1069]=-2. + 5.*m_z2k;

  Box2Abb[1070]=2. + 2.*Box2Abb[1069]*m_z12 + 5.*m_z12_2 + 4.*Box2Abb[61]*m_z2k;

  Box2Abb[1071]=-2. + 7.*m_z2k;

  Box2Abb[1072]=8. + 2.*Box2Abb[1071]*m_z12 + m_z12_2 + 4.*Box2Abb[1065]*m_z2k;

  Box2Abb[1073]=-4. + Box2Abb[1072]*m_z12 + 8.*m_z2k;

  Box2Abb[1074]=4. + 22.*m_z2k;

  Box2Abb[1075]=-2. + Box2Abb[1074]*m_z12 + m_z12_2 + 8.*Box2Abb[230]*m_z2k;

  Box2Abb[1076]=Box2Abb[1075]*m_z12 + 4.*m_z2k;

  Box2Abb[1077]=-Box2Abb[1073]*m_x_3 + 2.*Box2Abb[1068]*m_x_4 - 4.*m_x_5*m_z12 + Box2Abb[1076]*m_x_2*m_z2k - Box2Abb[1070]*m_x*m_z12*m_z2k_2 + m_z12_3*m_z2k_3;

  Box2Abb[1078]=Box2Abb[444]*m_x - m_x_2 - Box2Abb[65]*m_z2k;

  Box2Abb[1079]=-Box2Abb[1013]*m_x + m_x_2*m_z12 - Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1080]=2. - 2.*m_x + 3.*m_z12;

  Box2Abb[1081]=2. + Box2Abb[1080]*m_x - 3.*m_z12;

  Box2Abb[1082]=2. - 8.*m_z12;

  Box2Abb[1083]=2. + Box2Abb[1082]*m_x + 6.*m_x_2 + Box2Abb[817]*m_z12;

  Box2Abb[1084]=-4. - 6.*m_x + 5.*m_z12;

  Box2Abb[1085]=Box2Abb[1081]*m_x + Box2Abb[1083]*m_z2k + Box2Abb[1084]*m_z2k_2 + 2.*m_z2k_3;

  Box2Abb[1086]=1. - 2.*m_z12 + 6.*m_z2k_2;

  Box2Abb[1087]=6. - 5.*m_z12;

  Box2Abb[1088]=-1. + m_z12 + m_z2k + 2.*Box2Abb[4]*m_z12*m_z2k + Box2Abb[1087]*m_z2k_2 - 6.*m_z2k_3;

  Box2Abb[1089]=Box2Abb[1088]*m_x + Box2Abb[1086]*m_x_2 + Box2Abb[199]*m_x_3 + 2.*pow(Box2Abb[5],2.)*Box2Abb[61]*m_z2k;

  Box2Abb[1090]=-2. + m_x;

  Box2Abb[1091]=-2.*Box2Abb[1090]*m_x + Box2Abb[1]*Box2Abb[606]*m_z12;

  Box2Abb[1092]=m_x - 3.*Box2Abb[210]*m_z12 + 2.*m_z12_2;

  Box2Abb[1093]=-2. + 2.*Box2Abb[1092]*m_x + 3.*m_z12;

  Box2Abb[1094]=2. + 6.*m_z12;

  Box2Abb[1095]=4. + Box2Abb[1094]*m_x - Box2Abb[703]*m_z12;

  Box2Abb[1096]=Box2Abb[1091]*m_x + Box2Abb[1093]*m_z2k + Box2Abb[1095]*m_z2k_2 - 2.*Box2Abb[9]*m_z2k_3;

  Box2Abb[1097]=3. + 4.*m_z12 + 4.*m_z2k;

  Box2Abb[1098]=15. + 2.*m_z12 + 16.*m_z2k;

  Box2Abb[1099]=-2.*Box2Abb[187] + Box2Abb[1098]*m_z12;

  Box2Abb[1100]=3. + Box2Abb[725]*m_z2k;

  Box2Abb[1101]=-2.*Box2Abb[1100]*Box2Abb[12] + 5.*m_z12 + 4.*Box2Abb[256]*m_z12_2;

  Box2Abb[1102]=2. - 8.*Box2Abb[61]*m_z2k;

  Box2Abb[1103]=-9. + 8.*m_z2k;

  Box2Abb[1104]=-2. + Box2Abb[1103]*m_z2k;

  Box2Abb[1105]=-Box2Abb[1104]*Box2Abb[61]*m_z12 + Box2Abb[1102]*m_z12_2 - 2.*Box2Abb[230]*pow(Box2Abb[61],2.)*m_z2k + 2.*m_z12_3*m_z2k;

  Box2Abb[1106]=Box2Abb[1105]*m_x - Box2Abb[1101]*m_x_2 + Box2Abb[1099]*m_x_3 - 2.*Box2Abb[1097]*m_x_4 + 4.*m_x_5 + 2.*Box2Abb[5]*Box2Abb[61]*m_z12_2*m_z2k;

  Box2Abb[1107]=1. + 6.*m_z12 + 4.*m_z2k;

  Box2Abb[1108]=-2.*Box2Abb[230]*pow(Box2Abb[61],2.) + Box2Abb[61]*m_z12 + 4.*m_z12_2*m_z2k;

  Box2Abb[1109]=13. + 2.*m_z12 + 28.*m_z2k;

  Box2Abb[1110]=-2.*Box2Abb[488] + Box2Abb[1109]*m_z12;

  Box2Abb[1111]=1. + 22.*m_z2k;

  Box2Abb[1112]=1. + 6.*Box2Abb[179]*m_z2k;

  Box2Abb[1113]=-1. + m_z2k + 4.*m_z2k_2;

  Box2Abb[1114]=-2.*Box2Abb[1113]*Box2Abb[61] + 2.*Box2Abb[1112]*m_z12 + Box2Abb[1111]*m_z12_2;

  Box2Abb[1115]=Box2Abb[1108]*Box2Abb[5] - Box2Abb[1114]*m_x + Box2Abb[1110]*m_x_2 - 2.*Box2Abb[1107]*m_x_3 + 4.*m_x_4;

  Box2Abb[1116]=-1. + 5.*m_z12;

  Box2Abb[1117]=1. + Box2Abb[1116]*m_z12;

  Box2Abb[1118]=4. + 12.*m_z12 + 14.*m_z12_2 - 3.*m_z12_3 + 2.*Box2Abb[157]*Box2Abb[783]*m_z2k + 12.*Box2Abb[1117]*m_z2k_2 + 40.*m_z12*m_z2k_3;

  Box2Abb[1119]=6. + 3.*m_z12 + 10.*m_z2k;

  Box2Abb[1120]=2. + Box2Abb[1119]*m_z12;

  Box2Abb[1121]=5. + 8.*m_z12;

  Box2Abb[1122]=12. + 14.*m_z12 - m_z12_2 + 4.*Box2Abb[1121]*m_z2k + 40.*m_z2k_2;

  Box2Abb[1123]=8. + Box2Abb[1122]*m_z12 + 12.*m_z2k;

  Box2Abb[1124]=2. + Box2Abb[424]*m_z2k;

  Box2Abb[1125]=-5. + 18.*m_z2k;

  Box2Abb[1126]=1. + Box2Abb[1125]*m_z2k;

  Box2Abb[1127]=-2.*Box2Abb[512]*pow(Box2Abb[61],3.) - 2.*Box2Abb[1124]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[1126]*Box2Abb[61]*m_z12_2 - 8.*m_z12_3*m_z2k_2;

  Box2Abb[1128]=19. + 34.*m_z2k;

  Box2Abb[1129]=-3. + Box2Abb[1128]*m_z2k;

  Box2Abb[1130]=-13. + 24.*m_z2k;

  Box2Abb[1131]=-14. + Box2Abb[1130]*m_z2k;

  Box2Abb[1132]=5. + Box2Abb[1131]*m_z2k;

  Box2Abb[1133]=-9. + 5.*m_z2k;

  Box2Abb[1134]=1. + Box2Abb[1133]*m_z2k;

  Box2Abb[1135]=3. + 2.*Box2Abb[1134]*m_z2k;

  Box2Abb[1136]=6. + 2.*Box2Abb[1135]*m_z2k;

  Box2Abb[1137]=4.*pow(Box2Abb[61],3.) + Box2Abb[1136]*m_z12 + 2.*Box2Abb[1132]*m_z12_2 + Box2Abb[1129]*m_z12_3;

  Box2Abb[1138]=Box2Abb[1137]*m_x_2 - Box2Abb[1118]*m_x_3 + Box2Abb[1123]*m_x_4 - 2.*Box2Abb[1120]*m_x_5 + Box2Abb[1127]*m_x*m_z12 + 4.*m_x_6*m_z12 - pow(Box2Abb[61],3.)*m_z12_3*m_z2k;

  Box2Abb[1139]=-2. + 2.*m_z12 + m_z2k + m_z2k_2;

  Box2Abb[1140]=9. + 23.*m_z12;

  Box2Abb[1141]=-21. + Box2Abb[1140]*m_z12;

  Box2Abb[1142]=3. + Box2Abb[1141]*m_z12;

  Box2Abb[1143]=5. + 14.*m_z12;

  Box2Abb[1144]=-3. + 2.*Box2Abb[9]*m_z12;

  Box2Abb[1145]=-4. + 29.*m_z12;

  Box2Abb[1146]=9. + Box2Abb[1145]*m_z12;

  Box2Abb[1147]=6. - 32.*m_z12 + 27.*m_z12_2 + 2.*Box2Abb[1142]*m_z2k + 2.*Box2Abb[1143]*Box2Abb[1144]*m_z2k_2 + 4.*Box2Abb[1146]*m_z2k_3 + 30.*Box2Abb[72]*m_z2k_4;

  Box2Abb[1148]=8. + 9.*m_z12 + 12.*m_z2k;

  Box2Abb[1149]=12. - Box2Abb[1148]*m_z12 + 24.*m_z2k;

  Box2Abb[1150]=33. + 50.*m_z2k;

  Box2Abb[1151]=4. + 5.*m_z2k;

  Box2Abb[1152]=-2. + 6.*Box2Abb[1151]*m_z2k;

  Box2Abb[1153]=17. + 30.*m_z2k;

  Box2Abb[1154]=8. + Box2Abb[1153]*m_z2k;

  Box2Abb[1155]=-2.*Box2Abb[1154] + Box2Abb[1152]*m_z12 + Box2Abb[1150]*m_z12_2;

  Box2Abb[1156]=3. + 8.*m_z2k_2;

  Box2Abb[1157]=95. + 109.*m_z2k;

  Box2Abb[1158]=45. + Box2Abb[1157]*m_z2k;

  Box2Abb[1159]=-5. + Box2Abb[450]*m_z2k;

  Box2Abb[1160]=30. - 8.*Box2Abb[1159]*m_z2k;

  Box2Abb[1161]=2.*Box2Abb[1156]*Box2Abb[510] + Box2Abb[1160]*m_z12 - Box2Abb[1158]*m_z12_2 - 18.*m_z12_3*m_z2k;

  Box2Abb[1162]=3. + Box2Abb[179]*m_z2k_2;

  Box2Abb[1163]=-3. + 8.*m_z2k;

  Box2Abb[1164]=-8. + Box2Abb[1163]*m_z2k;

  Box2Abb[1165]=5. + Box2Abb[1164]*m_z2k;

  Box2Abb[1166]=-3. + 10.*m_z2k;

  Box2Abb[1167]=-5. + Box2Abb[1166]*m_z2k;

  Box2Abb[1168]=16. + Box2Abb[1167]*m_z2k;

  Box2Abb[1169]=2.*Box2Abb[1162]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[1168]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[1165]*m_z12_3 - 2.*Box2Abb[230]*pow(Box2Abb[61],3.)*m_z2k + 4.*Box2Abb[833]*m_z12_4*m_z2k;

  Box2Abb[1170]=19. + 7.*Box2Abb[234]*m_z2k;

  Box2Abb[1171]=43. + 32.*m_z2k - 59.*m_z2k_3;

  Box2Abb[1172]=-6. + Box2Abb[1171]*m_z2k;

  Box2Abb[1173]=-5. + 6.*m_z2k;

  Box2Abb[1174]=7. + 2.*Box2Abb[1173]*m_z2k;

  Box2Abb[1175]=-4. + Box2Abb[1174]*m_z2k;

  Box2Abb[1176]=2. + Box2Abb[1175]*m_z2k;

  Box2Abb[1177]=4. - 3.*Box2Abb[179]*m_z2k;

  Box2Abb[1178]=5. + 2.*Box2Abb[1177]*m_z2k;

  Box2Abb[1179]=-10. + Box2Abb[1178]*m_z2k;

  Box2Abb[1180]=5. + Box2Abb[1179]*m_z2k;

  Box2Abb[1181]=2.*Box2Abb[1176]*Box2Abb[61] + 2.*Box2Abb[1180]*m_z12 + Box2Abb[1172]*m_z12_2 - 2.*Box2Abb[1170]*m_z12_3*m_z2k - 8.*m_z12_4*m_z2k_2;

  Box2Abb[1182]=-Box2Abb[1181]*m_x_2 - Box2Abb[1147]*m_x_3 - Box2Abb[1161]*m_x_4 - Box2Abb[1155]*m_x_5 - Box2Abb[1149]*m_x_6 - 2.*Box2Abb[72]*m_x_7 - Box2Abb[1169]*m_x*m_z2k - Box2Abb[1139]*Box2Abb[61]*Box2Abb[65]*m_z12_2*m_z2k_2;

  Box2Abb[1183]=2. - 5.*m_z12 - 10.*m_z2k;

  Box2Abb[1184]=8. + Box2Abb[1183]*m_z12 + 20.*m_z2k;

  Box2Abb[1185]=20. + Box2Abb[872]*m_z12;

  Box2Abb[1186]=-8. + 9.*m_z12;

  Box2Abb[1187]=-2. + Box2Abb[1186]*m_z12;

  Box2Abb[1188]=-Box2Abb[1185]*m_z12 + 3.*Box2Abb[1187]*m_z2k + 20.*Box2Abb[72]*m_z2k_2;

  Box2Abb[1189]=2. + Box2Abb[420]*m_z2k;

  Box2Abb[1190]=19. - 20.*m_z2k;

  Box2Abb[1191]=1. + Box2Abb[1190]*m_z2k_2;

  Box2Abb[1192]=-8. + 15.*m_z2k;

  Box2Abb[1193]=3. + Box2Abb[1192]*m_z2k;

  Box2Abb[1194]=2.*Box2Abb[230]*pow(Box2Abb[61],4.) - 2.*Box2Abb[1189]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[1193]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[1191]*m_z12_3 - 8.*m_z12_4*m_z2k_2;

  Box2Abb[1195]=3. - 12.*m_z2k;

  Box2Abb[1196]=4. + m_z2k + 26.*m_z2k_2;

  Box2Abb[1197]=-2. + 5.*Box2Abb[797]*m_z2k;

  Box2Abb[1198]=19. - 10.*Box2Abb[841]*m_z2k;

  Box2Abb[1199]=20. + 2.*Box2Abb[1198]*m_z2k;

  Box2Abb[1200]=-10. + Box2Abb[1199]*m_z12 - 2.*Box2Abb[1196]*m_z12_2 + Box2Abb[1195]*m_z12_3 + 2.*Box2Abb[1197]*m_z2k;

  Box2Abb[1201]=-2. + Box2Abb[1166]*m_z2k;

  Box2Abb[1202]=13. + 32.*m_z2k;

  Box2Abb[1203]=-3. + Box2Abb[1202]*m_z2k;

  Box2Abb[1204]=-23. + 5.*m_z2k;

  Box2Abb[1205]=-9. + Box2Abb[1204]*m_z2k;

  Box2Abb[1206]=4. + Box2Abb[1205]*m_z2k;

  Box2Abb[1207]=-21. + 22.*m_z2k;

  Box2Abb[1208]=-22. + Box2Abb[1207]*m_z2k;

  Box2Abb[1209]=6. + 2.*Box2Abb[1208]*m_z2k;

  Box2Abb[1210]=-2.*Box2Abb[1201]*pow(Box2Abb[61],2.) + 2.*Box2Abb[1206]*Box2Abb[61]*m_z12 + Box2Abb[1209]*m_z12_2 + Box2Abb[1203]*m_z12_3;

  Box2Abb[1211]=Box2Abb[1194]*m_x + Box2Abb[1210]*m_x_2 + Box2Abb[1200]*m_x_3 + Box2Abb[1188]*m_x_4 + Box2Abb[1184]*m_x_5 + 2.*Box2Abb[72]*m_x_6 + Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k;

  Box2Abb[1244]=(2.*m_cL*m_z12)/m_s;

  Box2Abb[1245]=(2.*m_cL)/m_s;

  Box2Abb[1246]=(2.*m_cL*m_x)/m_s;

  Box2Abb[1247]=-m_cL*m_z12;

  Box2Abb[1248]=-m_cL;

  Box2Abb[1249]=-m_cL*m_s*m_z12;

  Box2Abb[1250]=-m_cL*m_s;

  Box2Abb[1251]=(-2.*m_cL)/m_s;

  Box2Abb[1252]=(2.*m_cL*m_z2k)/m_s;

  Box2Abb[1253]=m_cL;

  Box2Abb[1286]=(2.*m_cR*m_z12)/m_s;

  Box2Abb[1287]=(2.*m_cR)/m_s;

  Box2Abb[1288]=(2.*m_cR*m_x)/m_s;

  Box2Abb[1289]=-m_cR*m_z12;

  Box2Abb[1290]=-m_cR;

  Box2Abb[1291]=-m_cR*m_s*m_z12;

  Box2Abb[1292]=-m_cR*m_s;

  Box2Abb[1293]=(-2.*m_cR)/m_s;

  Box2Abb[1294]=(2.*m_cR*m_z2k)/m_s;

  Box2Abb[1295]=m_cR;

  Box2Abb[1296]=-3. + 3.*m_z12 + m_z2k;

  Box2Abb[1297]=-4. + 2.*Box2Abb[399]*m_z12 - 3.*m_z2k;

  Box2Abb[1298]=1. - 2.*m_z12 + 6.*m_z2k - 9.*m_z12*m_z2k;

  Box2Abb[1299]=1. + Box2Abb[1298]*m_z12;

  Box2Abb[1300]=Box2Abb[1299]*m_x_2 + Box2Abb[53]*m_x_3*m_z12 + Box2Abb[1297]*m_x*m_z12*m_z2k - Box2Abb[1296]*m_z12_2*m_z2k_2;

  Box2Abb[1301]=-11. + 2.*m_z12;

  Box2Abb[1302]=-2. + 2.*pow(Box2Abb[72],2.)*m_z12 + 6.*m_z2k + Box2Abb[1301]*m_z12*m_z2k - 4.*m_z2k_2;

  Box2Abb[1303]=-9. + 4.*m_z12 + 8.*m_z2k;

  Box2Abb[1304]=14. + 3.*Box2Abb[1303]*m_z12 - 36.*m_z2k;

  Box2Abb[1305]=-4. + 11.*m_z12;

  Box2Abb[1306]=-5. + m_z12 + 4.*m_z12_2 - 2.*m_z12_3 - Box2Abb[1305]*Box2Abb[72]*m_z2k - 12.*Box2Abb[72]*m_z2k_2;

  Box2Abb[1307]=2. + Box2Abb[1306]*m_z12 - 4.*m_z2k;

  Box2Abb[1308]=Box2Abb[1307]*m_x_2 + Box2Abb[1304]*m_x_3*m_z12 + 4.*Box2Abb[505]*m_x_4*m_z12 + Box2Abb[1302]*m_x*m_z12*m_z2k + Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[1309]=1. + m_z12 - 8.*m_z12_2 + 6.*m_z12_3;

  Box2Abb[1310]=2. + m_z12 - 4.*m_z12_2;

  Box2Abb[1311]=-2.*pow(Box2Abb[4],2.)*m_z12 - 2.*Box2Abb[1309]*m_z2k + 3.*Box2Abb[1310]*m_z2k_2 - 4.*m_z2k_3;

  Box2Abb[1312]=-7. + 4.*m_z12_2 + 6.*Box2Abb[460]*m_z2k + 25.*m_z12*m_z2k;

  Box2Abb[1313]=2. + m_z12 + Box2Abb[1312]*m_z12_2 - 4.*m_z2k;

  Box2Abb[1314]=-9. + 6.*m_z12 + 8.*m_z2k;

  Box2Abb[1315]=6. + Box2Abb[1314]*m_z12 - 4.*m_z2k;

  Box2Abb[1316]=8. - 3.*Box2Abb[1315]*m_z12;

  Box2Abb[1317]=Box2Abb[1313]*m_x_2 + Box2Abb[1316]*m_x_3 + Box2Abb[1311]*m_x*m_z12 + 4.*Box2Abb[155]*m_x_4*m_z12 + Box2Abb[5]*Box2Abb[668]*m_z12_3*m_z2k;

  Box2Abb[1318]=-5. + 2.*m_z12;

  Box2Abb[1319]=3. + Box2Abb[1318]*m_z12 - 3.*m_z2k;

  Box2Abb[1320]=13. - 5.*m_z12 - 3.*m_z2k;

  Box2Abb[1321]=-9. + Box2Abb[1320]*m_z12 + 6.*m_z2k;

  Box2Abb[1322]=1. + Box2Abb[1321]*m_z12;

  Box2Abb[1323]=Box2Abb[1322]*m_x_2 + Box2Abb[1319]*Box2Abb[5]*m_x*m_z12 + Box2Abb[62]*m_x_3*m_z12 + pow(Box2Abb[5],2.)*m_z12_2*m_z2k;

  Box2Abb[1324]=-4. + 5.*m_z12;

  Box2Abb[1325]=4. + m_z12 + 7.*m_z2k;

  Box2Abb[1326]=-5. + Box2Abb[1325]*m_z12 - 4.*m_z2k;

  Box2Abb[1327]=5. - 2.*m_z12 - 11.*m_z2k;

  Box2Abb[1328]=-7. + Box2Abb[1327]*m_z12 + 8.*m_z2k;

  Box2Abb[1329]=4. + Box2Abb[1328]*m_z12;

  Box2Abb[1330]=Box2Abb[1329]*m_x_2 + Box2Abb[1324]*m_x_3*m_z12 + Box2Abb[1326]*m_x*m_z12*m_z2k - Box2Abb[1296]*m_z12_2*m_z2k_2;

  Box2Abb[1331]=2. + Box2Abb[62]*m_z12;

  Box2Abb[1332]=-2. + m_z12 + m_z12_3;

  Box2Abb[1333]=Box2Abb[1331]*pow(Box2Abb[4],2.) + Box2Abb[1332]*m_z2k + Box2Abb[70]*m_z12*m_z2k_2;

  Box2Abb[1334]=18. - 7.*m_z12 - 5.*m_z2k;

  Box2Abb[1335]=-17. + Box2Abb[1334]*m_z12 + 2.*m_z2k;

  Box2Abb[1336]=6. + Box2Abb[1335]*m_z12;

  Box2Abb[1337]=Box2Abb[1333]*m_x + Box2Abb[1336]*m_x_2 + Box2Abb[155]*m_x_3*m_z12 + pow(Box2Abb[5],2.)*Box2Abb[72]*m_z12*m_z2k;

  Box2Abb[1338]=-1. - 2.*Box2Abb[538]*m_x;

  Box2Abb[1339]=m_x + m_z12 - 2.*Box2Abb[210]*m_x*m_z12;

  Box2Abb[1340]=Box2Abb[1]*Box2Abb[212]*m_x_2 + Box2Abb[1338]*m_x*m_z2k + Box2Abb[1339]*m_z2k_2 - 2.*m_z12*m_z2k_3 + m_z12*m_z2k_4;

  Box2Abb[1341]=m_x - 2.*Box2Abb[12]*m_x_2 + m_x_3 + Box2Abb[137]*m_x*m_z2k + Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1342]=1. - Box2Abb[234]*m_x + 2.*m_x_2 - 3.*m_z2k + 2.*Box2Abb[137]*m_z2k;

  Box2Abb[1343]=-1. + m_z12 + 8.*m_z2k;

  Box2Abb[1344]=-3. + m_z12 + 2.*m_z2k;

  Box2Abb[1345]=1. + Box2Abb[1344]*m_z2k;

  Box2Abb[1346]=-7. + 12.*m_z2k;

  Box2Abb[1347]=-4. + 3.*Box2Abb[12]*m_z12 + Box2Abb[1346]*m_z2k;

  Box2Abb[1348]=11. - 8.*m_z2k;

  Box2Abb[1349]=-8. + Box2Abb[1348]*m_z2k;

  Box2Abb[1350]=-2. - 3.*Box2Abb[179]*m_z2k;

  Box2Abb[1351]=3. + Box2Abb[1350]*m_z12 + Box2Abb[1349]*m_z2k;

  Box2Abb[1352]=Box2Abb[1351]*m_x + Box2Abb[1347]*m_x_2 - Box2Abb[1343]*m_x_3 + 2.*m_x_4 + Box2Abb[1345]*Box2Abb[61]*m_z2k;

  Box2Abb[1353]=m_x - Box2Abb[391]*m_x_2 + m_x_3*m_z12 + Box2Abb[11]*Box2Abb[4]*m_x*m_z2k - pow(Box2Abb[61],2.)*m_z12*m_z2k + 3.*m_x*m_z12*m_z2k_2;

  Box2Abb[1354]=-3. + m_z12 - 5.*m_z12*m_z2k;

  Box2Abb[1355]=-1. + 3.*Box2Abb[72]*m_z12;

  Box2Abb[1356]=-1. + m_z12 + Box2Abb[1036]*m_z2k + Box2Abb[1355]*m_z2k_2 + 5.*m_z12*m_z2k_3;

  Box2Abb[1357]=-2. + m_z12 + Box2Abb[72]*m_z2k + 10.*m_z2k_2;

  Box2Abb[1358]=2. + Box2Abb[1357]*m_z12 + 5.*m_z2k;

  Box2Abb[1359]=-4. + m_z2k;

  Box2Abb[1360]=11. + 2.*Box2Abb[1069]*m_z2k;

  Box2Abb[1361]=-2. + Box2Abb[1360]*m_z2k;

  Box2Abb[1362]=1. + Box2Abb[1361]*m_z12 + Box2Abb[470]*Box2Abb[61]*m_z12_2 + Box2Abb[1359]*m_z2k;

  Box2Abb[1363]=-Box2Abb[1362]*m_x_2 + Box2Abb[1358]*m_x_3 + Box2Abb[1354]*m_x_4 + m_x_5*m_z12 + Box2Abb[1356]*m_x*m_z2k - Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[1364]=-1. + m_z12 + 2.*m_z12*m_z2k;

  Box2Abb[1365]=-Box2Abb[1364]*m_x + m_x_2*m_z12 + Box2Abb[5]*m_z12*m_z2k;

  Box2Abb[1366]=3. - Box2Abb[411]*m_z12 + m_z12_2 + m_z2k + m_z2k_2;

  Box2Abb[1367]=-5. + 2.*m_z12 + 3.*m_z2k;

  Box2Abb[1368]=3. + Box2Abb[1367]*m_z12;

  Box2Abb[1369]=-3. + m_z12 + Box2Abb[201]*m_z2k + 3.*m_z2k_2;

  Box2Abb[1370]=2. + Box2Abb[1369]*m_z12 + 3.*m_z2k;

  Box2Abb[1371]=Box2Abb[1370]*m_x_2 - Box2Abb[1368]*m_x_3 + m_x_4*m_z12 - Box2Abb[1366]*m_x*m_z12*m_z2k - Box2Abb[4]*m_z12_2*m_z2k_2;

  Box2Abb[1372]=-24. + 15.*m_z12 + 20.*m_z2k;

  Box2Abb[1373]=12. + Box2Abb[1372]*m_z12;

  Box2Abb[1374]=-2. + 11.*m_z2k;

  Box2Abb[1375]=3. + 14.*Box2Abb[61]*m_z2k;

  Box2Abb[1376]=-2. + Box2Abb[1375]*m_z12 + Box2Abb[1374]*m_z12_2 + m_z12_3 + 6.*m_z2k + 4.*Box2Abb[179]*m_z2k_2;

  Box2Abb[1377]=-5. + 13.*m_z2k;

  Box2Abb[1378]=3. + 4.*Box2Abb[514]*m_z2k;

  Box2Abb[1379]=-8. + 17.*m_z2k;

  Box2Abb[1380]=3. + 2.*Box2Abb[1379]*m_z2k;

  Box2Abb[1381]=2. + Box2Abb[1380]*m_z12 + Box2Abb[1377]*m_z12_2 + 2.*m_z12_3 + 2.*Box2Abb[1378]*m_z2k;

  Box2Abb[1382]=-2. + Box2Abb[1381]*m_z12;

  Box2Abb[1383]=-7. + 6.*m_z2k;

  Box2Abb[1384]=-34. + 50.*m_z2k;

  Box2Abb[1385]=24. + Box2Abb[1384]*m_z12 + 15.*m_z12_2 + 8.*Box2Abb[1383]*m_z2k;

  Box2Abb[1386]=Box2Abb[1385]*m_z12 + 24.*m_z2k;

  Box2Abb[1387]=-4. + Box2Abb[1386]*m_z12;

  Box2Abb[1388]=Box2Abb[1387]*m_x_3 - Box2Abb[1382]*m_x_2*m_z12 - 2.*Box2Abb[1373]*m_x_4*m_z12 + 12.*m_x_5*m_z12_2 + Box2Abb[1376]*m_x*m_z12_2*m_z2k - Box2Abb[5]*m_z12_4*m_z2k_2;

  Box2Abb[1389]=-2. + m_z12 + 4.*m_z12_2;

  Box2Abb[1390]=3. + 2.*m_z12;

  Box2Abb[1391]=-3. + Box2Abb[1390]*m_z12;

  Box2Abb[1392]=-2.*pow(Box2Abb[4],3.)*m_z12 - Box2Abb[1389]*Box2Abb[4]*m_z2k - 2.*Box2Abb[1391]*m_z2k_2 + 2.*Box2Abb[201]*m_z2k_3 + 4.*m_z2k_4;

  Box2Abb[1393]=-8. + 5.*m_z12 + 2.*m_z2k;

  Box2Abb[1394]=4. + Box2Abb[1393]*m_z12;

  Box2Abb[1395]=6. + Box2Abb[485]*m_z12;

  Box2Abb[1396]=-12. + 5.*Box2Abb[1395]*m_z12;

  Box2Abb[1397]=2.*pow(Box2Abb[4],3.) + Box2Abb[1396]*m_z2k + 6.*Box2Abb[53]*m_z12*m_z2k_2 + 24.*m_z2k_3 - 12.*m_z2k_4;

  Box2Abb[1398]=Box2Abb[1397]*m_z12 + 2.*m_z2k;

  Box2Abb[1399]=36. + 23.*Box2Abb[72]*m_z12 - 88.*m_z2k + 64.*m_z12*m_z2k + 8.*m_z2k_2;

  Box2Abb[1400]=-16. + Box2Abb[1399]*m_z12 + 48.*m_z2k;

  Box2Abb[1401]=4. + Box2Abb[1400]*m_z12;

  Box2Abb[1402]=1. - 44.*m_z2k;

  Box2Abb[1403]=2. - 3.*m_z2k;

  Box2Abb[1404]=15. + 70.*m_z2k - 36.*m_z2k_2;

  Box2Abb[1405]=9. + 2.*m_z2k;

  Box2Abb[1406]=3. + Box2Abb[1405]*m_z2k;

  Box2Abb[1407]=6. + 2.*Box2Abb[1406]*Box2Abb[460]*m_z12 + Box2Abb[1404]*m_z12_2 + Box2Abb[1402]*m_z12_3 - 4.*m_z12_4 + 8.*Box2Abb[1403]*m_z2k;

  Box2Abb[1408]=Box2Abb[1407]*m_z12 - 4.*m_z2k;

  Box2Abb[1409]=Box2Abb[1408]*m_x_3 + Box2Abb[1401]*m_x_4 + Box2Abb[1398]*m_x_2*m_z12 - 6.*Box2Abb[1394]*m_x_5*m_z12 + 4.*m_x_6*m_z12_2 + Box2Abb[1392]*m_x*m_z12_2*m_z2k + Box2Abb[5]*m_z12_4*m_z2k_3;

  Box2Abb[1410]=pow(Box2Abb[4],3.) + 3.*pow(Box2Abb[4],2.)*m_z2k + 3.*m_z12*m_z2k_2;

  Box2Abb[1411]=1. + Box2Abb[389]*m_z12;

  Box2Abb[1412]=Box2Abb[1410]*m_x - 3.*Box2Abb[1411]*m_x_2 + m_x_3*m_z12 - m_z12*m_z2k_3;

  Box2Abb[1413]=-1. + m_z12 + 2.*m_z2k;

  Box2Abb[1414]=m_x + Box2Abb[1413]*m_x*m_z12 - 3.*m_x_2*m_z12 + Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1415]=m_s - m_s*m_z12;

  Box2Abb[1416]=2. + Box2Abb[841]*m_z12 + m_z12_2 + Box2Abb[230]*m_z2k;

  Box2Abb[1417]=-5. + 3.*m_z12 + 6.*m_z2k;

  Box2Abb[1418]=2. + Box2Abb[1417]*m_z12;

  Box2Abb[1419]=-1. + 8.*m_z2k;

  Box2Abb[1420]=-2. + Box2Abb[1419]*m_z12 + m_z12_2 + 12.*Box2Abb[61]*m_z2k;

  Box2Abb[1421]=2. + Box2Abb[1420]*m_z12 + 4.*m_z2k;

  Box2Abb[1422]=Box2Abb[1421]*m_x_2 - 2.*Box2Abb[1418]*m_x_3 + 4.*m_x_4*m_z12 - 2.*Box2Abb[1416]*m_x*m_z12*m_z2k - Box2Abb[4]*m_z12_2*m_z2k_2;

  Box2Abb[1423]=6. + Box2Abb[538]*m_z12;

  Box2Abb[1424]=1. - 2.*m_z12;

  Box2Abb[1425]=4. - Box2Abb[1423]*m_z12 + 3.*m_z2k + 3.*Box2Abb[1424]*m_z12*m_z2k - 2.*Box2Abb[538]*m_z2k_2;

  Box2Abb[1426]=-3. + 3.*Box2Abb[12]*m_z12 - Box2Abb[222]*m_z2k;

  Box2Abb[1427]=3. - 2.*m_z12 - 14.*m_z2k;

  Box2Abb[1428]=-3. + Box2Abb[1427]*m_z12 + 18.*m_z2k;

  Box2Abb[1429]=2. + Box2Abb[1428]*m_z12;

  Box2Abb[1430]=4. + Box2Abb[538]*m_z12 + Box2Abb[220]*m_z12*m_z2k + 6.*Box2Abb[62]*m_z2k_2;

  Box2Abb[1431]=-2.*Box2Abb[12] + Box2Abb[1430]*m_z12;

  Box2Abb[1432]=Box2Abb[1431]*m_x_3 + Box2Abb[1429]*m_x_4 + Box2Abb[710]*m_x_5*m_z12 + Box2Abb[1425]*m_x_2*m_z12*m_z2k + Box2Abb[1426]*m_x*m_z12_2*m_z2k_2 - Box2Abb[4]*m_z12_3*m_z2k_3;

  Box2Abb[1433]=-3. + Box2Abb[166]*m_x + Box2Abb[1046]*m_z12;

  Box2Abb[1434]=9. + 26.*m_x;

  Box2Abb[1435]=-3. + Box2Abb[1434]*m_x;

  Box2Abb[1436]=4. + 3.*m_x;

  Box2Abb[1437]=7. + 22.*m_x;

  Box2Abb[1438]=2. + 4.*m_x + 2.*Box2Abb[1435]*m_z12 - 12.*Box2Abb[1436]*m_x*m_z12_2 + Box2Abb[1437]*m_z12_3 - 3.*m_z12_4;

  Box2Abb[1439]=1. + 5.*m_x;

  Box2Abb[1440]=1. + 6.*Box2Abb[1439]*m_x;

  Box2Abb[1441]=9. + 14.*m_x;

  Box2Abb[1442]=m_z12 + 2.*m_x*m_z12;

  Box2Abb[1443]=9.*pow(Box2Abb[1442],2.) - 4.*m_x - 2.*Box2Abb[1440]*m_z12 - Box2Abb[1441]*m_z12_3 + 2.*m_z12_4;

  Box2Abb[1444]=7. - 3.*m_z12;

  Box2Abb[1445]=3. + Box2Abb[520]*m_z12;

  Box2Abb[1446]=2.*Box2Abb[1445]*m_x + 4.*Box2Abb[1444]*m_x_2 + Box2Abb[4]*m_z12_2;

  Box2Abb[1447]=-4.*m_x + m_z12_2;

  Box2Abb[1448]=Box2Abb[1433]*Box2Abb[51]*m_x_3*m_z12 + Box2Abb[1438]*m_x_2*m_z2k + Box2Abb[1443]*m_x*m_z2k_2 + Box2Abb[1446]*m_z12*m_z2k_3 + Box2Abb[1447]*m_z12*m_z2k_4;

  Box2Abb[1449]=-9. + 4.*m_z12;

  Box2Abb[1450]=3. + Box2Abb[1449]*m_z12;

  Box2Abb[1451]=3. + Box2Abb[72]*m_z12;

  Box2Abb[1452]=2.*pow(Box2Abb[4],3.)*m_z12 + 2.*m_z2k + Box2Abb[1450]*m_z12*m_z2k - 2.*Box2Abb[1451]*m_z2k_2 - 4.*Box2Abb[4]*m_z2k_3;

  Box2Abb[1453]=-12. + 11.*m_z12;

  Box2Abb[1454]=-2. + Box2Abb[1453]*m_z12_2;

  Box2Abb[1455]=2. + Box2Abb[1087]*m_z12_2;

  Box2Abb[1456]=-2.*pow(Box2Abb[4],3.)*m_z12 - Box2Abb[1454]*Box2Abb[4]*m_z2k + 2.*Box2Abb[1455]*m_z2k_2 + 4.*Box2Abb[4]*m_z12*m_z2k_3;

  Box2Abb[1457]=-24. + 32.*m_z12 - 15.*m_z12_2 - 20.*Box2Abb[4]*m_z2k;

  Box2Abb[1458]=8. + Box2Abb[1457]*m_z12;

  Box2Abb[1459]=5. - 2.*m_z2k;

  Box2Abb[1460]=-5. + 26.*m_z2k;

  Box2Abb[1461]=-1. + 4.*Box2Abb[1359]*m_z2k;

  Box2Abb[1462]=4. + 3.*Box2Abb[1461]*m_z12 + Box2Abb[1460]*m_z12_2 + 4.*m_z12_3 + 6.*Box2Abb[1459]*m_z2k;

  Box2Abb[1463]=Box2Abb[1462]*m_z12 - 12.*m_z2k;

  Box2Abb[1464]=Box2Abb[1456]*m_x_2 + Box2Abb[1463]*m_x_3 + Box2Abb[1458]*m_x_4 + 8.*Box2Abb[4]*m_x_5*m_z12 + Box2Abb[1452]*m_x*m_z12*m_z2k + Box2Abb[5]*m_z12_3*m_z2k_3;

  Box2Abb[1465]=pow(Box2Abb[4],3.) + 3.*pow(Box2Abb[4],2.)*m_z2k + 3.*Box2Abb[72]*m_z2k_2;

  Box2Abb[1466]=12. - 5.*m_z12 - 9.*m_z2k;

  Box2Abb[1467]=-9. + Box2Abb[1466]*m_z12 + 12.*m_z2k;

  Box2Abb[1468]=2. + Box2Abb[1467]*m_z12;

  Box2Abb[1469]=Box2Abb[1468]*m_x_2 + Box2Abb[1465]*m_x*m_z12 + Box2Abb[710]*m_x_3*m_z12 + m_z12_2*m_z2k_3;

  Box2Abb[1470]=-9. + 3.*m_z12 - 5.*m_z2k;

  Box2Abb[1471]=6. + Box2Abb[1470]*m_z12;

  Box2Abb[1472]=-3.*pow(Box2Abb[4],4.) - 6.*Box2Abb[4]*Box2Abb[538]*m_z12*m_z2k + 10.*m_z12_2*m_z2k_2;

  Box2Abb[1473]=Box2Abb[201]*pow(Box2Abb[4],2.) + 6.*Box2Abb[4]*m_z2k_2 + 2.*m_z2k_3;

  Box2Abb[1474]=6. + Box2Abb[201]*m_z12;

  Box2Abb[1475]=Box2Abb[1474]*pow(Box2Abb[4],2.) - 12.*Box2Abb[384]*m_z2k_2 - 10.*m_z12*m_z2k_3;

  Box2Abb[1476]=2. + Box2Abb[72]*m_z12_2;

  Box2Abb[1477]=Box2Abb[1476]*pow(Box2Abb[4],2.) + 6.*pow(Box2Abb[4],4.)*m_z2k - 36.*Box2Abb[4]*m_z12*m_z2k_2 - 20.*m_z12_2*m_z2k_3;

  Box2Abb[1478]=Box2Abb[1477]*m_x_3 + 2.*Box2Abb[1472]*m_x_4 + 2.*Box2Abb[1471]*m_x_5*m_z12 + 2.*m_x_6*m_z12_2 - Box2Abb[1475]*m_x_2*m_z12*m_z2k - Box2Abb[1473]*m_x*m_z12_2*m_z2k_2 - pow(Box2Abb[4],2.)*m_z12_3*m_z2k_3;

  Box2Abb[1479]=-3. + 7.*m_z12;

  Box2Abb[1480]=1. + Box2Abb[4]*m_z12;

  Box2Abb[1481]=-9. + 7.*m_z12;

  Box2Abb[1482]=6. + Box2Abb[1481]*m_z12;

  Box2Abb[1483]=-2. + Box2Abb[1482]*m_z12;

  Box2Abb[1484]=Box2Abb[1483]*pow(Box2Abb[4],3.)*m_x - Box2Abb[1480]*pow(Box2Abb[4],4.)*m_z12 - 2.*Box2Abb[1479]*pow(Box2Abb[4],3.)*m_x_2*m_z12 + 2.*Box2Abb[192]*Box2Abb[4]*m_x_3*m_z12_2 + 2.*m_x_4*m_z12_3;

  Box2Abb[1485]=-1. + 3.*m_z12;

  Box2Abb[1486]=6.*Box2Abb[1485]*pow(Box2Abb[4],3.)*m_x_2 + pow(Box2Abb[4],4.)*m_z12 - 4.*Box2Abb[4]*Box2Abb[485]*m_x_3*m_z12 - 4.*pow(Box2Abb[4],3.)*m_x*m_z12_2 - 10.*m_x_4*m_z12_2;

  Box2Abb[1487]=pow(Box2Abb[4],3.)*Box2Abb[538] + 24.*Box2Abb[4]*m_x_2 - 20.*m_x_3*m_z12;

  Box2Abb[1488]=pow(Box2Abb[4],3.) + m_x + Box2Abb[699]*m_x*m_z12 + 5.*m_x_2*m_z12;

  Box2Abb[1489]=4. + 5.*m_x - 3.*m_z12;

  Box2Abb[1490]=-1. + Box2Abb[1489]*m_z12;

  Box2Abb[1491]=Box2Abb[1484]*m_x + Box2Abb[1486]*m_z12*m_z2k - Box2Abb[1487]*m_z12_2*m_z2k_2 - 4.*Box2Abb[1488]*m_z12_2*m_z2k_3 + 2.*Box2Abb[1490]*m_z12_2*m_z2k_4 - 2.*m_z12_3*m_z2k_5;

  Box2Abb[1492]=-2. + 6.*m_x - m_z12;

  Box2Abb[1493]=2. + Box2Abb[1492]*m_z12;

  Box2Abb[1494]=1. + 6.*m_x - 2.*m_z12;

  Box2Abb[1495]=6.*m_x - m_z12;

  Box2Abb[1496]=Box2Abb[1493]*m_x_2 - 2.*Box2Abb[1494]*m_x*m_z12*m_z2k + Box2Abb[1495]*m_z12*m_z2k_2;

  Box2Abb[1497]=m_z12 + 2.*Box2Abb[9]*m_z2k + m_z2k_2;

  Box2Abb[1498]=6. - 3.*m_z12 + m_z12_2 + 6.*Box2Abb[9]*m_z2k + 12.*m_z2k_2;

  Box2Abb[1499]=-2. + Box2Abb[1498]*m_z12 + 6.*m_z2k;

  Box2Abb[1500]=m_z12 + 6.*m_z2k;

  Box2Abb[1501]=4. + Box2Abb[1500]*Box2Abb[194]*m_z12 + 24.*m_z2k;

  Box2Abb[1502]=-2. + Box2Abb[1501]*m_z12;

  Box2Abb[1503]=Box2Abb[1502]*m_x_4 - 6.*Box2Abb[379]*m_x_5*m_z12 + 6.*m_x_6*m_z12_2 - 2.*Box2Abb[1499]*m_x_3*m_z12*m_z2k + 6.*Box2Abb[1497]*m_x_2*m_z12_2*m_z2k_2 - 2.*Box2Abb[399]*m_x*m_z12_3*m_z2k_3 + m_z12_4*m_z2k_4;

  Box2Abb[1504]=-9. + 10.*m_z12;

  Box2Abb[1505]=Box2Abb[1481]*Box2Abb[4]*m_z12 + 4.*Box2Abb[1504]*Box2Abb[4]*m_z2k + 30.*m_z12*m_z2k_2;

  Box2Abb[1506]=-13. + 7.*m_z12 + 9.*m_z2k;

  Box2Abb[1507]=6. + Box2Abb[1506]*m_z12;

  Box2Abb[1508]=Box2Abb[512]*Box2Abb[61] + 4.*pow(Box2Abb[61],2.)*m_z12 + 5.*Box2Abb[61]*m_z12_2 + 2.*m_z12_3;

  Box2Abb[1509]=-3. + Box2Abb[62]*m_z12;

  Box2Abb[1510]=4. + Box2Abb[1509]*m_z12 + 6.*m_z2k + 8.*Box2Abb[538]*m_z12*m_z2k - 4.*Box2Abb[242]*m_z2k_2;

  Box2Abb[1511]=Box2Abb[1510]*m_z12 - 2.*m_z2k;

  Box2Abb[1512]=7. - Box2Abb[201]*Box2Abb[538]*m_z12 + 8.*Box2Abb[4]*m_z12*m_z2k + 14.*Box2Abb[4]*m_z2k_2 + 6.*m_z2k_3;

  Box2Abb[1513]=-1. + Box2Abb[1512]*m_z12;

  Box2Abb[1514]=4. + Box2Abb[538]*m_z12;

  Box2Abb[1515]=-5. + Box2Abb[1514]*m_z12;

  Box2Abb[1516]=-8. + 3.*m_z12;

  Box2Abb[1517]=6. + Box2Abb[1516]*m_z12;

  Box2Abb[1518]=-9. + 8.*m_z12;

  Box2Abb[1519]=5. + Box2Abb[1515]*m_z12 - 8.*m_z2k + 6.*Box2Abb[1517]*m_z12*m_z2k + 4.*Box2Abb[1518]*Box2Abb[4]*m_z2k_2 + 20.*m_z12*m_z2k_3;

  Box2Abb[1520]=2. - Box2Abb[1519]*m_z12 - 2.*m_z2k;

  Box2Abb[1521]=Box2Abb[1520]*m_x_3 + Box2Abb[1505]*m_x_4*m_z12 - 2.*Box2Abb[1507]*m_x_5*m_z12 + 4.*m_x_6*m_z12_2 + Box2Abb[1511]*Box2Abb[4]*m_x_2*m_z2k + Box2Abb[1513]*m_x*m_z12*m_z2k_2 - Box2Abb[1508]*m_z12_2*m_z2k_3;

  Box2Abb[1522]=1. - 3.*m_x + m_z12 + m_x_2*m_z12 + Box2Abb[1]*m_z12_2;

  Box2Abb[1523]=2. + 3.*m_x - 3.*m_z12;

  Box2Abb[1524]=m_x + m_z12 + Box2Abb[1523]*m_x*m_z12;

  Box2Abb[1525]=2. + 3.*m_x;

  Box2Abb[1526]=Box2Abb[1522]*m_x - Box2Abb[1524]*m_z2k + Box2Abb[1525]*m_z12*m_z2k_2 - m_z12*m_z2k_3;

  Box2Abb[1527]=15. + m_z12;

  Box2Abb[1528]=28. + 5.*m_z12;

  Box2Abb[1529]=-16. + Box2Abb[1528]*m_z12;

  Box2Abb[1530]=-10. + Box2Abb[1529]*m_z12;

  Box2Abb[1531]=5. + Box2Abb[1530]*m_z12;

  Box2Abb[1532]=4. + 17.*m_z12_2;

  Box2Abb[1533]=-12. + 19.*m_z12;

  Box2Abb[1534]=-15. + Box2Abb[1527]*m_z12 - 8.*m_z2k + 16.*Box2Abb[693]*m_z12*m_z2k + Box2Abb[1531]*m_z2k_2 + 4.*Box2Abb[1532]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[1533]*m_z12*m_z2k_4;

  Box2Abb[1535]=-3. + 8.*Box2Abb[201]*Box2Abb[4]*m_z12;

  Box2Abb[1536]=11. + 2.*m_z12;

  Box2Abb[1537]=-26. + Box2Abb[1536]*m_z12;

  Box2Abb[1538]=7. + Box2Abb[1537]*m_z12;

  Box2Abb[1539]=4. - 11.*m_z12;

  Box2Abb[1540]=-6. + Box2Abb[1539]*m_z12;

  Box2Abb[1541]=12. + Box2Abb[1540]*m_z12;

  Box2Abb[1542]=-11. + Box2Abb[1541]*m_z12;

  Box2Abb[1543]=23. - 13.*m_z12;

  Box2Abb[1544]=-6. + Box2Abb[1543]*m_z12;

  Box2Abb[1545]=1. + Box2Abb[1544]*m_z12;

  Box2Abb[1546]=8. - 15.*m_z12;

  Box2Abb[1547]=3. - 3.*m_z12 + Box2Abb[1535]*m_z2k - Box2Abb[1538]*Box2Abb[4]*m_z2k_2 + Box2Abb[1542]*m_z2k_3 + 4.*Box2Abb[1545]*m_z2k_4 + 3.*Box2Abb[1546]*m_z12*m_z2k_5;

  Box2Abb[1548]=7. + 2.*Box2Abb[1359]*m_z2k;

  Box2Abb[1549]=5. + Box2Abb[841]*m_z2k;

  Box2Abb[1550]=Box2Abb[1549]*Box2Abb[61] + Box2Abb[1548]*m_z12 + Box2Abb[179]*m_z12_2;

  Box2Abb[1551]=16. + 31.*m_z2k;

  Box2Abb[1552]=8. - Box2Abb[1551]*m_z12 + 24.*m_z2k;

  Box2Abb[1553]=4. + Box2Abb[1552]*m_z12;

  Box2Abb[1554]=46. + 81.*m_z2k;

  Box2Abb[1555]=18. + Box2Abb[1554]*m_z2k;

  Box2Abb[1556]=5. + Box2Abb[1555]*m_z12 - 12.*Box2Abb[450]*m_z2k + 10.*m_z12_2*m_z2k;

  Box2Abb[1557]=-17. + Box2Abb[1556]*m_z12 - 16.*m_z2k;

  Box2Abb[1558]=8. + 21.*m_z2k;

  Box2Abb[1559]=23. + 24.*m_z2k;

  Box2Abb[1560]=-1. + 4.*Box2Abb[510]*m_z2k;

  Box2Abb[1561]=-21. + 4.*Box2Abb[1560]*m_z2k;

  Box2Abb[1562]=18. + 115.*m_z2k;

  Box2Abb[1563]=17. + Box2Abb[1562]*m_z2k;

  Box2Abb[1564]=8. + Box2Abb[1563]*m_z2k;

  Box2Abb[1565]=25. + Box2Abb[1561]*m_z12 - Box2Abb[1564]*m_z12_2 + Box2Abb[1559]*m_z2k - 2.*Box2Abb[1558]*m_z12_3*m_z2k;

  Box2Abb[1566]=8. - m_z2k + 4.*m_z2k_3;

  Box2Abb[1567]=6. + 7.*m_z2k;

  Box2Abb[1568]=-5. + Box2Abb[1567]*m_z2k;

  Box2Abb[1569]=-14. + 9.*m_z2k;

  Box2Abb[1570]=-7. + Box2Abb[1569]*m_z2k;

  Box2Abb[1571]=9. + Box2Abb[1570]*m_z2k;

  Box2Abb[1572]=-1. + Box2Abb[1571]*m_z2k;

  Box2Abb[1573]=-31. + 11.*m_z2k;

  Box2Abb[1574]=7. + Box2Abb[1573]*m_z2k;

  Box2Abb[1575]=17. + Box2Abb[1574]*m_z2k;

  Box2Abb[1576]=-10. + Box2Abb[1575]*m_z2k;

  Box2Abb[1577]=-Box2Abb[1566]*pow(Box2Abb[61],2.) + Box2Abb[1576]*Box2Abb[61]*m_z12 + 2.*Box2Abb[1572]*m_z12_2 + Box2Abb[1568]*m_z12_3*m_z2k;

  Box2Abb[1578]=Box2Abb[1547]*m_x_2 + Box2Abb[1534]*m_x_3 + Box2Abb[1565]*m_x_4 + Box2Abb[1557]*m_x_5 + Box2Abb[1553]*m_x_6 + Box2Abb[1324]*m_x_7*m_z12 + Box2Abb[1577]*m_x*m_z12*m_z2k - Box2Abb[1550]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k_2;

  Box2Abb[1579]=2. + Box2Abb[389]*m_z12 - 2.*m_z2k;

  Box2Abb[1580]=3. + 7.*m_z2k;

  Box2Abb[1581]=-1. + 9.*m_z2k_2;

  Box2Abb[1582]=-11. + 9.*m_z2k;

  Box2Abb[1583]=-3. + Box2Abb[1582]*m_z2k;

  Box2Abb[1584]=-36. + 11.*m_z2k;

  Box2Abb[1585]=-4. + Box2Abb[1584]*m_z2k;

  Box2Abb[1586]=2.*pow(Box2Abb[61],3.) - 2.*Box2Abb[1581]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[1585]*Box2Abb[61]*m_z12_2*m_z2k + 2.*Box2Abb[1583]*m_z12_3*m_z2k + Box2Abb[1580]*m_z12_4*m_z2k;

  Box2Abb[1587]=10. + 31.*m_z2k;

  Box2Abb[1588]=12. - Box2Abb[1587]*m_z12 + 38.*m_z2k;

  Box2Abb[1589]=2. + Box2Abb[1588]*m_z12;

  Box2Abb[1590]=6. + m_z2k + 42.*m_z2k_2;

  Box2Abb[1591]=52. - 115.*m_z2k;

  Box2Abb[1592]=-4. + Box2Abb[1591]*m_z2k;

  Box2Abb[1593]=5. + Box2Abb[1592]*m_z2k;

  Box2Abb[1594]=3. + m_z2k + 25.*m_z2k_2;

  Box2Abb[1595]=1. + 2.*Box2Abb[1594]*m_z2k;

  Box2Abb[1596]=2. + 3.*Box2Abb[1595]*m_z12 + Box2Abb[1593]*m_z12_2 - Box2Abb[1590]*m_z12_3 + m_z2k + 4.*m_z2k_2;

  Box2Abb[1597]=14. + 81.*m_z2k;

  Box2Abb[1598]=5. + Box2Abb[1597]*m_z2k;

  Box2Abb[1599]=-6. + Box2Abb[1598]*m_z12 + 2.*Box2Abb[510]*m_z12_2 - 34.*Box2Abb[256]*m_z2k;

  Box2Abb[1600]=-7. + Box2Abb[1599]*m_z12 - 6.*m_z2k;

  Box2Abb[1601]=1. + Box2Abb[831]*m_z2k;

  Box2Abb[1602]=-15. + 17.*m_z2k;

  Box2Abb[1603]=3. + Box2Abb[1602]*m_z2k;

  Box2Abb[1604]=-38. + 65.*m_z2k;

  Box2Abb[1605]=15. + Box2Abb[1604]*m_z2k;

  Box2Abb[1606]=17. + Box2Abb[1605]*m_z2k;

  Box2Abb[1607]=-148. + 95.*m_z2k;

  Box2Abb[1608]=28. + Box2Abb[1607]*m_z2k;

  Box2Abb[1609]=46. + Box2Abb[1608]*m_z2k;

  Box2Abb[1610]=-10. + Box2Abb[1609]*m_z2k;

  Box2Abb[1611]=Box2Abb[1601]*Box2Abb[222] + Box2Abb[1610]*m_z12_2 + 2.*Box2Abb[1603]*Box2Abb[444]*m_z12_3 - 2.*Box2Abb[1606]*m_z12*m_z2k + 5.*m_z12_4*m_z2k_2;

  Box2Abb[1612]=9. - 11.*m_z2k;

  Box2Abb[1613]=1. - 2.*m_z2k;

  Box2Abb[1614]=-1. + 13.*m_z2k;

  Box2Abb[1615]=2. + Box2Abb[1614]*m_z2k;

  Box2Abb[1616]=-19. + 33.*m_z2k;

  Box2Abb[1617]=6. + Box2Abb[1616]*m_z2k;

  Box2Abb[1618]=-29. + 2.*Box2Abb[1617]*m_z2k;

  Box2Abb[1619]=3. + Box2Abb[1618]*m_z2k;

  Box2Abb[1620]=142. - 45.*m_z2k;

  Box2Abb[1621]=-86. + Box2Abb[1620]*m_z2k;

  Box2Abb[1622]=45. + Box2Abb[1621]*m_z2k;

  Box2Abb[1623]=-29. + Box2Abb[1622]*m_z2k;

  Box2Abb[1624]=5. + Box2Abb[1623]*m_z2k;

  Box2Abb[1625]=Box2Abb[1619]*Box2Abb[61]*m_z12 + Box2Abb[1624]*m_z12_2 - pow(Box2Abb[1613],2.)*Box2Abb[1615]*m_z12_3 - 3.*pow(Box2Abb[61],2.)*Box2Abb[833]*m_z2k + Box2Abb[1612]*m_z12_4*m_z2k_2;

  Box2Abb[1626]=Box2Abb[1625]*m_x_2 + Box2Abb[1611]*m_x_3 + Box2Abb[1596]*m_x_4 + Box2Abb[1600]*m_x_5 + Box2Abb[1589]*m_x_6 + Box2Abb[710]*m_x_7*m_z12 + Box2Abb[1586]*Box2Abb[61]*m_x*m_z2k - Box2Abb[1579]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[1627]=2. + Box2Abb[72]*m_z12;

  Box2Abb[1628]=2.*Box2Abb[4]*m_z12 + Box2Abb[1627]*m_z2k + Box2Abb[72]*m_z2k_2;

  Box2Abb[1629]=-7. + m_z12;

  Box2Abb[1630]=32. + Box2Abb[1629]*m_z12;

  Box2Abb[1631]=-23. + Box2Abb[1630]*m_z12;

  Box2Abb[1632]=-22. + 3.*m_z12;

  Box2Abb[1633]=-62. + Box2Abb[1632]*Box2Abb[201]*m_z12;

  Box2Abb[1634]=17. + Box2Abb[1633]*m_z12;

  Box2Abb[1635]=-10. + 7.*m_z12;

  Box2Abb[1636]=9. + Box2Abb[1635]*m_z12;

  Box2Abb[1637]=2. + 5.*m_z12;

  Box2Abb[1638]=1. + 18.*Box2Abb[4]*m_z12_2 + 2.*m_z2k + 2.*Box2Abb[1631]*m_z12*m_z2k + Box2Abb[1634]*m_z2k_2 + 4.*Box2Abb[1636]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[1637]*m_z12*m_z2k_4;

  Box2Abb[1639]=-3. + 5.*Box2Abb[328]*m_z12;

  Box2Abb[1640]=-23. + 9.*m_z12 - 6.*m_z12_2;

  Box2Abb[1641]=16. + Box2Abb[1640]*m_z12;

  Box2Abb[1642]=-28. + 9.*m_z12;

  Box2Abb[1643]=43. + Box2Abb[1642]*m_z12;

  Box2Abb[1644]=-18. + Box2Abb[1643]*m_z12;

  Box2Abb[1645]=36. - 5.*m_z12;

  Box2Abb[1646]=-82. + Box2Abb[1645]*m_z12;

  Box2Abb[1647]=70. + Box2Abb[1646]*m_z12;

  Box2Abb[1648]=-31. + Box2Abb[1647]*m_z12;

  Box2Abb[1649]=11. - 6.*m_z12;

  Box2Abb[1650]=-2. + Box2Abb[1649]*m_z12;

  Box2Abb[1651]=7. + Box2Abb[1650]*m_z12;

  Box2Abb[1652]=-2. + Box2Abb[1639]*m_z12 + m_z2k + Box2Abb[1641]*m_z12*m_z2k + Box2Abb[1644]*Box2Abb[4]*m_z2k_2 + Box2Abb[1648]*m_z2k_3 + 2.*Box2Abb[1651]*m_z2k_4 - 3.*Box2Abb[669]*m_z12*m_z2k_5;

  Box2Abb[1653]=8. + Box2Abb[893]*m_z12 - 10.*m_z2k;

  Box2Abb[1654]=6. - Box2Abb[1653]*m_z12;

  Box2Abb[1655]=6. + m_z2k - m_z2k_2;

  Box2Abb[1656]=-1. + Box2Abb[854]*m_z2k;

  Box2Abb[1657]=1. + Box2Abb[222]*Box2Abb[680]*m_z2k;

  Box2Abb[1658]=4. + 3.*m_z2k;

  Box2Abb[1659]=6. + Box2Abb[1359]*Box2Abb[1658]*m_z2k;

  Box2Abb[1660]=4. + Box2Abb[1659]*m_z2k;

  Box2Abb[1661]=-2.*Box2Abb[12]*Box2Abb[1656]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[1660]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[1657]*m_z12_3 + 2.*pow(Box2Abb[61],3.)*m_z2k + Box2Abb[1655]*m_z12_4*m_z2k;

  Box2Abb[1662]=27. + 44.*m_z2k;

  Box2Abb[1663]=3. - 22.*m_z2k;

  Box2Abb[1664]=-14. + Box2Abb[1663]*m_z2k;

  Box2Abb[1665]=-7. + 45.*m_z2k;

  Box2Abb[1666]=11. + Box2Abb[1665]*m_z2k;

  Box2Abb[1667]=-47. + 5.*m_z2k;

  Box2Abb[1668]=-13. + Box2Abb[1667]*m_z2k;

  Box2Abb[1669]=-17. + 2.*Box2Abb[1668]*m_z2k;

  Box2Abb[1670]=16. + Box2Abb[1669]*m_z12 - Box2Abb[1666]*Box2Abb[61]*m_z12_2 + Box2Abb[1664]*m_z12_3 + Box2Abb[1662]*m_z2k;

  Box2Abb[1671]=23. - 9.*m_z2k;

  Box2Abb[1672]=4. + 6.*m_z2k;

  Box2Abb[1673]=-10. + 39.*m_z2k;

  Box2Abb[1674]=-5. + Box2Abb[1673]*m_z2k;

  Box2Abb[1675]=28. + Box2Abb[1674]*m_z12 + Box2Abb[1672]*m_z12_2 + 2.*Box2Abb[1671]*m_z2k;

  Box2Abb[1676]=-21. + Box2Abb[1675]*m_z12 - 26.*m_z2k;

  Box2Abb[1677]=-Box2Abb[1661]*Box2Abb[61]*m_x + Box2Abb[1652]*m_x_2 + Box2Abb[1638]*m_x_3 + Box2Abb[1670]*m_x_4 + Box2Abb[1676]*m_x_5 + Box2Abb[1654]*m_x_6 + Box2Abb[155]*m_x_7*m_z12 + Box2Abb[1628]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[1678]=21. + 2.*m_z12;

  Box2Abb[1679]=-45. + Box2Abb[1678]*m_z12;

  Box2Abb[1680]=15. + Box2Abb[1679]*m_z12;

  Box2Abb[1681]=12. + Box2Abb[1318]*m_z12;

  Box2Abb[1682]=-6. + 11.*m_z12;

  Box2Abb[1683]=16. - 9.*m_z12;

  Box2Abb[1684]=-3. - Box2Abb[1680]*m_z12 + 5.*m_z2k + Box2Abb[1681]*m_z12*m_z2k - 2.*Box2Abb[1682]*Box2Abb[816]*m_z2k_2 + 5.*Box2Abb[1683]*m_z12*m_z2k_3;

  Box2Abb[1685]=5. + 6.*m_z12;

  Box2Abb[1686]=-41. + Box2Abb[1685]*m_z12;

  Box2Abb[1687]=36. + Box2Abb[1686]*m_z12;

  Box2Abb[1688]=3. + Box2Abb[538]*m_z12;

  Box2Abb[1689]=-7. + 3.*Box2Abb[1688]*m_z12;

  Box2Abb[1690]=-21. + m_z12;

  Box2Abb[1691]=37. + Box2Abb[1690]*m_z12;

  Box2Abb[1692]=-23. + Box2Abb[1691]*m_z12;

  Box2Abb[1693]=-30. + 7.*m_z12;

  Box2Abb[1694]=4. + Box2Abb[1693]*m_z12;

  Box2Abb[1695]=-12. + 5.*m_z12;

  Box2Abb[1696]=-7. + Box2Abb[1687]*m_z12 + 2.*Box2Abb[1689]*m_z12*m_z2k + Box2Abb[1485]*Box2Abb[1692]*m_z2k_2 + 4.*Box2Abb[1694]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[1695]*m_z12*m_z2k_4;

  Box2Abb[1697]=4. - 17.*m_z2k;

  Box2Abb[1698]=-4. + Box2Abb[1697]*m_z12 + 24.*m_z2k;

  Box2Abb[1699]=4. + Box2Abb[1698]*m_z12;

  Box2Abb[1700]=-9. + 5.*m_z12;

  Box2Abb[1701]=-20. + 13.*m_z12;

  Box2Abb[1702]=13. + 2.*Box2Abb[1700]*m_z12 + 36.*m_z2k + 6.*Box2Abb[1629]*m_z12*m_z2k + 3.*Box2Abb[1701]*m_z2k_2;

  Box2Abb[1703]=-11. + Box2Abb[1702]*m_z12 - 16.*m_z2k;

  Box2Abb[1704]=-7. - 4.*Box2Abb[179]*m_z2k;

  Box2Abb[1705]=10. + 3.*m_z2k;

  Box2Abb[1706]=5. + Box2Abb[1705]*Box2Abb[61]*m_z2k;

  Box2Abb[1707]=3. - Box2Abb[1706]*m_z12 + Box2Abb[1189]*m_z12_2 + Box2Abb[1704]*m_z2k;

  Box2Abb[1708]=16. - 5.*m_z2k;

  Box2Abb[1709]=3. + Box2Abb[1708]*m_z2k;

  Box2Abb[1710]=-6. + Box2Abb[1709]*m_z2k;

  Box2Abb[1711]=7. - 9.*m_z2k + 6.*m_z2k_2;

  Box2Abb[1712]=-11. + 4.*Box2Abb[1711]*m_z2k;

  Box2Abb[1713]=-23. + m_z2k;

  Box2Abb[1714]=31. + Box2Abb[1713]*m_z2k;

  Box2Abb[1715]=-40. + 3.*Box2Abb[1714]*m_z2k;

  Box2Abb[1716]=-1. + Box2Abb[1715]*m_z2k;

  Box2Abb[1717]=-73. - 12.*Box2Abb[420]*m_z2k;

  Box2Abb[1718]=8. + Box2Abb[1717]*m_z2k;

  Box2Abb[1719]=13. + Box2Abb[1718]*m_z2k;

  Box2Abb[1720]=pow(Box2Abb[61],3.)*Box2Abb[725] + Box2Abb[1712]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[1716]*Box2Abb[61]*m_z12_2 + Box2Abb[1719]*m_z12_3 + Box2Abb[1710]*m_z12_4;

  Box2Abb[1721]=Box2Abb[1720]*m_x_2 + Box2Abb[1696]*m_x_3 + Box2Abb[1684]*m_x_4 + Box2Abb[1703]*m_x_5 + Box2Abb[1699]*m_x_6 + Box2Abb[1707]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12 + Box2Abb[166]*m_x_7*m_z12 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_2*m_z2k;

  Box2Abb[1722]=-1. + Box2Abb[166]*Box2Abb[538]*m_z12;

  Box2Abb[1723]=-26. + 7.*m_z12;

  Box2Abb[1724]=32. + Box2Abb[1723]*m_z12;

  Box2Abb[1725]=-10. + Box2Abb[1724]*m_z12;

  Box2Abb[1726]=2. + Box2Abb[1725]*m_z12;

  Box2Abb[1727]=-13. + 6.*m_z12;

  Box2Abb[1728]=4. + Box2Abb[1727]*m_z12;

  Box2Abb[1729]=2. + Box2Abb[1445]*m_z12 - 4.*m_z2k + Box2Abb[1722]*m_z12*m_z2k + Box2Abb[1726]*m_z2k_2 + 3.*Box2Abb[1728]*m_z12*m_z2k_3 + Box2Abb[1305]*m_z12*m_z2k_4;

  Box2Abb[1730]=2. + 2.*Box2Abb[179]*m_z12 + m_z12_2 + Box2Abb[179]*m_z2k;

  Box2Abb[1731]=1. + Box2Abb[1730]*m_z12 - m_z2k;

  Box2Abb[1732]=5. - 6.*m_z12;

  Box2Abb[1733]=-27. + Box2Abb[1732]*m_z12;

  Box2Abb[1734]=14. + 9.*m_z12;

  Box2Abb[1735]=59. - 42.*m_z12;

  Box2Abb[1736]=-109. + Box2Abb[1735]*m_z12;

  Box2Abb[1737]=72. + Box2Abb[1736]*m_z12;

  Box2Abb[1738]=16. - 23.*m_z12;

  Box2Abb[1739]=17. + Box2Abb[1733]*m_z12 + 55.*m_z2k - 3.*Box2Abb[1734]*m_z12*m_z2k + Box2Abb[1737]*m_z2k_2 + 5.*Box2Abb[1738]*m_z12*m_z2k_3;

  Box2Abb[1740]=7. + Box2Abb[1739]*m_z12;

  Box2Abb[1741]=8. + 3.*Box2Abb[72]*m_z12;

  Box2Abb[1742]=88. - 21.*m_z12;

  Box2Abb[1743]=-41. + Box2Abb[1742]*m_z12;

  Box2Abb[1744]=82. + 5.*Box2Abb[520]*m_z12;

  Box2Abb[1745]=-30. + Box2Abb[1744]*m_z12;

  Box2Abb[1746]=-15. + Box2Abb[1745]*m_z12;

  Box2Abb[1747]=-22. + 17.*m_z12;

  Box2Abb[1748]=12. + Box2Abb[1747]*m_z12;

  Box2Abb[1749]=-1. + 2.*Box2Abb[1741]*m_z12 - 34.*m_z2k + Box2Abb[1743]*m_z12*m_z2k + Box2Abb[1746]*m_z2k_2 + 4.*Box2Abb[1748]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[1533]*m_z12*m_z2k_4;

  Box2Abb[1750]=-8. + Box2Abb[1749]*m_z12 + 16.*m_z2k;

  Box2Abb[1751]=11. + 31.*m_z2k;

  Box2Abb[1752]=5. + Box2Abb[1751]*m_z12 - 24.*m_z2k;

  Box2Abb[1753]=12. - Box2Abb[1752]*m_z12;

  Box2Abb[1754]=7. + Box2Abb[1597]*m_z2k;

  Box2Abb[1755]=28. + Box2Abb[1754]*m_z12 + 2.*Box2Abb[510]*m_z12_2 + 20.*Box2Abb[1403]*m_z2k;

  Box2Abb[1756]=-31. + Box2Abb[1755]*m_z12 - 48.*m_z2k;

  Box2Abb[1757]=-3. + m_z2k_2 + 12.*m_z2k_3;

  Box2Abb[1758]=-31. + 26.*m_z2k;

  Box2Abb[1759]=20. + Box2Abb[1758]*m_z2k;

  Box2Abb[1760]=-6. + Box2Abb[1759]*m_z2k;

  Box2Abb[1761]=1. + Box2Abb[1760]*m_z2k;

  Box2Abb[1762]=139. - 45.*m_z2k;

  Box2Abb[1763]=-122. + Box2Abb[1762]*m_z2k;

  Box2Abb[1764]=58. + Box2Abb[1763]*m_z2k;

  Box2Abb[1765]=-54. + Box2Abb[1764]*m_z2k;

  Box2Abb[1766]=6. + Box2Abb[1765]*m_z2k;

  Box2Abb[1767]=-79. + 24.*m_z2k;

  Box2Abb[1768]=70. + Box2Abb[1767]*m_z2k;

  Box2Abb[1769]=-27. + Box2Abb[1768]*m_z2k;

  Box2Abb[1770]=44. + Box2Abb[1769]*m_z2k;

  Box2Abb[1771]=-8. + Box2Abb[1770]*m_z2k;

  Box2Abb[1772]=pow(Box2Abb[61],2.) + Box2Abb[1757]*Box2Abb[61]*m_z12 + Box2Abb[1771]*m_z12_2 + Box2Abb[1766]*m_z12_3 - 2.*Box2Abb[1761]*m_z12_4 + Box2Abb[1612]*m_z12_5*m_z2k_2;

  Box2Abb[1773]=Box2Abb[1772]*m_x_2 + Box2Abb[1750]*m_x_3 + Box2Abb[1740]*m_x_4 + Box2Abb[1756]*m_x_5*m_z12 + Box2Abb[1753]*m_x_6*m_z12 + Box2Abb[1324]*m_x_7*m_z12_2 + Box2Abb[1729]*Box2Abb[61]*m_x*m_z12*m_z2k - Box2Abb[1731]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_2;

  Box2Abb[1774]=10. - 17.*m_z12;

  Box2Abb[1775]=12. - 2.*Box2Abb[500]*m_z12 + Box2Abb[1774]*m_z12*m_z2k;

  Box2Abb[1776]=7. - 8.*m_z2k + 10.*m_z2k_3;

  Box2Abb[1777]=-23. + 2.*Box2Abb[411]*m_z2k;

  Box2Abb[1778]=7. + Box2Abb[1777]*m_z2k;

  Box2Abb[1779]=10. - 3.*m_z2k;

  Box2Abb[1780]=16. + Box2Abb[1779]*m_z2k;

  Box2Abb[1781]=-28. + Box2Abb[1780]*m_z2k;

  Box2Abb[1782]=10. + Box2Abb[1781]*m_z2k;

  Box2Abb[1783]=2.*pow(Box2Abb[61],2.) + Box2Abb[1776]*Box2Abb[61]*m_z12 + Box2Abb[1782]*m_z12_2 - Box2Abb[1778]*m_z12_3 + Box2Abb[1189]*m_z12_4;

  Box2Abb[1784]=-23. + 72.*m_z2k;

  Box2Abb[1785]=27. + 4.*m_z2k + 22.*m_z2k_2;

  Box2Abb[1786]=74. - 45.*m_z2k;

  Box2Abb[1787]=-45. + Box2Abb[1786]*m_z2k;

  Box2Abb[1788]=44. + Box2Abb[1787]*m_z2k;

  Box2Abb[1789]=-72. + 5.*m_z2k;

  Box2Abb[1790]=26. + Box2Abb[1789]*m_z2k;

  Box2Abb[1791]=-7. + 2.*Box2Abb[1790]*m_z2k;

  Box2Abb[1792]=-32. + Box2Abb[1791]*m_z12 + Box2Abb[1788]*m_z12_2 - Box2Abb[1785]*m_z12_3 - 2.*m_z12_4 + Box2Abb[1784]*m_z2k;

  Box2Abb[1793]=20. + Box2Abb[1792]*m_z12 + 6.*m_z2k;

  Box2Abb[1794]=2. + Box2Abb[187]*m_z2k;

  Box2Abb[1795]=-19. + 28.*Box2Abb[179]*m_z2k;

  Box2Abb[1796]=11. + Box2Abb[1795]*m_z2k;

  Box2Abb[1797]=-21. + 16.*m_z2k;

  Box2Abb[1798]=10. + Box2Abb[1797]*m_z2k;

  Box2Abb[1799]=3. + Box2Abb[1798]*m_z2k;

  Box2Abb[1800]=58. + 5.*m_z2k;

  Box2Abb[1801]=-83. + Box2Abb[1800]*m_z2k;

  Box2Abb[1802]=7. + 2.*Box2Abb[1801]*m_z2k;

  Box2Abb[1803]=43. + Box2Abb[1802]*m_z2k;

  Box2Abb[1804]=-96. + 25.*m_z2k;

  Box2Abb[1805]=174. + Box2Abb[1804]*m_z2k;

  Box2Abb[1806]=34. + Box2Abb[1805]*m_z2k;

  Box2Abb[1807]=-46. + Box2Abb[1806]*m_z2k;

  Box2Abb[1808]=-4. - 3.*Box2Abb[1799]*m_z12 + Box2Abb[1803]*m_z12_2 + Box2Abb[1807]*m_z12_3 + Box2Abb[1796]*m_z12_4 + 3.*Box2Abb[1794]*m_z12_5 + 10.*m_z2k - 6.*m_z2k_2;

  Box2Abb[1809]=37. - 9.*m_z2k;

  Box2Abb[1810]=-16. + 39.*m_z2k;

  Box2Abb[1811]=-7. + Box2Abb[1810]*m_z2k;

  Box2Abb[1812]=8. + Box2Abb[1811]*m_z12 + 6.*Box2Abb[187]*m_z12_2 + 2.*Box2Abb[1809]*m_z2k;

  Box2Abb[1813]=-5. + Box2Abb[1812]*m_z12 - 48.*m_z2k;

  Box2Abb[1814]=-2. + Box2Abb[1813]*m_z12;

  Box2Abb[1815]=-13. + 12.*m_z2k;

  Box2Abb[1816]=8. + Box2Abb[1815]*m_z2k;

  Box2Abb[1817]=16. - 3.*m_z2k;

  Box2Abb[1818]=-91. + 4.*Box2Abb[1817]*m_z2k;

  Box2Abb[1819]=12. + Box2Abb[1818]*m_z2k;

  Box2Abb[1820]=11. + Box2Abb[1819]*m_z2k;

  Box2Abb[1821]=16. + 9.*m_z2k;

  Box2Abb[1822]=-41. + Box2Abb[1821]*m_z2k;

  Box2Abb[1823]=59. + 2.*Box2Abb[1822]*m_z2k;

  Box2Abb[1824]=-9. + Box2Abb[1823]*m_z2k;

  Box2Abb[1825]=-34. + 3.*m_z2k;

  Box2Abb[1826]=150. + Box2Abb[1825]*m_z2k;

  Box2Abb[1827]=-176. + Box2Abb[1826]*m_z2k;

  Box2Abb[1828]=63. + Box2Abb[1827]*m_z2k;

  Box2Abb[1829]=2. + Box2Abb[1828]*m_z2k;

  Box2Abb[1830]=2.*pow(Box2Abb[61],3.) + Box2Abb[1816]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[1824]*Box2Abb[61]*m_z12_2 - Box2Abb[1829]*m_z12_3 + Box2Abb[1820]*m_z12_4 + Box2Abb[1710]*m_z12_5;

  Box2Abb[1831]=Box2Abb[1830]*m_x_2 + Box2Abb[1808]*m_x_3 + Box2Abb[1793]*m_x_4 + Box2Abb[1814]*m_x_5 + Box2Abb[1783]*pow(Box2Abb[61],2.)*m_x*m_z12 + Box2Abb[1775]*m_x_6*m_z12 + Box2Abb[155]*m_x_7*m_z12_2 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*Box2Abb[72]*m_z12_2*m_z2k;

  Box2Abb[1832]=-2. + m_z12 + 2.*m_z12*m_z2k;

  Box2Abb[1833]=-Box2Abb[1832]*m_x + m_x_2*m_z12 + Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1834]=1. + m_z12 + 4.*m_z12*m_z2k;

  Box2Abb[1835]=2. + 9.*m_z12;

  Box2Abb[1836]=11. + 5.*m_z12;

  Box2Abb[1837]=-4. + Box2Abb[1836]*m_z12;

  Box2Abb[1838]=-3. + Box2Abb[1837]*m_z12;

  Box2Abb[1839]=-1. + 2.*m_z12 + 6.*m_z12_2;

  Box2Abb[1840]=5. - 5.*m_z12 + Box2Abb[1835]*Box2Abb[4]*m_z2k + Box2Abb[1838]*m_z2k_2 + 2.*Box2Abb[1839]*m_z2k_3 - 5.*m_z12*m_z2k_4;

  Box2Abb[1841]=5. + 7.*m_z12;

  Box2Abb[1842]=2. + Box2Abb[1841]*m_z12;

  Box2Abb[1843]=-1. + m_z12 + 4.*m_z12_2;

  Box2Abb[1844]=3. + Box2Abb[1843]*m_z12;

  Box2Abb[1845]=-7. + 4.*m_z12;

  Box2Abb[1846]=1. + Box2Abb[1845]*m_z12;

  Box2Abb[1847]=-1. + m_z12 + pow(Box2Abb[4],2.)*m_z2k - Box2Abb[1842]*Box2Abb[4]*m_z2k_2 - Box2Abb[1844]*m_z2k_3 + Box2Abb[1846]*m_z2k_4 + 4.*m_z12*m_z2k_5;

  Box2Abb[1848]=1. + m_z2k + 4.*m_z2k_2;

  Box2Abb[1849]=Box2Abb[471]*pow(Box2Abb[61],2.) + Box2Abb[1848]*Box2Abb[61]*m_z12 + Box2Abb[187]*m_z12_2*m_z2k;

  Box2Abb[1850]=-4. + m_z2k + 8.*m_z12*m_z2k + 5.*m_z2k_2;

  Box2Abb[1851]=5. + Box2Abb[1850]*m_z12 + 2.*m_z2k;

  Box2Abb[1852]=17. + 20.*m_z2k;

  Box2Abb[1853]=11. - Box2Abb[1852]*m_z12;

  Box2Abb[1854]=8. + Box2Abb[1853]*m_z2k;

  Box2Abb[1855]=-8. + Box2Abb[1854]*m_z12 + m_z2k;

  Box2Abb[1856]=Box2Abb[1847]*m_x + Box2Abb[1840]*m_x_2 + Box2Abb[1855]*m_x_3 + Box2Abb[1851]*m_x_4 - Box2Abb[1834]*m_x_5 + m_x_6*m_z12 - Box2Abb[1849]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[1857]=2. + 15.*m_z12 - 28.*m_z2k;

  Box2Abb[1858]=-8. + Box2Abb[1857]*m_z12;

  Box2Abb[1859]=-73. + 14.*Box2Abb[669]*m_z12;

  Box2Abb[1860]=8. + 5.*m_z12;

  Box2Abb[1861]=74. - Box2Abb[1860]*m_z12;

  Box2Abb[1862]=54. - 151.*m_z12;

  Box2Abb[1863]=-16. + Box2Abb[1862]*m_z12;

  Box2Abb[1864]=4.*Box2Abb[841] + Box2Abb[1859]*m_z12 + Box2Abb[1861]*m_z12*m_z2k + Box2Abb[1863]*m_z2k_2 + 60.*m_z12*m_z2k_3;

  Box2Abb[1865]=2. + Box2Abb[514]*m_z12 - 6.*m_z2k + 4.*m_z2k_2;

  Box2Abb[1866]=-9. + 4.*m_z2k;

  Box2Abb[1867]=-62. + 26.*m_z2k;

  Box2Abb[1868]=5. + Box2Abb[1867]*m_z12 - 4.*m_z12_2 + 4.*Box2Abb[1866]*m_z2k;

  Box2Abb[1869]=34. + Box2Abb[1868]*m_z12 + 20.*m_z2k;

  Box2Abb[1870]=8. + Box2Abb[179]*m_z2k;

  Box2Abb[1871]=-6. + m_z2k + 28.*m_z2k_2;

  Box2Abb[1872]=-20. + m_z2k - 6.*m_z2k_2 + 70.*m_z2k_3;

  Box2Abb[1873]=-7. + 15.*m_z2k;

  Box2Abb[1874]=55. + 4.*Box2Abb[1873]*m_z2k;

  Box2Abb[1875]=95. + 2.*Box2Abb[1874]*m_z2k;

  Box2Abb[1876]=-24. + 79.*m_z12 + 2.*Box2Abb[1872]*m_z12_2 + 3.*Box2Abb[1871]*m_z12_3 + 8.*Box2Abb[1870]*m_z2k - Box2Abb[1875]*m_z12*m_z2k;

  Box2Abb[1877]=-5. + 8.*m_z2k;

  Box2Abb[1878]=2. + Box2Abb[183]*Box2Abb[1877]*m_z2k;

  Box2Abb[1879]=3. + m_z2k + 8.*m_z2k_2;

  Box2Abb[1880]=2. + Box2Abb[1879]*m_z2k;

  Box2Abb[1881]=45. - 19.*m_z2k;

  Box2Abb[1882]=-31. + Box2Abb[1881]*m_z2k;

  Box2Abb[1883]=4. + Box2Abb[1882]*m_z2k;

  Box2Abb[1884]=2. + 2.*Box2Abb[1883]*m_z2k_2;

  Box2Abb[1885]=2.*Box2Abb[230]*pow(Box2Abb[61],4.) - Box2Abb[1878]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[1884]*m_z12_2 + Box2Abb[1880]*Box2Abb[61]*m_z12_3 + 12.*m_z12_4*m_z2k_3;

  Box2Abb[1886]=10. + m_z2k + m_z2k_2 - 86.*m_z2k_3;

  Box2Abb[1887]=24. + m_z2k;

  Box2Abb[1888]=46. + Box2Abb[1887]*m_z2k;

  Box2Abb[1889]=-12. + Box2Abb[1888]*m_z2k;

  Box2Abb[1890]=1. + Box2Abb[1889]*m_z2k;

  Box2Abb[1891]=-35. + 46.*m_z2k;

  Box2Abb[1892]=19. + Box2Abb[1891]*m_z2k;

  Box2Abb[1893]=-41. + 2.*Box2Abb[1892]*m_z2k;

  Box2Abb[1894]=23. + Box2Abb[1893]*m_z2k;

  Box2Abb[1895]=-4.*Box2Abb[222]*Box2Abb[230]*pow(Box2Abb[61],2.) + Box2Abb[1894]*Box2Abb[61]*m_z12 + Box2Abb[1890]*m_z12_2 + Box2Abb[1886]*m_z12_3 - 12.*m_z12_4*m_z2k_2;

  Box2Abb[1896]=Box2Abb[1885]*m_x + Box2Abb[1895]*m_x_2 + Box2Abb[1876]*m_x_3 + Box2Abb[1864]*m_x_4 + Box2Abb[1869]*m_x_5 + Box2Abb[1858]*m_x_6 + 8.*m_x_7*m_z12 + Box2Abb[1865]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[1897]=38. + 31.*m_z12 + 12.*m_z2k;

  Box2Abb[1898]=-16. + Box2Abb[1897]*m_z12;

  Box2Abb[1899]=2. + 3.*m_z2k_2;

  Box2Abb[1900]=2. - 2.*m_z2k + 8.*m_z2k_2;

  Box2Abb[1901]=Box2Abb[1899]*Box2Abb[61]*m_z12 + Box2Abb[1900]*m_z12_2 - 2.*Box2Abb[230]*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[1902]=77. + 34.*m_z2k;

  Box2Abb[1903]=132. + 109.*m_z2k;

  Box2Abb[1904]=23. + Box2Abb[1903]*m_z12 + 8.*m_z12_2 + 2.*Box2Abb[1902]*m_z2k;

  Box2Abb[1905]=-62. + Box2Abb[1904]*m_z12 - 60.*m_z2k;

  Box2Abb[1906]=32. + 43.*m_z2k;

  Box2Abb[1907]=33. + 38.*m_z2k;

  Box2Abb[1908]=14. + Box2Abb[1907]*m_z2k;

  Box2Abb[1909]=354. + 175.*m_z2k;

  Box2Abb[1910]=194. + Box2Abb[1909]*m_z2k;

  Box2Abb[1911]=29. + 26.*m_z2k;

  Box2Abb[1912]=-19. + 6.*Box2Abb[1911]*m_z2k;

  Box2Abb[1913]=-135. + Box2Abb[1912]*m_z2k;

  Box2Abb[1914]=-2.*Box2Abb[1908] + Box2Abb[1913]*m_z12 + Box2Abb[1910]*m_z12_2 + Box2Abb[1906]*m_z12_3;

  Box2Abb[1915]=5. + Box2Abb[1163]*m_z2k_2;

  Box2Abb[1916]=5. + 6.*m_z2k;

  Box2Abb[1917]=-15. + 2.*Box2Abb[1916]*m_z2k;

  Box2Abb[1918]=10. + Box2Abb[1917]*m_z2k;

  Box2Abb[1919]=2. + Box2Abb[1918]*m_z2k;

  Box2Abb[1920]=-68. + 45.*m_z2k;

  Box2Abb[1921]=29. + Box2Abb[1920]*m_z2k;

  Box2Abb[1922]=-8. + Box2Abb[1921]*m_z2k;

  Box2Abb[1923]=-4. + Box2Abb[1922]*m_z2k;

  Box2Abb[1924]=-63. + 71.*m_z2k;

  Box2Abb[1925]=18. + Box2Abb[1924]*m_z2k;

  Box2Abb[1926]=10. + Box2Abb[1925]*m_z2k;

  Box2Abb[1927]=-2. + Box2Abb[1926]*m_z2k;

  Box2Abb[1928]=-Box2Abb[1919]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[1923]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[1927]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[230]*pow(Box2Abb[61],4.)*m_z2k + 2.*Box2Abb[1915]*m_z12_4*m_z2k;

  Box2Abb[1929]=-11. + 6.*m_z2k;

  Box2Abb[1930]=7. + Box2Abb[1929]*m_z2k;

  Box2Abb[1931]=135. + 71.*m_z2k;

  Box2Abb[1932]=50. + Box2Abb[1931]*m_z2k;

  Box2Abb[1933]=334. + 229.*m_z2k;

  Box2Abb[1934]=232. + Box2Abb[1933]*m_z2k;

  Box2Abb[1935]=106. + Box2Abb[1934]*m_z2k;

  Box2Abb[1936]=-31. + 90.*m_z2k;

  Box2Abb[1937]=-58. + Box2Abb[1936]*m_z2k;

  Box2Abb[1938]=-83. + Box2Abb[1937]*m_z2k;

  Box2Abb[1939]=-187. + 2.*Box2Abb[1938]*m_z2k;

  Box2Abb[1940]=48. + Box2Abb[1939]*m_z12 + Box2Abb[1935]*m_z12_2 + Box2Abb[1932]*m_z12_3 - 4.*Box2Abb[1930]*m_z2k + 4.*m_z12_4*m_z2k;

  Box2Abb[1941]=7. + 4.*m_z2k;

  Box2Abb[1942]=-11. + Box2Abb[1065]*Box2Abb[12]*m_z2k;

  Box2Abb[1943]=9. + Box2Abb[1942]*m_z2k;

  Box2Abb[1944]=87. + 47.*m_z2k;

  Box2Abb[1945]=127. + 2.*Box2Abb[1944]*m_z2k;

  Box2Abb[1946]=38. + Box2Abb[1945]*m_z2k;

  Box2Abb[1947]=-4. + 29.*m_z2k;

  Box2Abb[1948]=-76. + 9.*Box2Abb[1947]*m_z2k;

  Box2Abb[1949]=-66. + Box2Abb[1948]*m_z2k;

  Box2Abb[1950]=-1. + Box2Abb[1949]*m_z2k;

  Box2Abb[1951]=-119. + 50.*m_z2k;

  Box2Abb[1952]=13. + Box2Abb[1951]*m_z2k;

  Box2Abb[1953]=36. + Box2Abb[1952]*m_z2k;

  Box2Abb[1954]=13. + 2.*Box2Abb[1953]*m_z2k;

  Box2Abb[1955]=-71. + Box2Abb[1954]*m_z2k;

  Box2Abb[1956]=4.*Box2Abb[1943] + Box2Abb[1955]*m_z12 + Box2Abb[1950]*m_z12_2 + Box2Abb[1946]*m_z12_3 + 2.*Box2Abb[1941]*m_z12_4*m_z2k;

  Box2Abb[1957]=13. + 6.*m_z2k;

  Box2Abb[1958]=9. + Box2Abb[1957]*m_z2k;

  Box2Abb[1959]=3. + Box2Abb[1201]*m_z2k;

  Box2Abb[1960]=14. - 59.*m_z2k;

  Box2Abb[1961]=13. + 2.*Box2Abb[1960]*m_z2k;

  Box2Abb[1962]=-27. + Box2Abb[1961]*m_z2k;

  Box2Abb[1963]=-14. + Box2Abb[1962]*m_z2k;

  Box2Abb[1964]=19. - 34.*m_z2k + 4.*m_z2k_2;

  Box2Abb[1965]=9. + Box2Abb[1964]*m_z2k;

  Box2Abb[1966]=-34. + 3.*Box2Abb[1965]*m_z2k;

  Box2Abb[1967]=2. + Box2Abb[1966]*m_z2k;

  Box2Abb[1968]=-312. + 175.*m_z2k;

  Box2Abb[1969]=108. + Box2Abb[1968]*m_z2k;

  Box2Abb[1970]=70. + Box2Abb[1969]*m_z2k;

  Box2Abb[1971]=59. - Box2Abb[1970]*m_z2k;

  Box2Abb[1972]=18. + Box2Abb[1971]*m_z2k;

  Box2Abb[1973]=-2.*Box2Abb[1959]*pow(Box2Abb[61],2.) - Box2Abb[1967]*Box2Abb[61]*m_z12 + Box2Abb[1972]*m_z12_2 + Box2Abb[1963]*m_z12_3 - 2.*Box2Abb[1958]*m_z12_4*m_z2k;

  Box2Abb[1974]=Box2Abb[1928]*m_x + Box2Abb[1973]*m_x_2 + Box2Abb[1956]*m_x_3 - Box2Abb[1940]*m_x_4 + Box2Abb[1914]*m_x_5 - Box2Abb[1905]*m_x_6 + Box2Abb[1898]*m_x_7 - Box2Abb[1901]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[1975]=2. + 7.*m_z12 - 2.*m_z2k;

  Box2Abb[1976]=-3. + 2.*Box2Abb[1975]*m_z12;

  Box2Abb[1977]=-4. + 2.*m_z12 + 3.*m_z2k;

  Box2Abb[1978]=1. + Box2Abb[1977]*m_z2k;

  Box2Abb[1979]=36. + 4.*m_z12 + 65.*m_z2k;

  Box2Abb[1980]=-31. + Box2Abb[1979]*m_z12 + 20.*m_z2k;

  Box2Abb[1981]=7. - Box2Abb[1980]*m_z12 + 14.*m_z2k;

  Box2Abb[1982]=7. + 11.*m_z2k;

  Box2Abb[1983]=101. + 120.*m_z2k;

  Box2Abb[1984]=22. + Box2Abb[1983]*m_z2k;

  Box2Abb[1985]=-27. + Box2Abb[1567]*m_z2k;

  Box2Abb[1986]=-31. + 2.*Box2Abb[1985]*m_z2k;

  Box2Abb[1987]=22. + 2.*Box2Abb[1986]*m_z12 + Box2Abb[1984]*m_z12_2 + 2.*Box2Abb[1982]*m_z12_3 - 5.*Box2Abb[510]*m_z2k;

  Box2Abb[1988]=1. + Box2Abb[841]*m_z2k;

  Box2Abb[1989]=4. + 5.*Box2Abb[61]*m_z2k;

  Box2Abb[1990]=-2. + Box2Abb[1989]*m_z2k;

  Box2Abb[1991]=-27. + 20.*m_z2k;

  Box2Abb[1992]=28. + Box2Abb[1991]*m_z2k;

  Box2Abb[1993]=-11. + Box2Abb[1992]*m_z2k;

  Box2Abb[1994]=pow(Box2Abb[61],4.) + 4.*Box2Abb[1990]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[1993]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[179]*Box2Abb[1988]*m_z12_3;

  Box2Abb[1995]=-9. + Box2Abb[450]*Box2Abb[460]*m_z2k;

  Box2Abb[1996]=16. + 19.*m_z2k;

  Box2Abb[1997]=9. + Box2Abb[1996]*m_z2k;

  Box2Abb[1998]=83. + 115.*m_z2k;

  Box2Abb[1999]=2. + Box2Abb[1998]*m_z2k;

  Box2Abb[2000]=-12. + Box2Abb[1999]*m_z2k;

  Box2Abb[2001]=4. - 7.*m_z2k;

  Box2Abb[2002]=107. + 10.*Box2Abb[2001]*m_z2k;

  Box2Abb[2003]=58. + Box2Abb[2002]*m_z2k;

  Box2Abb[2004]=27. + Box2Abb[2003]*m_z2k;

  Box2Abb[2005]=2.*Box2Abb[1995] + Box2Abb[2004]*m_z12 - Box2Abb[2000]*m_z12_2 - 2.*Box2Abb[1997]*m_z12_3;

  Box2Abb[2006]=-16. + 5.*m_z2k;

  Box2Abb[2007]=10. - Box2Abb[179]*Box2Abb[2006]*m_z2k;

  Box2Abb[2008]=11. + 12.*m_z2k;

  Box2Abb[2009]=8. + Box2Abb[2008]*m_z2k;

  Box2Abb[2010]=5. + Box2Abb[2009]*m_z2k;

  Box2Abb[2011]=-5. + 3.*m_z2k;

  Box2Abb[2012]=13. + 7.*Box2Abb[2011]*m_z2k;

  Box2Abb[2013]=7. + 2.*Box2Abb[2012]*m_z2k;

  Box2Abb[2014]=-6. + Box2Abb[2013]*m_z2k;

  Box2Abb[2015]=-9. + 35.*m_z2k;

  Box2Abb[2016]=-14. + Box2Abb[2015]*m_z2k;

  Box2Abb[2017]=-9. + Box2Abb[2016]*m_z2k;

  Box2Abb[2018]=-8. + Box2Abb[2017]*m_z2k;

  Box2Abb[2019]=5. + 2.*Box2Abb[2018]*m_z12_2 + 2.*Box2Abb[2010]*m_z12_3 + Box2Abb[2007]*m_z2k + 2.*Box2Abb[2014]*m_z12*m_z2k;

  Box2Abb[2020]=5. + 2.*m_z2k;

  Box2Abb[2021]=-4. + Box2Abb[2020]*m_z2k;

  Box2Abb[2022]=1. + Box2Abb[2021]*m_z2k;

  Box2Abb[2023]=-8. + Box2Abb[230]*m_z2k;

  Box2Abb[2024]=5. + Box2Abb[2023]*m_z2k;

  Box2Abb[2025]=1. + Box2Abb[2024]*m_z2k;

  Box2Abb[2026]=-25. + 14.*m_z2k;

  Box2Abb[2027]=83. + 4.*Box2Abb[2026]*m_z2k;

  Box2Abb[2028]=-3. + Box2Abb[2027]*m_z2k;

  Box2Abb[2029]=-17. + Box2Abb[2028]*m_z2k;

  Box2Abb[2030]=-1. + Box2Abb[2029]*m_z2k;

  Box2Abb[2031]=-82. + 39.*m_z2k;

  Box2Abb[2032]=60. + Box2Abb[2031]*m_z2k;

  Box2Abb[2033]=22. + Box2Abb[2032]*m_z2k;

  Box2Abb[2034]=23. - Box2Abb[2033]*m_z2k;

  Box2Abb[2035]=4. + Box2Abb[2034]*m_z2k;

  Box2Abb[2036]=-Box2Abb[2022]*pow(Box2Abb[61],2.) - Box2Abb[2030]*Box2Abb[61]*m_z12 + Box2Abb[2035]*m_z12_2 - 2.*Box2Abb[2025]*m_z12_3;

  Box2Abb[2037]=Box2Abb[2036]*m_x_2 + Box2Abb[2019]*m_x_3 + Box2Abb[2005]*m_x_4 + Box2Abb[1987]*m_x_5 + Box2Abb[1981]*m_x_6 + Box2Abb[1976]*m_x_7 + m_x_8*m_z12 + Box2Abb[1994]*Box2Abb[61]*m_x*m_z2k - Box2Abb[1978]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[2038]=9. + m_z12;

  Box2Abb[2039]=7. + 6.*m_z12;

  Box2Abb[2040]=-26. + Box2Abb[2039]*m_z12;

  Box2Abb[2041]=-1. + Box2Abb[2040]*m_z12;

  Box2Abb[2042]=5. + 22.*m_z12;

  Box2Abb[2043]=-4. + Box2Abb[341]*m_z12 + 12.*m_z2k - 2.*Box2Abb[2038]*m_z12*m_z2k + Box2Abb[2041]*m_z2k_2 + 2.*Box2Abb[2042]*m_z12*m_z2k_3;

  Box2Abb[2044]=-9. + m_z12;

  Box2Abb[2045]=1. + Box2Abb[2044]*m_z12;

  Box2Abb[2046]=7. - 12.*m_z12;

  Box2Abb[2047]=22. + Box2Abb[2046]*m_z12;

  Box2Abb[2048]=1. + Box2Abb[2047]*m_z12;

  Box2Abb[2049]=15. + 26.*m_z12;

  Box2Abb[2050]=1. + Box2Abb[4]*m_z12 - 4.*m_z2k + Box2Abb[533]*m_z12*m_z2k + 2.*Box2Abb[2045]*m_z2k_2 + Box2Abb[2048]*m_z2k_3 - Box2Abb[2049]*m_z12*m_z2k_4;

  Box2Abb[2051]=-9. - 7.*m_z12 - 15.*m_z2k + 2.*m_z12*m_z2k;

  Box2Abb[2052]=1. + Box2Abb[2051]*m_z12;

  Box2Abb[2053]=5. - 13.*m_z12;

  Box2Abb[2054]=2. + 9.*m_z12 + 22.*m_z2k + 2.*Box2Abb[2053]*m_z2k_2;

  Box2Abb[2055]=2. + Box2Abb[2054]*m_z12 - m_z2k;

  Box2Abb[2056]=-1. + m_z2k + 5.*m_z2k_2;

  Box2Abb[2057]=-1. + m_z2k - 2.*m_z2k_2 + 2.*m_z2k_4;

  Box2Abb[2058]=Box2Abb[2056]*pow(Box2Abb[61],2.) + Box2Abb[2057]*m_z12 + 6.*m_z12_2*m_z2k_3;

  Box2Abb[2059]=Box2Abb[2050]*m_x_2 + Box2Abb[2043]*m_x_3 + Box2Abb[2055]*m_x_4 + Box2Abb[2052]*m_x_5 + Box2Abb[1058]*m_x_6*m_z12 + Box2Abb[2058]*m_x*m_z12*m_z2k + Box2Abb[230]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_2;

  Box2Abb[2060]=6. + 17.*m_z12;

  Box2Abb[2061]=2. + Box2Abb[841]*m_z12 - 2.*m_z2k;

  Box2Abb[2062]=49. + 4.*m_z12 + 73.*m_z2k;

  Box2Abb[2063]=-37. + Box2Abb[2062]*m_z12 + 18.*m_z2k;

  Box2Abb[2064]=2. + Box2Abb[2063]*m_z12;

  Box2Abb[2065]=-17. + 2.*m_z2k;

  Box2Abb[2066]=9. + Box2Abb[2065]*m_z2k;

  Box2Abb[2067]=-23. + Box2Abb[1982]*m_z2k;

  Box2Abb[2068]=3. + Box2Abb[2067]*m_z2k;

  Box2Abb[2069]=4. + 17.*m_z2k;

  Box2Abb[2070]=-39. + Box2Abb[2069]*m_z2k;

  Box2Abb[2071]=10. + Box2Abb[2070]*m_z2k;

  Box2Abb[2072]=2.*pow(Box2Abb[61],3.) + Box2Abb[2066]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2071]*Box2Abb[61]*m_z12_2 + Box2Abb[2068]*m_z12_3;

  Box2Abb[2073]=7. + 8.*m_z2k;

  Box2Abb[2074]=28. + 31.*m_z2k;

  Box2Abb[2075]=18. + Box2Abb[2074]*m_z2k;

  Box2Abb[2076]=-16. + 9.*Box2Abb[559]*m_z2k;

  Box2Abb[2077]=6. + Box2Abb[2076]*m_z2k;

  Box2Abb[2078]=32. + 25.*m_z2k;

  Box2Abb[2079]=13. + Box2Abb[2078]*m_z2k;

  Box2Abb[2080]=26. + Box2Abb[2079]*m_z2k;

  Box2Abb[2081]=28. - 2.*Box2Abb[2080]*m_z12 + Box2Abb[2077]*m_z12_2 + Box2Abb[2075]*m_z12_3 + 2.*Box2Abb[2073]*m_z2k;

  Box2Abb[2082]=-113. + 2.*m_z2k;

  Box2Abb[2083]=106. + 109.*m_z2k;

  Box2Abb[2084]=44. + Box2Abb[2083]*m_z2k;

  Box2Abb[2085]=-84. + Box2Abb[2084]*m_z12 + 7.*Box2Abb[150]*m_z12_2 + Box2Abb[2082]*m_z2k;

  Box2Abb[2086]=34. + Box2Abb[2085]*m_z12 + 10.*m_z2k;

  Box2Abb[2087]=5. + 3.*m_z2k_2;

  Box2Abb[2088]=-17. + 4.*m_z2k;

  Box2Abb[2089]=-2. + Box2Abb[2088]*m_z2k;

  Box2Abb[2090]=53. - 35.*m_z2k;

  Box2Abb[2091]=15. + Box2Abb[2090]*m_z2k;

  Box2Abb[2092]=16. + Box2Abb[2091]*m_z2k;

  Box2Abb[2093]=-6. + Box2Abb[2092]*m_z2k;

  Box2Abb[2094]=56. + 45.*m_z2k;

  Box2Abb[2095]=50. + Box2Abb[2094]*m_z2k;

  Box2Abb[2096]=36. + Box2Abb[2095]*m_z2k;

  Box2Abb[2097]=9. + Box2Abb[2096]*m_z2k;

  Box2Abb[2098]=10. + 2.*Box2Abb[2093]*m_z12 - Box2Abb[2097]*m_z12_2 + 2.*Box2Abb[12]*Box2Abb[2087]*m_z12_3 + 2.*Box2Abb[2089]*m_z2k;

  Box2Abb[2099]=5. + m_z2k;

  Box2Abb[2100]=-1. + Box2Abb[2099]*m_z2k;

  Box2Abb[2101]=5. + 9.*m_z2k;

  Box2Abb[2102]=15. + Box2Abb[2101]*m_z2k;

  Box2Abb[2103]=-3. + Box2Abb[2102]*m_z2k;

  Box2Abb[2104]=-1. + Box2Abb[2103]*m_z2k;

  Box2Abb[2105]=-93. + 34.*m_z2k;

  Box2Abb[2106]=21. + Box2Abb[2105]*m_z2k;

  Box2Abb[2107]=29. + Box2Abb[2106]*m_z2k;

  Box2Abb[2108]=-1. + Box2Abb[2107]*m_z2k;

  Box2Abb[2109]=5. + 53.*m_z2k;

  Box2Abb[2110]=-66. + Box2Abb[2109]*m_z2k;

  Box2Abb[2111]=-24. + Box2Abb[2110]*m_z2k;

  Box2Abb[2112]=25. + Box2Abb[2111]*m_z2k;

  Box2Abb[2113]=3. + Box2Abb[2112]*m_z2k;

  Box2Abb[2114]=2.*Box2Abb[2100]*pow(Box2Abb[61],2.) + Box2Abb[2108]*Box2Abb[61]*m_z12 + Box2Abb[2113]*m_z12_2 + 2.*Box2Abb[2104]*m_z12_3;

  Box2Abb[2115]=Box2Abb[2114]*m_x_2 + Box2Abb[2098]*m_x_3 - Box2Abb[2081]*m_x_4 + Box2Abb[2086]*m_x_5 - Box2Abb[2064]*m_x_6 + Box2Abb[2060]*m_x_7*m_z12 - Box2Abb[2072]*Box2Abb[61]*m_x*m_z2k + Box2Abb[2061]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[2116]=-1. + m_z2k - 2.*m_z12*m_z2k;

  Box2Abb[2117]=-1. + m_z12 + 6.*m_z2k - 2.*Box2Abb[12]*m_z12*m_z2k;

  Box2Abb[2118]=Box2Abb[2117]*m_x_2 + Box2Abb[1424]*m_x_3 + m_x_4*m_z12 + Box2Abb[2116]*m_x*m_z2k + pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[2119]=-2. + 5.*m_z12 - 20.*m_z2k;

  Box2Abb[2120]=-4. + Box2Abb[2119]*m_z12;

  Box2Abb[2121]=-2. - Box2Abb[179]*m_z12 + m_z2k + m_z2k_2;

  Box2Abb[2122]=-17. + 3.*m_z2k;

  Box2Abb[2123]=-5. + Box2Abb[2122]*m_z12 - 4.*m_z2k + 40.*m_z2k_2;

  Box2Abb[2124]=6.*Box2Abb[833] + Box2Abb[2123]*m_z12;

  Box2Abb[2125]=21. - 38.*m_z2k_2;

  Box2Abb[2126]=9. + 4.*m_z2k + 6.*m_z2k_2;

  Box2Abb[2127]=15. + 8.*Box2Abb[509]*m_z2k;

  Box2Abb[2128]=-2. + Box2Abb[2127]*m_z2k;

  Box2Abb[2129]=-2.*Box2Abb[2126] + Box2Abb[2128]*m_z12 + Box2Abb[2125]*m_z12_2 + 3.*m_z12_3*m_z2k;

  Box2Abb[2130]=-7. + 11.*m_z2k;

  Box2Abb[2131]=7. - 2.*m_z2k + 4.*m_z2k_2;

  Box2Abb[2132]=2. + Box2Abb[2131]*m_z2k;

  Box2Abb[2133]=-2. + 3.*Box2Abb[854]*m_z2k;

  Box2Abb[2134]=2. + Box2Abb[2133]*m_z2k;

  Box2Abb[2135]=Box2Abb[2132]*pow(Box2Abb[61],2.) + Box2Abb[2134]*Box2Abb[61]*m_z12 + Box2Abb[12]*Box2Abb[2130]*m_z12_2*m_z2k;

  Box2Abb[2136]=-8. + 7.*m_z2k;

  Box2Abb[2137]=-1. - 8.*m_z2k + 46.*m_z2k_2;

  Box2Abb[2138]=-11. + Box2Abb[2137]*m_z2k;

  Box2Abb[2139]=5. + 4.*Box2Abb[493]*m_z2k;

  Box2Abb[2140]=4. + Box2Abb[2139]*m_z2k;

  Box2Abb[2141]=7. + Box2Abb[2140]*m_z2k;

  Box2Abb[2142]=2.*Box2Abb[179]*Box2Abb[444]*Box2Abb[61] + Box2Abb[2141]*m_z12 + Box2Abb[2138]*m_z12_2 + Box2Abb[2136]*m_z12_3*m_z2k;

  Box2Abb[2143]=Box2Abb[2142]*m_x_2 + Box2Abb[2129]*m_x_3 + Box2Abb[2124]*m_x_4 + Box2Abb[2120]*m_x_5 - Box2Abb[2135]*m_x*m_z12 + 4.*m_x_6*m_z12 - Box2Abb[2121]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k;

  Box2Abb[2144]=-2. + m_z12 + 2.*m_z2k;

  Box2Abb[2145]=-2. + 4.*m_z2k;

  Box2Abb[2146]=3. + 16.*m_z2k;

  Box2Abb[2147]=m_z12 + Box2Abb[2145]*m_z12_2 - m_z2k - Box2Abb[2146]*m_z12*m_z2k;

  Box2Abb[2148]=-1. + m_z2k + 8.*m_z2k_2;

  Box2Abb[2149]=1. + m_z2k + 8.*m_z2k_2;

  Box2Abb[2150]=1. + 4.*Box2Abb[59]*m_z2k;

  Box2Abb[2151]=-2.*pow(Box2Abb[61],3.) + 2.*Box2Abb[2148]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2150]*Box2Abb[61]*m_z12_2 + Box2Abb[2149]*m_z12_3;

  Box2Abb[2152]=7. + 2.*m_z2k;

  Box2Abb[2153]=3. + 8.*m_z2k;

  Box2Abb[2154]=1. + Box2Abb[61]*m_z2k;

  Box2Abb[2155]=-3. + 20.*Box2Abb[444]*m_z2k;

  Box2Abb[2156]=5. - 4.*m_z2k;

  Box2Abb[2157]=23. + 4.*Box2Abb[2156]*m_z2k;

  Box2Abb[2158]=-4. + Box2Abb[2157]*m_z2k;

  Box2Abb[2159]=-2. - 2.*Box2Abb[2154]*Box2Abb[2155]*m_z12 + Box2Abb[2158]*m_z12_2 + 2.*Box2Abb[2152]*m_z2k + Box2Abb[2153]*m_z12_3*m_z2k;

  Box2Abb[2160]=-2. + Box2Abb[280]*m_z2k;

  Box2Abb[2161]=12. + m_z2k + 35.*m_z2k_2;

  Box2Abb[2162]=-3. + Box2Abb[2161]*m_z2k;

  Box2Abb[2163]=1. + 2.*Box2Abb[2162]*m_z12 - 3.*Box2Abb[2160]*m_z12_2 + 4.*Box2Abb[187]*m_z2k - m_z12_3*m_z2k;

  Box2Abb[2164]=-1. + 4.*Box2Abb[222]*m_z2k;

  Box2Abb[2165]=-11. + 25.*m_z2k;

  Box2Abb[2166]=3. + Box2Abb[2165]*m_z2k;

  Box2Abb[2167]=-11. + Box2Abb[2166]*m_z2k;

  Box2Abb[2168]=1. + Box2Abb[2167]*m_z2k;

  Box2Abb[2169]=-24. + 43.*m_z2k;

  Box2Abb[2170]=18. + Box2Abb[2169]*m_z2k;

  Box2Abb[2171]=-6. + Box2Abb[2170]*m_z2k;

  Box2Abb[2172]=1. + Box2Abb[2171]*m_z2k;

  Box2Abb[2173]=-Box2Abb[2164]*pow(Box2Abb[61],2.) + 2.*Box2Abb[2168]*Box2Abb[61]*m_z12 + Box2Abb[2172]*m_z12_2 - Box2Abb[1580]*m_z12_3*m_z2k;

  Box2Abb[2174]=Box2Abb[2173]*m_x_2 + Box2Abb[2159]*m_x_3 + Box2Abb[2163]*m_x_4 + 2.*Box2Abb[2147]*m_x_5 + Box2Abb[1500]*m_x_6*m_z12 - Box2Abb[2151]*Box2Abb[61]*m_x*m_z2k + Box2Abb[2144]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[2175]=-2. + 3.*m_z12 + 2.*m_z2k;

  Box2Abb[2176]=14. + m_z12;

  Box2Abb[2177]=-16. + Box2Abb[2176]*m_z12;

  Box2Abb[2178]=20. + Box2Abb[2177]*m_z12;

  Box2Abb[2179]=52. + m_z12;

  Box2Abb[2180]=-26. + Box2Abb[2179]*m_z12;

  Box2Abb[2181]=3. - 22.*m_z12 + 26.*m_z12_2 + Box2Abb[2178]*m_z2k - Box2Abb[2180]*m_z2k_2 + 18.*m_z12*m_z2k_3;

  Box2Abb[2182]=28. + m_z12;

  Box2Abb[2183]=-14. + Box2Abb[2182]*m_z12;

  Box2Abb[2184]=-24. + Box2Abb[2183]*m_z12;

  Box2Abb[2185]=-98. + m_z12;

  Box2Abb[2186]=44. + Box2Abb[2185]*m_z12;

  Box2Abb[2187]=8. + 6.*Box2Abb[817]*m_z12 + m_z2k + Box2Abb[252]*Box2Abb[9]*m_z12*m_z2k - Box2Abb[2184]*m_z2k_2 + Box2Abb[2186]*m_z2k_3 + 10.*m_z12*m_z2k_4;

  Box2Abb[2188]=14. + 5.*m_z2k;

  Box2Abb[2189]=-6. + Box2Abb[2188]*m_z12 + 10.*Box2Abb[61]*m_z2k;

  Box2Abb[2190]=Box2Abb[2189]*m_z12 + 6.*m_z2k;

  Box2Abb[2191]=-5. + m_z2k + 5.*m_z2k_2;

  Box2Abb[2192]=1. + Box2Abb[2191]*m_z2k;

  Box2Abb[2193]=7. + m_z2k;

  Box2Abb[2194]=-4. + Box2Abb[2193]*m_z2k;

  Box2Abb[2195]=2. + Box2Abb[2194]*m_z2k;

  Box2Abb[2196]=8. + 11.*m_z2k;

  Box2Abb[2197]=-12. + Box2Abb[2196]*m_z2k;

  Box2Abb[2198]=4. + Box2Abb[2197]*m_z2k;

  Box2Abb[2199]=2.*Box2Abb[2192]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2198]*Box2Abb[61]*m_z12_2 + Box2Abb[2195]*m_z12_3 - 2.*pow(Box2Abb[61],3.)*m_z2k;

  Box2Abb[2200]=-2. + 6.*m_z2k - 26.*m_z2k_3 + 30.*m_z2k_4 + m_z2k_5 - 9.*m_z2k_6;

  Box2Abb[2201]=-9. + 2.*m_z2k;

  Box2Abb[2202]=-3. + 2.*Box2Abb[2201]*m_z2k;

  Box2Abb[2203]=7. + Box2Abb[2202]*m_z2k;

  Box2Abb[2204]=-3. + 2.*Box2Abb[1071]*m_z2k;

  Box2Abb[2205]=2. + Box2Abb[2204]*m_z2k;

  Box2Abb[2206]=47. + m_z2k;

  Box2Abb[2207]=-2. + Box2Abb[2206]*m_z2k;

  Box2Abb[2208]=-10. + Box2Abb[2207]*m_z2k;

  Box2Abb[2209]=2. + Box2Abb[2208]*m_z2k;

  Box2Abb[2210]=Box2Abb[2205]*pow(Box2Abb[61],2.) + 2.*Box2Abb[2200]*m_z12 - Box2Abb[2209]*Box2Abb[61]*m_z12_2 + Box2Abb[2203]*m_z12_3*m_z2k;

  Box2Abb[2211]=5. - 9.*m_z2k;

  Box2Abb[2212]=-5. + 4.*Box2Abb[2211]*m_z2k;

  Box2Abb[2213]=12. + Box2Abb[2212]*m_z2k;

  Box2Abb[2214]=68. - 9.*m_z2k;

  Box2Abb[2215]=15. + Box2Abb[2214]*m_z2k;

  Box2Abb[2216]=8. + Box2Abb[2215]*m_z2k;

  Box2Abb[2217]=-11. + Box2Abb[2216]*m_z2k;

  Box2Abb[2218]=36. + 5.*m_z2k;

  Box2Abb[2219]=-39. + Box2Abb[2218]*m_z2k;

  Box2Abb[2220]=-11. + Box2Abb[2219]*m_z2k;

  Box2Abb[2221]=-9. + Box2Abb[2220]*m_z2k;

  Box2Abb[2222]=9. + Box2Abb[2221]*m_z2k;

  Box2Abb[2223]=-7. + 2.*Box2Abb[2222]*m_z12 + Box2Abb[2217]*m_z12_2 + Box2Abb[2213]*m_z2k - Box2Abb[1405]*m_z12_3*m_z2k;

  Box2Abb[2224]=Box2Abb[2210]*m_x_2 + Box2Abb[2223]*m_x_3 + Box2Abb[2187]*m_x_4 - Box2Abb[2181]*m_x_5 + Box2Abb[2190]*m_x_6 - Box2Abb[755]*m_x_7*m_z12 + Box2Abb[2199]*Box2Abb[61]*m_x*m_z2k - Box2Abb[2175]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_3;

  Box2Abb[2225]=m_x + Box2Abb[1]*m_x*m_z12;

  Box2Abb[2226]=-4. + m_x + 12.*m_x_2;

  Box2Abb[2227]=-1. + 3.*m_x + 8.*m_x_2 + Box2Abb[1]*Box2Abb[2226]*m_z12 - 3.*pow(Box2Abb[1],2.)*m_z12_2;

  Box2Abb[2228]=5. - 2.*m_x;

  Box2Abb[2229]=10. + Box2Abb[2228]*m_x;

  Box2Abb[2230]=9. + 2.*Box2Abb[2229]*m_x;

  Box2Abb[2231]=-9. + Box2Abb[2230]*m_x;

  Box2Abb[2232]=-15. + 6.*m_x + 4.*m_x_2;

  Box2Abb[2233]=-4. + Box2Abb[2232]*m_x;

  Box2Abb[2234]=-32. + Box2Abb[2233]*m_x;

  Box2Abb[2235]=32. + Box2Abb[2234]*m_x;

  Box2Abb[2236]=-11. + 31.*m_x;

  Box2Abb[2237]=5. + Box2Abb[2236]*m_x;

  Box2Abb[2238]=-33. + Box2Abb[2237]*m_x;

  Box2Abb[2239]=-2. + Box2Abb[2238]*m_x;

  Box2Abb[2240]=1. + 10.*Box2Abb[210]*m_x;

  Box2Abb[2241]=Box2Abb[2231]*m_x + m_z12 + Box2Abb[2235]*m_x*m_z12 - Box2Abb[1]*Box2Abb[2239]*m_z12_2 + pow(Box2Abb[1],2.)*Box2Abb[2240]*m_z12_3;

  Box2Abb[2242]=-1. + 8.*m_x;

  Box2Abb[2243]=2. + Box2Abb[2242]*m_x;

  Box2Abb[2244]=25. + 2.*Box2Abb[2243]*m_x;

  Box2Abb[2245]=7. + 4.*m_x;

  Box2Abb[2246]=34. + 3.*Box2Abb[2245]*m_x;

  Box2Abb[2247]=22. + Box2Abb[2246]*m_x;

  Box2Abb[2248]=77. + 2.*Box2Abb[2247]*m_x;

  Box2Abb[2249]=6. + Box2Abb[2248]*m_x;

  Box2Abb[2250]=-13. + 17.*m_x;

  Box2Abb[2251]=32. + Box2Abb[2250]*m_x;

  Box2Abb[2252]=5. + Box2Abb[210]*Box2Abb[2251]*m_x;

  Box2Abb[2253]=-1. + 3.*m_x;

  Box2Abb[2254]=8. + Box2Abb[2253]*m_x;

  Box2Abb[2255]=2. + Box2Abb[2254]*m_x;

  Box2Abb[2256]=-Box2Abb[2244]*m_x + Box2Abb[2249]*m_z12 - 2.*Box2Abb[2252]*m_z12_2 + 2.*Box2Abb[2255]*m_z12_3;

  Box2Abb[2257]=31. + 26.*m_x + 24.*m_x_2;

  Box2Abb[2258]=9. + 5.*m_x;

  Box2Abb[2259]=12. + Box2Abb[2258]*m_x;

  Box2Abb[2260]=9. + Box2Abb[2259]*m_x;

  Box2Abb[2261]=3. + 2.*Box2Abb[2260]*m_x;

  Box2Abb[2262]=61. + 31.*m_x;

  Box2Abb[2263]=82. + Box2Abb[2262]*m_x;

  Box2Abb[2264]=28. + Box2Abb[2263]*m_x;

  Box2Abb[2265]=11. + 10.*m_x;

  Box2Abb[2266]=Box2Abb[2257]*m_x - 6.*Box2Abb[2261]*m_z12 + Box2Abb[2264]*m_z12_2 - Box2Abb[2265]*m_z12_3;

  Box2Abb[2267]=9. + 8.*m_x;

  Box2Abb[2268]=33. + 20.*m_x;

  Box2Abb[2269]=28. + Box2Abb[2268]*m_x;

  Box2Abb[2270]=8. + Box2Abb[2269]*m_x;

  Box2Abb[2271]=15. + 8.*m_x;

  Box2Abb[2272]=44. + 5.*Box2Abb[2271]*m_x;

  Box2Abb[2273]=2. + m_x;

  Box2Abb[2274]=-2.*Box2Abb[2267]*m_x + 4.*Box2Abb[2270]*m_z12 - Box2Abb[2272]*m_z12_2 + 7.*Box2Abb[2273]*m_z12_3;

  Box2Abb[2275]=11. + 26.*m_x + 20.*m_x_2;

  Box2Abb[2276]=34. + 33.*m_x;

  Box2Abb[2277]=4.*m_x - 3.*Box2Abb[2275]*m_z12 + Box2Abb[2276]*m_z12_2 - 6.*m_z12_3;

  Box2Abb[2278]=9. + 12.*m_x - 5.*m_z12;

  Box2Abb[2279]=pow(Box2Abb[1],3.)*pow(Box2Abb[2225],2.) - Box2Abb[1]*Box2Abb[212]*Box2Abb[2227]*m_x*m_z2k - Box2Abb[2241]*m_z2k_2 + Box2Abb[2256]*m_z2k_3 + Box2Abb[2266]*m_z2k_4 + Box2Abb[2274]*m_z2k_5 + Box2Abb[2277]*m_z2k_6 + 2.*Box2Abb[2278]*m_z12*m_z2k_7 - 4.*m_z12*m_z2k_8;

  Box2Abb[2280]=1. + m_z12 - m_z2k;

  Box2Abb[2281]=15. + Box2Abb[1485]*Box2Abb[533]*m_z12;

  Box2Abb[2282]=64. + 9.*m_z12;

  Box2Abb[2283]=18. + Box2Abb[2282]*m_z12;

  Box2Abb[2284]=-8. + Box2Abb[2283]*m_z12;

  Box2Abb[2285]=120. + 23.*m_z12;

  Box2Abb[2286]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box2Abb[2281]*m_z2k + Box2Abb[2284]*m_z2k_2 + Box2Abb[2285]*m_z12*m_z2k_3;

  Box2Abb[2287]=10. + 3.*m_z12;

  Box2Abb[2288]=-10. + Box2Abb[2287]*m_z12;

  Box2Abb[2289]=12. + Box2Abb[2288]*m_z12;

  Box2Abb[2290]=19. + 10.*m_z12;

  Box2Abb[2291]=32. + Box2Abb[2290]*m_z12;

  Box2Abb[2292]=-15. + Box2Abb[2291]*m_z12;

  Box2Abb[2293]=9. - 4.*Box2Abb[734]*m_z12;

  Box2Abb[2294]=2. + Box2Abb[2293]*m_z12;

  Box2Abb[2295]=-80. + 23.*m_z12;

  Box2Abb[2296]=pow(Box2Abb[4],2.) - Box2Abb[2289]*m_z2k - Box2Abb[2292]*m_z2k_2 + 2.*Box2Abb[2294]*m_z2k_3 + Box2Abb[2295]*m_z12*m_z2k_4;

  Box2Abb[2297]=-4. + 9.*m_z2k;

  Box2Abb[2298]=2. + Box2Abb[2297]*m_z12 + 20.*m_z2k;

  Box2Abb[2299]=20. + Box2Abb[2182]*m_z12;

  Box2Abb[2300]=80. + 33.*m_z12;

  Box2Abb[2301]=6. - 6.*m_z12 + Box2Abb[2299]*m_z2k + Box2Abb[2300]*m_z2k_2;

  Box2Abb[2302]=-1. + Box2Abb[2301]*m_z12 - 4.*m_z2k;

  Box2Abb[2303]=-10. + 9.*m_z2k;

  Box2Abb[2304]=-1. + 9.*m_z2k;

  Box2Abb[2305]=2. + Box2Abb[2304]*m_z2k;

  Box2Abb[2306]=2.*pow(Box2Abb[61],2.) + Box2Abb[2305]*m_z12_2 + Box2Abb[2303]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[2307]=-1. + 2.*m_z2k_2;

  Box2Abb[2308]=64. - 33.*m_z2k;

  Box2Abb[2309]=-19. + Box2Abb[2308]*m_z2k;

  Box2Abb[2310]=12. + Box2Abb[2309]*m_z2k;

  Box2Abb[2311]=7. - 4.*m_z2k;

  Box2Abb[2312]=3. + 2.*Box2Abb[2311]*m_z2k;

  Box2Abb[2313]=1. + Box2Abb[2312]*m_z2k;

  Box2Abb[2314]=pow(Box2Abb[61],2.) + Box2Abb[2313]*m_z12_3 + 10.*Box2Abb[2307]*Box2Abb[61]*m_z12*m_z2k + Box2Abb[2310]*m_z12_2*m_z2k;

  Box2Abb[2315]=Box2Abb[2296]*m_x_3 + Box2Abb[2286]*m_x_4 - Box2Abb[2302]*m_x_5 + Box2Abb[2298]*m_x_6*m_z12 + m_x_7*m_z12_2 + Box2Abb[2314]*m_x_2*m_z2k + Box2Abb[2306]*Box2Abb[61]*m_x*m_z12*m_z2k_2 - Box2Abb[2280]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_3;

  Box2Abb[2316]=2. + 7.*m_x;

  Box2Abb[2317]=1. - m_x + 10.*m_x_2 + Box2Abb[1]*Box2Abb[2316]*Box2Abb[768]*m_z12 - 3.*pow(Box2Abb[1],2.)*m_z12_2;

  Box2Abb[2318]=5. - 26.*m_x;

  Box2Abb[2319]=-7. + Box2Abb[2318]*m_x;

  Box2Abb[2320]=17. - 8.*m_x;

  Box2Abb[2321]=2. + Box2Abb[2320]*m_x;

  Box2Abb[2322]=-11. + 2.*Box2Abb[2321]*m_x;

  Box2Abb[2323]=4. + Box2Abb[2322]*m_x;

  Box2Abb[2324]=-59. + 40.*m_x;

  Box2Abb[2325]=7. + Box2Abb[2324]*m_x;

  Box2Abb[2326]=3. + m_x + Box2Abb[2325]*m_x_2;

  Box2Abb[2327]=8. + Box2Abb[2326]*m_x;

  Box2Abb[2328]=1. + 2.*m_x + 4.*m_x_2;

  Box2Abb[2329]=2. + Box2Abb[2319]*m_x - 7.*m_z12 + Box2Abb[2323]*m_x*m_z12 + Box2Abb[2327]*m_z12_2 - 3.*pow(Box2Abb[1],2.)*Box2Abb[2328]*m_z12_3;

  Box2Abb[2330]=38. - 15.*m_z12;

  Box2Abb[2331]=-7. + 5.*m_z12;

  Box2Abb[2332]=5. + 2.*Box2Abb[2331]*m_z12;

  Box2Abb[2333]=8. + Box2Abb[702]*m_z12;

  Box2Abb[2334]=19. + 7.*m_z12;

  Box2Abb[2335]=23. - Box2Abb[2334]*m_z12;

  Box2Abb[2336]=1. + Box2Abb[2335]*m_z12;

  Box2Abb[2337]=-37. + 9.*m_z12;

  Box2Abb[2338]=74. + Box2Abb[2337]*m_z12;

  Box2Abb[2339]=-4. + Box2Abb[2338]*m_z12;

  Box2Abb[2340]=2.*Box2Abb[2332]*Box2Abb[4]*m_x + Box2Abb[2336]*m_x_2 + 4.*Box2Abb[2333]*m_x_3 + Box2Abb[2339]*m_x_4 - 2.*pow(Box2Abb[4],2.)*m_z12 + 2.*Box2Abb[2330]*m_x_5*m_z12;

  Box2Abb[2341]=17. + 12.*m_x;

  Box2Abb[2342]=18. + Box2Abb[2341]*m_x;

  Box2Abb[2343]=99. + 70.*m_x;

  Box2Abb[2344]=70. + Box2Abb[2343]*m_x;

  Box2Abb[2345]=61. + 2.*Box2Abb[2344]*m_x;

  Box2Abb[2346]=97. - 30.*m_x;

  Box2Abb[2347]=116. + Box2Abb[2346]*m_x;

  Box2Abb[2348]=85. + Box2Abb[2347]*m_x;

  Box2Abb[2349]=-20. + Box2Abb[2348]*m_x;

  Box2Abb[2350]=-4. + 3.*m_x;

  Box2Abb[2351]=-33. + 4.*Box2Abb[2350]*m_x;

  Box2Abb[2352]=8. + Box2Abb[2351]*m_x;

  Box2Abb[2353]=Box2Abb[2342]*m_x + 12.*m_z12 - Box2Abb[2345]*m_x*m_z12 + Box2Abb[2349]*m_z12_2 + Box2Abb[2352]*m_z12_3;

  Box2Abb[2354]=7. + 6.*m_x;

  Box2Abb[2355]=2. + 5.*m_x;

  Box2Abb[2356]=-7. + 2.*Box2Abb[1525]*Box2Abb[2355]*m_x;

  Box2Abb[2357]=-49. + 54.*m_x;

  Box2Abb[2358]=-50. + Box2Abb[2357]*m_x;

  Box2Abb[2359]=31. + Box2Abb[2358]*m_x;

  Box2Abb[2360]=8. - 15.*m_x;

  Box2Abb[2361]=-7. + Box2Abb[2360]*m_x;

  Box2Abb[2362]=-2.*Box2Abb[2354]*m_x + 4.*Box2Abb[2356]*m_z12 + Box2Abb[2359]*m_z12_2 + Box2Abb[2361]*m_z12_3;

  Box2Abb[2363]=13. + 2.*m_z12;

  Box2Abb[2364]=32. - Box2Abb[2363]*m_z12;

  Box2Abb[2365]=2. + 15.*m_z12;

  Box2Abb[2366]=4. + Box2Abb[2365]*m_z12;

  Box2Abb[2367]=Box2Abb[2366]*m_x + Box2Abb[2364]*m_z12 - 8.*Box2Abb[1058]*m_x_2*m_z12;

  Box2Abb[2368]=-18. - 2.*Box2Abb[1637]*m_x + Box2Abb[485]*m_z12;

  Box2Abb[2369]=pow(Box2Abb[1],3.)*pow(Box2Abb[212],2.)*m_x_3 - Box2Abb[1]*Box2Abb[212]*Box2Abb[2317]*m_x_2*m_z2k + Box2Abb[2329]*m_x*m_z2k_2 + Box2Abb[2340]*m_z2k_3 + Box2Abb[2353]*m_z2k_4 + Box2Abb[2362]*m_z2k_5 + Box2Abb[2367]*m_z2k_6 + Box2Abb[2368]*m_z12*m_z2k_7 + Box2Abb[252]*m_z12*m_z2k_8;

  Box2Abb[2370]=2. + m_z12 - 5.*m_z2k + 2.*m_z12*m_z2k + 3.*m_z2k_2;

  Box2Abb[2371]=10. + m_z12 + 12.*m_z2k;

  Box2Abb[2372]=-2. + Box2Abb[2371]*m_z12;

  Box2Abb[2373]=-Box2Abb[2372]*m_x_2 + 2.*Box2Abb[2370]*m_x*m_z12 + 6.*m_x_3*m_z12 - pow(Box2Abb[61],2.)*m_z12_2;

  Box2Abb[2374]=18. - 7.*m_z12;

  Box2Abb[2375]=-7. + Box2Abb[2374]*m_z12;

  Box2Abb[2376]=-9. + 2.*m_z12;

  Box2Abb[2377]=4. + Box2Abb[2376]*Box2Abb[710]*m_z12;

  Box2Abb[2378]=-51. + 32.*m_z12;

  Box2Abb[2379]=12. + Box2Abb[2378]*m_z12;

  Box2Abb[2380]=4. + Box2Abb[2379]*m_z12;

  Box2Abb[2381]=-4. + Box2Abb[669]*m_z12;

  Box2Abb[2382]=-31. + 5.*Box2Abb[2381]*m_z12;

  Box2Abb[2383]=27. + Box2Abb[2382]*m_z12;

  Box2Abb[2384]=-5. + Box2Abb[2383]*m_z12;

  Box2Abb[2385]=-59. + 64.*m_z12;

  Box2Abb[2386]=144. + Box2Abb[2385]*m_z12;

  Box2Abb[2387]=-117. + Box2Abb[2386]*m_z12;

  Box2Abb[2388]=8. + Box2Abb[2387]*m_z12;

  Box2Abb[2389]=-108. + 43.*m_z12;

  Box2Abb[2390]=108. + Box2Abb[2389]*m_z12;

  Box2Abb[2391]=-2. + Box2Abb[2390]*m_z12;

  Box2Abb[2392]=-4. + Box2Abb[2375]*m_z12 + 8.*m_z2k - Box2Abb[2377]*m_z12*m_z2k + 2.*Box2Abb[2380]*m_z12*m_z2k_2 + 2.*Box2Abb[2384]*m_z2k_3 + Box2Abb[2388]*m_z2k_4 + Box2Abb[2391]*m_z2k_5 + 14.*Box2Abb[538]*m_z12*m_z2k_6;

  Box2Abb[2393]=-11. + m_z12 - 17.*m_z2k;

  Box2Abb[2394]=27. + 2.*Box2Abb[2393]*m_z12 + 42.*m_z2k;

  Box2Abb[2395]=2. + Box2Abb[2394]*m_z12;

  Box2Abb[2396]=6. + 5.*m_z2k;

  Box2Abb[2397]=-11. + 17.*m_z2k;

  Box2Abb[2398]=44. + 80.*m_z2k + 98.*m_z2k_2;

  Box2Abb[2399]=-41. + Box2Abb[2398]*m_z12 + Box2Abb[2397]*m_z12_2 - 18.*Box2Abb[1567]*m_z2k;

  Box2Abb[2400]=-2.*Box2Abb[2396] + Box2Abb[2399]*m_z12;

  Box2Abb[2401]=3. + Box2Abb[187]*m_z2k_2;

  Box2Abb[2402]=-8. + 5.*Box2Abb[833]*m_z2k;

  Box2Abb[2403]=1. + Box2Abb[2402]*m_z2k;

  Box2Abb[2404]=-10. + Box2Abb[2101]*m_z2k;

  Box2Abb[2405]=4. + Box2Abb[2404]*m_z2k;

  Box2Abb[2406]=Box2Abb[2401]*pow(Box2Abb[61],3.) + Box2Abb[2405]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2403]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[187]*m_z12_3*m_z2k_2;

  Box2Abb[2407]=33. + 53.*m_z2k;

  Box2Abb[2408]=25. - 2.*Box2Abb[2407]*m_z2k;

  Box2Abb[2409]=9. + 14.*m_z2k;

  Box2Abb[2410]=31. + 5.*Box2Abb[2409]*m_z2k;

  Box2Abb[2411]=22. + 3.*Box2Abb[2410]*m_z2k;

  Box2Abb[2412]=36. + 77.*m_z2k;

  Box2Abb[2413]=9. + Box2Abb[2412]*m_z2k;

  Box2Abb[2414]=59. + 2.*Box2Abb[2413]*m_z2k;

  Box2Abb[2415]=30. + Box2Abb[2411]*m_z12 - Box2Abb[2414]*m_z12_2 + Box2Abb[2408]*m_z12_3 + 4.*Box2Abb[887]*m_z2k + m_z12_4*m_z2k;

  Box2Abb[2416]=15. + 41.*m_z2k + 105.*m_z2k_3;

  Box2Abb[2417]=1. + Box2Abb[2416]*m_z2k;

  Box2Abb[2418]=5. + 2.*Box2Abb[559]*m_z2k;

  Box2Abb[2419]=18. + Box2Abb[2418]*m_z2k;

  Box2Abb[2420]=192. + 205.*m_z2k;

  Box2Abb[2421]=106. + Box2Abb[2420]*m_z2k;

  Box2Abb[2422]=-30. + Box2Abb[2421]*m_z2k;

  Box2Abb[2423]=-9. + 5.*Box2Abb[230]*m_z2k;

  Box2Abb[2424]=-95. + 14.*Box2Abb[2423]*m_z2k;

  Box2Abb[2425]=61. + Box2Abb[2424]*m_z2k;

  Box2Abb[2426]=-2.*Box2Abb[2419] - 2.*Box2Abb[2417]*m_z12 + Box2Abb[2425]*m_z12_2 + Box2Abb[2422]*m_z12_3 + Box2Abb[1460]*m_z12_4*m_z2k;

  Box2Abb[2427]=-32. - Box2Abb[1359]*Box2Abb[280]*m_z2k;

  Box2Abb[2428]=5. + Box2Abb[2427]*m_z2k;

  Box2Abb[2429]=4. - 9.*m_z2k + 6.*m_z2k_2;

  Box2Abb[2430]=-5. + Box2Abb[2429]*m_z2k;

  Box2Abb[2431]=-2. + Box2Abb[2430]*m_z2k;

  Box2Abb[2432]=13. + m_z2k;

  Box2Abb[2433]=-9. + Box2Abb[2432]*m_z2k;

  Box2Abb[2434]=-9. + 2.*Box2Abb[2433]*m_z2k;

  Box2Abb[2435]=3. + Box2Abb[2434]*m_z2k;

  Box2Abb[2436]=-35. + 6.*Box2Abb[488]*m_z2k;

  Box2Abb[2437]=-45. + Box2Abb[2436]*m_z2k;

  Box2Abb[2438]=21. + Box2Abb[2437]*m_z2k;

  Box2Abb[2439]=-1. + Box2Abb[2438]*m_z2k;

  Box2Abb[2440]=Box2Abb[2431]*pow(Box2Abb[61],3.) + Box2Abb[2435]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[2439]*Box2Abb[61]*m_z12_2 + Box2Abb[2428]*m_z12_3*m_z2k - 2.*Box2Abb[1941]*m_z12_4*m_z2k_3;

  Box2Abb[2441]=-10. + 77.*m_z2k;

  Box2Abb[2442]=3. + Box2Abb[1151]*m_z2k;

  Box2Abb[2443]=89. + 170.*m_z2k;

  Box2Abb[2444]=29. + Box2Abb[2443]*m_z2k;

  Box2Abb[2445]=96. + Box2Abb[2444]*m_z2k;

  Box2Abb[2446]=-20. + Box2Abb[2445]*m_z2k;

  Box2Abb[2447]=17. - 7.*m_z2k;

  Box2Abb[2448]=22. + 5.*Box2Abb[2447]*m_z2k;

  Box2Abb[2449]=53. + Box2Abb[2448]*m_z2k;

  Box2Abb[2450]=53. + Box2Abb[2449]*m_z2k;

  Box2Abb[2451]=-22. + Box2Abb[2450]*m_z2k;

  Box2Abb[2452]=-15. + 14.*m_z2k;

  Box2Abb[2453]=98. + 9.*Box2Abb[2452]*m_z2k;

  Box2Abb[2454]=4. + Box2Abb[2453]*m_z2k;

  Box2Abb[2455]=8. + Box2Abb[2454]*m_z2k;

  Box2Abb[2456]=5. + Box2Abb[2455]*m_z2k;

  Box2Abb[2457]=20. + Box2Abb[2456]*m_z12 + 2.*Box2Abb[2451]*m_z12_2 - Box2Abb[2446]*m_z12_3 + 2.*Box2Abb[179]*Box2Abb[2442]*m_z2k - Box2Abb[12]*Box2Abb[2441]*m_z12_4*m_z2k;

  Box2Abb[2458]=Box2Abb[2392]*m_x_3 + Box2Abb[2457]*m_x_4 + Box2Abb[2426]*m_x_5 + Box2Abb[2415]*m_x_6 + Box2Abb[2400]*m_x_7 + Box2Abb[2395]*m_x_8 + Box2Abb[2440]*m_x_2*m_z12 + Box2Abb[710]*m_x_9*m_z12 - Box2Abb[2406]*Box2Abb[61]*m_x*m_z12_2*m_z2k + Box2Abb[4]*pow(Box2Abb[61],5.)*m_z12_3*m_z2k_2;

  Box2Abb[2459]=-36. + 11.*m_z12 - 84.*m_z2k;

  Box2Abb[2460]=68. + Box2Abb[2459]*m_z12 + 116.*m_z2k;

  Box2Abb[2461]=-78. + 5.*Box2Abb[703]*m_z12;

  Box2Abb[2462]=96. + Box2Abb[2461]*m_z12;

  Box2Abb[2463]=18. + Box2Abb[72]*m_z12;

  Box2Abb[2464]=-14. + Box2Abb[2463]*m_z12;

  Box2Abb[2465]=119. - 57.*m_z12;

  Box2Abb[2466]=162. + Box2Abb[2465]*m_z12;

  Box2Abb[2467]=-172. + Box2Abb[2466]*m_z12;

  Box2Abb[2468]=10. + Box2Abb[2467]*m_z12;

  Box2Abb[2469]=130. - 183.*m_z12;

  Box2Abb[2470]=48. + Box2Abb[2469]*m_z12;

  Box2Abb[2471]=232. + Box2Abb[2470]*m_z12;

  Box2Abb[2472]=12. + Box2Abb[2471]*m_z12;

  Box2Abb[2473]=-78. + 53.*m_z12;

  Box2Abb[2474]=58. + Box2Abb[2473]*m_z12;

  Box2Abb[2475]=4. + Box2Abb[2474]*m_z12;

  Box2Abb[2476]=17. - 9.*m_z12;

  Box2Abb[2477]=-34. + Box2Abb[2462]*m_z12 + 52.*m_z2k + 10.*Box2Abb[2464]*m_z12*m_z2k + Box2Abb[2468]*m_z2k_2 + Box2Abb[2472]*m_z2k_3 - 10.*Box2Abb[2475]*m_z2k_4 + 28.*Box2Abb[2476]*m_z12*m_z2k_5;

  Box2Abb[2478]=-2.*Box2Abb[230]*pow(Box2Abb[61],3.) - 3.*pow(Box2Abb[61],2.)*Box2Abb[831]*m_z12 + 2.*Box2Abb[1988]*Box2Abb[61]*m_z12_2 + Box2Abb[12]*Box2Abb[444]*m_z12_3;

  Box2Abb[2479]=25. + m_z12;

  Box2Abb[2480]=69. - 2.*Box2Abb[2479]*m_z12;

  Box2Abb[2481]=16. + m_z12;

  Box2Abb[2482]=-13. + 9.*m_z12;

  Box2Abb[2483]=50. + 137.*m_z2k;

  Box2Abb[2484]=-2.*Box2Abb[2483] + Box2Abb[2480]*m_z12 + 6.*Box2Abb[2481]*m_z12*m_z2k + 28.*Box2Abb[2482]*m_z2k_2;

  Box2Abb[2485]=Box2Abb[2484]*m_z12 + 4.*m_z2k;

  Box2Abb[2486]=10. + 11.*m_z2k;

  Box2Abb[2487]=1. + 32.*m_z2k;

  Box2Abb[2488]=85. - 6.*Box2Abb[2487]*m_z2k;

  Box2Abb[2489]=6. - 35.*m_z2k;

  Box2Abb[2490]=5. + Box2Abb[2489]*m_z2k;

  Box2Abb[2491]=-129. + 12.*Box2Abb[2490]*m_z2k;

  Box2Abb[2492]=81. + 161.*m_z2k;

  Box2Abb[2493]=51. + Box2Abb[2492]*m_z2k;

  Box2Abb[2494]=53. + 2.*Box2Abb[2493]*m_z2k;

  Box2Abb[2495]=2.*Box2Abb[2494]*m_z12 + Box2Abb[2491]*m_z12_2 + Box2Abb[2488]*m_z12_3 + Box2Abb[2486]*m_z12_4 - 2.*Box2Abb[282]*m_z2k;

  Box2Abb[2496]=-3. + Box2Abb[844]*m_z2k;

  Box2Abb[2497]=1. + 2.*Box2Abb[2496]*m_z2k;

  Box2Abb[2498]=-9. + m_z2k;

  Box2Abb[2499]=3. + Box2Abb[2498]*m_z2k;

  Box2Abb[2500]=-1. + 4.*Box2Abb[2499]*m_z2k;

  Box2Abb[2501]=1. + Box2Abb[2500]*m_z2k;

  Box2Abb[2502]=25. - 78.*m_z2k + 56.*m_z2k_2;

  Box2Abb[2503]=-16. + Box2Abb[2502]*m_z2k;

  Box2Abb[2504]=5. + Box2Abb[2503]*m_z2k;

  Box2Abb[2505]=-15. + 67.*m_z2k;

  Box2Abb[2506]=17. + Box2Abb[2505]*m_z2k;

  Box2Abb[2507]=-3. + Box2Abb[2506]*m_z2k;

  Box2Abb[2508]=2. + Box2Abb[2507]*m_z2k;

  Box2Abb[2509]=-2.*Box2Abb[230]*pow(Box2Abb[61],5.) + 2.*Box2Abb[2497]*pow(Box2Abb[61],4.)*m_z12 - 3.*Box2Abb[2501]*pow(Box2Abb[61],3.)*m_z12_2 - Box2Abb[2504]*pow(Box2Abb[61],2.)*m_z12_3 - Box2Abb[2508]*Box2Abb[61]*m_z12_4 - 24.*m_z12_5*m_z2k_4;

  Box2Abb[2510]=15. + 4.*Box2Abb[1151]*m_z2k;

  Box2Abb[2511]=-25. + 46.*m_z2k;

  Box2Abb[2512]=-20. + Box2Abb[2511]*m_z2k;

  Box2Abb[2513]=67. + 243.*m_z2k;

  Box2Abb[2514]=11. + Box2Abb[2513]*m_z2k;

  Box2Abb[2515]=-30. + Box2Abb[2514]*m_z2k;

  Box2Abb[2516]=-9. + 7.*m_z2k;

  Box2Abb[2517]=-129. + 20.*Box2Abb[2516]*m_z2k;

  Box2Abb[2518]=-53. + Box2Abb[2517]*m_z2k;

  Box2Abb[2519]=50. + Box2Abb[2518]*m_z2k;

  Box2Abb[2520]=-11. + 70.*m_z2k;

  Box2Abb[2521]=64. + 5.*Box2Abb[2520]*m_z2k;

  Box2Abb[2522]=13. + Box2Abb[2521]*m_z2k;

  Box2Abb[2523]=64. + Box2Abb[2522]*m_z2k;

  Box2Abb[2524]=22. - 2.*Box2Abb[2523]*m_z12 + 3.*Box2Abb[2519]*m_z12_2 + 2.*Box2Abb[2515]*m_z12_3 + Box2Abb[2512]*m_z12_4 + 2.*Box2Abb[2510]*m_z2k;

  Box2Abb[2525]=-7. + 2.*Box2Abb[559]*m_z2k;

  Box2Abb[2526]=-6. + 97.*m_z2k;

  Box2Abb[2527]=-5. + Box2Abb[2526]*m_z2k;

  Box2Abb[2528]=5. + Box2Abb[2527]*m_z2k;

  Box2Abb[2529]=-5. + Box2Abb[2528]*m_z2k;

  Box2Abb[2530]=-65. + 98.*m_z2k;

  Box2Abb[2531]=6. + Box2Abb[2530]*m_z2k;

  Box2Abb[2532]=-37. + Box2Abb[2531]*m_z2k;

  Box2Abb[2533]=14. + Box2Abb[2532]*m_z2k;

  Box2Abb[2534]=-5. + m_z2k;

  Box2Abb[2535]=41. + 28.*Box2Abb[2534]*m_z2k;

  Box2Abb[2536]=-21. + Box2Abb[2535]*m_z2k;

  Box2Abb[2537]=23. + Box2Abb[2536]*m_z2k;

  Box2Abb[2538]=-3. + Box2Abb[2537]*m_z2k;

  Box2Abb[2539]=-195. + 137.*m_z2k;

  Box2Abb[2540]=-6. + Box2Abb[2539]*m_z2k;

  Box2Abb[2541]=-16. + Box2Abb[2540]*m_z2k;

  Box2Abb[2542]=-15. + Box2Abb[2541]*m_z2k;

  Box2Abb[2543]=7. + Box2Abb[2542]*m_z2k;

  Box2Abb[2544]=2.*Box2Abb[2525]*pow(Box2Abb[61],3.) - 2.*Box2Abb[2533]*pow(Box2Abb[61],2.)*m_z12 + 3.*Box2Abb[2538]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[2543]*m_z12_3 + 2.*Box2Abb[2529]*m_z12_4 + 24.*m_z12_5*m_z2k_3;

  Box2Abb[2545]=Box2Abb[2509]*m_x_2 + Box2Abb[2544]*m_x_3 + Box2Abb[2477]*m_x_4 + Box2Abb[2524]*m_x_5 + Box2Abb[2495]*m_x_6 + Box2Abb[2485]*m_x_7 + Box2Abb[2460]*m_x_8*m_z12 + 4.*Box2Abb[166]*m_x_9*m_z12 + Box2Abb[2478]*pow(Box2Abb[61],3.)*m_x*m_z12*m_z2k - Box2Abb[5]*pow(Box2Abb[61],5.)*m_z12_3*m_z2k_2;

  Box2Abb[2546]=119. + 22.*m_z12;

  Box2Abb[2547]=-81. + Box2Abb[2546]*m_z12;

  Box2Abb[2548]=-88. + Box2Abb[2547]*m_z12;

  Box2Abb[2549]=232. + 31.*m_z12;

  Box2Abb[2550]=81. + Box2Abb[2549]*m_z12;

  Box2Abb[2551]=-208. + Box2Abb[2550]*m_z12;

  Box2Abb[2552]=-5. + Box2Abb[735]*m_z12;

  Box2Abb[2553]=72. + 13.*Box2Abb[2552]*m_z12;

  Box2Abb[2554]=80. + Box2Abb[2548]*m_z12 + 138.*m_z2k + Box2Abb[2551]*m_z12*m_z2k + 2.*Box2Abb[2553]*m_z2k_2 - 336.*Box2Abb[4]*m_z12*m_z2k_3;

  Box2Abb[2555]=7. + 5.*m_z12;

  Box2Abb[2556]=-76. + 5.*Box2Abb[2555]*m_z12;

  Box2Abb[2557]=5. + Box2Abb[2556]*m_z12;

  Box2Abb[2558]=99. + 2.*m_z12;

  Box2Abb[2559]=299. + Box2Abb[2558]*m_z12;

  Box2Abb[2560]=-238. + Box2Abb[2559]*m_z12;

  Box2Abb[2561]=-98. + Box2Abb[2560]*m_z12;

  Box2Abb[2562]=424. + 51.*m_z12;

  Box2Abb[2563]=93. + Box2Abb[2562]*m_z12;

  Box2Abb[2564]=-212. + Box2Abb[2563]*m_z12;

  Box2Abb[2565]=130. + Box2Abb[2564]*m_z12;

  Box2Abb[2566]=418. - 59.*m_z12;

  Box2Abb[2567]=-257. + Box2Abb[2566]*m_z12;

  Box2Abb[2568]=110. + Box2Abb[2567]*m_z12;

  Box2Abb[2569]=50. + 2.*Box2Abb[2557]*m_z12 + 110.*m_z2k + Box2Abb[2561]*m_z12*m_z2k + Box2Abb[2565]*m_z2k_2 + 2.*Box2Abb[2568]*m_z2k_3 - 336.*Box2Abb[4]*m_z12*m_z2k_4;

  Box2Abb[2570]=-7. + 12.*m_z12;

  Box2Abb[2571]=-82. + 5.*Box2Abb[2570]*m_z12;

  Box2Abb[2572]=76. + Box2Abb[2571]*m_z12;

  Box2Abb[2573]=63. + 5.*Box2Abb[533]*m_z12;

  Box2Abb[2574]=-146. + Box2Abb[2573]*m_z12;

  Box2Abb[2575]=44. + Box2Abb[2574]*m_z12;

  Box2Abb[2576]=84. + 5.*m_z12;

  Box2Abb[2577]=335. + 2.*Box2Abb[2576]*m_z12;

  Box2Abb[2578]=-323. + Box2Abb[2577]*m_z12;

  Box2Abb[2579]=30. + Box2Abb[2578]*m_z12;

  Box2Abb[2580]=20. + Box2Abb[2579]*m_z12;

  Box2Abb[2581]=696. + m_z12;

  Box2Abb[2582]=-515. + Box2Abb[2581]*m_z12;

  Box2Abb[2583]=280. + Box2Abb[2582]*m_z12;

  Box2Abb[2584]=-60. + Box2Abb[2583]*m_z12;

  Box2Abb[2585]=65. - 11.*m_z12;

  Box2Abb[2586]=-345. + 8.*Box2Abb[2585]*m_z12;

  Box2Abb[2587]=100. + Box2Abb[2586]*m_z12;

  Box2Abb[2588]=-10. + Box2Abb[2572]*m_z12 + 46.*m_z2k + 2.*Box2Abb[2575]*m_z12*m_z2k + Box2Abb[2580]*m_z2k_2 + Box2Abb[2584]*m_z2k_3 + 2.*Box2Abb[2587]*m_z2k_4 - 168.*Box2Abb[4]*m_z12*m_z2k_5;

  Box2Abb[2589]=-8. + 15.*m_z12 - 60.*m_z2k;

  Box2Abb[2590]=16. + Box2Abb[2589]*m_z12 + 60.*m_z2k;

  Box2Abb[2591]=8. + Box2Abb[2590]*m_z12;

  Box2Abb[2592]=5. + 32.*m_z2k;

  Box2Abb[2593]=70. + 53.*m_z2k;

  Box2Abb[2594]=7. - 48.*m_z2k;

  Box2Abb[2595]=-7. + 4.*Box2Abb[2594]*m_z2k;

  Box2Abb[2596]=-44. + Box2Abb[2595]*m_z12 + Box2Abb[2593]*m_z12_2 + 4.*m_z12_3 + 6.*Box2Abb[2592]*m_z2k;

  Box2Abb[2597]=40. + Box2Abb[2596]*m_z12 + 52.*m_z2k;

  Box2Abb[2598]=3. + 4.*Box2Abb[61]*m_z2k;

  Box2Abb[2599]=2. + Box2Abb[2598]*m_z2k;

  Box2Abb[2600]=-3. + 5.*Box2Abb[12]*m_z2k;

  Box2Abb[2601]=1. + Box2Abb[2600]*m_z2k;

  Box2Abb[2602]=12. + 17.*m_z2k;

  Box2Abb[2603]=-9. + Box2Abb[2602]*m_z2k;

  Box2Abb[2604]=6. + Box2Abb[2603]*m_z2k;

  Box2Abb[2605]=21. + 23.*m_z2k;

  Box2Abb[2606]=-16. + Box2Abb[2605]*m_z2k;

  Box2Abb[2607]=6. + Box2Abb[2606]*m_z2k;

  Box2Abb[2608]=Box2Abb[2599]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[2604]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[2607]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[2601]*m_z12_4 - 2.*Box2Abb[230]*pow(Box2Abb[61],4.)*m_z2k;

  Box2Abb[2609]=1. + Box2Abb[797]*m_z2k;

  Box2Abb[2610]=1. + 3.*Box2Abb[12]*Box2Abb[2609]*m_z2k;

  Box2Abb[2611]=14. + Box2Abb[2396]*m_z2k;

  Box2Abb[2612]=-10. + Box2Abb[2611]*m_z2k;

  Box2Abb[2613]=5. + Box2Abb[2612]*m_z2k;

  Box2Abb[2614]=-38. + 66.*m_z2k + 69.*m_z2k_2;

  Box2Abb[2615]=-30. + Box2Abb[2614]*m_z2k;

  Box2Abb[2616]=27. + Box2Abb[2615]*m_z2k;

  Box2Abb[2617]=-2. + Box2Abb[2616]*m_z2k;

  Box2Abb[2618]=-165. + 84.*m_z2k + 86.*m_z2k_2;

  Box2Abb[2619]=16. + Box2Abb[2618]*m_z2k;

  Box2Abb[2620]=25. + Box2Abb[2619]*m_z2k;

  Box2Abb[2621]=-6. + Box2Abb[2620]*m_z2k;

  Box2Abb[2622]=-109. + 4.*Box2Abb[59]*m_z2k;

  Box2Abb[2623]=31. + Box2Abb[2622]*m_z2k;

  Box2Abb[2624]=4. + Box2Abb[2623]*m_z2k;

  Box2Abb[2625]=-6. + Box2Abb[2624]*m_z2k;

  Box2Abb[2626]=-2.*Box2Abb[2610]*pow(Box2Abb[61],4.)*m_z12 + Box2Abb[2625]*pow(Box2Abb[61],3.)*m_z12_2 + Box2Abb[2621]*pow(Box2Abb[61],2.)*m_z12_3 + Box2Abb[2617]*Box2Abb[61]*m_z12_4 + 2.*Box2Abb[230]*pow(Box2Abb[61],5.)*m_z2k + 2.*Box2Abb[2613]*m_z12_5*m_z2k;

  Box2Abb[2627]=5. + Box2Abb[2099]*m_z2k;

  Box2Abb[2628]=235. - 27.*m_z2k;

  Box2Abb[2629]=160. + Box2Abb[2628]*m_z2k;

  Box2Abb[2630]=10. + Box2Abb[2629]*m_z2k;

  Box2Abb[2631]=40. + Box2Abb[2630]*m_z2k;

  Box2Abb[2632]=2. + 9.*Box2Abb[514]*m_z2k;

  Box2Abb[2633]=-7. + 2.*Box2Abb[2632]*m_z2k;

  Box2Abb[2634]=7. + Box2Abb[2633]*m_z2k;

  Box2Abb[2635]=310. - 187.*m_z2k;

  Box2Abb[2636]=-106. + Box2Abb[2635]*m_z2k;

  Box2Abb[2637]=12. + Box2Abb[2636]*m_z2k;

  Box2Abb[2638]=9. + Box2Abb[2637]*m_z2k;

  Box2Abb[2639]=28. + 2.*Box2Abb[2638]*m_z2k;

  Box2Abb[2640]=-37. + m_z2k;

  Box2Abb[2641]=133. + 20.*Box2Abb[2640]*m_z2k;

  Box2Abb[2642]=136. + Box2Abb[2641]*m_z2k;

  Box2Abb[2643]=-43. + Box2Abb[2642]*m_z2k;

  Box2Abb[2644]=70. + Box2Abb[2643]*m_z2k;

  Box2Abb[2645]=-125. + 68.*m_z2k;

  Box2Abb[2646]=92. + 3.*Box2Abb[2645]*m_z2k;

  Box2Abb[2647]=-4. + 3.*Box2Abb[2646]*m_z2k;

  Box2Abb[2648]=-88. + Box2Abb[2647]*m_z2k;

  Box2Abb[2649]=17. + Box2Abb[2648]*m_z2k;

  Box2Abb[2650]=2.*Box2Abb[2634]*Box2Abb[61] + Box2Abb[2639]*m_z12 + Box2Abb[2649]*m_z12_2 - Box2Abb[2644]*m_z12_3 + Box2Abb[2631]*m_z12_4 + 4.*Box2Abb[2627]*m_z12_5*m_z2k;

  Box2Abb[2651]=-8. + m_z2k;

  Box2Abb[2652]=-5. + Box2Abb[2651]*m_z2k_2;

  Box2Abb[2653]=1. + m_z2k - 7.*m_z2k_2 + 16.*m_z2k_3;

  Box2Abb[2654]=43. + 24.*m_z2k;

  Box2Abb[2655]=-66. + Box2Abb[2654]*m_z2k;

  Box2Abb[2656]=3. + Box2Abb[2655]*m_z2k;

  Box2Abb[2657]=2. + Box2Abb[2656]*m_z2k;

  Box2Abb[2658]=2. + Box2Abb[2657]*m_z2k;

  Box2Abb[2659]=194. + 41.*m_z2k;

  Box2Abb[2660]=6. + Box2Abb[2659]*m_z2k;

  Box2Abb[2661]=66. + Box2Abb[2660]*m_z2k;

  Box2Abb[2662]=-45. + Box2Abb[2661]*m_z2k;

  Box2Abb[2663]=14. + Box2Abb[2662]*m_z2k;

  Box2Abb[2664]=17. + 6.*m_z2k;

  Box2Abb[2665]=-581. + 8.*Box2Abb[2664]*m_z2k;

  Box2Abb[2666]=300. + Box2Abb[2665]*m_z2k;

  Box2Abb[2667]=-64. + Box2Abb[2666]*m_z2k;

  Box2Abb[2668]=20. + Box2Abb[2667]*m_z2k;

  Box2Abb[2669]=-27. + Box2Abb[2668]*m_z2k;

  Box2Abb[2670]=140. + 59.*m_z2k;

  Box2Abb[2671]=-595. + 2.*Box2Abb[2670]*m_z2k;

  Box2Abb[2672]=212. + Box2Abb[2671]*m_z2k;

  Box2Abb[2673]=-128. + Box2Abb[2672]*m_z2k;

  Box2Abb[2674]=68. + Box2Abb[2673]*m_z2k;

  Box2Abb[2675]=-35. + Box2Abb[2674]*m_z2k;

  Box2Abb[2676]=2.*Box2Abb[2653]*pow(Box2Abb[61],3.) - 2.*Box2Abb[2658]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2669]*Box2Abb[61]*m_z12_2 + Box2Abb[2675]*m_z12_3 + Box2Abb[2663]*m_z12_4 - 4.*Box2Abb[2652]*m_z12_5*m_z2k;

  Box2Abb[2677]=-Box2Abb[2626]*m_x_2 + Box2Abb[2676]*m_x_3 - Box2Abb[2650]*m_x_4 + Box2Abb[2588]*m_x_5 - Box2Abb[2569]*m_x_6 + Box2Abb[2554]*m_x_7 - Box2Abb[2597]*m_x_8 + Box2Abb[2591]*m_x_9 + 8.*Box2Abb[4]*m_x_10*m_z12 + Box2Abb[2608]*pow(Box2Abb[61],2.)*m_x*m_z12*m_z2k + Box2Abb[5]*pow(Box2Abb[61],5.)*m_z12_3*m_z2k_3;

  Box2Abb[2678]=26. + m_z12;

  Box2Abb[2679]=6. + Box2Abb[2678]*m_z12;

  Box2Abb[2680]=-56. + Box2Abb[2679]*m_z12;

  Box2Abb[2681]=29. + Box2Abb[2680]*m_z12;

  Box2Abb[2682]=91. + 34.*m_z12;

  Box2Abb[2683]=-39. + Box2Abb[2682]*m_z12;

  Box2Abb[2684]=-26. + Box2Abb[2683]*m_z12;

  Box2Abb[2685]=240. + 17.*m_z12;

  Box2Abb[2686]=-9. + Box2Abb[2685]*m_z12;

  Box2Abb[2687]=30. + Box2Abb[2686]*m_z12;

  Box2Abb[2688]=10. + Box2Abb[2681]*m_z12 + 20.*m_z2k + Box2Abb[2684]*m_z12*m_z2k + Box2Abb[2687]*m_z2k_2 + 84.*Box2Abb[505]*m_z12*m_z2k_3;

  Box2Abb[2689]=-14. + Box2Abb[734]*m_z12;

  Box2Abb[2690]=3. - 5.*Box2Abb[2689]*m_z12;

  Box2Abb[2691]=-55. + Box2Abb[2690]*m_z12;

  Box2Abb[2692]=12. - 5.*Box2Abb[735]*m_z12;

  Box2Abb[2693]=100. + Box2Abb[2692]*m_z12;

  Box2Abb[2694]=-87. + Box2Abb[2693]*m_z12;

  Box2Abb[2695]=18. + 5.*m_z12;

  Box2Abb[2696]=-27. + 2.*Box2Abb[2695]*m_z12;

  Box2Abb[2697]=7. + Box2Abb[2696]*m_z12;

  Box2Abb[2698]=-614. + 115.*m_z12;

  Box2Abb[2699]=285. + Box2Abb[2698]*m_z12;

  Box2Abb[2700]=-40. + Box2Abb[2699]*m_z12;

  Box2Abb[2701]=20. + Box2Abb[2691]*m_z12 + Box2Abb[2694]*m_z12*m_z2k - 5.*Box2Abb[2697]*m_z12*m_z2k_2 + Box2Abb[2700]*m_z2k_3 + 42.*Box2Abb[1635]*m_z12*m_z2k_4;

  Box2Abb[2702]=3. + m_z12 - m_z12_2;

  Box2Abb[2703]=-11. + 2.*Box2Abb[2702]*m_z12;

  Box2Abb[2704]=29. + 5.*Box2Abb[2703]*m_z12;

  Box2Abb[2705]=-17. + 4.*m_z12;

  Box2Abb[2706]=102. + 5.*Box2Abb[2705]*m_z12;

  Box2Abb[2707]=-53. + Box2Abb[2706]*m_z12;

  Box2Abb[2708]=2. + Box2Abb[62]*Box2Abb[735]*m_z12;

  Box2Abb[2709]=-73. + 26.*m_z12;

  Box2Abb[2710]=116. + Box2Abb[2709]*m_z12;

  Box2Abb[2711]=-91. + Box2Abb[2710]*m_z12;

  Box2Abb[2712]=10. + Box2Abb[2711]*m_z12;

  Box2Abb[2713]=825. - 510.*m_z12 + 86.*m_z12_2;

  Box2Abb[2714]=-595. + Box2Abb[2713]*m_z12;

  Box2Abb[2715]=40. + Box2Abb[2714]*m_z12;

  Box2Abb[2716]=-558. + 185.*m_z12;

  Box2Abb[2717]=495. + Box2Abb[2716]*m_z12;

  Box2Abb[2718]=-12. + Box2Abb[2717]*m_z12;

  Box2Abb[2719]=-4. + Box2Abb[2704]*m_z12 + 12.*m_z2k + Box2Abb[2707]*m_z12*m_z2k + 2.*Box2Abb[2708]*m_z2k_2 - 4.*Box2Abb[2712]*m_z2k_3 + Box2Abb[2715]*m_z2k_4 + Box2Abb[2718]*m_z2k_5 + 84.*Box2Abb[72]*m_z12*m_z2k_6;

  Box2Abb[2720]=-9. + Box2Abb[70]*m_z12;

  Box2Abb[2721]=39. + 5.*Box2Abb[2720]*m_z12;

  Box2Abb[2722]=-7. + 2.*Box2Abb[2721]*m_z12;

  Box2Abb[2723]=23. - 14.*m_z12 + 10.*m_z12_3;

  Box2Abb[2724]=-12. + Box2Abb[2723]*m_z12;

  Box2Abb[2725]=24. + 5.*m_z12;

  Box2Abb[2726]=-27. + Box2Abb[2725]*m_z12;

  Box2Abb[2727]=14. + Box2Abb[2726]*m_z12;

  Box2Abb[2728]=-12. + Box2Abb[2727]*m_z12;

  Box2Abb[2729]=55. - 2.*m_z12;

  Box2Abb[2730]=-515. + 7.*Box2Abb[2729]*m_z12;

  Box2Abb[2731]=340. + Box2Abb[2730]*m_z12;

  Box2Abb[2732]=-40. + Box2Abb[2731]*m_z12;

  Box2Abb[2733]=160. - 47.*m_z12;

  Box2Abb[2734]=-111. + Box2Abb[2733]*m_z12;

  Box2Abb[2735]=6. + Box2Abb[2734]*m_z12;

  Box2Abb[2736]=8. - 5.*m_z12;

  Box2Abb[2737]=2.*Box2Abb[1877] + Box2Abb[2722]*m_z12 + Box2Abb[2724]*m_z12*m_z2k + 2.*Box2Abb[2728]*m_z12*m_z2k_2 + Box2Abb[2732]*m_z2k_3 + 5.*Box2Abb[2735]*m_z2k_4 + 42.*Box2Abb[2736]*m_z12*m_z2k_5;

  Box2Abb[2738]=-2. + 10.*m_z12 - 39.*m_z2k;

  Box2Abb[2739]=15. + Box2Abb[2738]*m_z12 + 48.*m_z2k;

  Box2Abb[2740]=2. + Box2Abb[2739]*m_z12;

  Box2Abb[2741]=5. + m_z12;

  Box2Abb[2742]=-19. + 6.*Box2Abb[2741]*m_z12;

  Box2Abb[2743]=34. + 37.*m_z12;

  Box2Abb[2744]=14. - 11.*m_z12;

  Box2Abb[2745]=-3. + Box2Abb[2742]*m_z12 + 51.*m_z2k + Box2Abb[2743]*m_z12*m_z2k + 12.*Box2Abb[2744]*m_z2k_2;

  Box2Abb[2746]=8. + Box2Abb[2745]*m_z12 + 12.*m_z2k;

  Box2Abb[2747]=1. - 3.*m_z2k + 6.*m_z2k_2;

  Box2Abb[2748]=-15. + m_z2k;

  Box2Abb[2749]=-8. + Box2Abb[2748]*Box2Abb[61]*m_z2k;

  Box2Abb[2750]=2. + Box2Abb[2749]*m_z2k;

  Box2Abb[2751]=6. + 5.*Box2Abb[179]*m_z2k;

  Box2Abb[2752]=-4. + Box2Abb[2751]*m_z2k;

  Box2Abb[2753]=1. + Box2Abb[2752]*m_z2k;

  Box2Abb[2754]=22. + 3.*m_z2k;

  Box2Abb[2755]=-39. + Box2Abb[2754]*m_z2k;

  Box2Abb[2756]=16. + Box2Abb[2755]*m_z2k;

  Box2Abb[2757]=-4. + Box2Abb[2756]*m_z2k;

  Box2Abb[2758]=37. + 7.*m_z2k;

  Box2Abb[2759]=-54. + Box2Abb[2758]*m_z2k;

  Box2Abb[2760]=24. + Box2Abb[2759]*m_z2k;

  Box2Abb[2761]=-6. + Box2Abb[2760]*m_z2k;

  Box2Abb[2762]=-Box2Abb[2747]*pow(Box2Abb[61],5.) - Box2Abb[2757]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[2761]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[2750]*Box2Abb[61]*m_z12_3 + Box2Abb[2753]*m_z12_4;

  Box2Abb[2763]=5. - 10.*m_z2k + 16.*m_z2k_3 - 13.*m_z2k_4;

  Box2Abb[2764]=53. - 75.*m_z2k + 48.*m_z2k_2;

  Box2Abb[2765]=-21. + Box2Abb[2764]*m_z2k;

  Box2Abb[2766]=7. + Box2Abb[2765]*m_z2k;

  Box2Abb[2767]=-46. + 29.*m_z2k;

  Box2Abb[2768]=2. + Box2Abb[2767]*Box2Abb[61]*m_z2k;

  Box2Abb[2769]=-15. + Box2Abb[2768]*m_z2k;

  Box2Abb[2770]=7. + Box2Abb[2769]*m_z2k;

  Box2Abb[2771]=-232. + 53.*m_z2k;

  Box2Abb[2772]=218. + Box2Abb[2771]*m_z2k;

  Box2Abb[2773]=-62. + Box2Abb[2772]*m_z2k;

  Box2Abb[2774]=10. + m_z2k + Box2Abb[2773]*m_z2k_2;

  Box2Abb[2775]=-40. + 3.*m_z2k;

  Box2Abb[2776]=217. + 4.*Box2Abb[2775]*m_z2k;

  Box2Abb[2777]=-110. + Box2Abb[2776]*m_z2k;

  Box2Abb[2778]=33. + Box2Abb[2777]*m_z2k;

  Box2Abb[2779]=-4. + Box2Abb[2778]*m_z2k;

  Box2Abb[2780]=2.*pow(Box2Abb[61],6.) + Box2Abb[2766]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[2779]*pow(Box2Abb[61],2.)*m_z12_2 - Box2Abb[2774]*Box2Abb[61]*m_z12_3 - 2.*Box2Abb[2770]*m_z12_4 + Box2Abb[2763]*m_z12_5;

  Box2Abb[2781]=Box2Abb[2780]*m_x_2 + Box2Abb[2719]*m_x_3 + Box2Abb[2737]*m_x_4 + Box2Abb[2701]*m_x_5 + Box2Abb[2688]*m_x_6 - Box2Abb[2746]*m_x_7 + Box2Abb[2740]*m_x_8 + Box2Abb[2762]*Box2Abb[61]*m_x*m_z12 + Box2Abb[710]*m_x_9*m_z12 + Box2Abb[5]*Box2Abb[519]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_3;

  Box2Abb[2782]=-3. + m_z12 - 6.*m_z2k;

  Box2Abb[2783]=-2. + Box2Abb[2782]*m_z12;

  Box2Abb[2784]=-1. + Box2Abb[1941]*m_z2k;

  Box2Abb[2785]=2.*Box2Abb[256]*pow(Box2Abb[61],2.) + Box2Abb[2784]*m_z12;

  Box2Abb[2786]=36. - Box2Abb[2182]*m_z12 + 54.*m_z2k + 12.*m_z12*m_z2k + 90.*m_z2k_2;

  Box2Abb[2787]=64. + Box2Abb[2786]*m_z12 + 48.*m_z2k;

  Box2Abb[2788]=-2. + Box2Abb[2787]*m_z12;

  Box2Abb[2789]=-2. + 3.*m_z2k + m_z2k_3;

  Box2Abb[2790]=-5. + 20.*m_z2k - 46.*m_z2k_3 + 31.*m_z2k_4;

  Box2Abb[2791]=-7. + 3.*m_z2k;

  Box2Abb[2792]=12. + Box2Abb[2791]*m_z2k;

  Box2Abb[2793]=1. + Box2Abb[2792]*m_z2k;

  Box2Abb[2794]=6.*Box2Abb[2789]*pow(Box2Abb[61],3.) + 4.*Box2Abb[2793]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2790]*m_z12_2 + 12.*m_z12_3*m_z2k_3;

  Box2Abb[2795]=15. + 2.*m_z2k - 6.*m_z2k_3 + 15.*m_z2k_4;

  Box2Abb[2796]=-5. + 31.*m_z2k;

  Box2Abb[2797]=-10. + Box2Abb[2796]*m_z2k;

  Box2Abb[2798]=29. + 20.*m_z2k + 42.*m_z2k_2;

  Box2Abb[2799]=-10. + Box2Abb[2798]*m_z2k;

  Box2Abb[2800]=-4. + 3.*m_z2k;

  Box2Abb[2801]=-13. + 4.*Box2Abb[2800]*m_z2k;

  Box2Abb[2802]=64. + 4.*Box2Abb[2801]*m_z2k;

  Box2Abb[2803]=-2.*Box2Abb[1359]*Box2Abb[61] + Box2Abb[2802]*m_z12 + 6.*Box2Abb[2795]*m_z12_2 + 4.*Box2Abb[2799]*m_z12_3 + Box2Abb[2797]*m_z12_4;

  Box2Abb[2804]=5. + 4.*m_z2k;

  Box2Abb[2805]=-25. + 51.*m_z2k;

  Box2Abb[2806]=13. + Box2Abb[756]*m_z2k;

  Box2Abb[2807]=6. + Box2Abb[286]*m_z2k;

  Box2Abb[2808]=10. - 8.*Box2Abb[2806]*m_z12 - 12.*Box2Abb[2807]*m_z12_2 - 2.*Box2Abb[12]*Box2Abb[2805]*m_z12_3 + Box2Abb[2804]*m_z12_4 + 4.*m_z2k;

  Box2Abb[2809]=-9. + m_z2k + 3.*m_z2k_2 - 3.*m_z2k_3 + 6.*m_z2k_4;

  Box2Abb[2810]=23. + 34.*m_z2k;

  Box2Abb[2811]=5. + Box2Abb[2810]*m_z2k;

  Box2Abb[2812]=-5. + Box2Abb[2811]*m_z2k;

  Box2Abb[2813]=-40. + 51.*m_z2k;

  Box2Abb[2814]=-54. + Box2Abb[2813]*m_z2k;

  Box2Abb[2815]=60. + Box2Abb[2814]*m_z2k;

  Box2Abb[2816]=-5. + Box2Abb[2815]*m_z2k;

  Box2Abb[2817]=2.*Box2Abb[470]*pow(Box2Abb[61],2.)*Box2Abb[841] + 3.*Box2Abb[2809]*Box2Abb[61]*m_z12 + Box2Abb[2816]*m_z12_2 + Box2Abb[2812]*m_z12_3;

  Box2Abb[2818]=Box2Abb[2803]*m_x_4 + Box2Abb[2808]*m_x_5 + Box2Abb[2788]*m_x_6 - 2.*Box2Abb[2817]*m_x_3*m_z12 + 6.*Box2Abb[2783]*m_x_7*m_z12 + Box2Abb[2794]*m_x_2*m_z12_2 + 6.*m_x_8*m_z12_2 + Box2Abb[2785]*pow(Box2Abb[61],3.)*m_x*m_z12_3 - pow(Box2Abb[61],5.)*m_z12_4*m_z2k;

  Box2Abb[2819]=3. + 7.*m_z12 - 15.*m_z2k;

  Box2Abb[2820]=-6. + Box2Abb[2819]*m_z12;

  Box2Abb[2821]=-26. + 5.*Box2Abb[11]*m_z12;

  Box2Abb[2822]=23. + Box2Abb[2821]*m_z12;

  Box2Abb[2823]=18. - Box2Abb[2822]*m_z12;

  Box2Abb[2824]=32. + 5.*Box2Abb[2376]*m_z12;

  Box2Abb[2825]=-5. + Box2Abb[2824]*m_z12;

  Box2Abb[2826]=9. + 2.*m_z12;

  Box2Abb[2827]=3. + 5.*Box2Abb[2826]*m_z12;

  Box2Abb[2828]=-14. + 5.*m_z12;

  Box2Abb[2829]=86. + 5.*Box2Abb[2828]*m_z12;

  Box2Abb[2830]=-54. + Box2Abb[2829]*m_z12;

  Box2Abb[2831]=3. + Box2Abb[2830]*m_z12;

  Box2Abb[2832]=-402. + 83.*m_z12;

  Box2Abb[2833]=433. + Box2Abb[2832]*m_z12;

  Box2Abb[2834]=-212. + Box2Abb[2833]*m_z12;

  Box2Abb[2835]=8. + Box2Abb[2834]*m_z12;

  Box2Abb[2836]=-88. + 73.*m_z12;

  Box2Abb[2837]=36. + Box2Abb[2836]*m_z12;

  Box2Abb[2838]=-6. + Box2Abb[2823]*m_z12 + 16.*m_z2k + 2.*Box2Abb[2825]*m_z12*m_z2k + 2.*Box2Abb[2827]*Box2Abb[4]*m_z2k_2 - 4.*Box2Abb[2831]*m_z2k_3 + Box2Abb[2835]*m_z2k_4 + 2.*Box2Abb[2837]*m_z12*m_z2k_5;

  Box2Abb[2839]=1. - 2.*m_z2k + 4.*m_z2k_2;

  Box2Abb[2840]=1. + Box2Abb[1071]*m_z2k;

  Box2Abb[2841]=Box2Abb[512]*pow(Box2Abb[61],3.) + 2.*Box2Abb[2839]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[2840]*Box2Abb[61]*m_z12_2 + 2.*m_z12_3*m_z2k_2;

  Box2Abb[2842]=54. + 58.*m_z2k;

  Box2Abb[2843]=-3. + Box2Abb[2842]*m_z12 + 7.*m_z12_2 + 32.*Box2Abb[1403]*m_z2k;

  Box2Abb[2844]=60. - Box2Abb[2843]*m_z12 + 72.*m_z2k;

  Box2Abb[2845]=32. + 39.*m_z2k;

  Box2Abb[2846]=41. + 90.*m_z2k;

  Box2Abb[2847]=82. + 31.*m_z2k;

  Box2Abb[2848]=32. + Box2Abb[2847]*m_z2k;

  Box2Abb[2849]=-61. + 42.*m_z2k;

  Box2Abb[2850]=35. + 4.*Box2Abb[2849]*m_z2k;

  Box2Abb[2851]=60. + Box2Abb[2850]*m_z2k;

  Box2Abb[2852]=2. - Box2Abb[2846]*Box2Abb[444]*m_z12 - Box2Abb[2851]*m_z12_2 + 2.*Box2Abb[2848]*m_z12_3 + Box2Abb[2845]*m_z12_4 + m_z12_5 - 2.*m_z2k;

  Box2Abb[2853]=-119. + 31.*m_z2k;

  Box2Abb[2854]=-33. + Box2Abb[2853]*m_z2k;

  Box2Abb[2855]=76. + 59.*m_z2k;

  Box2Abb[2856]=55. + Box2Abb[2855]*m_z2k;

  Box2Abb[2857]=37. + 60.*m_z2k;

  Box2Abb[2858]=11. + Box2Abb[2857]*m_z2k;

  Box2Abb[2859]=-19. + 7.*m_z2k;

  Box2Abb[2860]=65. + 8.*Box2Abb[2859]*m_z2k;

  Box2Abb[2861]=148. + 3.*Box2Abb[2860]*m_z2k;

  Box2Abb[2862]=53. + Box2Abb[2861]*m_z2k;

  Box2Abb[2863]=18. + Box2Abb[2858]*Box2Abb[831]*m_z12 + Box2Abb[2862]*m_z12_2 - Box2Abb[2856]*m_z12_4 - 5.*Box2Abb[12]*m_z12_5 + 4.*m_z2k + 2.*Box2Abb[2854]*m_z12_3*m_z2k + 8.*m_z2k_2;

  Box2Abb[2864]=26. + m_z2k;

  Box2Abb[2865]=-32. + Box2Abb[2864]*m_z2k;

  Box2Abb[2866]=18. + Box2Abb[2865]*m_z2k;

  Box2Abb[2867]=-3. + Box2Abb[2866]*m_z2k;

  Box2Abb[2868]=-17. + 4.*Box2Abb[514]*m_z2k;

  Box2Abb[2869]=8. + Box2Abb[2868]*m_z2k;

  Box2Abb[2870]=-3. + Box2Abb[2869]*m_z2k;

  Box2Abb[2871]=6. + 11.*m_z2k;

  Box2Abb[2872]=-20. + Box2Abb[2871]*m_z2k;

  Box2Abb[2873]=11. + Box2Abb[2872]*m_z2k;

  Box2Abb[2874]=-2. + Box2Abb[2873]*m_z2k;

  Box2Abb[2875]=pow(Box2Abb[61],4.) - Box2Abb[2870]*pow(Box2Abb[61],3.)*m_z12 - 2.*Box2Abb[2874]*pow(Box2Abb[61],2.)*m_z12_2 - Box2Abb[2867]*Box2Abb[61]*m_z12_3 + Box2Abb[2753]*m_z12_4;

  Box2Abb[2876]=58. - 7.*m_z2k;

  Box2Abb[2877]=15. + Box2Abb[2876]*m_z2k;

  Box2Abb[2878]=40. + Box2Abb[2877]*m_z2k;

  Box2Abb[2879]=-2. + m_z2k + 3.*m_z2k_2;

  Box2Abb[2880]=1. + Box2Abb[2879]*m_z2k;

  Box2Abb[2881]=14. + 3.*Box2Abb[470]*Box2Abb[854]*m_z2k;

  Box2Abb[2882]=5. + 2.*Box2Abb[2881]*m_z2k;

  Box2Abb[2883]=-31. + 19.*m_z2k;

  Box2Abb[2884]=-3. + Box2Abb[2883]*Box2Abb[461]*m_z2k;

  Box2Abb[2885]=25. + Box2Abb[2884]*m_z2k;

  Box2Abb[2886]=-447. + 440.*m_z2k - 84.*m_z2k_2;

  Box2Abb[2887]=4. + Box2Abb[2886]*m_z2k;

  Box2Abb[2888]=17. + Box2Abb[2887]*m_z2k;

  Box2Abb[2889]=14. + Box2Abb[2888]*m_z2k;

  Box2Abb[2890]=-4.*Box2Abb[2880] - 2.*Box2Abb[2882]*m_z12 + Box2Abb[2889]*m_z12_2 - 2.*Box2Abb[2885]*m_z12_3 + Box2Abb[2878]*m_z12_4 + 10.*Box2Abb[890]*m_z12_5;

  Box2Abb[2891]=-5. + 10.*m_z2k - 16.*m_z2k_3 + 13.*m_z2k_4;

  Box2Abb[2892]=-8. + 3.*m_z2k;

  Box2Abb[2893]=15. + 4.*Box2Abb[2892]*m_z2k;

  Box2Abb[2894]=-3. + Box2Abb[2893]*m_z2k;

  Box2Abb[2895]=-35. + 4.*m_z2k + 8.*m_z2k_2;

  Box2Abb[2896]=86. + 3.*Box2Abb[2895]*m_z2k;

  Box2Abb[2897]=-41. + Box2Abb[2896]*m_z2k;

  Box2Abb[2898]=-75. + 11.*m_z2k;

  Box2Abb[2899]=89. + Box2Abb[2898]*m_z2k;

  Box2Abb[2900]=-41. + Box2Abb[2899]*m_z2k;

  Box2Abb[2901]=8. + Box2Abb[2900]*m_z2k;

  Box2Abb[2902]=2. + Box2Abb[2901]*m_z2k;

  Box2Abb[2903]=-28. + 11.*m_z2k;

  Box2Abb[2904]=94. + 5.*Box2Abb[2903]*m_z2k;

  Box2Abb[2905]=4. + Box2Abb[2904]*m_z2k;

  Box2Abb[2906]=-25. + Box2Abb[2905]*m_z2k;

  Box2Abb[2907]=8. + Box2Abb[2906]*m_z2k;

  Box2Abb[2908]=2.*pow(Box2Abb[61],5.) + Box2Abb[2894]*pow(Box2Abb[61],3.)*m_z12 + 2.*Box2Abb[2902]*Box2Abb[61]*m_z12_3 + Box2Abb[2907]*m_z12_4 + Box2Abb[2891]*m_z12_5 - Box2Abb[2897]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k;

  Box2Abb[2909]=-Box2Abb[2908]*m_x_2 + Box2Abb[2838]*m_x_3 + Box2Abb[2890]*m_x_4 + Box2Abb[2863]*m_x_5 + Box2Abb[2852]*m_x_6 + Box2Abb[2875]*Box2Abb[61]*m_x*m_z12 + Box2Abb[2844]*m_x_7*m_z12 + 2.*Box2Abb[2820]*m_x_8*m_z12 + 4.*m_x_9*m_z12_2 + Box2Abb[2841]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k;

  Box2Abb[2981]=(-2.*m_cL*m_z12)/m_s;

  Box2Abb[2982]=(4.*m_cL)/m_s_2;

  Box2Abb[2983]=(-4.*m_cL)/m_s_2;

  Box2Abb[2984]=(4.*m_cL*m_z12)/m_s_2;

  Box2Abb[2985]=(-4.*m_cL*m_x)/m_s_2;

  Box2Abb[2986]=(-4.*m_cL*m_x)/(m_s_2*m_z2k);

  Box2Abb[2987]=-m_cL*m_z12_2;

  Box2Abb[2988]=(4.*m_cL*m_z12)/m_s;

  Box2Abb[2989]=(4.*m_cL*m_x*m_z12)/m_s;

  Box2Abb[2990]=(4.*m_cL*m_z12_2)/m_s;

  Box2Abb[2991]=(-4.*m_cL*m_x*m_z12)/m_s;

  Box2Abb[2992]=-m_cL*m_s*m_z12_2;

  Box2Abb[2993]=4.*m_cL*m_z12;

  Box2Abb[2994]=2.*m_cL;

  Box2Abb[2995]=4.*m_cL*m_x*m_z12;

  Box2Abb[2996]=4.*m_cL*m_z12_2;

  Box2Abb[2997]=-4.*m_cL*m_x*m_z12;

  Box2Abb[2998]=-4.*m_cL;

  Box2Abb[2999]=-2.*m_cL;

  Box2Abb[3000]=-4.*m_cL*m_x;

  Box2Abb[3001]=-4.*m_cL*m_z12;

  Box2Abb[3002]=(-4.*m_cL)/m_s;

  Box2Abb[3003]=(-4.*m_cL*m_x)/m_s;

  Box2Abb[3004]=(8.*m_cL*m_x)/m_s_2;

  Box2Abb[3005]=(4.*m_cL*m_x)/m_s_2;

  Box2Abb[3006]=(-4.*m_cL*m_z2k)/m_s_2;

  Box2Abb[3007]=(4.*m_cL*m_x)/(m_s_2*m_z2k);

  Box2Abb[3008]=(4.*m_cL)/m_s;

  Box2Abb[3009]=(4.*m_cL*m_x)/m_s;

  Box2Abb[3081]=(-2.*m_cR*m_z12)/m_s;

  Box2Abb[3082]=(4.*m_cR)/m_s_2;

  Box2Abb[3083]=(-4.*m_cR)/m_s_2;

  Box2Abb[3084]=(4.*m_cR*m_z12)/m_s_2;

  Box2Abb[3085]=(-4.*m_cR*m_x)/m_s_2;

  Box2Abb[3086]=(-4.*m_cR*m_x)/(m_s_2*m_z2k);

  Box2Abb[3087]=-m_cR*m_z12_2;

  Box2Abb[3088]=(4.*m_cR*m_z12)/m_s;

  Box2Abb[3089]=(4.*m_cR*m_x*m_z12)/m_s;

  Box2Abb[3090]=(4.*m_cR*m_z12_2)/m_s;

  Box2Abb[3091]=(-4.*m_cR*m_x*m_z12)/m_s;

  Box2Abb[3092]=-m_cR*m_s*m_z12_2;

  Box2Abb[3093]=4.*m_cR*m_z12;

  Box2Abb[3094]=2.*m_cR;

  Box2Abb[3095]=4.*m_cR*m_x*m_z12;

  Box2Abb[3096]=4.*m_cR*m_z12_2;

  Box2Abb[3097]=-4.*m_cR*m_x*m_z12;

  Box2Abb[3098]=-4.*m_cR;

  Box2Abb[3099]=-2.*m_cR;

  Box2Abb[3100]=-4.*m_cR*m_x;

  Box2Abb[3101]=-4.*m_cR*m_z12;

  Box2Abb[3102]=(-4.*m_cR)/m_s;

  Box2Abb[3103]=(-4.*m_cR*m_x)/m_s;

  Box2Abb[3104]=(8.*m_cR*m_x)/m_s_2;

  Box2Abb[3105]=(4.*m_cR*m_x)/m_s_2;

  Box2Abb[3106]=(-4.*m_cR*m_z2k)/m_s_2;

  Box2Abb[3107]=(4.*m_cR*m_x)/(m_s_2*m_z2k);

  Box2Abb[3108]=(4.*m_cR)/m_s;

  Box2Abb[3109]=(4.*m_cR*m_x)/m_s;

  Box2Abb[3110]=Box2Abb[768]*m_cL + Box2Abb[2]*m_cR;

  Box2Abb[3111]=Box2Abb[72]*m_cL + 2.*Box2Abb[2]*m_cR;

  Box2Abb[3112]=Box2Abb[3110]*Box2Abb[72]*m_x + Box2Abb[3111]*m_z12*m_z2k + 2.*Box2Abb[50]*m_z12*m_z2k_2;

  Box2Abb[3113]=-3. + 2.*Box2Abb[2741]*m_z12;

  Box2Abb[3114]=10. + 17.*m_z12;

  Box2Abb[3115]=-1. + m_z12 + Box2Abb[3113]*m_z2k + Box2Abb[3114]*m_z2k_2;

  Box2Abb[3116]=2. + 7.*m_z2k;

  Box2Abb[3117]=2. + 9.*m_z2k;

  Box2Abb[3118]=-2. + Box2Abb[3116]*m_z12 + Box2Abb[3117]*m_z2k;

  Box2Abb[3119]=3. + Box2Abb[3116]*m_z12 + 20.*m_z2k;

  Box2Abb[3120]=-2. + Box2Abb[3119]*m_z12;

  Box2Abb[3121]=Box2Abb[3120]*m_x_3 - Box2Abb[3115]*m_x_2*m_z12 + Box2Abb[872]*m_x_4*m_z12 + Box2Abb[3118]*m_x*m_z12_2*m_z2k - Box2Abb[397]*m_z12_3*m_z2k_2;

  Box2Abb[3122]=10. - 11.*m_z12;

  Box2Abb[3123]=2.*pow(Box2Abb[4],2.) + 5.*Box2Abb[4]*m_z2k + m_z2k_2;

  Box2Abb[3124]=-10. + 13.*m_z12;

  Box2Abb[3125]=2. + 2.*Box2Abb[72]*m_z12 + m_z2k + Box2Abb[1453]*m_z12*m_z2k + Box2Abb[3124]*m_z2k_2;

  Box2Abb[3126]=-1. + 4.*m_z12 + 23.*m_z2k;

  Box2Abb[3127]=-5.*Box2Abb[183] + Box2Abb[3126]*m_z12;

  Box2Abb[3128]=2. + Box2Abb[3127]*m_z12;

  Box2Abb[3129]=Box2Abb[3128]*m_x_2 - Box2Abb[3125]*m_x*m_z12 + Box2Abb[3122]*m_x_3*m_z12 + Box2Abb[3123]*m_z12_2*m_z2k;

  Box2Abb[3130]=Box2Abb[3121]*m_cL + Box2Abb[3129]*m_cR*m_x;

  Box2Abb[3131]=-10. + 3.*m_z12;

  Box2Abb[3132]=10. + 11.*m_z12;

  Box2Abb[3133]=-2. + 2.*m_z12 - 9.*m_z2k + 3.*Box2Abb[203]*m_z12*m_z2k + Box2Abb[3132]*m_z2k_2;

  Box2Abb[3134]=4. + 9.*m_z2k;

  Box2Abb[3135]=Box2Abb[61]*Box2Abb[719] + Box2Abb[3134]*m_z12;

  Box2Abb[3136]=9. - 2.*m_z12 + m_z2k;

  Box2Abb[3137]=-3. + Box2Abb[3136]*m_z12 + 20.*m_z2k;

  Box2Abb[3138]=-2. + Box2Abb[3137]*m_z12;

  Box2Abb[3139]=Box2Abb[3138]*m_x_3 - Box2Abb[3133]*m_x_2*m_z12 + Box2Abb[3131]*m_x_4*m_z12 + Box2Abb[3135]*m_x*m_z12_2*m_z2k - 2.*Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[3140]=10. - 13.*m_z12;

  Box2Abb[3141]=-7. + 19.*m_z12;

  Box2Abb[3142]=-10. + 19.*m_z12;

  Box2Abb[3143]=Box2Abb[11]*pow(Box2Abb[4],2.) + Box2Abb[3141]*Box2Abb[4]*m_z2k + Box2Abb[3142]*m_z2k_2;

  Box2Abb[3144]=-14. + 11.*m_z12 + 29.*m_z2k;

  Box2Abb[3145]=1. + Box2Abb[3144]*m_z12 - 20.*m_z2k;

  Box2Abb[3146]=2. + Box2Abb[3145]*m_z12;

  Box2Abb[3147]=Box2Abb[3146]*m_x_2 - Box2Abb[3143]*m_x*m_z12 + Box2Abb[3140]*m_x_3*m_z12 + 3.*pow(Box2Abb[5],2.)*m_z12_2*m_z2k;

  Box2Abb[3148]=Box2Abb[3139]*m_cL + Box2Abb[3147]*m_cR*m_x;

  Box2Abb[3149]=-10. + m_z12 + 4.*m_z12_2;

  Box2Abb[3150]=-34. + 33.*m_z12;

  Box2Abb[3151]=5. + Box2Abb[3149]*m_z12 + 7.*m_z2k + Box2Abb[3150]*m_z12*m_z2k + Box2Abb[3142]*m_z2k_2;

  Box2Abb[3152]=2. + 13.*m_z2k;

  Box2Abb[3153]=-2. + Box2Abb[3152]*m_z12 + Box2Abb[2791]*m_z2k;

  Box2Abb[3154]=-18. + 17.*m_z12 + 29.*m_z2k;

  Box2Abb[3155]=1. + Box2Abb[3154]*m_z12 - 20.*m_z2k;

  Box2Abb[3156]=2. + Box2Abb[3155]*m_z12;

  Box2Abb[3157]=Box2Abb[3156]*m_x_3 - Box2Abb[3151]*m_x_2*m_z12 + Box2Abb[3140]*m_x_4*m_z12 + Box2Abb[3153]*Box2Abb[5]*m_x*m_z12_2 - 2.*pow(Box2Abb[5],2.)*m_z12_3*m_z2k;

  Box2Abb[3158]=Box2Abb[3139]*m_cL + Box2Abb[3157]*m_cR;

  Box2Abb[3159]=-1. + m_z12 + 5.*m_z2k;

  Box2Abb[3160]=pow(Box2Abb[4],2.)*Box2Abb[62] + 3.*Box2Abb[4]*Box2Abb[815]*m_z2k - 5.*Box2Abb[70]*m_z2k_2;

  Box2Abb[3161]=9. - 18.*m_z12 + 7.*m_z12_2 + 5.*Box2Abb[201]*m_z2k;

  Box2Abb[3162]=2. + Box2Abb[3161]*m_z12;

  Box2Abb[3163]=-Box2Abb[3162]*m_x_2 + Box2Abb[3160]*m_x*m_z12 + 5.*Box2Abb[72]*m_x_3*m_z12 + Box2Abb[3159]*Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[3164]=1. + Box2Abb[62]*m_z12 - 8.*m_z2k + 11.*m_z12*m_z2k + 5.*m_z2k_2;

  Box2Abb[3165]=-15. - 13.*Box2Abb[72]*m_z12;

  Box2Abb[3166]=26. - 21.*m_z12;

  Box2Abb[3167]=2. + Box2Abb[3165]*m_z12 - 13.*m_z2k + 2.*Box2Abb[3166]*m_z12*m_z2k + 5.*Box2Abb[879]*m_z2k_2;

  Box2Abb[3168]=-33. + 25.*m_z12 + 35.*m_z2k;

  Box2Abb[3169]=7. + Box2Abb[3168]*m_z12 - 20.*m_z2k;

  Box2Abb[3170]=2. + Box2Abb[3169]*m_z12;

  Box2Abb[3171]=Box2Abb[3170]*m_x_3 + Box2Abb[3167]*m_x_2*m_z12 + 5.*Box2Abb[699]*m_x_4*m_z12 + Box2Abb[3164]*Box2Abb[5]*m_x*m_z12_2 - pow(Box2Abb[5],2.)*m_z12_3*m_z2k;

  Box2Abb[3172]=Box2Abb[3171]*m_cR + Box2Abb[3163]*m_cL*m_x;

  Box2Abb[3173]=2.*Box2Abb[187] + 2.*Box2Abb[538]*m_z12 + Box2Abb[53]*m_z12*m_z2k + 4.*Box2Abb[693]*m_z2k_2;

  Box2Abb[3174]=-3. + m_z12 + 3.*m_z2k;

  Box2Abb[3175]=2. + Box2Abb[3174]*m_z12 - 2.*m_z2k;

  Box2Abb[3176]=-14. + 9.*m_z12 - 2.*m_z12_2 + 12.*m_z2k;

  Box2Abb[3177]=8. + Box2Abb[3176]*m_z12;

  Box2Abb[3178]=-4. + 3.*m_z12 + 12.*m_z2k;

  Box2Abb[3179]=-3. + m_z12 + Box2Abb[3178]*m_z12*m_z2k;

  Box2Abb[3180]=2. + Box2Abb[3179]*m_z12 + 4.*m_z2k;

  Box2Abb[3181]=-Box2Abb[3180]*m_x_2 + Box2Abb[3177]*m_x_3 + 4.*Box2Abb[72]*m_x_4*m_z12 + Box2Abb[3173]*m_x*m_z12*m_z2k - Box2Abb[3175]*m_z12_2*m_z2k_2;

  Box2Abb[3182]=2. + 2.*Box2Abb[537]*m_z12 - 6.*m_z2k + 3.*Box2Abb[274]*m_z12*m_z2k + 4.*Box2Abb[11]*m_z2k_2;

  Box2Abb[3183]=-11. + 12.*m_z12 + 48.*m_z2k;

  Box2Abb[3184]=10. + Box2Abb[3183]*m_z12 - 12.*m_z2k;

  Box2Abb[3185]=-8. + Box2Abb[3184]*m_z12;

  Box2Abb[3186]=7. + 2.*Box2Abb[72]*m_z12 - 2.*m_z2k + 23.*m_z12*m_z2k + 36.*m_z2k_2;

  Box2Abb[3187]=7. - Box2Abb[3186]*m_z12 + 8.*m_z2k;

  Box2Abb[3188]=-2. + Box2Abb[3187]*m_z12 + 4.*m_z2k;

  Box2Abb[3189]=Box2Abb[3188]*m_x_2 + Box2Abb[3185]*m_x_3 + 4.*Box2Abb[879]*m_x_4*m_z12 + Box2Abb[3182]*m_x*m_z12*m_z2k - 3.*Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[3190]=Box2Abb[3181]*m_cL + Box2Abb[3189]*m_cR;

  Box2Abb[3191]=-13. + 10.*m_z12;

  Box2Abb[3192]=-2. + 2.*Box2Abb[1627]*m_z12 + 6.*m_z2k + Box2Abb[3191]*m_z12*m_z2k + 4.*Box2Abb[693]*m_z2k_2;

  Box2Abb[3193]=27. - 8.*m_z12;

  Box2Abb[3194]=-26. + Box2Abb[3193]*m_z12 + 12.*m_z2k;

  Box2Abb[3195]=8. + Box2Abb[3194]*m_z12;

  Box2Abb[3196]=11. + 2.*Box2Abb[201]*m_z12 - 2.*m_z2k - 5.*m_z12*m_z2k - 12.*m_z2k_2;

  Box2Abb[3197]=-7. + Box2Abb[3196]*m_z12 + 8.*m_z2k;

  Box2Abb[3198]=2. + Box2Abb[3197]*m_z12 - 4.*m_z2k;

  Box2Abb[3199]=Box2Abb[3198]*m_x_2 + Box2Abb[3195]*m_x_3 + 4.*Box2Abb[72]*m_x_4*m_z12 + Box2Abb[3192]*m_x*m_z12*m_z2k - Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[3200]=3. + Box2Abb[62]*m_z12;

  Box2Abb[3201]=-2. + Box2Abb[3200]*m_z12 - 6.*m_z2k + 3.*Box2Abb[703]*m_z12*m_z2k + 4.*Box2Abb[11]*m_z2k_2;

  Box2Abb[3202]=-45. + 34.*m_z12 + 48.*m_z2k;

  Box2Abb[3203]=22. + Box2Abb[3202]*m_z12 - 12.*m_z2k;

  Box2Abb[3204]=-8. + Box2Abb[3203]*m_z12;

  Box2Abb[3205]=-29. + 49.*m_z2k;

  Box2Abb[3206]=25. + Box2Abb[3205]*m_z12 + 15.*m_z12_2 + 36.*Box2Abb[61]*m_z2k;

  Box2Abb[3207]=17. - Box2Abb[3206]*m_z12;

  Box2Abb[3208]=-6. + Box2Abb[3207]*m_z12 + 4.*m_z2k;

  Box2Abb[3209]=Box2Abb[3208]*m_x_2 + Box2Abb[3204]*m_x_3 + Box2Abb[3201]*Box2Abb[5]*m_x*m_z12 + 4.*Box2Abb[879]*m_x_4*m_z12 - pow(Box2Abb[5],2.)*Box2Abb[70]*m_z12_2*m_z2k;

  Box2Abb[3210]=Box2Abb[3199]*m_cL + Box2Abb[3209]*m_cR;

  Box2Abb[3211]=Box2Abb[707]*m_cL - Box2Abb[429]*m_cR;

  Box2Abb[3212]=Box2Abb[435]*m_cL + Box2Abb[429]*m_cR;

  Box2Abb[3213]=Box2Abb[1340]*m_cR - Box2Abb[429]*pow(Box2Abb[8],2.)*m_cL*m_z2k;

  Box2Abb[3214]=m_cL - m_cR;

  Box2Abb[3215]=1. + m_x - m_z12;

  Box2Abb[3216]=-2.*m_x + Box2Abb[1]*m_z12;

  Box2Abb[3217]=2.*Box2Abb[210]*m_x + m_z12 - 3.*m_x*m_z12 + Box2Abb[1]*m_z12_2;

  Box2Abb[3218]=Box2Abb[3217]*m_cL + Box2Abb[3215]*Box2Abb[3216]*m_cR;

  Box2Abb[3219]=m_x + 3.*Box2Abb[210]*m_z12 - 3.*m_z12_2;

  Box2Abb[3220]=m_x + Box2Abb[2273]*m_z12 - m_z12_2;

  Box2Abb[3221]=-Box2Abb[3219]*m_cL + Box2Abb[3220]*m_cR;

  Box2Abb[3222]=2.*m_cL - m_cR;

  Box2Abb[3223]=pow(Box2Abb[1],2.)*Box2Abb[3214]*Box2Abb[4]*m_x + Box2Abb[3218]*m_z2k + Box2Abb[3221]*m_z2k_2 + Box2Abb[3222]*m_z12*m_z2k_3;

  Box2Abb[3224]=-1. + m_x + 2.*m_z12 - 3.*m_x*m_z12;

  Box2Abb[3225]=Box2Abb[1]*Box2Abb[4]*m_cL + Box2Abb[3224]*m_cR;

  Box2Abb[3226]=Box2Abb[72]*m_x - m_z12;

  Box2Abb[3227]=-2. + 6.*m_z12;

  Box2Abb[3228]=-2. + Box2Abb[3227]*m_x + Box2Abb[1046]*m_z12;

  Box2Abb[3229]=-Box2Abb[3215]*Box2Abb[3226]*m_cL + Box2Abb[3228]*m_cR*m_x;

  Box2Abb[3230]=Box2Abb[3220]*m_cL + Box2Abb[1485]*m_cR*m_x;

  Box2Abb[3231]=Box2Abb[1]*Box2Abb[3225]*m_x + Box2Abb[3229]*m_z2k - Box2Abb[3230]*m_z2k_2 + m_cL*m_z12*m_z2k_3;

  Box2Abb[3232]=1. + Box2Abb[1485]*m_x - 2.*m_z12;

  Box2Abb[3233]=m_x - 3.*m_x*m_z12 + m_z12_2;

  Box2Abb[3234]=2. + 2.*Box2Abb[3233] - 5.*m_z12;

  Box2Abb[3235]=Box2Abb[1]*Box2Abb[3232] + Box2Abb[3234]*m_z2k + Box2Abb[1485]*m_z2k_2;

  Box2Abb[3236]=1. - 4.*m_z12;

  Box2Abb[3237]=2. + 2.*m_z12 + 17.*m_z2k;

  Box2Abb[3238]=-2. + Box2Abb[3237]*m_z12 - 4.*m_z2k;

  Box2Abb[3239]=5. + 7.*m_z2k;

  Box2Abb[3240]=5. - 28.*m_z2k;

  Box2Abb[3241]=4. + Box2Abb[3240]*m_z2k;

  Box2Abb[3242]=1. + Box2Abb[3241]*m_z12 - Box2Abb[3239]*m_z12_2 + 2.*m_z2k + 6.*m_z2k_2;

  Box2Abb[3243]=-2. + Box2Abb[509]*m_z2k;

  Box2Abb[3244]=-7. + 8.*m_z2k;

  Box2Abb[3245]=3. + Box2Abb[3244]*m_z2k;

  Box2Abb[3246]=-Box2Abb[3245]*Box2Abb[61]*m_z12 + Box2Abb[3243]*m_z12_2 + pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[3247]=3. + Box2Abb[2297]*m_z2k;

  Box2Abb[3248]=-19. + 22.*m_z2k;

  Box2Abb[3249]=11. + Box2Abb[3248]*m_z2k;

  Box2Abb[3250]=-4. + Box2Abb[3249]*m_z2k;

  Box2Abb[3251]=Box2Abb[3250]*m_z12 + Box2Abb[3247]*m_z12_2 + 2.*Box2Abb[687]*m_z2k;

  Box2Abb[3252]=Box2Abb[3251]*m_x_2 + Box2Abb[3242]*m_x_3 + Box2Abb[3238]*m_x_4 + Box2Abb[3236]*m_x_5 + Box2Abb[3246]*m_x*m_z2k + Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[3253]=Box2Abb[3252]*m_cR + Box2Abb[3235]*pow(Box2Abb[8],2.)*m_cL*m_x;

  Box2Abb[3254]=-1. + Box2Abb[234]*m_x - 2.*m_x_2 + 3.*m_z2k - 2.*Box2Abb[137]*m_z2k;

  Box2Abb[3255]=Box2Abb[1341]*m_cL + Box2Abb[3254]*m_cR*m_x;

  Box2Abb[3256]=7. - 12.*m_z2k;

  Box2Abb[3257]=4. - 3.*Box2Abb[12]*m_z12 + Box2Abb[3256]*m_z2k;

  Box2Abb[3258]=2. + 3.*Box2Abb[179]*m_z2k;

  Box2Abb[3259]=-11. + 8.*m_z2k;

  Box2Abb[3260]=8. + Box2Abb[3259]*m_z2k;

  Box2Abb[3261]=-3. + Box2Abb[3258]*m_z12 + Box2Abb[3260]*m_z2k;

  Box2Abb[3262]=Box2Abb[3261]*m_x + Box2Abb[3257]*m_x_2 + Box2Abb[1343]*m_x_3 - 2.*m_x_4 - Box2Abb[1345]*Box2Abb[61]*m_z2k;

  Box2Abb[3263]=Box2Abb[1342]*pow(Box2Abb[8],2.)*m_cL + Box2Abb[3262]*m_cR;

  Box2Abb[3264]=1. + Box2Abb[3215]*m_z12;

  Box2Abb[3265]=-2. - 4.*m_x + m_z12;

  Box2Abb[3266]=-pow(Box2Abb[72],2.)*m_x_2 + 2.*Box2Abb[3264]*m_x*m_z2k + Box2Abb[3265]*m_z12*m_z2k_2 + 2.*m_z12*m_z2k_3;

  Box2Abb[3267]=-2. + 5.*m_z12 + 2.*m_z2k;

  Box2Abb[3268]=4. - 2.*Box2Abb[184]*m_z12 + 5.*m_z12_2;

  Box2Abb[3269]=12. + 2.*Box2Abb[520]*m_z12 - 11.*m_z12*m_z2k - 6.*m_z2k_2;

  Box2Abb[3270]=2.*Box2Abb[61] + Box2Abb[3269]*m_z12;

  Box2Abb[3271]=8. + 3.*m_z2k;

  Box2Abb[3272]=4. + m_z12 - m_z12_2 + 2.*Box2Abb[3271]*m_z2k - m_z12*m_z2k;

  Box2Abb[3273]=-4. + Box2Abb[3272]*m_z12 - 6.*m_z2k;

  Box2Abb[3274]=Box2Abb[3273]*m_x_2 + Box2Abb[3268]*m_x_3 + Box2Abb[3270]*m_x*m_z2k + Box2Abb[3267]*Box2Abb[5]*m_z12*m_z2k_2;

  Box2Abb[3275]=Box2Abb[3266]*Box2Abb[8]*m_cL + Box2Abb[3274]*m_cR;

  Box2Abb[3276]=3. + m_x - m_z12;

  Box2Abb[3277]=-3. + Box2Abb[3276]*m_z12;

  Box2Abb[3278]=-3. + 2.*m_x + m_z12;

  Box2Abb[3279]=m_x - m_z12;

  Box2Abb[3280]=Box2Abb[3277]*m_x_2 - Box2Abb[3278]*m_x*m_z12*m_z2k + Box2Abb[3279]*m_z12*m_z2k_2;

  Box2Abb[3281]=pow(Box2Abb[4],3.) + 3.*pow(Box2Abb[4],2.)*m_z2k + 2.*Box2Abb[4]*m_z2k_2 + m_z2k_3;

  Box2Abb[3282]=1. + 5.*m_z12;

  Box2Abb[3283]=-pow(Box2Abb[4],3.) - Box2Abb[3282]*pow(Box2Abb[4],2.)*m_z2k - 3.*Box2Abb[1485]*Box2Abb[4]*m_z2k_2 - 4.*m_z12*m_z2k_3;

  Box2Abb[3284]=1. + 3.*Box2Abb[831]*m_z12_2 + 2.*m_z12_3 + 6.*m_z2k + 6.*Box2Abb[841]*m_z12*m_z2k;

  Box2Abb[3285]=-8. + 5.*m_z12 + 4.*m_z2k;

  Box2Abb[3286]=3. + Box2Abb[3285]*m_z12;

  Box2Abb[3287]=Box2Abb[3283]*m_x + Box2Abb[3284]*m_x_2 - Box2Abb[3286]*m_x_3 + m_x_4*m_z12 + Box2Abb[3281]*m_z12*m_z2k;

  Box2Abb[3288]=-Box2Abb[3280]*pow(Box2Abb[8],2.)*m_cL + Box2Abb[3287]*m_cR*m_x;

  Box2Abb[3289]=-3. + 2.*m_z12 + m_z2k;

  Box2Abb[3290]=1. + Box2Abb[3289]*m_z12;

  Box2Abb[3291]=9. + Box2Abb[1695]*m_z12 - 12.*m_z2k + 9.*m_z12*m_z2k + 3.*m_z2k_2;

  Box2Abb[3292]=2. - Box2Abb[3291]*m_z12 - 3.*m_z2k;

  Box2Abb[3293]=Box2Abb[3292]*m_x + 3.*Box2Abb[3290]*m_x_2 + pow(Box2Abb[5],3.)*m_z12 - m_x_3*m_z12;

  Box2Abb[3294]=Box2Abb[148]*Box2Abb[3280]*m_cL + Box2Abb[3293]*m_cR*m_x;

  Box2Abb[3295]=2. + Box2Abb[162]*m_z12;

  Box2Abb[3296]=-Box2Abb[3295]*m_x_2 + 2.*m_x_3*m_z12 + Box2Abb[9]*m_x*m_z2k + Box2Abb[61]*m_z12*m_z2k_2;

  Box2Abb[3297]=15. - 2.*m_z12;

  Box2Abb[3298]=-15. + Box2Abb[3297]*m_z12;

  Box2Abb[3299]=2. + Box2Abb[3298]*m_z12 - 6.*m_z2k + 6.*Box2Abb[70]*m_z12*m_z2k + 4.*m_z2k_2;

  Box2Abb[3300]=4. + Box2Abb[4]*Box2Abb[62]*m_z12 + 30.*m_z2k + 6.*Box2Abb[166]*m_z12*m_z2k + 12.*Box2Abb[1485]*m_z2k_2;

  Box2Abb[3301]=-4.*Box2Abb[256] + Box2Abb[3300]*m_z12;

  Box2Abb[3302]=24. - 11.*m_z12 - 36.*m_z2k;

  Box2Abb[3303]=-6. + 5.*m_z2k;

  Box2Abb[3304]=4.*Box2Abb[3303] + Box2Abb[3302]*m_z12;

  Box2Abb[3305]=8. + Box2Abb[3304]*m_z12;

  Box2Abb[3306]=2. + m_z2k + m_z2k_2;

  Box2Abb[3307]=-18. + 12.*Box2Abb[3306]*m_z12 + Box2Abb[801]*m_z12_2 + 3.*m_z12_3 + 4.*m_z2k_2;

  Box2Abb[3308]=2. + Box2Abb[3307]*m_z12 - 4.*m_z2k;

  Box2Abb[3309]=Box2Abb[3301]*m_x_3 + Box2Abb[3305]*m_x_4 + 4.*Box2Abb[155]*m_x_5*m_z12 - Box2Abb[3308]*m_x_2*m_z2k + Box2Abb[3299]*m_x*m_z12*m_z2k_2 - 3.*Box2Abb[5]*m_z12_3*m_z2k_3;

  Box2Abb[3310]=-2.*Box2Abb[3296]*Box2Abb[382]*Box2Abb[8]*m_cL + Box2Abb[3309]*m_cR;

  Box2Abb[3311]=5. + 2.*Box2Abb[457]*m_x;

  Box2Abb[3312]=5. + 7.*m_x;

  Box2Abb[3313]=-2. - 4.*m_x + Box2Abb[3311]*m_z12 - Box2Abb[3312]*m_z12_2 + 2.*m_z12_3;

  Box2Abb[3314]=2. + m_z12_3;

  Box2Abb[3315]=-12. + m_z12;

  Box2Abb[3316]=15. + Box2Abb[3315]*m_z12;

  Box2Abb[3317]=-6. + Box2Abb[3316]*m_z12;

  Box2Abb[3318]=Box2Abb[3314]*Box2Abb[4] + 2.*Box2Abb[3317]*m_x - 4.*Box2Abb[815]*m_x_2*m_z12;

  Box2Abb[3319]=1. - 6.*m_x_2;

  Box2Abb[3320]=3. + 14.*m_x;

  Box2Abb[3321]=-3.*pow(Box2Abb[1442],2.) + 4.*m_x + 2.*Box2Abb[3319]*m_z12 + Box2Abb[3320]*m_z12_3 - 2.*m_z12_4;

  Box2Abb[3322]=3. + Box2Abb[710]*m_z12;

  Box2Abb[3323]=-2.*Box2Abb[3322]*m_x + 4.*Box2Abb[1116]*m_x_2 + Box2Abb[4]*m_z12_2;

  Box2Abb[3324]=Box2Abb[1025]*m_x + m_z12_2;

  Box2Abb[3325]=Box2Abb[3313]*Box2Abb[72]*m_x_3 + Box2Abb[3318]*m_x_2*m_z2k + Box2Abb[3321]*m_x*m_z2k_2 + Box2Abb[3323]*m_z12*m_z2k_3 + Box2Abb[3324]*m_z12*m_z2k_4;

  Box2Abb[3326]=-pow(Box2Abb[4],3.) + Box2Abb[384]*m_z2k + m_z12*m_z2k_2;

  Box2Abb[3327]=-6. + 4.*m_z12 + 5.*m_z2k;

  Box2Abb[3328]=2. + Box2Abb[3327]*m_z12;

  Box2Abb[3329]=-8. + 5.*m_z12;

  Box2Abb[3330]=pow(Box2Abb[4],2.) + Box2Abb[3329]*m_z2k + 3.*m_z2k_2;

  Box2Abb[3331]=Box2Abb[3330]*m_z12 + 3.*m_z2k;

  Box2Abb[3332]=Box2Abb[3331]*m_x_2 - Box2Abb[3328]*m_x_3 + 2.*m_x_4*m_z12 + Box2Abb[3326]*m_x*m_z2k - pow(Box2Abb[5],2.)*m_z12*m_z2k_2;

  Box2Abb[3333]=Box2Abb[3325]*m_cL - 2.*Box2Abb[3332]*Box2Abb[388]*m_cR;

  Box2Abb[3334]=-Box2Abb[1627]*m_x + 2.*m_x_2*m_z12 - Box2Abb[2144]*m_z12*m_z2k;

  Box2Abb[3335]=2. - 7.*m_z12;

  Box2Abb[3336]=2. + m_z12 - 2.*m_z2k;

  Box2Abb[3337]=1. + m_z12 - 5.*m_z12_2;

  Box2Abb[3338]=2. - 2.*m_z12 + m_z12_2 - m_z12_3 + 2.*Box2Abb[3337]*m_z2k - 2.*Box2Abb[70]*m_z2k_2;

  Box2Abb[3339]=-2. + 7.*m_z12 + 18.*m_z2k;

  Box2Abb[3340]=2. + Box2Abb[3339]*m_z12;

  Box2Abb[3341]=-4. + Box2Abb[3340]*m_z12;

  Box2Abb[3342]=Box2Abb[3341]*m_x_2 + Box2Abb[3338]*m_x*m_z12 + 2.*Box2Abb[3335]*m_x_3*m_z12 + Box2Abb[3336]*Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[3343]=Box2Abb[3334]*Box2Abb[382]*m_cL + Box2Abb[3342]*m_cR;

  Box2Abb[3344]=-1. + m_z12 + 6.*m_z2k;

  Box2Abb[3345]=1. + Box2Abb[3344]*m_z12;

  Box2Abb[3346]=-Box2Abb[3345]*m_x_2 + 3.*m_x_3*m_z12 + Box2Abb[399]*m_x*m_z12*m_z2k - m_z12_2*m_z2k_2;

  Box2Abb[3347]=-6. + 7.*m_z12;

  Box2Abb[3348]=1. + 5.*Box2Abb[693]*m_z12;

  Box2Abb[3349]=pow(Box2Abb[4],2.)*Box2Abb[693]*Box2Abb[70] + 2.*Box2Abb[3348]*Box2Abb[4]*m_z2k + 6.*Box2Abb[53]*m_z12*m_z2k_2;

  Box2Abb[3350]=pow(Box2Abb[4],3.) + 3.*pow(Box2Abb[4],2.)*m_z2k + 4.*Box2Abb[4]*m_z2k_2 + m_z2k_3;

  Box2Abb[3351]=-1. + 7.*m_z12;

  Box2Abb[3352]=-pow(Box2Abb[4],3.) - Box2Abb[3351]*pow(Box2Abb[4],2.)*m_z2k - 3.*Box2Abb[1116]*Box2Abb[4]*m_z2k_2 + 2.*Box2Abb[341]*m_z2k_3;

  Box2Abb[3353]=16. - 9.*m_z12 - 22.*m_z2k;

  Box2Abb[3354]=-9. + Box2Abb[3353]*m_z12 + 18.*m_z2k;

  Box2Abb[3355]=2. + Box2Abb[3354]*m_z12;

  Box2Abb[3356]=Box2Abb[3349]*m_x_2 + Box2Abb[3355]*m_x_3 + Box2Abb[3352]*m_x*m_z12 + Box2Abb[3347]*m_x_4*m_z12 + Box2Abb[3350]*m_z12_2*m_z2k;

  Box2Abb[3357]=-Box2Abb[3346]*Box2Abb[382]*Box2Abb[8]*m_cL + Box2Abb[3356]*m_cR*m_x;

  Box2Abb[3358]=Box2Abb[272]*pow(Box2Abb[4],2.) + 3.*Box2Abb[220]*Box2Abb[4]*m_z2k + 3.*Box2Abb[155]*m_z2k_2;

  Box2Abb[3359]=27. - 14.*m_z12 - 15.*m_z2k;

  Box2Abb[3360]=3.*Box2Abb[725] + Box2Abb[3359]*m_z12;

  Box2Abb[3361]=2. + Box2Abb[3360]*m_z12;

  Box2Abb[3362]=Box2Abb[3361]*m_x_2 + Box2Abb[3358]*m_x*m_z12 + Box2Abb[3347]*m_x_3*m_z12 - pow(Box2Abb[5],3.)*m_z12_2;

  Box2Abb[3363]=Box2Abb[148]*Box2Abb[3346]*Box2Abb[382]*m_cL - Box2Abb[3362]*Box2Abb[8]*m_cR*m_x;

  Box2Abb[3364]=3. + Box2Abb[199]*m_z2k;

  Box2Abb[3365]=3. + Box2Abb[1042]*m_z12;

  Box2Abb[3366]=3. + 3.*Box2Abb[61]*m_z12 + m_z12_2 + 4.*Box2Abb[841]*m_z2k;

  Box2Abb[3367]=-4. + Box2Abb[3366]*m_z12;

  Box2Abb[3368]=6. - Box2Abb[411]*m_z12 + m_z12_2;

  Box2Abb[3369]=-4. + Box2Abb[3368]*m_z12 + 12.*m_z2k;

  Box2Abb[3370]=2. + Box2Abb[3369]*m_z12;

  Box2Abb[3371]=Box2Abb[3370]*m_x_3 - 4.*Box2Abb[3365]*m_x_4*m_z12 + 2.*m_x_5*m_z12_2 + Box2Abb[3367]*m_x_2*m_z12*m_z2k + Box2Abb[3364]*m_x*m_z12_2*m_z2k_2 + Box2Abb[61]*m_z12_3*m_z2k_3;

  Box2Abb[3372]=12. + 4.*Box2Abb[2651]*m_z12 + 15.*m_z12_2;

  Box2Abb[3373]=-34. + 15.*m_z12;

  Box2Abb[3374]=12. + Box2Abb[3373]*m_z12;

  Box2Abb[3375]=4. + 3.*Box2Abb[201]*m_z12;

  Box2Abb[3376]=Box2Abb[3374]*m_z12 + 6.*Box2Abb[3375]*m_z2k;

  Box2Abb[3377]=4. + Box2Abb[3376]*m_z12;

  Box2Abb[3378]=3. + 2.*Box2Abb[179]*m_z2k;

  Box2Abb[3379]=1. - 6.*Box2Abb[61]*m_z2k;

  Box2Abb[3380]=2. + 3.*Box2Abb[3379]*m_z12 - 3.*Box2Abb[3116]*m_z12_2 + m_z12_3 - 2.*Box2Abb[3378]*m_z2k;

  Box2Abb[3381]=-3. + 6.*Box2Abb[280]*m_z2k;

  Box2Abb[3382]=2. + Box2Abb[3381]*m_z12 + Box2Abb[2211]*m_z12_2 - 2.*m_z12_3 - 30.*m_z2k + 8.*m_z2k_3;

  Box2Abb[3383]=-2. + Box2Abb[3382]*m_z12;

  Box2Abb[3384]=Box2Abb[3377]*m_x_3 + Box2Abb[3383]*m_x_2*m_z12 - 2.*Box2Abb[3372]*m_x_4*m_z12 + 4.*m_x_5*m_z12_2 + Box2Abb[3380]*m_x*m_z12_2*m_z2k + 3.*Box2Abb[5]*m_z12_4*m_z2k_2;

  Box2Abb[3385]=2.*Box2Abb[3371]*m_cL - Box2Abb[3384]*m_cR;

  Box2Abb[3386]=12. + 4.*Box2Abb[1359]*m_z12 + 7.*m_z12_2;

  Box2Abb[3387]=-10. + m_z2k;

  Box2Abb[3388]=9. + 2.*Box2Abb[2534]*m_z2k;

  Box2Abb[3389]=-2. + Box2Abb[3388]*m_z12 + Box2Abb[3387]*m_z12_2 + 3.*m_z12_3 + 6.*m_z2k + 4.*Box2Abb[179]*m_z2k_2;

  Box2Abb[3390]=-16. + 7.*m_z2k;

  Box2Abb[3391]=9. + 11.*m_z2k;

  Box2Abb[3392]=9. + 4.*Box2Abb[1359]*m_z2k;

  Box2Abb[3393]=10. - 15.*m_z12 + Box2Abb[3391]*m_z12_2 - 2.*m_z12_3 + 2.*Box2Abb[3392]*m_z2k + 2.*Box2Abb[3390]*m_z12*m_z2k;

  Box2Abb[3394]=-2. + Box2Abb[3393]*m_z12;

  Box2Abb[3395]=-42. + 11.*m_z12 + 2.*m_z2k;

  Box2Abb[3396]=60. + Box2Abb[3395]*m_z12 - 8.*m_z2k;

  Box2Abb[3397]=8.*Box2Abb[2800] + Box2Abb[3396]*m_z12;

  Box2Abb[3398]=4. + Box2Abb[3397]*m_z12;

  Box2Abb[3399]=Box2Abb[3398]*m_x_3 + Box2Abb[3394]*m_x_2*m_z12 - 2.*Box2Abb[3386]*m_x_4*m_z12 + 4.*m_x_5*m_z12_2 - Box2Abb[3389]*m_x*m_z12_2*m_z2k - Box2Abb[5]*m_z12_4*m_z2k_2;

  Box2Abb[3400]=pow(Box2Abb[4],3.) - 3.*Box2Abb[4]*Box2Abb[693]*m_z2k + 2.*Box2Abb[341]*m_z2k_2 - 2.*m_z2k_3;

  Box2Abb[3401]=4. + Box2Abb[1453]*m_z12;

  Box2Abb[3402]=-8. + m_z12;

  Box2Abb[3403]=4. + Box2Abb[3402]*m_z12;

  Box2Abb[3404]=-pow(Box2Abb[4],3.)*m_z12_2 - Box2Abb[3401]*pow(Box2Abb[4],2.)*m_z2k - Box2Abb[3403]*Box2Abb[4]*m_z2k_2 + 6.*Box2Abb[155]*m_z12*m_z2k_3 + 6.*m_z12*m_z2k_4;

  Box2Abb[3405]=-6. + 4.*m_z12 + m_z2k;

  Box2Abb[3406]=2. + Box2Abb[3405]*m_z12;

  Box2Abb[3407]=-89. + 40.*m_z12;

  Box2Abb[3408]=63. + Box2Abb[3407]*m_z12;

  Box2Abb[3409]=-16. + Box2Abb[3408]*m_z12;

  Box2Abb[3410]=3. + Box2Abb[1700]*m_z12;

  Box2Abb[3411]=3.*Box2Abb[155]*pow(Box2Abb[4],2.)*m_z12 + Box2Abb[3409]*m_z2k + 4.*Box2Abb[3410]*m_z2k_2 - 4.*m_z12*m_z2k_3;

  Box2Abb[3412]=Box2Abb[3411]*m_z12 + 2.*m_z2k;

  Box2Abb[3413]=-19. + m_z2k;

  Box2Abb[3414]=-58. + 51.*m_z2k;

  Box2Abb[3415]=42. + Box2Abb[3414]*m_z12 + 26.*m_z12_2 + 4.*Box2Abb[3413]*m_z2k;

  Box2Abb[3416]=-12. + Box2Abb[3415]*m_z12 + 24.*m_z2k;

  Box2Abb[3417]=2. + Box2Abb[3416]*m_z12;

  Box2Abb[3418]=Box2Abb[3412]*m_x_3 - Box2Abb[3417]*m_x_4 + Box2Abb[3404]*m_x_2*m_z12 + 6.*Box2Abb[3406]*m_x_5*m_z12 - 2.*m_x_6*m_z12_2 + Box2Abb[3400]*Box2Abb[5]*m_x*m_z12_2*m_z2k + pow(Box2Abb[5],3.)*m_z12_3*m_z2k_2;

  Box2Abb[3419]=Box2Abb[3399]*Box2Abb[8]*m_cL + 2.*Box2Abb[3418]*m_cR;

  Box2Abb[3420]=3. + 4.*m_x - 3.*m_z12;

  Box2Abb[3421]=-1. + m_x - 2.*m_x_2 + m_z12 - m_x*m_z12 + Box2Abb[3420]*m_z2k - 2.*m_z2k_2;

  Box2Abb[3422]=-1. + 3.*m_z12 + 2.*m_z2k;

  Box2Abb[3423]=m_x - Box2Abb[234]*m_x_2 + 2.*m_x_3 + Box2Abb[3422]*m_x*m_z2k + Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[3424]=Box2Abb[3423]*m_cL + Box2Abb[3421]*m_cR*m_x;

  Box2Abb[3425]=-1. + Box2Abb[201]*m_x + m_z12;

  Box2Abb[3426]=23. - 24.*m_x;

  Box2Abb[3427]=-11. + Box2Abb[3426]*m_x;

  Box2Abb[3428]=6. + Box2Abb[3427]*m_x;

  Box2Abb[3429]=-3. + m_x + m_x_2;

  Box2Abb[3430]=2. + m_x + 3.*Box2Abb[3429]*m_x_2;

  Box2Abb[3431]=Box2Abb[3428]*m_x + Box2Abb[3430]*m_z12 - 2.*pow(Box2Abb[1],2.)*Box2Abb[210]*m_z12_2;

  Box2Abb[3432]=11. + m_z12;

  Box2Abb[3433]=20. + m_z12;

  Box2Abb[3434]=14. + Box2Abb[520]*m_z12;

  Box2Abb[3435]=-21. + 2.*m_z12;

  Box2Abb[3436]=-38. + Box2Abb[3435]*m_z12;

  Box2Abb[3437]=13. + m_z12;

  Box2Abb[3438]=-29. + Box2Abb[3437]*m_z12;

  Box2Abb[3439]=6. + Box2Abb[3438]*m_z12;

  Box2Abb[3440]=Box2Abb[3439]*m_x_2 - Box2Abb[3432]*Box2Abb[72]*Box2Abb[9]*m_x_3 + Box2Abb[3436]*m_x_4 + 3.*Box2Abb[3433]*m_x_5 + Box2Abb[3434]*m_x*m_z12 - Box2Abb[4]*m_z12_2;

  Box2Abb[3441]=9. + 40.*m_x;

  Box2Abb[3442]=11. + Box2Abb[3441]*m_x;

  Box2Abb[3443]=5. + 2.*Box2Abb[3442]*m_x;

  Box2Abb[3444]=-24. + 25.*m_x;

  Box2Abb[3445]=-9. + Box2Abb[3444]*m_x;

  Box2Abb[3446]=6. + Box2Abb[3445]*m_x;

  Box2Abb[3447]=2. + Box2Abb[2354]*m_x;

  Box2Abb[3448]=-4. + Box2Abb[3447]*m_x;

  Box2Abb[3449]=8. + m_x;

  Box2Abb[3450]=-5. + Box2Abb[3449]*m_x;

  Box2Abb[3451]=Box2Abb[3443]*m_x + Box2Abb[3446]*m_x*m_z12 - 2.*Box2Abb[3448]*m_z12_2 + Box2Abb[3450]*m_z12_3;

  Box2Abb[3452]=4. + 3.*m_z12;

  Box2Abb[3453]=9. - 28.*m_z12;

  Box2Abb[3454]=37. + Box2Abb[3453]*m_z12;

  Box2Abb[3455]=3. - Box2Abb[1629]*m_z12;

  Box2Abb[3456]=14. + m_z12 - 14.*m_z12_2 + 5.*m_z12_3;

  Box2Abb[3457]=Box2Abb[3456]*m_x + Box2Abb[3454]*m_x_2 + 15.*Box2Abb[3452]*m_x_3 + Box2Abb[3455]*m_z12;

  Box2Abb[3458]=8. + 13.*m_z12;

  Box2Abb[3459]=13. + 24.*m_z12 - 22.*m_z12_2;

  Box2Abb[3460]=Box2Abb[3459]*m_x + 3.*Box2Abb[3458]*m_x_2 + 3.*Box2Abb[1451]*m_z12;

  Box2Abb[3461]=4. + 17.*m_z12;

  Box2Abb[3462]=-Box2Abb[3461]*m_x + 3.*Box2Abb[62]*m_z12;

  Box2Abb[3463]=pow(Box2Abb[1],4.)*Box2Abb[3425]*m_x_2 - Box2Abb[1]*Box2Abb[3431]*m_x*m_z2k - Box2Abb[3440]*m_z2k_2 + Box2Abb[3451]*m_z2k_3 - Box2Abb[3457]*m_z2k_4 + Box2Abb[3460]*m_z2k_5 + Box2Abb[3462]*m_z2k_6 + 3.*m_z12*m_z2k_7;

  Box2Abb[3464]=4. - 5.*m_z12;

  Box2Abb[3465]=20. - 27.*m_z12;

  Box2Abb[3466]=19. - 18.*m_z12 + 47.*m_z2k - 10.*Box2Abb[2741]*m_z12*m_z2k + 3.*Box2Abb[3465]*m_z2k_2;

  Box2Abb[3467]=7. - 4.*m_z12;

  Box2Abb[3468]=-13. + Box2Abb[70]*m_z12;

  Box2Abb[3469]=2. + 11.*m_z12;

  Box2Abb[3470]=19. + Box2Abb[3469]*m_z12;

  Box2Abb[3471]=-18. + Box2Abb[3470]*m_z12;

  Box2Abb[3472]=37. - 76.*m_z12 + 52.*m_z12_2;

  Box2Abb[3473]=-8. + 15.*m_z12;

  Box2Abb[3474]=-6. + 2.*Box2Abb[3467]*m_z12 + 11.*m_z2k + 2.*Box2Abb[3468]*m_z12*m_z2k + Box2Abb[3471]*m_z2k_2 + Box2Abb[3472]*m_z2k_3 + 3.*Box2Abb[3473]*m_z2k_4;

  Box2Abb[3475]=-7. + 2.*m_z12;

  Box2Abb[3476]=-17. + 5.*Box2Abb[669]*m_z12;

  Box2Abb[3477]=-16. + Box2Abb[3476]*m_z12;

  Box2Abb[3478]=9. - 22.*m_z12 + 34.*m_z12_2;

  Box2Abb[3479]=-1. + m_z12 + Box2Abb[3475]*m_z2k + Box2Abb[3477]*m_z2k_2 + 2.*Box2Abb[3478]*m_z2k_3 + 5.*Box2Abb[1533]*m_z2k_4;

  Box2Abb[3480]=5. + 8.*m_z2k;

  Box2Abb[3481]=-3.*Box2Abb[3480] + Box2Abb[1551]*m_z12;

  Box2Abb[3482]=-2. + Box2Abb[831]*m_z2k;

  Box2Abb[3483]=5. - Box2Abb[1567]*m_z2k;

  Box2Abb[3484]=-16. + 11.*m_z2k;

  Box2Abb[3485]=4. + Box2Abb[3484]*m_z2k_2;

  Box2Abb[3486]=-4. + m_z2k + 11.*m_z2k_2 - 9.*m_z2k_3;

  Box2Abb[3487]=1. + Box2Abb[3486]*m_z2k;

  Box2Abb[3488]=Box2Abb[3482]*pow(Box2Abb[61],3.) - Box2Abb[3485]*pow(Box2Abb[61],2.)*m_z12 + 2.*Box2Abb[3487]*m_z12_2 + Box2Abb[3483]*m_z12_3*m_z2k;

  Box2Abb[3489]=19. + 40.*m_z2k;

  Box2Abb[3490]=15. + Box2Abb[3489]*m_z2k;

  Box2Abb[3491]=34. + 115.*m_z2k;

  Box2Abb[3492]=11. + Box2Abb[3491]*m_z2k;

  Box2Abb[3493]=8. + Box2Abb[3492]*m_z2k;

  Box2Abb[3494]=-9. + Box2Abb[3493]*m_z12 - 2.*Box2Abb[3490]*m_z2k + 2.*Box2Abb[1558]*m_z12_2*m_z2k;

  Box2Abb[3495]=-Box2Abb[3479]*m_x_3 + Box2Abb[3494]*m_x_4 + Box2Abb[3466]*m_x_5 + Box2Abb[3481]*m_x_6 + Box2Abb[3464]*m_x_7 + Box2Abb[3488]*m_x*m_z2k + Box2Abb[3474]*m_x_2*m_z2k + Box2Abb[179]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[3496]=Box2Abb[3463]*m_cL + Box2Abb[3495]*m_cR;

  Box2Abb[3497]=-6. + Box2Abb[3452]*m_z12;

  Box2Abb[3498]=9. + 7.*Box2Abb[72]*m_z12;

  Box2Abb[3499]=1. - m_z12 + Box2Abb[3497]*m_z2k + Box2Abb[3498]*m_z2k_2 + Box2Abb[1305]*m_z2k_3;

  Box2Abb[3500]=9. + 31.*m_z2k;

  Box2Abb[3501]=-11. + Box2Abb[3500]*m_z12 - 24.*m_z2k;

  Box2Abb[3502]=9. + 20.*m_z2k;

  Box2Abb[3503]=16. + 81.*m_z2k;

  Box2Abb[3504]=-6. + 3.*m_z12 + 2.*Box2Abb[510]*m_z12_2 - 3.*Box2Abb[3502]*m_z2k + Box2Abb[3503]*m_z12*m_z2k;

  Box2Abb[3505]=15. + 22.*m_z2k - 68.*m_z2k_2;

  Box2Abb[3506]=-6. + Box2Abb[3505]*m_z2k;

  Box2Abb[3507]=-29. + 30.*m_z2k;

  Box2Abb[3508]=11. + Box2Abb[3507]*m_z2k;

  Box2Abb[3509]=-1. + m_z2k + Box2Abb[3508]*m_z2k_2;

  Box2Abb[3510]=-9. + 95.*m_z2k;

  Box2Abb[3511]=22. + Box2Abb[3510]*Box2Abb[61]*m_z2k;

  Box2Abb[3512]=8. - Box2Abb[3511]*m_z2k;

  Box2Abb[3513]=2.*Box2Abb[3509] + Box2Abb[3512]*m_z12 + Box2Abb[3506]*m_z12_2 - 5.*m_z12_3*m_z2k_2;

  Box2Abb[3514]=-3. + m_z2k - 40.*m_z2k_2;

  Box2Abb[3515]=3. + m_z2k + 21.*m_z2k_2;

  Box2Abb[3516]=-31. + 115.*m_z2k;

  Box2Abb[3517]=-9. + Box2Abb[3516]*m_z2k;

  Box2Abb[3518]=-5. + Box2Abb[3517]*m_z2k;

  Box2Abb[3519]=1. + Box2Abb[3518]*m_z12 + 2.*Box2Abb[3515]*m_z12_2 + 2.*Box2Abb[3514]*m_z2k;

  Box2Abb[3520]=-9. + 11.*m_z2k;

  Box2Abb[3521]=8. + 3.*Box2Abb[1163]*m_z2k;

  Box2Abb[3522]=-2. + Box2Abb[3521]*m_z2k;

  Box2Abb[3523]=-11. + 45.*m_z2k;

  Box2Abb[3524]=-14. + Box2Abb[3523]*Box2Abb[61]*m_z2k;

  Box2Abb[3525]=4. + Box2Abb[3524]*m_z2k;

  Box2Abb[3526]=-25. + 26.*m_z2k;

  Box2Abb[3527]=5. + Box2Abb[3526]*m_z2k;

  Box2Abb[3528]=-3. + Box2Abb[3527]*m_z2k;

  Box2Abb[3529]=1. + Box2Abb[3528]*m_z2k;

  Box2Abb[3530]=-Box2Abb[3522]*pow(Box2Abb[61],2.) + Box2Abb[3525]*Box2Abb[61]*m_z12 + 2.*Box2Abb[3529]*m_z12_2 + Box2Abb[3520]*m_z12_3*m_z2k_2;

  Box2Abb[3531]=Box2Abb[3530]*m_x_2 + Box2Abb[3513]*m_x_3 + Box2Abb[3519]*m_x_4 - Box2Abb[3504]*m_x_5 + Box2Abb[3501]*m_x_6 + Box2Abb[3464]*m_x_7 - Box2Abb[3499]*Box2Abb[5]*Box2Abb[61]*m_x*m_z2k + pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[3532]=11. + 2.*m_z12 - 3.*Box2Abb[3402]*m_z2k;

  Box2Abb[3533]=-11. + m_z12;

  Box2Abb[3534]=5. + 14.*m_z12 + 27.*m_z2k - 2.*Box2Abb[3533]*m_z12*m_z2k + 3.*Box2Abb[3433]*m_z2k_2;

  Box2Abb[3535]=41. + 10.*m_z12;

  Box2Abb[3536]=2. + Box2Abb[3535]*m_z12;

  Box2Abb[3537]=-28. + m_z12;

  Box2Abb[3538]=1. + Box2Abb[3537]*m_z12;

  Box2Abb[3539]=16. + 5.*m_z12;

  Box2Abb[3540]=-9. + 20.*m_z12 + Box2Abb[3536]*m_z2k - 2.*Box2Abb[3538]*m_z2k_2 + 5.*Box2Abb[3539]*m_z2k_3;

  Box2Abb[3541]=13. + 3.*Box2Abb[9]*m_z12;

  Box2Abb[3542]=-17. + 2.*Box2Abb[3541]*m_z12;

  Box2Abb[3543]=-23. + m_z12_2;

  Box2Abb[3544]=-57. + 14.*Box2Abb[693]*m_z12;

  Box2Abb[3545]=-2. + 2.*m_z12 + 6.*m_z2k + 22.*Box2Abb[4]*m_z12*m_z2k + Box2Abb[3542]*m_z2k_2 + Box2Abb[3543]*Box2Abb[72]*m_z2k_3 + Box2Abb[3544]*m_z2k_4 + 3.*Box2Abb[3458]*m_z2k_5;

  Box2Abb[3546]=Box2Abb[150]*m_z12 + 3.*Box2Abb[61]*m_z2k;

  Box2Abb[3547]=-4. - 8.*m_z2k + 17.*m_z2k_3;

  Box2Abb[3548]=9. + Box2Abb[3303]*m_z2k;

  Box2Abb[3549]=-4. + 11.*m_z2k;

  Box2Abb[3550]=2. + Box2Abb[3549]*m_z2k;

  Box2Abb[3551]=-2. + Box2Abb[3550]*m_z2k;

  Box2Abb[3552]=Box2Abb[3547]*pow(Box2Abb[61],2.)*m_z12 + 2.*Box2Abb[3551]*Box2Abb[61]*m_z12_2 + pow(Box2Abb[61],3.)*Box2Abb[725]*m_z2k + Box2Abb[3548]*m_z12_3*m_z2k;

  Box2Abb[3553]=11. + 6.*m_z2k;

  Box2Abb[3554]=15. + Box2Abb[3553]*m_z2k;

  Box2Abb[3555]=8. + Box2Abb[3507]*m_z2k;

  Box2Abb[3556]=1. + 2.*Box2Abb[3555]*m_z2k;

  Box2Abb[3557]=44. + 45.*m_z2k;

  Box2Abb[3558]=13. + Box2Abb[3557]*m_z2k;

  Box2Abb[3559]=-2. + Box2Abb[3558]*m_z2k;

  Box2Abb[3560]=11. + Box2Abb[3559]*m_z2k;

  Box2Abb[3561]=9. - Box2Abb[3560]*m_z12 - Box2Abb[3556]*m_z2k - 2.*Box2Abb[3554]*m_z12_2*m_z2k + m_z12_3*m_z2k_2;

  Box2Abb[3562]=Box2Abb[3545]*m_x_2 + Box2Abb[3561]*m_x_3 + Box2Abb[3540]*m_x_4 - Box2Abb[3534]*m_x_5 + Box2Abb[3532]*m_x_6 + Box2Abb[201]*m_x_7 - Box2Abb[3552]*m_x*m_z2k + Box2Abb[3546]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[3563]=Box2Abb[3562]*m_cL + Box2Abb[3531]*m_cR;

  Box2Abb[3564]=-11. + 9.*m_z12;

  Box2Abb[3565]=16. + Box2Abb[1301]*m_z12;

  Box2Abb[3566]=-2. + Box2Abb[3565]*m_z12;

  Box2Abb[3567]=-26. + 5.*m_z12;

  Box2Abb[3568]=13. + Box2Abb[3567]*m_z12;

  Box2Abb[3569]=-10. + Box2Abb[3568]*m_z12;

  Box2Abb[3570]=29. - 72.*m_z12 + 34.*m_z12_2;

  Box2Abb[3571]=6. + 2.*Box2Abb[3564]*m_z12 + Box2Abb[3566]*m_z2k + Box2Abb[3569]*m_z2k_2 + 2.*Box2Abb[3570]*m_z2k_3 + 5.*Box2Abb[1533]*m_z2k_4;

  Box2Abb[3572]=7. + 31.*m_z2k;

  Box2Abb[3573]=-11. + Box2Abb[3572]*m_z12 - 24.*m_z2k;

  Box2Abb[3574]=4. + 81.*m_z2k;

  Box2Abb[3575]=3. - Box2Abb[3574]*m_z2k;

  Box2Abb[3576]=4. + Box2Abb[3575]*m_z12 - 2.*Box2Abb[450]*m_z12_2 + 3.*Box2Abb[3502]*m_z2k;

  Box2Abb[3577]=-6. + m_z2k + 7.*m_z2k_2;

  Box2Abb[3578]=-15. + 11.*m_z2k;

  Box2Abb[3579]=-1. + Box2Abb[3578]*m_z2k;

  Box2Abb[3580]=2. + Box2Abb[3579]*m_z2k;

  Box2Abb[3581]=-Box2Abb[3482]*pow(Box2Abb[61],2.) + Box2Abb[3580]*Box2Abb[61]*m_z12 + Box2Abb[3577]*m_z12_2*m_z2k;

  Box2Abb[3582]=1. + m_z2k - 40.*m_z2k_2;

  Box2Abb[3583]=1. + 7.*m_z2k;

  Box2Abb[3584]=7. + 3.*Box2Abb[3583]*m_z2k;

  Box2Abb[3585]=-61. + 115.*m_z2k;

  Box2Abb[3586]=-23. + Box2Abb[3585]*m_z2k;

  Box2Abb[3587]=-15. + Box2Abb[3586]*m_z2k;

  Box2Abb[3588]=5. + Box2Abb[3587]*m_z12 + 2.*Box2Abb[3584]*m_z12_2 + 2.*Box2Abb[3582]*m_z2k;

  Box2Abb[3589]=-2. + Box2Abb[1163]*m_z2k_2;

  Box2Abb[3590]=6. + Box2Abb[2130]*m_z2k;

  Box2Abb[3591]=-33. + 26.*m_z2k;

  Box2Abb[3592]=5. + Box2Abb[3591]*m_z2k;

  Box2Abb[3593]=-7. + Box2Abb[3592]*m_z2k;

  Box2Abb[3594]=5. + Box2Abb[3593]*m_z2k;

  Box2Abb[3595]=-86. + 45.*m_z2k;

  Box2Abb[3596]=17. + Box2Abb[3595]*m_z2k;

  Box2Abb[3597]=-10. + Box2Abb[3596]*m_z2k;

  Box2Abb[3598]=16. + Box2Abb[3597]*m_z2k;

  Box2Abb[3599]=-3.*Box2Abb[3589]*pow(Box2Abb[61],2.) + Box2Abb[3598]*Box2Abb[61]*m_z12 + 2.*Box2Abb[3594]*m_z12_2 + Box2Abb[3590]*m_z12_3*m_z2k;

  Box2Abb[3600]=-Box2Abb[3581]*Box2Abb[5]*Box2Abb[61]*m_x + Box2Abb[3599]*m_x_2 - Box2Abb[3571]*m_x_3 + Box2Abb[3588]*m_x_4 + Box2Abb[3576]*m_x_5 + Box2Abb[3573]*m_x_6 + Box2Abb[3464]*m_x_7 + Box2Abb[179]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[3601]=Box2Abb[3562]*m_cL + Box2Abb[3600]*m_cR;

  Box2Abb[3602]=-24. + 31.*m_z12;

  Box2Abb[3603]=-7. + Box2Abb[3602]*m_z2k;

  Box2Abb[3604]=23. + 2.*m_z12;

  Box2Abb[3605]=-22. + Box2Abb[3604]*m_z12;

  Box2Abb[3606]=1. + Box2Abb[538]*m_z12;

  Box2Abb[3607]=-16. + 23.*m_z12;

  Box2Abb[3608]=2. + Box2Abb[3605]*m_z12 + 10.*m_z2k + 8.*m_z12_2*m_z2k + 42.*Box2Abb[3606]*m_z2k_2 + 5.*Box2Abb[3607]*m_z2k_3;

  Box2Abb[3609]=11. + 6.*m_z12;

  Box2Abb[3610]=-26. + Box2Abb[3609]*m_z12;

  Box2Abb[3611]=10. + m_z12 + 8.*m_z12_2;

  Box2Abb[3612]=-54. + 5.*m_z12;

  Box2Abb[3613]=101. + Box2Abb[3612]*m_z12;

  Box2Abb[3614]=-40. + Box2Abb[3613]*m_z12;

  Box2Abb[3615]=49. + 34.*Box2Abb[538]*m_z12;

  Box2Abb[3616]=10. + Box2Abb[3610]*m_z12 - 8.*m_z2k + Box2Abb[3611]*m_z12*m_z2k + Box2Abb[3614]*m_z2k_2 + 2.*Box2Abb[3615]*m_z2k_3 + 5.*Box2Abb[1533]*m_z2k_4;

  Box2Abb[3617]=7. + 60.*m_z2k;

  Box2Abb[3618]=7. + 30.*m_z2k - 81.*m_z2k_2;

  Box2Abb[3619]=-5. + Box2Abb[3618]*m_z12 - 10.*Box2Abb[12]*m_z12_2 + Box2Abb[3617]*m_z2k;

  Box2Abb[3620]=2. + 7.*m_z2k_2;

  Box2Abb[3621]=-6. + 11.*m_z2k;

  Box2Abb[3622]=3. + Box2Abb[3621]*m_z2k;

  Box2Abb[3623]=-pow(Box2Abb[61],2.)*Box2Abb[831] + Box2Abb[3622]*Box2Abb[61]*m_z12 + Box2Abb[3620]*m_z12_2;

  Box2Abb[3624]=-5. + 24.*m_z2k;

  Box2Abb[3625]=3. + Box2Abb[3624]*m_z2k;

  Box2Abb[3626]=9. - 50.*m_z2k + 52.*m_z2k_2;

  Box2Abb[3627]=7. + Box2Abb[3626]*m_z2k;

  Box2Abb[3628]=-12. + 11.*m_z2k;

  Box2Abb[3629]=3. + Box2Abb[3628]*m_z2k;

  Box2Abb[3630]=6. + Box2Abb[3629]*m_z2k;

  Box2Abb[3631]=-22. + 15.*m_z2k;

  Box2Abb[3632]=4. + Box2Abb[3631]*m_z2k;

  Box2Abb[3633]=-2. + 3.*Box2Abb[3632]*m_z2k;

  Box2Abb[3634]=-Box2Abb[3625]*pow(Box2Abb[61],3.) + Box2Abb[3633]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3627]*Box2Abb[61]*m_z12_2 + Box2Abb[3630]*m_z12_3;

  Box2Abb[3635]=-Box2Abb[3623]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x + Box2Abb[3634]*m_x_2 - Box2Abb[3616]*m_x_3 + Box2Abb[3608]*m_x_4 + Box2Abb[3619]*m_x_5 + Box2Abb[3603]*m_x_6 + Box2Abb[3464]*m_x_7 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12*m_z2k;

  Box2Abb[3636]=7. - 3.*m_z2k;

  Box2Abb[3637]=7. + Box2Abb[3636]*m_z12 + 24.*m_z2k;

  Box2Abb[3638]=-22. + m_z12;

  Box2Abb[3639]=6. + Box2Abb[3435]*m_z12 - 7.*m_z2k + 2.*Box2Abb[3638]*m_z12*m_z2k - 3.*Box2Abb[3433]*m_z2k_2;

  Box2Abb[3640]=26. - 27.*m_z12;

  Box2Abb[3641]=-2. + Box2Abb[3640]*m_z12;

  Box2Abb[3642]=-50. + m_z12;

  Box2Abb[3643]=29. + Box2Abb[3642]*m_z12;

  Box2Abb[3644]=-34. + Box2Abb[3643]*m_z12;

  Box2Abb[3645]=-49. + 32.*m_z12 + 6.*m_z12_2;

  Box2Abb[3646]=6.*pow(Box2Abb[4],2.) + Box2Abb[3641]*m_z2k + Box2Abb[3644]*m_z2k_2 - 2.*Box2Abb[3645]*m_z2k_3 - 15.*Box2Abb[3452]*m_z2k_4;

  Box2Abb[3647]=3. + 17.*m_z2k;

  Box2Abb[3648]=-5. + Box2Abb[3647]*m_z2k;

  Box2Abb[3649]=pow(Box2Abb[61],2.)*Box2Abb[725] + Box2Abb[3648]*Box2Abb[61]*m_z12 + Box2Abb[510]*m_z12_2*m_z2k;

  Box2Abb[3650]=-21. + 40.*m_z2k;

  Box2Abb[3651]=-7. + Box2Abb[3650]*m_z2k;

  Box2Abb[3652]=91. + 25.*m_z2k;

  Box2Abb[3653]=37. + Box2Abb[3652]*m_z2k;

  Box2Abb[3654]=21. + Box2Abb[3653]*m_z2k;

  Box2Abb[3655]=-13. + Box2Abb[3654]*m_z12 - 2.*Box2Abb[2499]*m_z12_2 + 2.*Box2Abb[3651]*m_z2k;

  Box2Abb[3656]=-1. + 7.*Box2Abb[12]*m_z2k;

  Box2Abb[3657]=-1. + Box2Abb[230]*Box2Abb[3656]*m_z2k;

  Box2Abb[3658]=-29. + 24.*m_z2k;

  Box2Abb[3659]=4. + Box2Abb[3658]*m_z2k;

  Box2Abb[3660]=-2. + Box2Abb[3659]*m_z2k;

  Box2Abb[3661]=20. + 39.*m_z2k;

  Box2Abb[3662]=-35. + Box2Abb[3661]*m_z2k;

  Box2Abb[3663]=6. + Box2Abb[3662]*m_z2k;

  Box2Abb[3664]=-4. + Box2Abb[3663]*m_z2k;

  Box2Abb[3665]=Box2Abb[3660]*pow(Box2Abb[61],2.) + Box2Abb[3664]*Box2Abb[61]*m_z12 + 2.*Box2Abb[3657]*m_z12_2 + Box2Abb[2432]*m_z12_3*m_z2k_2;

  Box2Abb[3666]=Box2Abb[3665]*m_x_2 + Box2Abb[3646]*m_x_3 + Box2Abb[3655]*m_x_4 + Box2Abb[3639]*m_x_5 + Box2Abb[3637]*m_x_6 + Box2Abb[201]*m_x_7 - Box2Abb[3649]*Box2Abb[5]*Box2Abb[61]*m_x*m_z2k + 3.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[3667]=Box2Abb[3666]*m_cL + Box2Abb[3635]*m_cR;

  Box2Abb[3668]=9. + 5.*Box2Abb[192]*m_z12;

  Box2Abb[3669]=5. + Box2Abb[3668]*m_z12;

  Box2Abb[3670]=-20. + m_z12;

  Box2Abb[3671]=24. + Box2Abb[3670]*m_z12;

  Box2Abb[3672]=-22. + Box2Abb[3671]*m_z12;

  Box2Abb[3673]=5. + Box2Abb[3672]*m_z12;

  Box2Abb[3674]=5. + 3.*m_z12;

  Box2Abb[3675]=-9. + Box2Abb[3674]*m_z12;

  Box2Abb[3676]=2. - 9.*m_z12;

  Box2Abb[3677]=-17. + 27.*m_z12 - 11.*m_z12_2 - 2.*Box2Abb[3669]*m_z2k + Box2Abb[3673]*m_z2k_2 - 4.*Box2Abb[3675]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[3676]*m_z12*m_z2k_4;

  Box2Abb[3678]=-4. + Box2Abb[1403]*m_z12 + 10.*m_z2k;

  Box2Abb[3679]=6. + Box2Abb[3678]*m_z12;

  Box2Abb[3680]=14. + Box2Abb[872]*m_z12;

  Box2Abb[3681]=31. - 14.*m_z12 + 2.*Box2Abb[3680]*m_z2k - 3.*Box2Abb[669]*m_z2k_2;

  Box2Abb[3682]=-23. + Box2Abb[3681]*m_z12 - 26.*m_z2k;

  Box2Abb[3683]=-5. + 6.*Box2Abb[61]*m_z2k;

  Box2Abb[3684]=2. + Box2Abb[1359]*m_z2k;

  Box2Abb[3685]=1. + Box2Abb[3684]*m_z2k;

  Box2Abb[3686]=3.*Box2Abb[3685]*m_z12 + Box2Abb[3683]*m_z12_2 + Box2Abb[150]*m_z12_3 - 2.*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[3687]=35. + 44.*m_z2k;

  Box2Abb[3688]=17. + 25.*m_z2k;

  Box2Abb[3689]=20. + Box2Abb[12]*Box2Abb[3688]*m_z2k;

  Box2Abb[3690]=-32. + 5.*m_z2k;

  Box2Abb[3691]=-24. + Box2Abb[3690]*m_z2k;

  Box2Abb[3692]=-47. + 2.*Box2Abb[3691]*m_z2k;

  Box2Abb[3693]=31. + Box2Abb[3692]*m_z12 + Box2Abb[3689]*m_z12_2 + Box2Abb[3687]*m_z2k - 2.*Box2Abb[2534]*m_z12_3*m_z2k;

  Box2Abb[3694]=-6. - 11.*m_z2k + 10.*m_z2k_3;

  Box2Abb[3695]=-18. + 11.*m_z2k;

  Box2Abb[3696]=12. + Box2Abb[3695]*m_z2k;

  Box2Abb[3697]=-11. + Box2Abb[3696]*m_z2k;

  Box2Abb[3698]=4. + 2.*Box2Abb[3697]*m_z2k;

  Box2Abb[3699]=-43. + 17.*m_z2k;

  Box2Abb[3700]=7. + Box2Abb[3699]*m_z2k;

  Box2Abb[3701]=-9. + Box2Abb[3700]*m_z2k;

  Box2Abb[3702]=10. + Box2Abb[3701]*m_z2k;

  Box2Abb[3703]=-Box2Abb[3694]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3702]*Box2Abb[61]*m_z12_2 + Box2Abb[3698]*m_z12_3 + 2.*pow(Box2Abb[61],3.)*m_z2k + Box2Abb[3548]*m_z12_4*m_z2k;

  Box2Abb[3704]=-9. + 14.*m_z2k;

  Box2Abb[3705]=-3. + Box2Abb[3704]*m_z2k_2;

  Box2Abb[3706]=-7. + 4.*Box2Abb[1071]*m_z2k;

  Box2Abb[3707]=22. + Box2Abb[3706]*m_z2k;

  Box2Abb[3708]=-97. + 39.*m_z2k;

  Box2Abb[3709]=87. + Box2Abb[3708]*m_z2k;

  Box2Abb[3710]=-56. + Box2Abb[3709]*m_z2k;

  Box2Abb[3711]=2. + Box2Abb[3710]*m_z2k;

  Box2Abb[3712]=-26. + Box2Abb[3117]*m_z2k;

  Box2Abb[3713]=39. + 2.*Box2Abb[3712]*m_z2k;

  Box2Abb[3714]=-38. + Box2Abb[3713]*m_z2k;

  Box2Abb[3715]=5. + Box2Abb[3714]*m_z2k;

  Box2Abb[3716]=Box2Abb[3705]*Box2Abb[61] - Box2Abb[3715]*m_z12 + Box2Abb[12]*Box2Abb[3711]*m_z12_2 + Box2Abb[3707]*m_z12_3*m_z2k + Box2Abb[184]*m_z12_4*m_z2k_2;

  Box2Abb[3717]=Box2Abb[3716]*m_x_2 + Box2Abb[3677]*m_x_3 + Box2Abb[3693]*m_x_4 + Box2Abb[3682]*m_x_5 + Box2Abb[3679]*m_x_6 + Box2Abb[72]*m_x_7*m_z12 - Box2Abb[3703]*m_x*m_z2k + Box2Abb[3686]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[3718]=10. + Box2Abb[192]*m_z12;

  Box2Abb[3719]=-35. + 9.*m_z12;

  Box2Abb[3720]=9. + Box2Abb[3719]*m_z12;

  Box2Abb[3721]=26. - 5.*m_z12;

  Box2Abb[3722]=-60. + Box2Abb[3721]*m_z12;

  Box2Abb[3723]=14. + Box2Abb[3722]*m_z12;

  Box2Abb[3724]=-11. + Box2Abb[3723]*m_z12;

  Box2Abb[3725]=-10. + 17.*m_z12;

  Box2Abb[3726]=9. + Box2Abb[3725]*m_z12;

  Box2Abb[3727]=2. + 19.*m_z12;

  Box2Abb[3728]=13. - 2.*Box2Abb[3718]*m_z12 + 10.*m_z2k + 2.*Box2Abb[3720]*m_z12*m_z2k + Box2Abb[3724]*m_z2k_2 - 4.*Box2Abb[3726]*Box2Abb[4]*m_z2k_3 - 5.*Box2Abb[3727]*m_z12*m_z2k_4;

  Box2Abb[3729]=2.*Box2Abb[61] + Box2Abb[389]*m_z12;

  Box2Abb[3730]=9. + 5.*m_z12;

  Box2Abb[3731]=-2. + 9.*m_z12;

  Box2Abb[3732]=34. + Box2Abb[1058]*m_z12 + 46.*m_z2k + 2.*Box2Abb[3730]*m_z12*m_z2k + 9.*Box2Abb[3731]*m_z2k_2;

  Box2Abb[3733]=-23. + Box2Abb[3732]*m_z12 - 26.*m_z2k;

  Box2Abb[3734]=8. + Box2Abb[1587]*m_z12 - 10.*m_z2k;

  Box2Abb[3735]=-6. + Box2Abb[3734]*m_z12;

  Box2Abb[3736]=1. + Box2Abb[450]*m_z2k;

  Box2Abb[3737]=-27. + 11.*m_z2k;

  Box2Abb[3738]=20. + Box2Abb[3737]*m_z2k;

  Box2Abb[3739]=-4. + Box2Abb[3738]*m_z2k_2;

  Box2Abb[3740]=-2.*pow(Box2Abb[61],3.) + 2.*Box2Abb[3736]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3739]*m_z12_2 + 2.*Box2Abb[1583]*m_z12_3*m_z2k + Box2Abb[1580]*m_z12_4*m_z2k;

  Box2Abb[3741]=47. - 5.*m_z2k;

  Box2Abb[3742]=27. + Box2Abb[3741]*m_z2k;

  Box2Abb[3743]=37. + 2.*Box2Abb[3742]*m_z2k;

  Box2Abb[3744]=-32. + 115.*m_z2k;

  Box2Abb[3745]=20. + Box2Abb[3744]*m_z2k;

  Box2Abb[3746]=-5. + Box2Abb[3745]*m_z2k;

  Box2Abb[3747]=-26. + Box2Abb[3743]*m_z12 + Box2Abb[3746]*m_z12_2 + Box2Abb[1590]*m_z12_3 - 11.*Box2Abb[234]*m_z2k;

  Box2Abb[3748]=-1. + 14.*m_z2k;

  Box2Abb[3749]=4. + Box2Abb[3748]*m_z2k;

  Box2Abb[3750]=11. + 9.*m_z2k;

  Box2Abb[3751]=-6. + Box2Abb[3750]*m_z2k;

  Box2Abb[3752]=29. + 2.*Box2Abb[3751]*m_z2k;

  Box2Abb[3753]=-7. + Box2Abb[3752]*m_z2k;

  Box2Abb[3754]=-34. + 15.*m_z2k;

  Box2Abb[3755]=86. + 3.*Box2Abb[3754]*m_z2k;

  Box2Abb[3756]=-37. + Box2Abb[3755]*m_z2k;

  Box2Abb[3757]=45. + Box2Abb[3756]*m_z2k;

  Box2Abb[3758]=-5. + Box2Abb[3757]*m_z2k;

  Box2Abb[3759]=-Box2Abb[3749]*pow(Box2Abb[61],2.) + Box2Abb[3753]*Box2Abb[61]*m_z12 + Box2Abb[3758]*m_z12_2 + pow(Box2Abb[1613],2.)*Box2Abb[1615]*m_z12_3 + Box2Abb[3520]*m_z12_4*m_z2k_2;

  Box2Abb[3760]=Box2Abb[3759]*m_x_2 + Box2Abb[3728]*m_x_3 + Box2Abb[3747]*m_x_4 - Box2Abb[3733]*m_x_5 + Box2Abb[3735]*m_x_6 + Box2Abb[879]*m_x_7*m_z12 - Box2Abb[3740]*Box2Abb[61]*m_x*m_z2k + Box2Abb[3729]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[3761]=Box2Abb[3717]*m_cL + Box2Abb[3760]*m_cR;

  Box2Abb[3762]=-1. + Box2Abb[1527]*Box2Abb[693]*m_z12;

  Box2Abb[3763]=63. + 16.*m_z12;

  Box2Abb[3764]=-18. + Box2Abb[3763]*m_z12;

  Box2Abb[3765]=-31. + 21.*m_z12;

  Box2Abb[3766]=62. + Box2Abb[3765]*m_z12;

  Box2Abb[3767]=-22. + Box2Abb[3766]*m_z12;

  Box2Abb[3768]=-2. + 23.*m_z12;

  Box2Abb[3769]=-3. + Box2Abb[3762]*m_z12 - 19.*m_z2k + Box2Abb[3764]*m_z12*m_z2k + 2.*Box2Abb[3767]*m_z2k_2 + 5.*Box2Abb[3768]*m_z12*m_z2k_3;

  Box2Abb[3770]=17. + 6.*m_z12;

  Box2Abb[3771]=-25. + Box2Abb[3770]*m_z12;

  Box2Abb[3772]=6. + Box2Abb[3771]*m_z12;

  Box2Abb[3773]=3. + 4.*m_z12;

  Box2Abb[3774]=20. + Box2Abb[3773]*m_z12;

  Box2Abb[3775]=-19. + Box2Abb[3774]*m_z12;

  Box2Abb[3776]=172. + 5.*Box2Abb[3402]*m_z12;

  Box2Abb[3777]=-130. + Box2Abb[3776]*m_z12;

  Box2Abb[3778]=29. + Box2Abb[3777]*m_z12;

  Box2Abb[3779]=-15. + 17.*m_z12;

  Box2Abb[3780]=9. + Box2Abb[3779]*m_z12;

  Box2Abb[3781]=-1. + Box2Abb[3772]*m_z12 + 8.*m_z2k + 2.*Box2Abb[3775]*m_z12*m_z2k + Box2Abb[3778]*m_z2k_2 + 4.*Box2Abb[3780]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[3727]*m_z12*m_z2k_4;

  Box2Abb[3782]=6. + 31.*m_z2k;

  Box2Abb[3783]=12. + Box2Abb[3782]*m_z12 - 10.*m_z2k;

  Box2Abb[3784]=-6. + Box2Abb[3783]*m_z12;

  Box2Abb[3785]=32. - 9.*m_z2k;

  Box2Abb[3786]=6. + 81.*m_z2k_2;

  Box2Abb[3787]=19. + Box2Abb[3786]*m_z12 + 2.*Box2Abb[2396]*m_z12_2 + 2.*Box2Abb[3785]*m_z2k;

  Box2Abb[3788]=-19. + Box2Abb[3787]*m_z12 - 26.*m_z2k;

  Box2Abb[3789]=-3. + 2.*Box2Abb[1151]*m_z2k;

  Box2Abb[3790]=16. + Box2Abb[3520]*m_z2k;

  Box2Abb[3791]=-3. + Box2Abb[3790]*m_z2k;

  Box2Abb[3792]=-2.*pow(Box2Abb[61],2.) + Box2Abb[3789]*Box2Abb[61]*m_z12 + Box2Abb[3791]*m_z12_2 + Box2Abb[3620]*m_z12_3;

  Box2Abb[3793]=3. + 14.*m_z2k;

  Box2Abb[3794]=81. - 88.*m_z2k + 52.*m_z2k_2;

  Box2Abb[3795]=-5. + Box2Abb[3794]*m_z2k_2;

  Box2Abb[3796]=20. + 9.*m_z2k;

  Box2Abb[3797]=-19. + Box2Abb[3796]*m_z2k;

  Box2Abb[3798]=-5. + 2.*Box2Abb[3797]*m_z2k;

  Box2Abb[3799]=-19. + 15.*m_z2k;

  Box2Abb[3800]=39. + Box2Abb[3799]*m_z2k;

  Box2Abb[3801]=-10. + Box2Abb[3800]*m_z2k;

  Box2Abb[3802]=-1. + 3.*Box2Abb[3801]*m_z2k;

  Box2Abb[3803]=-Box2Abb[3793]*pow(Box2Abb[61],3.) + Box2Abb[3798]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3802]*Box2Abb[61]*m_z12_2 + Box2Abb[3795]*m_z12_3 + Box2Abb[3630]*m_z12_4;

  Box2Abb[3804]=-Box2Abb[3792]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x + Box2Abb[3803]*m_x_2 - Box2Abb[3781]*m_x_3 + Box2Abb[3769]*m_x_4 - Box2Abb[3788]*m_x_5 + Box2Abb[3784]*m_x_6 + Box2Abb[879]*m_x_7*m_z12 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*Box2Abb[70]*m_z12*m_z2k;

  Box2Abb[3805]=-31. + 15.*m_z12;

  Box2Abb[3806]=19. + Box2Abb[3805]*m_z12;

  Box2Abb[3807]=76. + Box2Abb[3642]*m_z12;

  Box2Abb[3808]=-62. + Box2Abb[3807]*m_z12;

  Box2Abb[3809]=23. + Box2Abb[3808]*m_z12;

  Box2Abb[3810]=-9. + Box2Abb[2287]*m_z12;

  Box2Abb[3811]=-1. + 2.*Box2Abb[166]*Box2Abb[4]*m_z12 - 2.*m_z2k - 2.*Box2Abb[3806]*m_z12*m_z2k + Box2Abb[3809]*m_z2k_2 - 4.*Box2Abb[3810]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[3676]*m_z12*m_z2k_4;

  Box2Abb[3812]=2. + 3.*Box2Abb[389]*m_z12 - 2.*m_z2k;

  Box2Abb[3813]=-8. - 3.*Box2Abb[179]*m_z12 + 10.*m_z2k;

  Box2Abb[3814]=6. + Box2Abb[3813]*m_z12;

  Box2Abb[3815]=-19. + 2.*m_z12;

  Box2Abb[3816]=-19. + m_z12;

  Box2Abb[3817]=30. + Box2Abb[3815]*m_z12 + 46.*m_z2k + 2.*Box2Abb[3816]*m_z12*m_z2k - 3.*Box2Abb[669]*m_z2k_2;

  Box2Abb[3818]=-19. + Box2Abb[3817]*m_z12 - 26.*m_z2k;

  Box2Abb[3819]=-13. + 11.*m_z2k;

  Box2Abb[3820]=-3. + Box2Abb[3819]*m_z2k;

  Box2Abb[3821]=4. + 2.*Box2Abb[3820]*m_z2k;

  Box2Abb[3822]=-32. + 17.*m_z2k;

  Box2Abb[3823]=-8. + Box2Abb[3822]*m_z2k;

  Box2Abb[3824]=12. + Box2Abb[3823]*m_z2k;

  Box2Abb[3825]=2.*pow(Box2Abb[61],3.) - 2.*Box2Abb[1159]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3824]*Box2Abb[61]*m_z12_2 + Box2Abb[3821]*m_z12_3 + Box2Abb[510]*m_z12_4*m_z2k;

  Box2Abb[3826]=21. + 44.*m_z2k;

  Box2Abb[3827]=19. - 2.*m_z2k;

  Box2Abb[3828]=-6. + Box2Abb[3827]*m_z2k;

  Box2Abb[3829]=-15. + Box2Abb[1667]*m_z2k;

  Box2Abb[3830]=-25. + 2.*Box2Abb[3829]*m_z2k;

  Box2Abb[3831]=72. + 25.*m_z2k;

  Box2Abb[3832]=4. + Box2Abb[3831]*m_z2k;

  Box2Abb[3833]=21. + Box2Abb[3832]*m_z2k;

  Box2Abb[3834]=14. + Box2Abb[3830]*m_z12 + Box2Abb[3833]*m_z12_2 + Box2Abb[3828]*m_z12_3 + Box2Abb[3826]*m_z2k;

  Box2Abb[3835]=-5. + 14.*m_z2k;

  Box2Abb[3836]=-37. + 4.*Box2Abb[3116]*m_z2k;

  Box2Abb[3837]=5. + Box2Abb[3836]*m_z2k;

  Box2Abb[3838]=-2. + Box2Abb[3837]*m_z2k;

  Box2Abb[3839]=-26. + Box2Abb[3750]*m_z2k;

  Box2Abb[3840]=17. + 2.*Box2Abb[3839]*m_z2k;

  Box2Abb[3841]=-3. + Box2Abb[3840]*m_z2k;

  Box2Abb[3842]=-19. + 39.*m_z2k;

  Box2Abb[3843]=-57. + Box2Abb[3842]*m_z2k;

  Box2Abb[3844]=12. + Box2Abb[3843]*m_z2k;

  Box2Abb[3845]=-5. + Box2Abb[3844]*m_z2k;

  Box2Abb[3846]=-Box2Abb[3841]*Box2Abb[61]*m_z12 + Box2Abb[3845]*Box2Abb[61]*m_z12_2 + Box2Abb[3838]*m_z12_3 + Box2Abb[3835]*pow(Box2Abb[61],2.)*m_z2k + Box2Abb[2432]*m_z12_4*m_z2k_2;

  Box2Abb[3847]=Box2Abb[3846]*m_x_2 + Box2Abb[3811]*m_x_3 + Box2Abb[3834]*m_x_4 + Box2Abb[3818]*m_x_5 + Box2Abb[3814]*m_x_6 + Box2Abb[72]*m_x_7*m_z12 - Box2Abb[3825]*Box2Abb[61]*m_x*m_z2k + Box2Abb[3812]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[3848]=Box2Abb[3847]*m_cL + Box2Abb[3804]*m_cR;

  Box2Abb[3849]=Box2Abb[431]*m_cR*m_x - Box2Abb[429]*Box2Abb[8]*m_cL*m_z2k;

  Box2Abb[3850]=18. - 11.*m_z12;

  Box2Abb[3851]=2. + Box2Abb[4]*m_z12;

  Box2Abb[3852]=-3. + Box2Abb[3351]*m_z12;

  Box2Abb[3853]=-9. + Box2Abb[3852]*m_z12;

  Box2Abb[3854]=-34. + 29.*m_z12;

  Box2Abb[3855]=13. + Box2Abb[3854]*m_z12;

  Box2Abb[3856]=-4. + 9.*m_z12;

  Box2Abb[3857]=-7. + Box2Abb[3850]*m_z12 + 10.*m_z2k - 8.*Box2Abb[3851]*m_z12*m_z2k + Box2Abb[3853]*m_z2k_2 + 2.*Box2Abb[3855]*m_z2k_3 + 5.*Box2Abb[3856]*m_z2k_4;

  Box2Abb[3858]=13. - 4.*Box2Abb[3134]*m_z12 + 5.*m_z12_2 + 20.*m_z2k;

  Box2Abb[3859]=13. + 20.*m_z2k;

  Box2Abb[3860]=26. + 75.*m_z2k;

  Box2Abb[3861]=21. + Box2Abb[3860]*m_z2k;

  Box2Abb[3862]=-11. + Box2Abb[3861]*m_z12 + Box2Abb[2088]*m_z12_2 - 2.*Box2Abb[3859]*m_z2k;

  Box2Abb[3863]=5. + 44.*m_z2k;

  Box2Abb[3864]=21. - Box2Abb[3863]*m_z2k;

  Box2Abb[3865]=3. - 10.*m_z2k;

  Box2Abb[3866]=11. + 8.*Box2Abb[3865]*m_z2k;

  Box2Abb[3867]=-26. + Box2Abb[3866]*m_z2k;

  Box2Abb[3868]=7. + Box2Abb[3867]*m_z12 + Box2Abb[3864]*m_z12_2 + 7.*m_z2k + 3.*m_z12_3*m_z2k + 40.*m_z2k_3;

  Box2Abb[3869]=-11. + m_z2k + 33.*m_z2k_2 - 25.*m_z2k_3;

  Box2Abb[3870]=2. + Box2Abb[3869]*m_z2k;

  Box2Abb[3871]=-1. + 4.*Box2Abb[2800]*m_z2k;

  Box2Abb[3872]=4. + Box2Abb[3871]*m_z2k;

  Box2Abb[3873]=Box2Abb[3482]*pow(Box2Abb[61],3.) - Box2Abb[3872]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3870]*m_z12_2 - Box2Abb[12]*Box2Abb[2130]*m_z12_3*m_z2k;

  Box2Abb[3874]=Box2Abb[3873]*m_x + Box2Abb[3857]*m_x_2 + Box2Abb[3868]*m_x_3 + Box2Abb[3862]*m_x_4 + Box2Abb[3858]*m_x_5 + Box2Abb[272]*m_x_6 + Box2Abb[179]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[3875]=-1. + 3.*m_z12 + m_z2k;

  Box2Abb[3876]=-3. + 2.*m_z12 + 3.*m_z2k;

  Box2Abb[3877]=4. + 11.*m_z12;

  Box2Abb[3878]=-13. + Box2Abb[3877]*m_z12;

  Box2Abb[3879]=-14. + 9.*m_z12 + Box2Abb[3878]*m_z2k + 6.*Box2Abb[201]*m_z2k_2;

  Box2Abb[3880]=13. - 9.*Box2Abb[12]*m_z12 + 16.*m_z2k;

  Box2Abb[3881]=-20. + Box2Abb[2101]*m_z2k;

  Box2Abb[3882]=11. - Box2Abb[2196]*m_z2k;

  Box2Abb[3883]=-pow(Box2Abb[61],2.)*Box2Abb[725] - Box2Abb[3881]*Box2Abb[61]*m_z12 + Box2Abb[3882]*m_z12_2 + 6.*m_z12_3*m_z2k;

  Box2Abb[3884]=-13. + 16.*m_z2k;

  Box2Abb[3885]=25. + 2.*Box2Abb[488]*m_z2k;

  Box2Abb[3886]=-3. + Box2Abb[3885]*m_z2k;

  Box2Abb[3887]=5. + Box2Abb[3886]*m_z12 - 11.*Box2Abb[187]*m_z12_2*m_z2k + Box2Abb[3884]*m_z2k_2;

  Box2Abb[3888]=Box2Abb[3887]*m_x_2 + Box2Abb[3879]*m_x_3 + Box2Abb[3880]*m_x_4 + Box2Abb[166]*m_x_5 + Box2Abb[3883]*m_x*m_z2k + Box2Abb[3875]*Box2Abb[3876]*Box2Abb[61]*m_z12*m_z2k_2;

  Box2Abb[3889]=-Box2Abb[3888]*Box2Abb[8]*m_cL + Box2Abb[3874]*m_cR;

  Box2Abb[3890]=-3. + 5.*m_z12 + 3.*m_z2k;

  Box2Abb[3891]=21. + Box2Abb[2038]*m_z12;

  Box2Abb[3892]=6. + Box2Abb[1449]*m_z12;

  Box2Abb[3893]=-3. + Box2Abb[3856]*m_z12 + m_z2k + Box2Abb[3891]*m_z12*m_z2k - 4.*Box2Abb[3892]*m_z2k_2 + 40.*m_z2k_3;

  Box2Abb[3894]=-4. + Box2Abb[2038]*m_z12;

  Box2Abb[3895]=19. - 5.*m_z12;

  Box2Abb[3896]=-11. + Box2Abb[3895]*m_z12;

  Box2Abb[3897]=27. + Box2Abb[3896]*m_z12;

  Box2Abb[3898]=21. + Box2Abb[1516]*m_z12;

  Box2Abb[3899]=1. + Box2Abb[166]*m_z12 - 6.*m_z2k + 2.*Box2Abb[3894]*m_z12*m_z2k + Box2Abb[3897]*m_z2k_2 - 2.*Box2Abb[3898]*m_z2k_3 + 5.*Box2Abb[3452]*m_z2k_4;

  Box2Abb[3900]=9. - 2.*m_z2k;

  Box2Abb[3901]=1. + m_z12 + Box2Abb[3900]*m_z12_2 + 10.*Box2Abb[183]*m_z2k + 3.*Box2Abb[738]*m_z12*m_z2k;

  Box2Abb[3902]=2. - 3.*m_z12 + 12.*m_z2k;

  Box2Abb[3903]=-9. + Box2Abb[3902]*m_z12 - 20.*m_z2k;

  Box2Abb[3904]=1. - Box2Abb[2651]*m_z2k;

  Box2Abb[3905]=-13. + 6.*Box2Abb[444]*m_z2k;

  Box2Abb[3906]=-7. + Box2Abb[763]*m_z2k;

  Box2Abb[3907]=pow(Box2Abb[61],3.)*Box2Abb[725] + Box2Abb[3905]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[3906]*Box2Abb[61]*m_z12_2 + Box2Abb[3904]*m_z12_3;

  Box2Abb[3908]=Box2Abb[3899]*m_x_2 - Box2Abb[3893]*m_x_3 + Box2Abb[3901]*m_x_4 + Box2Abb[3903]*m_x_5 + Box2Abb[505]*m_x_6 - Box2Abb[3907]*m_x*m_z2k + Box2Abb[3890]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[3909]=m_z12_2 - 4.*m_z2k + 7.*m_z12*m_z2k;

  Box2Abb[3910]=-4. + 7.*m_z2k;

  Box2Abb[3911]=9. + 16.*m_z2k;

  Box2Abb[3912]=7. + 29.*m_z2k;

  Box2Abb[3913]=2. - Box2Abb[3912]*m_z2k;

  Box2Abb[3914]=Box2Abb[3913]*m_z12 + Box2Abb[3910]*m_z12_2 + Box2Abb[3911]*m_z2k;

  Box2Abb[3915]=20. + 3.*m_z2k;

  Box2Abb[3916]=-6. + Box2Abb[3915]*m_z2k;

  Box2Abb[3917]=8. - 23.*m_z2k;

  Box2Abb[3918]=-9. + Box2Abb[3917]*m_z2k;

  Box2Abb[3919]=3. + Box2Abb[3918]*m_z2k;

  Box2Abb[3920]=-1. + 2.*Box2Abb[3919]*m_z12 + Box2Abb[3916]*m_z12_2 + m_z12_3*m_z2k + m_z2k_2 + 24.*m_z2k_3;

  Box2Abb[3921]=-4. + 17.*m_z2k - 21.*m_z2k_3;

  Box2Abb[3922]=7. + 16.*m_z2k;

  Box2Abb[3923]=-2. + Box2Abb[3922]*m_z2k;

  Box2Abb[3924]=-16. + m_z2k + 29.*m_z2k_2 - 17.*m_z2k_3;

  Box2Abb[3925]=3. + Box2Abb[3924]*m_z2k;

  Box2Abb[3926]=Box2Abb[3923]*pow(Box2Abb[61],2.) + 2.*Box2Abb[3925]*m_z12 + Box2Abb[3921]*m_z12_2 + Box2Abb[1580]*m_z12_3*m_z2k;

  Box2Abb[3927]=Box2Abb[3499]*Box2Abb[5]*Box2Abb[61]*m_x + Box2Abb[3926]*m_x_2 - Box2Abb[3920]*m_x_3 + Box2Abb[3914]*m_x_4 + Box2Abb[3909]*m_x_5 - pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[3928]=Box2Abb[3927]*Box2Abb[8]*m_cR + Box2Abb[3908]*m_cL*m_z2k;

  Box2Abb[3929]=3.*m_z12_2 - 4.*m_z2k + 7.*m_z12*m_z2k;

  Box2Abb[3930]=5. - 4.*m_z12;

  Box2Abb[3931]=3. + 5.*m_z12;

  Box2Abb[3932]=12. + Box2Abb[3931]*m_z12;

  Box2Abb[3933]=-47. + 19.*m_z12;

  Box2Abb[3934]=12. + Box2Abb[3933]*m_z12;

  Box2Abb[3935]=-8. + 6.*Box2Abb[3930]*m_z12 + 6.*m_z2k - Box2Abb[3932]*m_z12*m_z2k + 2.*Box2Abb[242]*Box2Abb[4]*Box2Abb[9]*m_z2k_2 - 2.*Box2Abb[3934]*m_z2k_3 + 40.*Box2Abb[1424]*m_z2k_4;

  Box2Abb[3936]=5. + 36.*m_z2k;

  Box2Abb[3937]=6. - Box2Abb[3936]*m_z2k;

  Box2Abb[3938]=Box2Abb[3937]*m_z12 - 2.*Box2Abb[2152]*m_z12_2 + Box2Abb[3502]*m_z2k;

  Box2Abb[3939]=13. + 5.*Box2Abb[12]*m_z2k;

  Box2Abb[3940]=-1. + 5.*Box2Abb[183]*m_z2k;

  Box2Abb[3941]=-19. + 75.*m_z2k;

  Box2Abb[3942]=-22. + Box2Abb[3941]*m_z2k_2;

  Box2Abb[3943]=3. + Box2Abb[3942]*m_z12 + 2.*Box2Abb[3939]*m_z12_2 - 2.*Box2Abb[3940]*m_z2k + m_z12_3*m_z2k;

  Box2Abb[3944]=-1. + 10.*m_z2k;

  Box2Abb[3945]=-7. + Box2Abb[3944]*m_z2k;

  Box2Abb[3946]=-2. + 3.*Box2Abb[725]*m_z2k;

  Box2Abb[3947]=2. + Box2Abb[3946]*m_z2k;

  Box2Abb[3948]=-Box2Abb[3482]*pow(Box2Abb[61],2.) + Box2Abb[3947]*Box2Abb[61]*m_z12 + Box2Abb[3945]*m_z12_2*m_z2k;

  Box2Abb[3949]=9. + Box2Abb[1065]*m_z2k;

  Box2Abb[3950]=-1. + Box2Abb[3944]*m_z2k;

  Box2Abb[3951]=-7. + 2.*Box2Abb[3950]*m_z2k;

  Box2Abb[3952]=-73. + 45.*m_z2k;

  Box2Abb[3953]=-5. + Box2Abb[3952]*m_z2k;

  Box2Abb[3954]=-3. + Box2Abb[3953]*m_z2k;

  Box2Abb[3955]=18. + Box2Abb[3954]*m_z2k;

  Box2Abb[3956]=-50. + 49.*m_z2k;

  Box2Abb[3957]=-6. + Box2Abb[3956]*m_z2k;

  Box2Abb[3958]=-12. + Box2Abb[3957]*m_z2k;

  Box2Abb[3959]=11. + Box2Abb[3958]*m_z2k;

  Box2Abb[3960]=-Box2Abb[3951]*pow(Box2Abb[61],2.) + Box2Abb[3955]*Box2Abb[61]*m_z12 + Box2Abb[3959]*m_z12_2 + Box2Abb[3949]*m_z12_3*m_z2k;

  Box2Abb[3961]=-Box2Abb[3948]*Box2Abb[5]*Box2Abb[61]*m_x + Box2Abb[3960]*m_x_2 + Box2Abb[3935]*m_x_3 + Box2Abb[3943]*m_x_4 + Box2Abb[3938]*m_x_5 + Box2Abb[3929]*m_x_6 + Box2Abb[179]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[3962]=Box2Abb[3961]*m_cR + Box2Abb[3908]*m_cL*m_z2k;

  Box2Abb[3963]=m_z12_2 - 4.*m_z2k + 3.*m_z12*m_z2k;

  Box2Abb[3964]=-8. + 4.*m_z12 + m_z12_3;

  Box2Abb[3965]=37. + 2.*m_z12;

  Box2Abb[3966]=-6. + Box2Abb[3965]*m_z12;

  Box2Abb[3967]=8. - 3.*m_z12;

  Box2Abb[3968]=-1. + 6.*m_z12 - 6.*m_z12_2 + Box2Abb[3964]*m_z2k + Box2Abb[3966]*m_z2k_2 + 5.*Box2Abb[3967]*m_z2k_3;

  Box2Abb[3969]=-13. + 3.*m_z12;

  Box2Abb[3970]=2. + Box2Abb[3969]*m_z12;

  Box2Abb[3971]=2. + Box2Abb[3970]*m_z12;

  Box2Abb[3972]=-1. + Box2Abb[242]*m_z12;

  Box2Abb[3973]=-33. + 5.*m_z12;

  Box2Abb[3974]=24. + Box2Abb[3973]*m_z12;

  Box2Abb[3975]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box2Abb[3971]*m_z2k + 4.*Box2Abb[3972]*m_z12*m_z2k_2 - 2.*Box2Abb[3974]*m_z2k_3 + 40.*m_z2k_4;

  Box2Abb[3976]=2. + 3.*Box2Abb[815]*m_z12;

  Box2Abb[3977]=10. + Box2Abb[2828]*m_z12;

  Box2Abb[3978]=11. + Box2Abb[3977]*m_z12;

  Box2Abb[3979]=25. + 2.*Box2Abb[2044]*m_z12;

  Box2Abb[3980]=-26. + Box2Abb[3979]*m_z12;

  Box2Abb[3981]=58. + Box2Abb[1723]*m_z12;

  Box2Abb[3982]=pow(Box2Abb[4],2.) - Box2Abb[3976]*Box2Abb[4]*m_z2k + Box2Abb[3978]*m_z2k_2 + 2.*Box2Abb[3980]*m_z2k_3 + Box2Abb[3981]*m_z2k_4 - 5.*Box2Abb[3452]*m_z2k_5;

  Box2Abb[3983]=5. - 12.*m_z2k;

  Box2Abb[3984]=2. + Box2Abb[3983]*m_z2k;

  Box2Abb[3985]=Box2Abb[3984]*m_z12 + 2.*Box2Abb[179]*m_z12_2 + 5.*Box2Abb[183]*m_z2k;

  Box2Abb[3986]=-11. + 4.*m_z2k;

  Box2Abb[3987]=1. + Box2Abb[3986]*m_z2k;

  Box2Abb[3988]=5. + m_z2k - 4.*m_z2k_2;

  Box2Abb[3989]=-2. + Box2Abb[3988]*m_z2k;

  Box2Abb[3990]=-pow(Box2Abb[61],2.)*Box2Abb[725] + 3.*Box2Abb[3989]*m_z12 + Box2Abb[3987]*m_z12_2;

  Box2Abb[3991]=Box2Abb[3982]*m_x_2 + Box2Abb[3975]*m_x_3 - Box2Abb[3968]*m_x_4 + Box2Abb[3985]*m_x_5 + Box2Abb[3963]*m_x_6 - Box2Abb[3990]*Box2Abb[5]*Box2Abb[61]*m_x*m_z2k - 3.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[3992]=13. + 6.*m_z12;

  Box2Abb[3993]=-18. + Box2Abb[3992]*m_z12;

  Box2Abb[3994]=33. + 10.*m_z12;

  Box2Abb[3995]=-16. + Box2Abb[3994]*m_z12;

  Box2Abb[3996]=9. + Box2Abb[3995]*m_z12;

  Box2Abb[3997]=1. - 9.*m_z12 + 6.*m_z12_2;

  Box2Abb[3998]=3. + 2.*Box2Abb[710]*m_z12 + 8.*m_z2k + 2.*Box2Abb[3993]*m_z12*m_z2k + Box2Abb[3996]*m_z2k_2 + 6.*Box2Abb[3997]*m_z2k_3 + 5.*Box2Abb[3473]*m_z2k_4;

  Box2Abb[3999]=5. + 12.*m_z2k;

  Box2Abb[4000]=2. - 7.*m_z2k_2;

  Box2Abb[4001]=Box2Abb[4000]*m_z12 - Box2Abb[3999]*m_z12_2 + 4.*m_z2k_2;

  Box2Abb[4002]=3. - 2.*m_z2k + 8.*m_z2k_2;

  Box2Abb[4003]=1. + Box2Abb[514]*m_z2k;

  Box2Abb[4004]=-pow(Box2Abb[61],2.)*Box2Abb[831] + 4.*Box2Abb[4003]*Box2Abb[61]*m_z12 + Box2Abb[4002]*m_z12_2;

  Box2Abb[4005]=10. + 7.*Box2Abb[2804]*m_z2k;

  Box2Abb[4006]=-5. + Box2Abb[2304]*m_z2k;

  Box2Abb[4007]=-2. + Box2Abb[4006]*m_z2k;

  Box2Abb[4008]=1. + 4.*Box2Abb[4007]*m_z12 + Box2Abb[4005]*m_z12_2 + 3.*m_z12_3*m_z2k - 5.*Box2Abb[183]*m_z2k_2;

  Box2Abb[4009]=12. - 9.*m_z2k + 5.*m_z2k_3;

  Box2Abb[4010]=-1. + 2.*m_z2k + 20.*m_z2k_2;

  Box2Abb[4011]=1. + Box2Abb[4010]*m_z2k;

  Box2Abb[4012]=-17. + 23.*m_z2k;

  Box2Abb[4013]=-7. + Box2Abb[4012]*m_z2k;

  Box2Abb[4014]=21. + 2.*Box2Abb[4013]*m_z2k;

  Box2Abb[4015]=-1. + Box2Abb[4014]*m_z2k;

  Box2Abb[4016]=-58. + 45.*m_z2k;

  Box2Abb[4017]=-6. + Box2Abb[4016]*m_z2k;

  Box2Abb[4018]=10. + Box2Abb[4017]*m_z2k;

  Box2Abb[4019]=-2. + Box2Abb[4018]*m_z2k;

  Box2Abb[4020]=-Box2Abb[4011]*pow(Box2Abb[61],3.) + Box2Abb[4019]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4015]*Box2Abb[61]*m_z12_2 + Box2Abb[4009]*m_z12_3*m_z2k;

  Box2Abb[4021]=18. + 5.*Box2Abb[187]*m_z2k;

  Box2Abb[4022]=-17. + m_z2k + 48.*m_z2k_2 - 40.*m_z2k_3;

  Box2Abb[4023]=5. + Box2Abb[4022]*m_z2k;

  Box2Abb[4024]=-17. + 47.*m_z2k;

  Box2Abb[4025]=7. + Box2Abb[4024]*m_z2k;

  Box2Abb[4026]=-12. + Box2Abb[4025]*m_z2k;

  Box2Abb[4027]=5. + Box2Abb[4026]*m_z2k;

  Box2Abb[4028]=47. + 16.*Box2Abb[1133]*m_z2k;

  Box2Abb[4029]=8. + Box2Abb[4028]*m_z2k;

  Box2Abb[4030]=-7. + Box2Abb[4029]*m_z2k;

  Box2Abb[4031]=-8. + Box2Abb[4030]*m_z2k;

  Box2Abb[4032]=3. + Box2Abb[4031]*m_z12 + Box2Abb[4027]*m_z12_2 + Box2Abb[4023]*m_z2k + Box2Abb[4021]*m_z12_3*m_z2k;

  Box2Abb[4033]=-Box2Abb[4020]*m_x_2 + Box2Abb[4032]*m_x_3 - Box2Abb[3998]*m_x_4 + Box2Abb[4008]*m_x_5 + Box2Abb[4001]*m_x_6 + m_x_7*m_z12_2 + Box2Abb[4004]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z2k - pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12*m_z2k_2;

  Box2Abb[4034]=Box2Abb[4033]*m_cR + Box2Abb[3991]*m_cL*m_z2k;

  Box2Abb[4035]=2.*pow(Box2Abb[61],2.) + 9.*Box2Abb[61]*m_z12 + 5.*m_z12_2;

  Box2Abb[4036]=12. + Box2Abb[1186]*m_z12;

  Box2Abb[4037]=16. - 9.*m_z12_2 + 6.*m_z2k - pow(Box2Abb[242],2.)*m_z12*m_z2k + 2.*Box2Abb[4036]*m_z2k_2;

  Box2Abb[4038]=15. + 2.*m_z12;

  Box2Abb[4039]=-16. + Box2Abb[4038]*m_z12;

  Box2Abb[4040]=18. - 5.*m_z12;

  Box2Abb[4041]=-11. + Box2Abb[4040]*m_z12;

  Box2Abb[4042]=20. + Box2Abb[4041]*m_z12;

  Box2Abb[4043]=6. + Box2Abb[158]*m_z12;

  Box2Abb[4044]=-2. + 3.*Box2Abb[4]*m_z12 - 6.*m_z2k + Box2Abb[4039]*m_z12*m_z2k + Box2Abb[4042]*m_z2k_2 - 2.*Box2Abb[4043]*m_z2k_3 + 10.*m_z12*m_z2k_4;

  Box2Abb[4045]=-2. - 3.*m_z12 + 8.*m_z2k;

  Box2Abb[4046]=6. + Box2Abb[4045]*m_z12;

  Box2Abb[4047]=6. - 5.*m_z2k;

  Box2Abb[4048]=7. - 3.*Box2Abb[841]*m_z12 + 2.*Box2Abb[4047]*m_z2k;

  Box2Abb[4049]=-20.*Box2Abb[12] + Box2Abb[4048]*m_z12;

  Box2Abb[4050]=-1. + Box2Abb[2651]*m_z2k;

  Box2Abb[4051]=-13. - 2.*m_z2k + 8.*m_z2k_2;

  Box2Abb[4052]=2.*pow(Box2Abb[61],3.) - Box2Abb[4051]*pow(Box2Abb[61],2.)*m_z12 - 3.*Box2Abb[1069]*Box2Abb[12]*Box2Abb[61]*m_z12_2 + Box2Abb[4050]*m_z12_3;

  Box2Abb[4053]=Box2Abb[4044]*m_x_2 + Box2Abb[4037]*m_x_3 + Box2Abb[4049]*m_x_4 + Box2Abb[4046]*m_x_5 - 2.*m_x_6*m_z12 + Box2Abb[4052]*m_x*m_z2k + Box2Abb[4035]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[4054]=m_z12 + Box2Abb[2145]*m_z12_2 - 3.*m_z2k + Box2Abb[2311]*m_z12*m_z2k;

  Box2Abb[4055]=1. + 4.*Box2Abb[460]*m_z2k;

  Box2Abb[4056]=-1. + Box2Abb[2804]*m_z2k;

  Box2Abb[4057]=2.*pow(Box2Abb[61],3.) - 2.*Box2Abb[4056]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4055]*Box2Abb[61]*m_z12_2 + Box2Abb[2149]*m_z12_3;

  Box2Abb[4058]=24. + 19.*m_z2k;

  Box2Abb[4059]=-6. + Box2Abb[4058]*m_z2k;

  Box2Abb[4060]=29. - 5.*m_z2k;

  Box2Abb[4061]=8. + Box2Abb[4060]*m_z2k;

  Box2Abb[4062]=3. + Box2Abb[4061]*m_z2k;

  Box2Abb[4063]=-1. + 2.*Box2Abb[4062]*m_z12 + Box2Abb[4059]*m_z12_2 - 4.*Box2Abb[2396]*m_z2k + m_z12_3*m_z2k;

  Box2Abb[4064]=23. + 36.*m_z2k;

  Box2Abb[4065]=-4. + Box2Abb[4064]*m_z2k;

  Box2Abb[4066]=5. + 6.*Box2Abb[444]*m_z2k;

  Box2Abb[4067]=1. + Box2Abb[4066]*m_z2k;

  Box2Abb[4068]=-1. + 40.*m_z2k;

  Box2Abb[4069]=-3. + Box2Abb[4068]*m_z2k;

  Box2Abb[4070]=6. + 2.*Box2Abb[4069]*m_z2k;

  Box2Abb[4071]=-2.*Box2Abb[4067] + Box2Abb[4070]*m_z12 + Box2Abb[4065]*m_z12_2 + Box2Abb[2153]*m_z12_3*m_z2k;

  Box2Abb[4072]=1. + 4.*Box2Abb[256]*m_z2k;

  Box2Abb[4073]=21. + 5.*m_z2k;

  Box2Abb[4074]=-1. + Box2Abb[4073]*m_z2k;

  Box2Abb[4075]=-1. + m_z2k + Box2Abb[4074]*m_z2k_2;

  Box2Abb[4076]=-40. + 19.*m_z2k;

  Box2Abb[4077]=-6. + Box2Abb[4076]*m_z2k;

  Box2Abb[4078]=-6. + Box2Abb[4077]*m_z2k;

  Box2Abb[4079]=1. + Box2Abb[4078]*m_z2k;

  Box2Abb[4080]=Box2Abb[4072]*pow(Box2Abb[61],2.) - 2.*Box2Abb[4075]*Box2Abb[61]*m_z12 + Box2Abb[4079]*m_z12_2 - Box2Abb[1580]*m_z12_3*m_z2k;

  Box2Abb[4081]=Box2Abb[4080]*m_x_2 + Box2Abb[4071]*m_x_3 - Box2Abb[4063]*m_x_4 + 2.*Box2Abb[4054]*m_x_5 + Box2Abb[56]*m_x_6*m_z12 - Box2Abb[4057]*Box2Abb[61]*m_x*m_z2k + Box2Abb[3336]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[4082]=Box2Abb[4081]*m_cR + Box2Abb[4053]*m_cL*m_z2k;

  Box2Abb[4083]=-1. + 2.*m_z12 - 4.*m_z2k;

  Box2Abb[4084]=Box2Abb[4083]*Box2Abb[61]*m_z12 - 3.*m_z2k;

  Box2Abb[4085]=-6. + m_z12_2;

  Box2Abb[4086]=-11. + 2.*Box2Abb[242]*m_z12;

  Box2Abb[4087]=2. + Box2Abb[4086]*m_z12;

  Box2Abb[4088]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box2Abb[1485]*Box2Abb[4085]*m_z2k + 2.*Box2Abb[4087]*m_z2k_2 - 8.*Box2Abb[4]*Box2Abb[538]*m_z2k_3;

  Box2Abb[4089]=-4. + Box2Abb[734]*m_z12;

  Box2Abb[4090]=26. + 7.*m_z12;

  Box2Abb[4091]=20. - Box2Abb[4090]*m_z12;

  Box2Abb[4092]=1. + 6.*Box2Abb[4]*m_z12 + 16.*m_z2k - Box2Abb[4089]*m_z12*m_z2k + Box2Abb[4091]*m_z2k_2 + 10.*m_z12*m_z2k_3;

  Box2Abb[4093]=11. - 4.*m_z2k;

  Box2Abb[4094]=-1. + Box2Abb[4093]*m_z2k;

  Box2Abb[4095]=-3. + m_z2k + 4.*m_z2k_2;

  Box2Abb[4096]=-5. + 4.*Box2Abb[222]*m_z2k;

  Box2Abb[4097]=-2.*pow(Box2Abb[61],3.) + 2.*Box2Abb[4095]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4096]*Box2Abb[61]*m_z12_2 + Box2Abb[4094]*m_z12_3;

  Box2Abb[4098]=1. + 4.*Box2Abb[470]*m_z2k;

  Box2Abb[4099]=-3. + Box2Abb[2804]*m_z2k;

  Box2Abb[4100]=6. - 12.*m_z2k + 5.*m_z2k_3;

  Box2Abb[4101]=1. + Box2Abb[4100]*m_z2k_2;

  Box2Abb[4102]=-24. + 11.*m_z2k;

  Box2Abb[4103]=2. + Box2Abb[4102]*m_z2k;

  Box2Abb[4104]=10. + Box2Abb[4103]*m_z2k;

  Box2Abb[4105]=1. + Box2Abb[4104]*m_z2k;

  Box2Abb[4106]=Box2Abb[4098]*pow(Box2Abb[61],2.) - 2.*Box2Abb[4101]*m_z12 + Box2Abb[4105]*m_z12_2 + Box2Abb[4099]*m_z12_3*m_z2k;

  Box2Abb[4107]=Box2Abb[4106]*m_x_2 + Box2Abb[4088]*m_x_3 + Box2Abb[4092]*m_x_4 + 2.*Box2Abb[4084]*m_x_5 + Box2Abb[56]*m_x_6*m_z12 + Box2Abb[4097]*Box2Abb[61]*m_x*m_z2k - Box2Abb[2175]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[4108]=-2. + 5.*m_z12 + 14.*m_z12*m_z2k + 2.*m_z2k_2;

  Box2Abb[4109]=-8. + Box2Abb[1527]*m_z12;

  Box2Abb[4110]=-18. + 37.*m_z12;

  Box2Abb[4111]=6. + Box2Abb[4110]*m_z12;

  Box2Abb[4112]=1. + 2.*Box2Abb[1324]*m_z12 + 3.*Box2Abb[4109]*m_z12*m_z2k + Box2Abb[4111]*m_z2k_2 + 10.*m_z12*m_z2k_3;

  Box2Abb[4113]=-1. + 6.*m_z12_2;

  Box2Abb[4114]=-14. + Box2Abb[505]*m_z12;

  Box2Abb[4115]=3. + Box2Abb[4114]*m_z12;

  Box2Abb[4116]=5. + 6.*m_z12 + 3.*m_z12_3;

  Box2Abb[4117]=121. + 8.*m_z12;

  Box2Abb[4118]=-200. + Box2Abb[4117]*m_z12;

  Box2Abb[4119]=26. + Box2Abb[4118]*m_z12;

  Box2Abb[4120]=22. + 5.*m_z12;

  Box2Abb[4121]=-6. + Box2Abb[4120]*m_z12;

  Box2Abb[4122]=-pow(Box2Abb[4],2.) - 2.*Box2Abb[4]*Box2Abb[4113]*m_z2k + Box2Abb[4115]*m_z2k_2 + 2.*Box2Abb[4116]*m_z2k_3 + Box2Abb[4119]*m_z2k_4 + 6.*Box2Abb[4121]*m_z2k_5 + 10.*m_z12*m_z2k_6;

  Box2Abb[4123]=1. - 3.*m_z2k + 2.*m_z2k_3;

  Box2Abb[4124]=-3. + 2.*Box2Abb[2152]*m_z2k;

  Box2Abb[4125]=-2.*pow(Box2Abb[61],2.) + 5.*Box2Abb[4123]*m_z12 + Box2Abb[4124]*m_z12_2;

  Box2Abb[4126]=10. + 13.*m_z2k;

  Box2Abb[4127]=5. + Box2Abb[4126]*m_z2k;

  Box2Abb[4128]=23. + 5.*Box2Abb[488]*m_z2k;

  Box2Abb[4129]=5. + Box2Abb[4128]*m_z2k;

  Box2Abb[4130]=-41. + 18.*Box2Abb[2534]*m_z2k;

  Box2Abb[4131]=-52. + Box2Abb[4130]*m_z2k;

  Box2Abb[4132]=-12. + Box2Abb[4131]*m_z2k;

  Box2Abb[4133]=3. + Box2Abb[4132]*m_z12 + 2.*Box2Abb[4129]*m_z12_2 + 2.*Box2Abb[4127]*m_z2k + 12.*Box2Abb[12]*m_z12_3*m_z2k;

  Box2Abb[4134]=-1. + 2.*Box2Abb[3116]*m_z2k;

  Box2Abb[4135]=26. + 9.*m_z2k;

  Box2Abb[4136]=8. + Box2Abb[4135]*m_z2k;

  Box2Abb[4137]=-3. + Box2Abb[4136]*m_z2k_2;

  Box2Abb[4138]=-25. + 9.*Box2Abb[222]*m_z2k;

  Box2Abb[4139]=-3. + Box2Abb[4138]*m_z2k;

  Box2Abb[4140]=-1. + 2.*Box2Abb[4139]*m_z2k;

  Box2Abb[4141]=86. + 25.*m_z2k;

  Box2Abb[4142]=-30. + Box2Abb[4141]*m_z2k;

  Box2Abb[4143]=-2. + Box2Abb[4142]*m_z2k;

  Box2Abb[4144]=-5. + Box2Abb[4143]*m_z2k;

  Box2Abb[4145]=Box2Abb[4134]*pow(Box2Abb[61],3.) - Box2Abb[4140]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[4144]*Box2Abb[61]*m_z12_2 - Box2Abb[4137]*m_z12_3;

  Box2Abb[4146]=9. + 4.*m_z2k;

  Box2Abb[4147]=9. + Box2Abb[4146]*m_z2k;

  Box2Abb[4148]=-5. + 11.*Box2Abb[444]*m_z2k;

  Box2Abb[4149]=11. + 2.*Box2Abb[4148]*m_z2k;

  Box2Abb[4150]=52. + 9.*m_z2k;

  Box2Abb[4151]=2. + Box2Abb[4150]*m_z2k;

  Box2Abb[4152]=-8. + Box2Abb[4151]*m_z2k;

  Box2Abb[4153]=5. - Box2Abb[4152]*m_z2k;

  Box2Abb[4154]=-84. + 5.*m_z2k;

  Box2Abb[4155]=29. + Box2Abb[4154]*m_z2k;

  Box2Abb[4156]=2. + Box2Abb[4155]*m_z2k;

  Box2Abb[4157]=-31. + 2.*Box2Abb[4156]*m_z2k;

  Box2Abb[4158]=-8. + Box2Abb[4157]*m_z2k;

  Box2Abb[4159]=3. + Box2Abb[4158]*m_z12 + Box2Abb[4153]*m_z12_2 + Box2Abb[4149]*m_z2k + 2.*Box2Abb[4147]*m_z12_3*m_z2k;

  Box2Abb[4160]=Box2Abb[4122]*m_x_3 + Box2Abb[4159]*m_x_4 - Box2Abb[4133]*m_x_5 + Box2Abb[4112]*m_x_6 - Box2Abb[4108]*m_x_7*m_z12 + m_x_8*m_z12_2 + Box2Abb[4145]*m_x_2*m_z2k + Box2Abb[4125]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z2k_2 - 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12*m_z2k_3;

  Box2Abb[4161]=Box2Abb[4160]*m_cR + Box2Abb[4107]*Box2Abb[8]*m_cL*m_z2k;

  Box2Abb[4162]=-4. + 7.*m_z12 + 6.*m_z12*m_z2k + 4.*Box2Abb[4]*m_z2k_2;

  Box2Abb[4163]=-2. + m_z2k + 8.*m_z2k_2;

  Box2Abb[4164]=2.*pow(Box2Abb[61],2.) - Box2Abb[4163]*Box2Abb[61]*m_z12 - 2.*Box2Abb[187]*m_z12_2*m_z2k;

  Box2Abb[4165]=Box2Abb[4164]*m_x + Box2Abb[4162]*m_x_2 - Box2Abb[783]*m_x_3 + 2.*Box2Abb[70]*m_x_4 + 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[4166]=7. + 2.*m_z12;

  Box2Abb[4167]=-6. + Box2Abb[4166]*m_z12;

  Box2Abb[4168]=2. + Box2Abb[1916]*m_z2k;

  Box2Abb[4169]=2.*Box2Abb[4168]*m_z12 - 3.*m_z12_2 + 4.*Box2Abb[841]*m_z2k;

  Box2Abb[4170]=1. - Box2Abb[2804]*m_z2k;

  Box2Abb[4171]=5. - 8.*m_z2k;

  Box2Abb[4172]=6. + Box2Abb[4171]*m_z2k;

  Box2Abb[4173]=-3. + Box2Abb[4172]*m_z2k;

  Box2Abb[4174]=2.*pow(Box2Abb[61],2.) + Box2Abb[4173]*m_z12 + Box2Abb[4170]*m_z12_2;

  Box2Abb[4175]=Box2Abb[4174]*m_x + Box2Abb[4169]*m_x_2 + Box2Abb[4167]*m_x_3 - 2.*Box2Abb[158]*m_x_4 + 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[4176]=Box2Abb[4165]*m_cL + Box2Abb[4175]*m_cR;

  Box2Abb[4177]=-6. + 13.*m_z12;

  Box2Abb[4178]=-7. - 2.*m_z12 + 4.*m_z12_2;

  Box2Abb[4179]=75. + 8.*m_z12;

  Box2Abb[4180]=-49. + Box2Abb[4179]*m_z12;

  Box2Abb[4181]=-2. + m_z12 - 11.*m_z12_2;

  Box2Abb[4182]=12. + 17.*Box2Abb[4181]*m_z12;

  Box2Abb[4183]=12. + 2.*Box2Abb[4178]*m_z12 + 8.*m_z2k - 2.*Box2Abb[4180]*m_z12_2*m_z2k + Box2Abb[4182]*m_z2k_2 - 5.*Box2Abb[734]*m_z12*m_z2k_3;

  Box2Abb[4184]=44. + 75.*m_z12;

  Box2Abb[4185]=-10. + 69.*m_z12;

  Box2Abb[4186]=-10. + 2.*Box2Abb[3432]*m_z12 - 13.*m_z2k + Box2Abb[4184]*m_z12*m_z2k + Box2Abb[4185]*m_z2k_2;

  Box2Abb[4187]=-8.*Box2Abb[12] + Box2Abb[4186]*m_z12;

  Box2Abb[4188]=2. + Box2Abb[2534]*m_z2k;

  Box2Abb[4189]=1. + Box2Abb[230]*m_z2k;

  Box2Abb[4190]=Box2Abb[4188]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4189]*Box2Abb[61]*m_z12 + 2.*Box2Abb[187]*m_z12_2*m_z2k;

  Box2Abb[4191]=33. + 2.*m_z12 + 53.*m_z2k;

  Box2Abb[4192]=19. - Box2Abb[4191]*m_z12 + 20.*m_z2k;

  Box2Abb[4193]=2. + Box2Abb[4192]*m_z12;

  Box2Abb[4194]=17. + 41.*m_z2k;

  Box2Abb[4195]=35. + 9.*Box2Abb[1567]*m_z2k;

  Box2Abb[4196]=-8. + Box2Abb[4195]*m_z2k;

  Box2Abb[4197]=-89. + 65.*m_z2k;

  Box2Abb[4198]=135. + Box2Abb[4197]*m_z2k;

  Box2Abb[4199]=-13. + Box2Abb[4198]*m_z2k;

  Box2Abb[4200]=-2. + Box2Abb[424]*m_z2k;

  Box2Abb[4201]=7. + 5.*Box2Abb[4200]*m_z2k;

  Box2Abb[4202]=5. + Box2Abb[4201]*m_z2k;

  Box2Abb[4203]=-8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4202]*m_z12 - Box2Abb[12]*Box2Abb[4199]*m_z12_2 + 2.*Box2Abb[4196]*m_z12_3 + 2.*Box2Abb[4194]*m_z12_4*m_z2k;

  Box2Abb[4204]=41. + 31.*m_z2k;

  Box2Abb[4205]=9. + Box2Abb[4204]*m_z2k;

  Box2Abb[4206]=9. + 44.*m_z2k;

  Box2Abb[4207]=6. + Box2Abb[4206]*m_z2k;

  Box2Abb[4208]=-3. + Box2Abb[4207]*m_z2k;

  Box2Abb[4209]=10. + Box2Abb[1065]*m_z2k;

  Box2Abb[4210]=1. + Box2Abb[4209]*m_z2k;

  Box2Abb[4211]=5. + 2.*Box2Abb[4210]*m_z2k;

  Box2Abb[4212]=-5. + 57.*m_z2k;

  Box2Abb[4213]=10. + Box2Abb[179]*Box2Abb[4212]*m_z2k;

  Box2Abb[4214]=37. + Box2Abb[4213]*m_z2k;

  Box2Abb[4215]=-15. + Box2Abb[4214]*m_z2k;

  Box2Abb[4216]=2.*pow(Box2Abb[61],4.) - Box2Abb[4208]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4215]*m_z12_2 + 2.*Box2Abb[4211]*m_z12_3 - 2.*Box2Abb[4205]*m_z12_4*m_z2k - 10.*m_z12_5*m_z2k_2;

  Box2Abb[4217]=-7. + Box2Abb[559]*m_z2k;

  Box2Abb[4218]=1. + Box2Abb[4217]*m_z2k;

  Box2Abb[4219]=-42. + 17.*m_z2k;

  Box2Abb[4220]=11. + Box2Abb[4219]*m_z2k;

  Box2Abb[4221]=2. + Box2Abb[4220]*m_z2k;

  Box2Abb[4222]=-4. + Box2Abb[4221]*m_z2k;

  Box2Abb[4223]=66. - 41.*m_z2k;

  Box2Abb[4224]=2. + Box2Abb[4223]*m_z2k;

  Box2Abb[4225]=-30. + Box2Abb[4224]*m_z2k;

  Box2Abb[4226]=5. + Box2Abb[4225]*m_z2k;

  Box2Abb[4227]=-2. + Box2Abb[4226]*m_z2k;

  Box2Abb[4228]=Box2Abb[230]*Box2Abb[450]*pow(Box2Abb[61],4.) - Box2Abb[4222]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4227]*m_z12_2 - 2.*Box2Abb[4218]*m_z12_3*m_z2k + 2.*Box2Abb[1941]*m_z12_4*m_z2k_2;

  Box2Abb[4229]=Box2Abb[4216]*m_x_2 + Box2Abb[4203]*m_x_3 + Box2Abb[4183]*m_x_4 + Box2Abb[4187]*m_x_5 + Box2Abb[4193]*m_x_6 + Box2Abb[4228]*m_x*m_z12 + Box2Abb[4177]*m_x_7*m_z12 + Box2Abb[4190]*Box2Abb[5]*Box2Abb[61]*m_z12_2*m_z2k;

  Box2Abb[4230]=17. + 2.*m_z12;

  Box2Abb[4231]=22. - Box2Abb[4230]*m_z12;

  Box2Abb[4232]=4. + Box2Abb[4231]*m_z12;

  Box2Abb[4233]=11. + 18.*m_z12;

  Box2Abb[4234]=-19. + Box2Abb[4233]*m_z12;

  Box2Abb[4235]=6. + Box2Abb[4234]*m_z12;

  Box2Abb[4236]=12. - 2.*Box2Abb[1629]*Box2Abb[4]*m_z12 + 8.*m_z2k + Box2Abb[4232]*m_z12*m_z2k + 2.*Box2Abb[4235]*m_z2k_2 + 5.*Box2Abb[1186]*m_z12*m_z2k_3;

  Box2Abb[4237]=6. + m_z12 + 3.*m_z2k;

  Box2Abb[4238]=-17. + Box2Abb[4237]*m_z12 - 20.*m_z2k;

  Box2Abb[4239]=-2. + Box2Abb[4238]*m_z12;

  Box2Abb[4240]=1. + Box2Abb[1069]*m_z2k;

  Box2Abb[4241]=-2. + Box2Abb[3117]*m_z2k;

  Box2Abb[4242]=7. + 34.*m_z2k;

  Box2Abb[4243]=4. + Box2Abb[4242]*Box2Abb[61]*m_z2k;

  Box2Abb[4244]=-1. + 4.*Box2Abb[61]*m_z2k;

  Box2Abb[4245]=1. + 2.*Box2Abb[4244]*m_z2k;

  Box2Abb[4246]=-Box2Abb[4241]*pow(Box2Abb[61],4.) - Box2Abb[4243]*pow(Box2Abb[61],2.)*m_z12 - 2.*Box2Abb[4245]*Box2Abb[61]*m_z12_2 + Box2Abb[4240]*m_z12_3*m_z2k;

  Box2Abb[4247]=7. + 10.*m_z2k;

  Box2Abb[4248]=4. + 7.*Box2Abb[150]*m_z2k;

  Box2Abb[4249]=-9. - Box2Abb[4248]*m_z12 + Box2Abb[285]*m_z12_2 - Box2Abb[4247]*m_z2k;

  Box2Abb[4250]=-8.*Box2Abb[12] + Box2Abb[4249]*m_z12;

  Box2Abb[4251]=2. + m_z2k + 39.*m_z2k_2 + 133.*m_z2k_3;

  Box2Abb[4252]=1. - 5.*m_z2k;

  Box2Abb[4253]=11. + 7.*Box2Abb[4252]*m_z2k;

  Box2Abb[4254]=21. + Box2Abb[4253]*m_z2k;

  Box2Abb[4255]=-8. + Box2Abb[4254]*m_z2k;

  Box2Abb[4256]=-44. + 15.*m_z2k;

  Box2Abb[4257]=-106. + Box2Abb[4256]*m_z2k;

  Box2Abb[4258]=-8. + Box2Abb[4257]*m_z2k;

  Box2Abb[4259]=9. + Box2Abb[4258]*m_z2k;

  Box2Abb[4260]=8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4255]*m_z12 + Box2Abb[4259]*m_z12_2 + Box2Abb[4251]*m_z12_3 + 4.*Box2Abb[61]*m_z12_4*m_z2k;

  Box2Abb[4261]=1. - 7.*m_z2k;

  Box2Abb[4262]=2. - 7.*m_z2k;

  Box2Abb[4263]=2. + 3.*Box2Abb[4262]*m_z2k;

  Box2Abb[4264]=-2. + Box2Abb[4263]*m_z2k;

  Box2Abb[4265]=-31. + 24.*m_z2k;

  Box2Abb[4266]=14. + Box2Abb[4265]*m_z2k;

  Box2Abb[4267]=4. + Box2Abb[4266]*m_z2k;

  Box2Abb[4268]=-1. + Box2Abb[4267]*m_z2k;

  Box2Abb[4269]=-46. + 33.*m_z2k;

  Box2Abb[4270]=-46. + Box2Abb[4269]*m_z2k;

  Box2Abb[4271]=108. + Box2Abb[4270]*m_z2k;

  Box2Abb[4272]=-51. + Box2Abb[4271]*m_z2k;

  Box2Abb[4273]=2. + Box2Abb[4272]*m_z2k_2;

  Box2Abb[4274]=Box2Abb[230]*Box2Abb[510]*pow(Box2Abb[61],4.) + Box2Abb[4273]*m_z12 + Box2Abb[12]*Box2Abb[4268]*m_z12_2 + 2.*Box2Abb[4264]*m_z12_3*m_z2k + Box2Abb[4261]*m_z12_4*m_z2k_2;

  Box2Abb[4275]=2. + 29.*m_z2k;

  Box2Abb[4276]=22. - 95.*m_z2k;

  Box2Abb[4277]=49. + Box2Abb[4276]*m_z2k;

  Box2Abb[4278]=-9. + Box2Abb[4277]*m_z2k;

  Box2Abb[4279]=-3. + Box2Abb[4278]*m_z2k;

  Box2Abb[4280]=-41. + 44.*m_z2k;

  Box2Abb[4281]=3. + Box2Abb[4280]*m_z2k;

  Box2Abb[4282]=21. + Box2Abb[4281]*m_z2k;

  Box2Abb[4283]=-3. + Box2Abb[4282]*m_z2k;

  Box2Abb[4284]=34. + 33.*m_z2k;

  Box2Abb[4285]=34. + Box2Abb[4284]*m_z2k;

  Box2Abb[4286]=-44. + Box2Abb[4285]*m_z2k;

  Box2Abb[4287]=15. + Box2Abb[4286]*m_z2k;

  Box2Abb[4288]=2. + Box2Abb[4287]*m_z2k;

  Box2Abb[4289]=-2.*pow(Box2Abb[61],4.) + Box2Abb[4283]*Box2Abb[61]*m_z12 + Box2Abb[4288]*m_z12_2 + Box2Abb[4279]*m_z12_3 - 2.*Box2Abb[4275]*m_z12_4*m_z2k_2 + m_z12_5*m_z2k_2;

  Box2Abb[4290]=-Box2Abb[4289]*m_x_3 - Box2Abb[4260]*m_x_4 + Box2Abb[4236]*m_x_5 + Box2Abb[4250]*m_x_6 - Box2Abb[4239]*m_x_7 + Box2Abb[4274]*m_x_2*m_z12 + 3.*Box2Abb[72]*m_x_8*m_z12 + Box2Abb[4246]*m_x*m_z12_2*m_z2k + Box2Abb[470]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],2.)*m_z12_3*m_z2k_2;

  Box2Abb[4291]=Box2Abb[4290]*m_cL - Box2Abb[4229]*m_cR*m_x;

  Box2Abb[4292]=-3. + 4.*m_z12 + 3.*m_z2k;

  Box2Abb[4293]=43. - 18.*Box2Abb[70]*m_z12;

  Box2Abb[4294]=8. + Box2Abb[4293]*m_z12;

  Box2Abb[4295]=7. - 4.*Box2Abb[3432]*m_z12;

  Box2Abb[4296]=7. + Box2Abb[4295]*m_z12;

  Box2Abb[4297]=-36. + 37.*m_z12;

  Box2Abb[4298]=-6. + Box2Abb[4297]*m_z12;

  Box2Abb[4299]=-35. + Box2Abb[4298]*m_z12;

  Box2Abb[4300]=4. + Box2Abb[4299]*m_z12;

  Box2Abb[4301]=62. + 33.*m_z12;

  Box2Abb[4302]=5. + Box2Abb[4301]*m_z12;

  Box2Abb[4303]=-4. + Box2Abb[4302]*m_z12;

  Box2Abb[4304]=14. - 39.*m_z12;

  Box2Abb[4305]=8.*Box2Abb[61] + Box2Abb[4294]*m_z12 + 2.*Box2Abb[4296]*m_z12*m_z2k + 2.*Box2Abb[4300]*m_z2k_2 + 2.*Box2Abb[4303]*m_z2k_3 + 5.*Box2Abb[4304]*m_z12*m_z2k_4;

  Box2Abb[4306]=14. + 14.*m_z12 - 15.*m_z2k;

  Box2Abb[4307]=5. + Box2Abb[4306]*m_z12 + 20.*m_z2k;

  Box2Abb[4308]=2. + Box2Abb[4307]*m_z12;

  Box2Abb[4309]=-1. + Box2Abb[260]*m_z2k;

  Box2Abb[4310]=3. + 47.*m_z2k;

  Box2Abb[4311]=-2. + Box2Abb[4310]*Box2Abb[61]*m_z2k;

  Box2Abb[4312]=-5. + 5.*m_z2k + m_z2k_3;

  Box2Abb[4313]=-1. + Box2Abb[4312]*m_z2k;

  Box2Abb[4314]=15. + 7.*Box2Abb[493]*m_z2k;

  Box2Abb[4315]=1. + Box2Abb[4314]*m_z2k;

  Box2Abb[4316]=Box2Abb[4309]*pow(Box2Abb[61],4.) - Box2Abb[4315]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[4311]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[4313]*m_z12_3 + 12.*m_z12_4*m_z2k_3;

  Box2Abb[4317]=58. - 13.*m_z2k;

  Box2Abb[4318]=-37. + 10.*m_z2k;

  Box2Abb[4319]=34. + 11.*m_z2k;

  Box2Abb[4320]=11. + Box2Abb[4319]*m_z2k;

  Box2Abb[4321]=-39. + 3.*Box2Abb[4320]*m_z12 + Box2Abb[4317]*m_z12_2 + 4.*m_z12_3 + Box2Abb[4318]*m_z2k;

  Box2Abb[4322]=8.*Box2Abb[12] + Box2Abb[4321]*m_z12;

  Box2Abb[4323]=-7. + m_z2k;

  Box2Abb[4324]=53. - 97.*m_z2k;

  Box2Abb[4325]=78. + Box2Abb[4324]*m_z2k;

  Box2Abb[4326]=47. + 20.*m_z2k;

  Box2Abb[4327]=18. + Box2Abb[4326]*m_z2k;

  Box2Abb[4328]=23. + Box2Abb[4327]*m_z2k;

  Box2Abb[4329]=122. + 145.*m_z2k;

  Box2Abb[4330]=103. + Box2Abb[4329]*m_z2k;

  Box2Abb[4331]=-20. + Box2Abb[4330]*m_z2k;

  Box2Abb[4332]=4.*Box2Abb[892] - 2.*Box2Abb[4328]*m_z12 + Box2Abb[4331]*m_z12_2 + Box2Abb[4325]*m_z12_3 - 2.*Box2Abb[4323]*m_z12_4;

  Box2Abb[4333]=17. + Box2Abb[3628]*Box2Abb[444]*m_z2k;

  Box2Abb[4334]=5. - 37.*m_z2k;

  Box2Abb[4335]=8. + Box2Abb[4334]*m_z2k;

  Box2Abb[4336]=5. + Box2Abb[4335]*m_z2k;

  Box2Abb[4337]=-1. + 44.*m_z2k;

  Box2Abb[4338]=-26. + Box2Abb[4337]*m_z2k;

  Box2Abb[4339]=-1. + Box2Abb[4338]*m_z2k;

  Box2Abb[4340]=-57. + 41.*m_z2k;

  Box2Abb[4341]=-7. + 3.*Box2Abb[4340]*m_z2k;

  Box2Abb[4342]=57. + Box2Abb[4341]*m_z2k;

  Box2Abb[4343]=12. + Box2Abb[4342]*m_z2k;

  Box2Abb[4344]=2.*pow(Box2Abb[61],4.) - Box2Abb[4339]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4343]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[4336]*m_z12_4 + 2.*Box2Abb[4333]*m_z12_3*m_z2k - 12.*m_z12_5*m_z2k_2;

  Box2Abb[4345]=Box2Abb[4344]*m_x_2 + Box2Abb[4305]*m_x_3 + Box2Abb[4332]*m_x_4 - Box2Abb[4322]*m_x_5 + Box2Abb[4308]*m_x_6 + Box2Abb[4316]*m_x*m_z12 + Box2Abb[3347]*m_x_7*m_z12 + Box2Abb[4292]*Box2Abb[5]*pow(Box2Abb[61],4.)*m_z12_2*m_z2k;

  Box2Abb[4346]=-11. + m_z12 - 20.*m_z2k + 13.*m_z12*m_z2k;

  Box2Abb[4347]=2. - Box2Abb[4346]*m_z12;

  Box2Abb[4348]=-28. + 33.*m_z12;

  Box2Abb[4349]=11. + Box2Abb[4348]*m_z12;

  Box2Abb[4350]=14. - 6.*Box2Abb[203]*m_z12 + Box2Abb[4349]*m_z2k - Box2Abb[2287]*m_z2k_2;

  Box2Abb[4351]=-8.*Box2Abb[12] + Box2Abb[4350]*m_z12;

  Box2Abb[4352]=-2. + m_z2k + 8.*m_z2k_2 - 7.*m_z2k_3;

  Box2Abb[4353]=-pow(Box2Abb[61],3.)*Box2Abb[719] + 2.*Box2Abb[4352]*m_z12 + 2.*Box2Abb[460]*m_z12_2*m_z2k;

  Box2Abb[4354]=4. + 15.*m_z2k;

  Box2Abb[4355]=22. + 35.*m_z2k;

  Box2Abb[4356]=-9. + Box2Abb[4355]*m_z2k_2;

  Box2Abb[4357]=72. - 25.*m_z2k;

  Box2Abb[4358]=46. + Box2Abb[4357]*m_z2k;

  Box2Abb[4359]=11. + Box2Abb[4358]*m_z2k_2;

  Box2Abb[4360]=22. - 9.*m_z2k + 6.*m_z2k_2;

  Box2Abb[4361]=12. + Box2Abb[4360]*m_z2k;

  Box2Abb[4362]=-8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4356]*Box2Abb[61]*m_z12 + Box2Abb[4359]*m_z12_2 - 2.*Box2Abb[4361]*m_z12_3 + 2.*Box2Abb[4354]*m_z12_4*m_z2k;

  Box2Abb[4363]=28. + 57.*m_z2k;

  Box2Abb[4364]=20. - Box2Abb[4363]*m_z2k;

  Box2Abb[4365]=4. + 5.*Box2Abb[2804]*m_z2k;

  Box2Abb[4366]=19. + Box2Abb[4365]*m_z2k;

  Box2Abb[4367]=68. + 5.*Box2Abb[3239]*m_z2k;

  Box2Abb[4368]=22. + Box2Abb[4367]*m_z2k;

  Box2Abb[4369]=4.*Box2Abb[892] - 2.*Box2Abb[4366]*m_z12 + Box2Abb[4368]*m_z12_2 + Box2Abb[4364]*m_z12_3 - 8.*m_z12_4*m_z2k;

  Box2Abb[4370]=-9. + 10.*m_z2k;

  Box2Abb[4371]=-2. + Box2Abb[4370]*m_z2k;

  Box2Abb[4372]=-30. + 13.*m_z2k;

  Box2Abb[4373]=3. + Box2Abb[4372]*m_z2k;

  Box2Abb[4374]=8. + Box2Abb[4373]*m_z2k;

  Box2Abb[4375]=41. + 23.*m_z2k;

  Box2Abb[4376]=-8. + Box2Abb[4375]*m_z2k;

  Box2Abb[4377]=4. + Box2Abb[4376]*Box2Abb[61]*m_z2k;

  Box2Abb[4378]=55. + 9.*m_z2k;

  Box2Abb[4379]=-33. + Box2Abb[4378]*m_z2k;

  Box2Abb[4380]=-15. + Box2Abb[4379]*m_z2k;

  Box2Abb[4381]=2. + Box2Abb[4380]*m_z2k;

  Box2Abb[4382]=-Box2Abb[4371]*pow(Box2Abb[61],4.) - Box2Abb[4377]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[4381]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[4374]*m_z12_3*m_z2k + 2.*m_z12_4*m_z2k_3;

  Box2Abb[4383]=6. + Box2Abb[1069]*m_z2k;

  Box2Abb[4384]=-15. + 44.*m_z2k;

  Box2Abb[4385]=-18. + Box2Abb[4384]*m_z2k;

  Box2Abb[4386]=-3. + Box2Abb[4385]*m_z2k;

  Box2Abb[4387]=-44. + 23.*m_z2k;

  Box2Abb[4388]=-31. + Box2Abb[4387]*m_z2k;

  Box2Abb[4389]=15. + Box2Abb[4388]*m_z2k;

  Box2Abb[4390]=6. + Box2Abb[4389]*m_z2k;

  Box2Abb[4391]=71. + 15.*m_z2k;

  Box2Abb[4392]=-22. + Box2Abb[4391]*m_z2k;

  Box2Abb[4393]=-112. + Box2Abb[4392]*m_z2k;

  Box2Abb[4394]=35. + Box2Abb[4393]*m_z2k;

  Box2Abb[4395]=17. + Box2Abb[4394]*m_z2k;

  Box2Abb[4396]=2.*pow(Box2Abb[61],4.) - Box2Abb[4386]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[4395]*m_z12_2 + 2.*Box2Abb[4390]*m_z12_3 + 2.*Box2Abb[4383]*m_z12_4*m_z2k - 4.*m_z12_5*m_z2k_2;

  Box2Abb[4397]=Box2Abb[4396]*m_x_3 + Box2Abb[4362]*m_x_4 + Box2Abb[4369]*m_x_5 + Box2Abb[4351]*m_x_6 + Box2Abb[4347]*m_x_7 - Box2Abb[4382]*m_x_2*m_z12 + Box2Abb[710]*m_x_8*m_z12 + Box2Abb[4353]*Box2Abb[5]*Box2Abb[61]*m_x*m_z12_2*m_z2k + 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12_3*m_z2k_2;

  Box2Abb[4398]=Box2Abb[4397]*m_cL - Box2Abb[4345]*m_cR*m_x;

  Box2Abb[4399]=6. - 7.*m_z12;

  Box2Abb[4400]=18. + 28.*m_z12 - 15.*m_z2k;

  Box2Abb[4401]=5. + Box2Abb[4400]*m_z12 + 20.*m_z2k;

  Box2Abb[4402]=2. + Box2Abb[4401]*m_z12;

  Box2Abb[4403]=15. + 8.*Box2Abb[1527]*m_z12;

  Box2Abb[4404]=118. + 65.*m_z12;

  Box2Abb[4405]=10. + 33.*m_z12;

  Box2Abb[4406]=-43. + Box2Abb[4403]*m_z12 - 37.*m_z2k + Box2Abb[4404]*m_z12*m_z2k + Box2Abb[4405]*m_z2k_2;

  Box2Abb[4407]=8.*Box2Abb[12] + Box2Abb[4406]*m_z12;

  Box2Abb[4408]=2. - 24.*m_z2k_2 + 22.*m_z2k_3;

  Box2Abb[4409]=-2. + Box2Abb[2791]*m_z2k;

  Box2Abb[4410]=Box2Abb[4409]*pow(Box2Abb[61],2.) + Box2Abb[4408]*m_z12 + 10.*Box2Abb[12]*m_z12_2*m_z2k;

  Box2Abb[4411]=3. + Box2Abb[187]*m_z2k;

  Box2Abb[4412]=-5. + Box2Abb[260]*m_z2k;

  Box2Abb[4413]=-13. + 7.*m_z2k;

  Box2Abb[4414]=-19. + 5.*Box2Abb[4413]*m_z2k;

  Box2Abb[4415]=1. + Box2Abb[4414]*m_z2k;

  Box2Abb[4416]=14. + 41.*m_z2k;

  Box2Abb[4417]=-13. + Box2Abb[4416]*m_z2k;

  Box2Abb[4418]=-7. + Box2Abb[4417]*m_z2k;

  Box2Abb[4419]=-50. + 121.*m_z2k;

  Box2Abb[4420]=-67. + Box2Abb[4419]*m_z2k;

  Box2Abb[4421]=-18. + Box2Abb[4420]*m_z2k;

  Box2Abb[4422]=-Box2Abb[4412]*pow(Box2Abb[61],4.) + Box2Abb[4415]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[4421]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[4418]*Box2Abb[61]*m_z12_3 + 6.*Box2Abb[4411]*m_z12_4*m_z2k;

  Box2Abb[4423]=32. + 34.*m_z2k;

  Box2Abb[4424]=47. + 17.*m_z2k;

  Box2Abb[4425]=176. + 5.*Box2Abb[4424]*m_z2k;

  Box2Abb[4426]=22. + Box2Abb[4326]*m_z2k;

  Box2Abb[4427]=19. + Box2Abb[4426]*m_z2k;

  Box2Abb[4428]=142. + 145.*m_z2k;

  Box2Abb[4429]=9. + Box2Abb[4428]*m_z2k;

  Box2Abb[4430]=-108. + Box2Abb[4429]*m_z2k;

  Box2Abb[4431]=4.*Box2Abb[892] - 2.*Box2Abb[4427]*m_z12 + Box2Abb[4430]*m_z12_2 + Box2Abb[4425]*m_z12_3 + Box2Abb[4423]*m_z12_4;

  Box2Abb[4432]=44. + 27.*m_z2k;

  Box2Abb[4433]=47. + Box2Abb[4432]*m_z2k;

  Box2Abb[4434]=38. + 2.*Box2Abb[4433]*m_z2k;

  Box2Abb[4435]=-34. + Box2Abb[4337]*m_z2k;

  Box2Abb[4436]=-25. + Box2Abb[4435]*m_z2k;

  Box2Abb[4437]=-51. + 107.*m_z2k;

  Box2Abb[4438]=-65. + Box2Abb[4437]*m_z2k;

  Box2Abb[4439]=-31. + Box2Abb[4438]*m_z2k;

  Box2Abb[4440]=-2. + Box2Abb[4439]*m_z2k;

  Box2Abb[4441]=-191. + 123.*m_z2k;

  Box2Abb[4442]=-111. + Box2Abb[4441]*m_z2k;

  Box2Abb[4443]=31. + Box2Abb[4442]*m_z2k;

  Box2Abb[4444]=58. + Box2Abb[4443]*m_z2k;

  Box2Abb[4445]=2.*pow(Box2Abb[61],4.) - Box2Abb[4436]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4444]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[4440]*m_z12_3 + Box2Abb[4434]*m_z12_4 + 2.*Box2Abb[2193]*m_z12_5*m_z2k;

  Box2Abb[4446]=53. + 13.*m_z2k;

  Box2Abb[4447]=25. + Box2Abb[4446]*m_z2k;

  Box2Abb[4448]=52. + 41.*m_z2k;

  Box2Abb[4449]=35. + Box2Abb[4448]*m_z2k;

  Box2Abb[4450]=47. + 2.*Box2Abb[4449]*m_z2k;

  Box2Abb[4451]=-7. + m_z2k + 7.*m_z2k_2;

  Box2Abb[4452]=9. - 5.*Box2Abb[4451]*m_z2k;

  Box2Abb[4453]=12. + Box2Abb[4452]*m_z2k;

  Box2Abb[4454]=-124. + 195.*m_z2k;

  Box2Abb[4455]=-144. + Box2Abb[4454]*m_z2k;

  Box2Abb[4456]=-108. + Box2Abb[4455]*m_z2k;

  Box2Abb[4457]=-153. + Box2Abb[4456]*m_z2k;

  Box2Abb[4458]=8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4453]*m_z12 + Box2Abb[4457]*m_z12_2 + 2.*Box2Abb[4450]*m_z12_3 + 2.*Box2Abb[4447]*m_z12_4 + 4.*m_z12_5*m_z2k;

  Box2Abb[4459]=-Box2Abb[4445]*m_x_3 + Box2Abb[4458]*m_x_4 - Box2Abb[4431]*m_x_5 + Box2Abb[4407]*m_x_6 - Box2Abb[4402]*m_x_7 + Box2Abb[4422]*m_x_2*m_z12 + Box2Abb[4399]*m_x_8*m_z12 - Box2Abb[4410]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12_2 + 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_3*m_z2k;

  Box2Abb[4460]=Box2Abb[4397]*m_cL + Box2Abb[4459]*m_cR;

  Box2Abb[4461]=31. + 7.*m_z12;

  Box2Abb[4462]=4. + Box2Abb[4461]*m_z12;

  Box2Abb[4463]=-31. + Box2Abb[4462]*m_z12;

  Box2Abb[4464]=6. + Box2Abb[4463]*m_z12;

  Box2Abb[4465]=157. + 6.*m_z12;

  Box2Abb[4466]=51. + Box2Abb[4465]*m_z12;

  Box2Abb[4467]=-36. + Box2Abb[4466]*m_z12;

  Box2Abb[4468]=8. + Box2Abb[4467]*m_z12;

  Box2Abb[4469]=-2. + Box2Abb[70]*m_z12;

  Box2Abb[4470]=12. + 55.*Box2Abb[4469]*m_z12;

  Box2Abb[4471]=8. - 37.*m_z12;

  Box2Abb[4472]=-2.*Box2Abb[4464] - Box2Abb[4468]*m_z2k - Box2Abb[4470]*m_z2k_2 + 5.*Box2Abb[4471]*m_z12*m_z2k_3;

  Box2Abb[4473]=11. - 5.*m_z2k;

  Box2Abb[4474]=-7. + 4.*Box2Abb[61]*m_z12 + 8.*m_z12_2 + Box2Abb[4473]*m_z2k;

  Box2Abb[4475]=1. + Box2Abb[4474]*m_z2k;

  Box2Abb[4476]=42. + 14.*m_z12 + 25.*m_z2k;

  Box2Abb[4477]=-3. + Box2Abb[4476]*m_z12 + 20.*m_z2k;

  Box2Abb[4478]=2. + Box2Abb[4477]*m_z12;

  Box2Abb[4479]=-61. + 10.*m_z2k;

  Box2Abb[4480]=158. + 105.*m_z2k;

  Box2Abb[4481]=69. + Box2Abb[4480]*m_z2k;

  Box2Abb[4482]=-59. + Box2Abb[4481]*m_z12 + 27.*Box2Abb[187]*m_z12_2 + 4.*m_z12_3 + Box2Abb[4479]*m_z2k;

  Box2Abb[4483]=8.*Box2Abb[12] + Box2Abb[4482]*m_z12;

  Box2Abb[4484]=3. + 5.*Box2Abb[460]*m_z2k;

  Box2Abb[4485]=49. + 5.*m_z2k;

  Box2Abb[4486]=-27. + Box2Abb[4485]*m_z2k;

  Box2Abb[4487]=7. + Box2Abb[4486]*m_z2k;

  Box2Abb[4488]=1. + Box2Abb[3484]*m_z2k;

  Box2Abb[4489]=-6. + 5.*Box2Abb[4488]*m_z2k;

  Box2Abb[4490]=-42. + 31.*m_z2k;

  Box2Abb[4491]=7. + Box2Abb[4490]*m_z2k;

  Box2Abb[4492]=3. + Box2Abb[4491]*m_z2k;

  Box2Abb[4493]=1. + Box2Abb[4492]*m_z2k;

  Box2Abb[4494]=-Box2Abb[4484]*pow(Box2Abb[61],4.) - Box2Abb[4487]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[4489]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[4493]*m_z12_3 + 4.*Box2Abb[411]*m_z12_4*m_z2k_2;

  Box2Abb[4495]=47. + 19.*m_z2k;

  Box2Abb[4496]=14. + Box2Abb[4495]*m_z2k;

  Box2Abb[4497]=5. + Box2Abb[4496]*m_z2k;

  Box2Abb[4498]=-25. + 44.*m_z2k;

  Box2Abb[4499]=-42. + Box2Abb[4498]*m_z2k;

  Box2Abb[4500]=7. + Box2Abb[4499]*m_z2k;

  Box2Abb[4501]=68. + 69.*m_z2k - 78.*m_z2k_2;

  Box2Abb[4502]=-33. + Box2Abb[4501]*m_z2k;

  Box2Abb[4503]=8. + Box2Abb[4502]*m_z2k;

  Box2Abb[4504]=-175. + 51.*m_z2k;

  Box2Abb[4505]=-51. + Box2Abb[4504]*m_z2k;

  Box2Abb[4506]=121. + Box2Abb[4505]*m_z2k;

  Box2Abb[4507]=-12. + Box2Abb[4506]*m_z2k;

  Box2Abb[4508]=-2.*pow(Box2Abb[61],4.) + Box2Abb[4500]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[4507]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[4503]*m_z12_3 - 2.*Box2Abb[4497]*m_z12_4 + 4.*m_z12_5*m_z2k_2;

  Box2Abb[4509]=14. - 5.*m_z2k;

  Box2Abb[4510]=9. + Box2Abb[4509]*m_z2k;

  Box2Abb[4511]=-8. - Box2Abb[12]*Box2Abb[2516]*Box2Abb[510]*m_z2k;

  Box2Abb[4512]=80. + 71.*m_z2k;

  Box2Abb[4513]=88. + Box2Abb[4512]*m_z2k;

  Box2Abb[4514]=6. + Box2Abb[4513]*m_z2k;

  Box2Abb[4515]=-172. + 155.*m_z2k;

  Box2Abb[4516]=-108. + Box2Abb[4515]*m_z2k;

  Box2Abb[4517]=-138. + Box2Abb[4516]*m_z2k;

  Box2Abb[4518]=-15. + Box2Abb[4517]*m_z2k;

  Box2Abb[4519]=8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4511]*m_z12 + Box2Abb[4518]*m_z12_2 + 2.*Box2Abb[4514]*m_z12_3 + 2.*Box2Abb[4510]*m_z12_4;

  Box2Abb[4520]=Box2Abb[4508]*m_x_2 + Box2Abb[4519]*m_x_3 + Box2Abb[4472]*m_x_4 + Box2Abb[4483]*m_x_5 - Box2Abb[4478]*m_x_6 + Box2Abb[4494]*m_x*m_z12 + Box2Abb[669]*m_x_7*m_z12 - Box2Abb[4475]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k;

  Box2Abb[4521]=19. + 9.*m_z12;

  Box2Abb[4522]=-35. + Box2Abb[4521]*m_z12;

  Box2Abb[4523]=-29. + 2.*Box2Abb[4522]*m_z12;

  Box2Abb[4524]=35. + Box2Abb[4523]*m_z12;

  Box2Abb[4525]=103. + 25.*m_z12;

  Box2Abb[4526]=-91. + Box2Abb[4525]*m_z12;

  Box2Abb[4527]=-91. + 2.*Box2Abb[4526]*m_z12;

  Box2Abb[4528]=32. + Box2Abb[4527]*m_z12;

  Box2Abb[4529]=-197. + 60.*Box2Abb[500]*m_z12;

  Box2Abb[4530]=33. + Box2Abb[4529]*m_z12;

  Box2Abb[4531]=-148. + Box2Abb[4530]*m_z12;

  Box2Abb[4532]=8. + Box2Abb[4531]*m_z12;

  Box2Abb[4533]=6. + 19.*m_z12;

  Box2Abb[4534]=-207. + 26.*Box2Abb[4533]*m_z12;

  Box2Abb[4535]=240. + Box2Abb[4534]*m_z12;

  Box2Abb[4536]=-16. + Box2Abb[4535]*m_z12;

  Box2Abb[4537]=-990. + 1121.*m_z12;

  Box2Abb[4538]=-55. + Box2Abb[4537]*m_z12;

  Box2Abb[4539]=30. + Box2Abb[4538]*m_z12;

  Box2Abb[4540]=-16. + 65.*m_z12;

  Box2Abb[4541]=2. + Box2Abb[4524]*m_z12 + 8.*m_z2k + Box2Abb[4528]*m_z12*m_z2k + Box2Abb[4532]*m_z2k_2 + Box2Abb[4536]*m_z2k_3 + Box2Abb[4539]*m_z2k_4 + 14.*Box2Abb[4540]*m_z12*m_z2k_5;

  Box2Abb[4542]=31. + 60.*m_z12 - 11.*m_z2k;

  Box2Abb[4543]=-5. + Box2Abb[4542]*m_z12 + 32.*m_z2k;

  Box2Abb[4544]=2. + Box2Abb[4543]*m_z12;

  Box2Abb[4545]=4. - 4.*m_z2k + 8.*m_z2k_2;

  Box2Abb[4546]=1. + Box2Abb[866]*m_z2k;

  Box2Abb[4547]=-8. + 9.*m_z2k;

  Box2Abb[4548]=3. + Box2Abb[4547]*m_z2k;

  Box2Abb[4549]=Box2Abb[4546]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4548]*Box2Abb[61]*m_z12 + Box2Abb[4545]*m_z12_2;

  Box2Abb[4550]=-81. + 56.*m_z2k;

  Box2Abb[4551]=89. + 157.*m_z2k;

  Box2Abb[4552]=-88. + 202.*m_z2k + 68.*m_z2k_2;

  Box2Abb[4553]=-54. + Box2Abb[4552]*m_z12 + 2.*Box2Abb[4551]*m_z12_2 + 31.*m_z12_3 + Box2Abb[4550]*m_z2k;

  Box2Abb[4554]=8. + Box2Abb[4553]*m_z12 + 12.*m_z2k;

  Box2Abb[4555]=97. + 186.*m_z2k;

  Box2Abb[4556]=2. + Box2Abb[1151]*m_z2k;

  Box2Abb[4557]=121. + 81.*Box2Abb[756]*m_z2k;

  Box2Abb[4558]=120. + 277.*m_z2k;

  Box2Abb[4559]=42. - Box2Abb[4558]*m_z2k;

  Box2Abb[4560]=100. + 91.*m_z2k;

  Box2Abb[4561]=-407. + 4.*Box2Abb[4560]*m_z2k;

  Box2Abb[4562]=-252. + Box2Abb[4561]*m_z2k;

  Box2Abb[4563]=6.*Box2Abb[4556] + Box2Abb[4559]*m_z12 + Box2Abb[4562]*m_z12_2 + Box2Abb[4557]*m_z12_3 + Box2Abb[4555]*m_z12_4 + 4.*m_z12_5;

  Box2Abb[4564]=14. + 27.*m_z2k;

  Box2Abb[4565]=2. + Box2Abb[450]*m_z2k;

  Box2Abb[4566]=347. + 417.*m_z2k;

  Box2Abb[4567]=104. + Box2Abb[4566]*m_z2k;

  Box2Abb[4568]=844. + 1053.*m_z2k;

  Box2Abb[4569]=95. + Box2Abb[4568]*m_z2k;

  Box2Abb[4570]=-50. + Box2Abb[4569]*m_z2k;

  Box2Abb[4571]=69. + 28.*m_z2k;

  Box2Abb[4572]=58. - 5.*Box2Abb[4571]*m_z2k;

  Box2Abb[4573]=46. + Box2Abb[4572]*m_z2k;

  Box2Abb[4574]=62. + Box2Abb[4573]*m_z2k;

  Box2Abb[4575]=1. + 35.*m_z2k;

  Box2Abb[4576]=-360. + 11.*Box2Abb[4575]*m_z2k;

  Box2Abb[4577]=-409. + 2.*Box2Abb[4576]*m_z2k;

  Box2Abb[4578]=-143. + Box2Abb[4577]*m_z2k;

  Box2Abb[4579]=8. + Box2Abb[4574]*m_z12 + Box2Abb[4578]*m_z12_2 + Box2Abb[4570]*m_z12_3 + Box2Abb[4567]*m_z12_4 + Box2Abb[4564]*m_z12_5 + 8.*Box2Abb[4565]*m_z2k;

  Box2Abb[4580]=-13. + 10.*m_z2k;

  Box2Abb[4581]=2. + Box2Abb[4580]*m_z2k;

  Box2Abb[4582]=-103. + 159.*m_z2k;

  Box2Abb[4583]=59. + Box2Abb[4582]*m_z2k;

  Box2Abb[4584]=-11. + Box2Abb[4583]*m_z2k;

  Box2Abb[4585]=26. + 7.*Box2Abb[2800]*m_z2k;

  Box2Abb[4586]=-10. + Box2Abb[4585]*m_z2k;

  Box2Abb[4587]=7. + Box2Abb[4586]*m_z2k;

  Box2Abb[4588]=-133. + 59.*m_z2k;

  Box2Abb[4589]=88. + Box2Abb[4588]*m_z2k;

  Box2Abb[4590]=-31. + Box2Abb[4589]*m_z2k;

  Box2Abb[4591]=1. + Box2Abb[4590]*m_z2k;

  Box2Abb[4592]=-158. + 115.*m_z2k;

  Box2Abb[4593]=108. + Box2Abb[4592]*m_z2k;

  Box2Abb[4594]=-50. + Box2Abb[4593]*m_z2k;

  Box2Abb[4595]=17. + Box2Abb[4594]*m_z2k;

  Box2Abb[4596]=Box2Abb[4591]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[4584]*pow(Box2Abb[61],3.)*m_z12_2 + Box2Abb[4595]*Box2Abb[61]*m_z12_3 + Box2Abb[4587]*m_z12_4 - Box2Abb[4581]*pow(Box2Abb[61],4.)*m_z2k;

  Box2Abb[4597]=1. + 3.*m_z2k_2;

  Box2Abb[4598]=58. + 63.*m_z2k;

  Box2Abb[4599]=40. + Box2Abb[4598]*m_z2k;

  Box2Abb[4600]=10. + Box2Abb[4599]*m_z2k;

  Box2Abb[4601]=16. + 79.*m_z2k;

  Box2Abb[4602]=2. + 5.*Box2Abb[4601]*m_z2k;

  Box2Abb[4603]=74. + Box2Abb[4602]*m_z2k;

  Box2Abb[4604]=-3. + Box2Abb[4603]*m_z2k;

  Box2Abb[4605]=-757. + 322.*m_z2k;

  Box2Abb[4606]=520. + Box2Abb[4605]*m_z2k;

  Box2Abb[4607]=-99. + Box2Abb[4606]*m_z2k;

  Box2Abb[4608]=-130. + Box2Abb[4607]*m_z2k;

  Box2Abb[4609]=54. + Box2Abb[4608]*m_z2k;

  Box2Abb[4610]=-884. + 905.*m_z2k;

  Box2Abb[4611]=213. + Box2Abb[4610]*m_z2k;

  Box2Abb[4612]=174. + Box2Abb[4611]*m_z2k;

  Box2Abb[4613]=-200. + Box2Abb[4612]*m_z2k;

  Box2Abb[4614]=-18. + Box2Abb[4613]*m_z2k;

  Box2Abb[4615]=-35. + 24.*m_z2k;

  Box2Abb[4616]=18. + 7.*Box2Abb[4615]*m_z2k;

  Box2Abb[4617]=140. + Box2Abb[4616]*m_z2k;

  Box2Abb[4618]=-124. + Box2Abb[4617]*m_z2k;

  Box2Abb[4619]=37. + Box2Abb[4618]*m_z2k;

  Box2Abb[4620]=10. - Box2Abb[4619]*m_z2k;

  Box2Abb[4621]=Box2Abb[4620]*m_z12 + Box2Abb[4614]*m_z12_3 + Box2Abb[4604]*m_z12_4 + Box2Abb[4600]*m_z12_5 + 4.*Box2Abb[4597]*pow(Box2Abb[61],2.)*m_z2k + 2.*Box2Abb[4609]*m_z12_2*m_z2k;

  Box2Abb[4622]=12. + m_z2k + 5.*m_z2k_2 + 20.*m_z2k_3;

  Box2Abb[4623]=1. + Box2Abb[4622]*m_z2k;

  Box2Abb[4624]=-81. + 64.*m_z2k;

  Box2Abb[4625]=22. + Box2Abb[4624]*m_z2k;

  Box2Abb[4626]=-13. + Box2Abb[4625]*m_z2k;

  Box2Abb[4627]=-275. + 254.*m_z2k;

  Box2Abb[4628]=166. + Box2Abb[4627]*m_z2k;

  Box2Abb[4629]=6. + Box2Abb[4628]*m_z2k;

  Box2Abb[4630]=-16. + Box2Abb[4629]*m_z2k;

  Box2Abb[4631]=-3. + Box2Abb[4630]*m_z2k;

  Box2Abb[4632]=-185. + 67.*m_z2k;

  Box2Abb[4633]=679. + 4.*Box2Abb[4632]*m_z2k;

  Box2Abb[4634]=-327. + Box2Abb[4633]*m_z2k;

  Box2Abb[4635]=121. + Box2Abb[4634]*m_z2k;

  Box2Abb[4636]=-31. + Box2Abb[4635]*m_z2k;

  Box2Abb[4637]=-2. + Box2Abb[4636]*m_z2k;

  Box2Abb[4638]=-1136. + 499.*m_z2k;

  Box2Abb[4639]=1003. + Box2Abb[4638]*m_z2k;

  Box2Abb[4640]=-490. + Box2Abb[4639]*m_z2k;

  Box2Abb[4641]=133. + Box2Abb[4640]*m_z2k;

  Box2Abb[4642]=-40. + Box2Abb[4641]*m_z2k;

  Box2Abb[4643]=-1. + Box2Abb[4642]*m_z2k;

  Box2Abb[4644]=Box2Abb[4637]*Box2Abb[61]*m_z12_2 + Box2Abb[4643]*m_z12_3 + Box2Abb[4631]*m_z12_4 + 2.*Box2Abb[4623]*m_z12_5 + 2.*pow(Box2Abb[61],4.)*m_z2k_2 - Box2Abb[4626]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[4645]=-Box2Abb[4644]*m_x_3 + Box2Abb[4621]*m_x_4 - Box2Abb[4541]*m_x_5 + Box2Abb[4579]*m_x_6 - Box2Abb[4563]*m_x_7 + Box2Abb[4554]*m_x_8 - Box2Abb[4544]*m_x_9 + Box2Abb[1087]*m_x_10*m_z12 + Box2Abb[4596]*m_x_2*m_z12*m_z2k - Box2Abb[4549]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12_2*m_z2k_2 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_3*m_z2k_3;

  Box2Abb[4646]=Box2Abb[4645]*m_cR - Box2Abb[4520]*pow(Box2Abb[8],2.)*m_cL*m_x;

  Box2Abb[4647]=2.*Box2Abb[9] - 3.*Box2Abb[70]*m_z2k + 4.*m_z2k_2;

  Box2Abb[4648]=46. + 15.*m_z12 + 12.*m_z2k;

  Box2Abb[4649]=-16. + Box2Abb[4648]*m_z12;

  Box2Abb[4650]=62. - 18.*m_z2k;

  Box2Abb[4651]=95. + Box2Abb[4650]*m_z12 + 4.*m_z12_2 + 28.*Box2Abb[2020]*m_z2k;

  Box2Abb[4652]=-70. + Box2Abb[4651]*m_z12 - 44.*m_z2k;

  Box2Abb[4653]=7. - 31.*m_z2k;

  Box2Abb[4654]=31. + 2.*Box2Abb[4653]*m_z2k;

  Box2Abb[4655]=40. + 2.*Box2Abb[4654]*m_z2k;

  Box2Abb[4656]=27. - 68.*m_z2k + 40.*m_z2k_2;

  Box2Abb[4657]=-1. + 2.*Box2Abb[4656]*m_z2k;

  Box2Abb[4658]=-32. + m_z12 + Box2Abb[4655]*m_z12_2 - 3.*Box2Abb[1871]*m_z12_3 + Box2Abb[4657]*m_z12*m_z2k + 8.*Box2Abb[2193]*m_z2k_2;

  Box2Abb[4659]=-23. + 4.*Box2Abb[2152]*m_z2k;

  Box2Abb[4660]=2. + Box2Abb[4659]*m_z2k;

  Box2Abb[4661]=-23. + 7.*m_z2k;

  Box2Abb[4662]=10. + Box2Abb[4661]*m_z2k;

  Box2Abb[4663]=-1. + Box2Abb[4662]*m_z2k;

  Box2Abb[4664]=2.*Box2Abb[230]*pow(Box2Abb[61],4.) - Box2Abb[4660]*pow(Box2Abb[61],3.)*m_z12 + 2.*Box2Abb[4663]*pow(Box2Abb[61],2.)*m_z12_2 - Box2Abb[1880]*Box2Abb[61]*m_z12_3 - 12.*m_z12_4*m_z2k_3;

  Box2Abb[4665]=-14. + 5.*m_z2k;

  Box2Abb[4666]=-28. + 127.*m_z2k;

  Box2Abb[4667]=-84. + Box2Abb[4666]*m_z2k;

  Box2Abb[4668]=37. + 50.*m_z2k;

  Box2Abb[4669]=59. + Box2Abb[4668]*m_z2k;

  Box2Abb[4670]=47. + 2.*Box2Abb[4669]*m_z2k;

  Box2Abb[4671]=76. - Box2Abb[4670]*m_z12 + Box2Abb[4667]*m_z12_2 + Box2Abb[4665]*m_z12_3 + 4.*Box2Abb[2153]*m_z2k;

  Box2Abb[4672]=-3. + Box2Abb[234]*m_z2k;

  Box2Abb[4673]=21. + 43.*m_z2k;

  Box2Abb[4674]=10. + Box2Abb[4673]*m_z2k;

  Box2Abb[4675]=-71. + 10.*m_z2k;

  Box2Abb[4676]=82. + Box2Abb[4675]*m_z2k;

  Box2Abb[4677]=-48. + m_z2k + 2.*Box2Abb[4676]*m_z2k_2;

  Box2Abb[4678]=5. + Box2Abb[4677]*m_z2k;

  Box2Abb[4679]=16. + 15.*m_z2k;

  Box2Abb[4680]=18. + Box2Abb[4679]*m_z2k;

  Box2Abb[4681]=-44. + Box2Abb[4680]*m_z2k;

  Box2Abb[4682]=-1. + Box2Abb[4681]*m_z2k;

  Box2Abb[4683]=-4.*Box2Abb[4672]*pow(Box2Abb[61],2.) - Box2Abb[4678]*m_z12 + Box2Abb[4682]*m_z12_2 + Box2Abb[230]*Box2Abb[4674]*m_z12_3 + 12.*m_z12_4*m_z2k_2;

  Box2Abb[4684]=Box2Abb[4664]*m_x + Box2Abb[4683]*m_x_2 + Box2Abb[4658]*m_x_3 + Box2Abb[4671]*m_x_4 + Box2Abb[4652]*m_x_5 - Box2Abb[4649]*m_x_6 + Box2Abb[4647]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[4685]=-6. + m_z12 - 6.*m_z12*m_z2k;

  Box2Abb[4686]=17. + 4.*m_z12;

  Box2Abb[4687]=-16. + Box2Abb[4686]*m_z12;

  Box2Abb[4688]=4. + m_z12 + 26.*m_z12_2;

  Box2Abb[4689]=-4.*Box2Abb[4323] - 10.*Box2Abb[9]*m_z12 + Box2Abb[4687]*m_z12*m_z2k + Box2Abb[4688]*m_z2k_2 - 10.*m_z12*m_z2k_3;

  Box2Abb[4690]=-2. + Box2Abb[230]*m_z2k;

  Box2Abb[4691]=Box2Abb[4690]*pow(Box2Abb[61],3.) + Box2Abb[258]*Box2Abb[444]*pow(Box2Abb[61],2.)*m_z12 + 2.*Box2Abb[12]*Box2Abb[230]*Box2Abb[61]*m_z12_2 + Box2Abb[231]*m_z12_3*m_z2k;

  Box2Abb[4692]=-3. + 17.*m_z2k;

  Box2Abb[4693]=-13. + Box2Abb[4692]*m_z12 + 2.*Box2Abb[841]*m_z2k;

  Box2Abb[4694]=23. + Box2Abb[4693]*m_z12 + 14.*m_z2k;

  Box2Abb[4695]=7. + 6.*m_z2k;

  Box2Abb[4696]=4. + Box2Abb[4695]*m_z2k;

  Box2Abb[4697]=16. + Box2Abb[4692]*m_z2k;

  Box2Abb[4698]=12. + Box2Abb[4697]*m_z2k;

  Box2Abb[4699]=8. - 5.*Box2Abb[179]*m_z2k;

  Box2Abb[4700]=4. + Box2Abb[4699]*m_z2k;

  Box2Abb[4701]=7. + 2.*Box2Abb[4700]*m_z2k;

  Box2Abb[4702]=10. + Box2Abb[4701]*m_z12 - Box2Abb[4698]*m_z12_2 - 2.*Box2Abb[4696]*m_z2k + 2.*Box2Abb[3116]*m_z12_3*m_z2k;

  Box2Abb[4703]=-15. + 8.*m_z2k;

  Box2Abb[4704]=2. + Box2Abb[4703]*m_z2k;

  Box2Abb[4705]=4. + Box2Abb[4704]*m_z2k;

  Box2Abb[4706]=-28. + m_z2k;

  Box2Abb[4707]=14. + Box2Abb[4706]*m_z2k;

  Box2Abb[4708]=9. + Box2Abb[4707]*m_z2k;

  Box2Abb[4709]=-1. + Box2Abb[4708]*m_z2k;

  Box2Abb[4710]=-13. - 2.*m_z2k + 6.*m_z2k_2;

  Box2Abb[4711]=2. + Box2Abb[4710]*m_z2k;

  Box2Abb[4712]=2. + Box2Abb[4711]*m_z2k;

  Box2Abb[4713]=Box2Abb[444]*pow(Box2Abb[61],4.) - Box2Abb[4712]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4709]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[4705]*m_z12_3*m_z2k + m_z12_4*m_z2k_3;

  Box2Abb[4714]=1. + 5.*m_z2k_2;

  Box2Abb[4715]=6.*m_z2k + 8.*m_z2k_3;

  Box2Abb[4716]=-23. + 4.*Box2Abb[2303]*m_z2k;

  Box2Abb[4717]=10. + Box2Abb[4716]*m_z2k;

  Box2Abb[4718]=6. + Box2Abb[4717]*m_z2k;

  Box2Abb[4719]=2. + Box2Abb[1405]*m_z2k;

  Box2Abb[4720]=-46. + Box2Abb[4719]*m_z2k;

  Box2Abb[4721]=14. + Box2Abb[4720]*m_z2k;

  Box2Abb[4722]=9. + Box2Abb[4721]*m_z2k;

  Box2Abb[4723]=2.*Box2Abb[4714]*pow(Box2Abb[61],2.) - Box2Abb[4722]*m_z12 + Box2Abb[4718]*m_z12_2 + Box2Abb[4715]*m_z12_3 - 2.*m_z12_4*m_z2k_2;

  Box2Abb[4724]=-Box2Abb[4713]*m_x_2 + Box2Abb[4723]*m_x_3 + Box2Abb[4702]*m_x_4 - Box2Abb[4689]*m_x_5 + Box2Abb[4694]*m_x_6 + Box2Abb[4685]*m_x_7 + 2.*m_x_8*m_z12 - Box2Abb[4691]*Box2Abb[61]*m_x*m_z12*m_z2k + pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_2;

  Box2Abb[4725]=2.*Box2Abb[4724]*m_cL + Box2Abb[4684]*m_cR*m_x;

  Box2Abb[4726]=-2.*Box2Abb[230]*pow(Box2Abb[61],2.) + Box2Abb[150]*Box2Abb[61]*m_z12 + 8.*m_z12_2*m_z2k;

  Box2Abb[4727]=10. + 5.*m_z12 + 4.*m_z2k;

  Box2Abb[4728]=-16. + 3.*Box2Abb[4727]*m_z12;

  Box2Abb[4729]=29. + 7.*m_z2k;

  Box2Abb[4730]=31. + 2.*Box2Abb[4729]*m_z12 + 4.*m_z12_2 + 92.*m_z2k + 56.*m_z2k_2;

  Box2Abb[4731]=-54. + Box2Abb[4730]*m_z12 - 44.*m_z2k;

  Box2Abb[4732]=-2. + m_z2k + 7.*Box2Abb[1163]*m_z2k_2;

  Box2Abb[4733]=-2. + 23.*m_z2k;

  Box2Abb[4734]=-3. + Box2Abb[4733]*Box2Abb[61]*m_z2k;

  Box2Abb[4735]=-7. + 4.*Box2Abb[833]*m_z2k;

  Box2Abb[4736]=6. + Box2Abb[4735]*m_z2k;

  Box2Abb[4737]=2.*Box2Abb[230]*pow(Box2Abb[61],4.) - Box2Abb[4736]*pow(Box2Abb[61],3.)*m_z12 + 2.*Box2Abb[4734]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[4732]*Box2Abb[61]*m_z12_3 + 4.*Box2Abb[411]*m_z12_4*m_z2k_2;

  Box2Abb[4738]=14. + 3.*m_z2k;

  Box2Abb[4739]=-7. + m_z2k - 8.*m_z2k_2;

  Box2Abb[4740]=96. + m_z2k;

  Box2Abb[4741]=68. + Box2Abb[4740]*m_z2k;

  Box2Abb[4742]=21. + 50.*m_z2k;

  Box2Abb[4743]=11. + Box2Abb[4742]*m_z2k;

  Box2Abb[4744]=-45. + 2.*Box2Abb[4743]*m_z2k;

  Box2Abb[4745]=4.*Box2Abb[4739] + Box2Abb[4744]*m_z12 + Box2Abb[4741]*m_z12_2 + Box2Abb[4738]*m_z12_3;

  Box2Abb[4746]=-17. + 20.*m_z2k;

  Box2Abb[4747]=-18. + Box2Abb[4746]*m_z2k;

  Box2Abb[4748]=43. + 50.*m_z2k + 34.*m_z2k_2;

  Box2Abb[4749]=8. + Box2Abb[4748]*m_z2k;

  Box2Abb[4750]=-4. + Box2Abb[2099]*m_z2k;

  Box2Abb[4751]=2. + Box2Abb[4750]*m_z2k;

  Box2Abb[4752]=5. + 52.*m_z2k - 40.*m_z2k_2;

  Box2Abb[4753]=29. + 2.*Box2Abb[4752]*m_z2k;

  Box2Abb[4754]=51. + Box2Abb[4753]*m_z2k;

  Box2Abb[4755]=-8.*Box2Abb[4751] + Box2Abb[4754]*m_z12 - 2.*Box2Abb[4749]*m_z12_2 + Box2Abb[4747]*m_z12_3;

  Box2Abb[4756]=85. + 26.*m_z2k;

  Box2Abb[4757]=13. + Box2Abb[4756]*m_z2k;

  Box2Abb[4758]=10. + Box2Abb[4757]*m_z2k;

  Box2Abb[4759]=-3. + Box2Abb[4318]*m_z2k;

  Box2Abb[4760]=23. + 2.*Box2Abb[4759]*m_z2k;

  Box2Abb[4761]=-1. + Box2Abb[4760]*m_z2k;

  Box2Abb[4762]=-88. + 113.*m_z2k;

  Box2Abb[4763]=-66. + Box2Abb[4762]*m_z2k;

  Box2Abb[4764]=20. + Box2Abb[4763]*m_z2k;

  Box2Abb[4765]=-15. + Box2Abb[4764]*m_z2k;

  Box2Abb[4766]=4.*Box2Abb[1601]*pow(Box2Abb[61],2.) + Box2Abb[4761]*Box2Abb[61]*m_z12 + Box2Abb[4765]*m_z12_2 + Box2Abb[4758]*m_z12_3 - 4.*m_z12_4*m_z2k_2;

  Box2Abb[4767]=-Box2Abb[4737]*m_x + Box2Abb[4766]*m_x_2 + Box2Abb[4755]*m_x_3 + Box2Abb[4745]*m_x_4 - Box2Abb[4731]*m_x_5 + Box2Abb[4728]*m_x_6 + Box2Abb[4726]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[4768]=50. - 117.*m_z12 + 7.*Box2Abb[2038]*m_z12_3;

  Box2Abb[4769]=349. + 26.*m_z12;

  Box2Abb[4770]=88. + Box2Abb[4769]*m_z12;

  Box2Abb[4771]=-365. + Box2Abb[4770]*m_z12;

  Box2Abb[4772]=64. + Box2Abb[4771]*m_z12;

  Box2Abb[4773]=77. + 48.*m_z12;

  Box2Abb[4774]=-585. + 8.*Box2Abb[4773]*m_z12;

  Box2Abb[4775]=78. + Box2Abb[4774]*m_z12;

  Box2Abb[4776]=-12. + 61.*m_z12;

  Box2Abb[4777]=-25. + 3.*Box2Abb[4776]*m_z12;

  Box2Abb[4778]=2.*Box2Abb[4768] + Box2Abb[4772]*m_z2k + Box2Abb[4775]*m_z2k_2 + 4.*Box2Abb[4777]*m_z2k_3 + 336.*m_z12*m_z2k_4;

  Box2Abb[4779]=42. + 73.*m_z12 + 12.*m_z2k;

  Box2Abb[4780]=-16. + Box2Abb[4779]*m_z12;

  Box2Abb[4781]=-3. + 19.*m_z2k_2;

  Box2Abb[4782]=-1. + Box2Abb[460]*m_z2k;

  Box2Abb[4783]=2.*Box2Abb[4782]*pow(Box2Abb[61],2.) + Box2Abb[4781]*Box2Abb[61]*m_z12 + 8.*Box2Abb[12]*m_z12_2*m_z2k;

  Box2Abb[4784]=12. + 5.*m_z2k;

  Box2Abb[4785]=232. + 362.*m_z2k;

  Box2Abb[4786]=-107. + Box2Abb[4785]*m_z12 + 34.*m_z12_2 + 16.*Box2Abb[4784]*m_z2k;

  Box2Abb[4787]=-50. + Box2Abb[4786]*m_z12 - 76.*m_z2k;

  Box2Abb[4788]=9. + 17.*m_z2k;

  Box2Abb[4789]=110. + 193.*m_z2k;

  Box2Abb[4790]=203. + 714.*m_z2k + 724.*m_z2k_2;

  Box2Abb[4791]=31. + 28.*m_z2k;

  Box2Abb[4792]=-233. + 4.*Box2Abb[4791]*m_z2k;

  Box2Abb[4793]=-333. + 2.*Box2Abb[4792]*m_z2k;

  Box2Abb[4794]=72. + Box2Abb[4793]*m_z12 + Box2Abb[4790]*m_z12_2 + Box2Abb[4789]*m_z12_3 + 4.*m_z12_4 - 8.*Box2Abb[4788]*m_z2k;

  Box2Abb[4795]=23. + 26.*m_z2k;

  Box2Abb[4796]=9. + Box2Abb[4795]*m_z2k;

  Box2Abb[4797]=-1. + 2.*Box2Abb[4547]*m_z2k;

  Box2Abb[4798]=7. + Box2Abb[4797]*m_z2k;

  Box2Abb[4799]=293. + 323.*m_z2k;

  Box2Abb[4800]=186. + Box2Abb[4799]*m_z2k;

  Box2Abb[4801]=56. + Box2Abb[4800]*m_z2k;

  Box2Abb[4802]=-327. - 226.*m_z2k + 358.*m_z2k_2;

  Box2Abb[4803]=-182. + Box2Abb[4802]*m_z2k;

  Box2Abb[4804]=-67. + Box2Abb[4803]*m_z2k;

  Box2Abb[4805]=-33. + 14.*m_z2k;

  Box2Abb[4806]=48. + 5.*Box2Abb[4805]*m_z2k;

  Box2Abb[4807]=45. + Box2Abb[4806]*m_z2k;

  Box2Abb[4808]=-17. + Box2Abb[4807]*m_z2k_2;

  Box2Abb[4809]=8.*Box2Abb[4798] + 4.*Box2Abb[4808]*m_z12 + Box2Abb[4804]*m_z12_2 + Box2Abb[4801]*m_z12_3 + 2.*Box2Abb[4796]*m_z12_4;

  Box2Abb[4810]=18. - 65.*m_z2k + 16.*m_z2k_3;

  Box2Abb[4811]=-1. + Box2Abb[4810]*m_z2k;

  Box2Abb[4812]=29. + 19.*m_z2k;

  Box2Abb[4813]=-6. + Box2Abb[4812]*m_z2k;

  Box2Abb[4814]=5. + 2.*Box2Abb[4813]*m_z2k;

  Box2Abb[4815]=10. + 7.*m_z2k;

  Box2Abb[4816]=-9. + Box2Abb[4815]*m_z2k;

  Box2Abb[4817]=11. + Box2Abb[4816]*m_z2k;

  Box2Abb[4818]=-3. + Box2Abb[4817]*m_z2k;

  Box2Abb[4819]=63. + 68.*m_z2k;

  Box2Abb[4820]=-95. + Box2Abb[4819]*m_z2k;

  Box2Abb[4821]=43. + Box2Abb[4820]*m_z2k;

  Box2Abb[4822]=-15. + Box2Abb[4821]*m_z2k;

  Box2Abb[4823]=Box2Abb[4811]*pow(Box2Abb[61],3.)*m_z12 + 2.*Box2Abb[4814]*pow(Box2Abb[61],3.)*m_z12_2 + Box2Abb[4822]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[4818]*m_z12_4 - 2.*Box2Abb[460]*pow(Box2Abb[61],4.)*m_z2k;

  Box2Abb[4824]=17. + 19.*Box2Abb[12]*m_z2k;

  Box2Abb[4825]=5. + Box2Abb[4824]*m_z2k;

  Box2Abb[4826]=-51. + 35.*m_z2k;

  Box2Abb[4827]=-58. + Box2Abb[4826]*m_z2k;

  Box2Abb[4828]=29. + Box2Abb[4827]*m_z2k;

  Box2Abb[4829]=2. + Box2Abb[4828]*m_z2k;

  Box2Abb[4830]=131. - 326.*m_z2k + 4.*m_z2k_2;

  Box2Abb[4831]=73. + Box2Abb[4830]*m_z2k;

  Box2Abb[4832]=-82. + Box2Abb[4831]*m_z2k;

  Box2Abb[4833]=-13. + Box2Abb[4832]*m_z2k;

  Box2Abb[4834]=-3. + 22.*m_z2k;

  Box2Abb[4835]=-56. + Box2Abb[4834]*m_z2k;

  Box2Abb[4836]=40. + Box2Abb[4835]*m_z2k;

  Box2Abb[4837]=-26. + Box2Abb[4836]*m_z2k;

  Box2Abb[4838]=7. + Box2Abb[4837]*m_z2k;

  Box2Abb[4839]=-38. + 7.*m_z2k;

  Box2Abb[4840]=963. + 16.*Box2Abb[4839]*m_z2k;

  Box2Abb[4841]=-254. + Box2Abb[4840]*m_z2k;

  Box2Abb[4842]=-80. + Box2Abb[4841]*m_z2k;

  Box2Abb[4843]=114. + Box2Abb[4842]*m_z2k;

  Box2Abb[4844]=-3. + Box2Abb[4843]*m_z2k;

  Box2Abb[4845]=2.*Box2Abb[4838] + Box2Abb[4844]*m_z12 + 2.*Box2Abb[4833]*m_z12_2 + 2.*Box2Abb[4829]*m_z12_3 + 2.*Box2Abb[4825]*m_z12_4;

  Box2Abb[4846]=-2. + 3.*Box2Abb[61]*m_z2k;

  Box2Abb[4847]=-5. + 7.*m_z2k;

  Box2Abb[4848]=-1. + Box2Abb[187]*Box2Abb[4847]*m_z2k;

  Box2Abb[4849]=97. + 14.*Box2Abb[3986]*m_z2k;

  Box2Abb[4850]=-89. + 3.*Box2Abb[4849]*m_z2k;

  Box2Abb[4851]=25. + Box2Abb[4850]*m_z2k;

  Box2Abb[4852]=3. + Box2Abb[4851]*m_z2k;

  Box2Abb[4853]=204. + 73.*m_z2k;

  Box2Abb[4854]=-124. + Box2Abb[4853]*m_z2k;

  Box2Abb[4855]=58. + Box2Abb[4854]*m_z2k;

  Box2Abb[4856]=15. + Box2Abb[4855]*m_z2k;

  Box2Abb[4857]=2. + Box2Abb[4856]*m_z2k;

  Box2Abb[4858]=133. + 54.*m_z2k;

  Box2Abb[4859]=-745. + 2.*Box2Abb[4858]*m_z2k;

  Box2Abb[4860]=478. + Box2Abb[4859]*m_z2k;

  Box2Abb[4861]=-204. + Box2Abb[4860]*m_z2k;

  Box2Abb[4862]=30. + Box2Abb[4861]*m_z2k;

  Box2Abb[4863]=3. + Box2Abb[4862]*m_z2k;

  Box2Abb[4864]=Box2Abb[4852]*Box2Abb[61]*m_z12 + Box2Abb[4863]*m_z12_2 + Box2Abb[4857]*m_z12_3 + 2.*Box2Abb[4848]*m_z12_4 - 8.*Box2Abb[4846]*pow(Box2Abb[61],2.)*m_z2k_2;

  Box2Abb[4865]=Box2Abb[4864]*m_x_3 + Box2Abb[4845]*m_x_4 - Box2Abb[4809]*m_x_5 + Box2Abb[4778]*m_x_6 - Box2Abb[4794]*m_x_7 + Box2Abb[4787]*m_x_8 - Box2Abb[4780]*m_x_9 - Box2Abb[4823]*m_x_2*m_z2k + Box2Abb[4783]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12*m_z2k_2 - 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_2*m_z2k_3;

  Box2Abb[4866]=Box2Abb[4865]*m_cR + Box2Abb[4767]*pow(Box2Abb[8],2.)*m_cL*m_x;

  Box2Abb[4867]=-4. + 2.*Box2Abb[184]*m_z12 + 5.*m_z12_2;

  Box2Abb[4868]=-9. + 2.*Box2Abb[2176]*m_z12;

  Box2Abb[4869]=-2.*Box2Abb[4695] + 3.*pow(Box2Abb[70],2.)*m_z12 + 2.*Box2Abb[11]*Box2Abb[2741]*m_z12*m_z2k + 2.*Box2Abb[4868]*m_z2k_2 + 20.*m_z12*m_z2k_3;

  Box2Abb[4870]=-12. + m_z2k;

  Box2Abb[4871]=1. + Box2Abb[4870]*m_z2k;

  Box2Abb[4872]=-1. + 5.*Box2Abb[12]*m_z2k;

  Box2Abb[4873]=2.*pow(Box2Abb[61],3.) - 2.*Box2Abb[4872]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4871]*Box2Abb[61]*m_z12_2 + Box2Abb[1848]*m_z12_3;

  Box2Abb[4874]=23. + 5.*m_z2k;

  Box2Abb[4875]=14. + 11.*m_z2k;

  Box2Abb[4876]=20. + Box2Abb[4875]*m_z12 + m_z12_2 + 2.*Box2Abb[4874]*m_z2k;

  Box2Abb[4877]=-2.*Box2Abb[679] + Box2Abb[4876]*m_z12;

  Box2Abb[4878]=-3. + 2.*Box2Abb[280]*m_z2k_2;

  Box2Abb[4879]=13. - 2.*m_z2k;

  Box2Abb[4880]=6. + Box2Abb[4879]*m_z2k;

  Box2Abb[4881]=2. + 2.*Box2Abb[4880]*m_z2k;

  Box2Abb[4882]=-2.*Box2Abb[1151]*pow(Box2Abb[61],2.) + 2.*Box2Abb[4878]*Box2Abb[61]*m_z12 + Box2Abb[4881]*m_z12_2 + 3.*Box2Abb[12]*m_z12_3;

  Box2Abb[4883]=-Box2Abb[4873]*Box2Abb[61]*m_x - Box2Abb[4882]*m_x_2 + Box2Abb[4869]*m_x_3 - Box2Abb[4877]*m_x_4 + Box2Abb[4867]*m_x_5 + Box2Abb[3336]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k;

  Box2Abb[4884]=2.*pow(Box2Abb[61],2.) + 7.*Box2Abb[61]*m_z12 + 4.*m_z12_2;

  Box2Abb[4885]=4. + m_z12 + 2.*m_z2k;

  Box2Abb[4886]=-4. + Box2Abb[4885]*m_z12;

  Box2Abb[4887]=5. + 3.*m_z2k + 9.*m_z2k_2;

  Box2Abb[4888]=1. + m_z2k + 4.*m_z2k_2 + 5.*m_z2k_3;

  Box2Abb[4889]=3. + 2.*Box2Abb[2020]*m_z2k;

  Box2Abb[4890]=-2.*Box2Abb[4887] + 4.*Box2Abb[4888]*m_z12 + Box2Abb[4889]*m_z12_2;

  Box2Abb[4891]=8. + Box2Abb[222]*m_z12 + 2.*Box2Abb[280]*m_z2k;

  Box2Abb[4892]=-2.*Box2Abb[1567] + Box2Abb[4891]*m_z12;

  Box2Abb[4893]=9. + 19.*m_z2k;

  Box2Abb[4894]=-6. + Box2Abb[4893]*m_z2k;

  Box2Abb[4895]=-2.*pow(Box2Abb[61],3.) + 10.*Box2Abb[12]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[4894]*Box2Abb[61]*m_z12_2 + 4.*Box2Abb[187]*m_z12_3*m_z2k;

  Box2Abb[4896]=-4. + Box2Abb[1069]*m_z2k;

  Box2Abb[4897]=-4. + Box2Abb[4896]*m_z2k;

  Box2Abb[4898]=9. + 8.*m_z2k;

  Box2Abb[4899]=15. + 2.*Box2Abb[4898]*m_z2k;

  Box2Abb[4900]=1. + Box2Abb[4899]*m_z2k;

  Box2Abb[4901]=-2. + Box2Abb[4900]*m_z12_2 + 4.*Box2Abb[4897]*m_z12*m_z2k + 2.*Box2Abb[4047]*m_z2k_2;

  Box2Abb[4902]=Box2Abb[4901]*m_x_2 - Box2Abb[4890]*m_x_3 + Box2Abb[4892]*m_x_4 - Box2Abb[4886]*m_x_5 - Box2Abb[4895]*m_x*m_z2k + Box2Abb[4884]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[4903]=Box2Abb[4902]*m_cL + Box2Abb[4883]*m_cR;

  Box2Abb[4904]=pow(Box2Abb[61],3.) + 6.*pow(Box2Abb[61],2.)*m_z12 + 6.*Box2Abb[61]*m_z12_2 + 2.*m_z12_3;

  Box2Abb[4905]=-16. + 9.*m_z12;

  Box2Abb[4906]=29. + 2.*Box2Abb[4905]*m_z12;

  Box2Abb[4907]=20. + 9.*m_z12 + 3.*m_z12_3;

  Box2Abb[4908]=-16. + m_z12 + 6.*m_z12_2;

  Box2Abb[4909]=6. + Box2Abb[4908]*m_z12;

  Box2Abb[4910]=9. + 2.*Box2Abb[4]*m_z12;

  Box2Abb[4911]=15. + Box2Abb[2828]*m_z12 + 2.*Box2Abb[4906]*m_z12*m_z2k + 2.*Box2Abb[4907]*m_z12*m_z2k_2 - 5.*Box2Abb[4909]*m_z2k_3 + 5.*Box2Abb[4910]*m_z2k_4 + 56.*m_z12*m_z2k_5;

  Box2Abb[4912]=-59. + 24.*m_z12;

  Box2Abb[4913]=54. + Box2Abb[4912]*m_z12;

  Box2Abb[4914]=41. + 9.*Box2Abb[72]*m_z12;

  Box2Abb[4915]=-27. + Box2Abb[4914]*m_z12;

  Box2Abb[4916]=12. + 5.*m_z12;

  Box2Abb[4917]=-46. + Box2Abb[4916]*m_z12;

  Box2Abb[4918]=20. + Box2Abb[4917]*m_z12;

  Box2Abb[4919]=-15. + Box2Abb[4918]*m_z12;

  Box2Abb[4920]=9. + m_z12 - 12.*m_z12_2;

  Box2Abb[4921]=-9. + Box2Abb[4920]*m_z12;

  Box2Abb[4922]=-18. + 11.*m_z12;

  Box2Abb[4923]=18. + Box2Abb[4922]*m_z12;

  Box2Abb[4924]=3. + Box2Abb[538]*m_z12 - 6.*m_z2k + Box2Abb[4913]*m_z12*m_z2k + 2.*Box2Abb[4915]*m_z12*m_z2k_2 - 2.*Box2Abb[4919]*m_z2k_3 + 5.*Box2Abb[4921]*m_z2k_4 + Box2Abb[4923]*m_z2k_5 + 28.*m_z12*m_z2k_6;

  Box2Abb[4925]=2. + m_z12 + 8.*m_z2k;

  Box2Abb[4926]=3. + Box2Abb[4925]*m_z12;

  Box2Abb[4927]=3. + 5.*Box2Abb[61]*m_z2k;

  Box2Abb[4928]=-18. + 5.*Box2Abb[12]*m_z2k;

  Box2Abb[4929]=Box2Abb[222]*pow(Box2Abb[61],5.) + Box2Abb[4928]*pow(Box2Abb[61],3.)*m_z12 + 6.*Box2Abb[222]*pow(Box2Abb[61],3.)*m_z12_2 - 2.*Box2Abb[4927]*Box2Abb[61]*m_z12_3 + 4.*Box2Abb[1613]*m_z12_4*m_z2k;

  Box2Abb[4930]=5. + 14.*m_z2k;

  Box2Abb[4931]=-5. + 5.*Box2Abb[12]*m_z12 + 2.*Box2Abb[4930]*m_z2k;

  Box2Abb[4932]=3.*Box2Abb[1916] + Box2Abb[4931]*m_z12;

  Box2Abb[4933]=9. + 28.*m_z2k;

  Box2Abb[4934]=-10. + Box2Abb[4933]*m_z2k_2;

  Box2Abb[4935]=30. + 2.*Box2Abb[4934]*m_z12 + 10.*Box2Abb[890]*m_z12_2 + 45.*Box2Abb[12]*m_z2k + 6.*m_z12_3*m_z2k;

  Box2Abb[4936]=9. + Box2Abb[1941]*m_z2k;

  Box2Abb[4937]=9. + Box2Abb[4936]*Box2Abb[61]*m_z2k;

  Box2Abb[4938]=-8. - Box2Abb[179]*Box2Abb[2396]*m_z2k;

  Box2Abb[4939]=1. + Box2Abb[4938]*m_z2k;

  Box2Abb[4940]=-41. + 9.*m_z2k + 5.*m_z2k_3;

  Box2Abb[4941]=36. + Box2Abb[4940]*m_z2k;

  Box2Abb[4942]=-9. + Box2Abb[4941]*m_z2k;

  Box2Abb[4943]=3.*pow(Box2Abb[61],5.) + 2.*Box2Abb[4937]*pow(Box2Abb[61],2.)*m_z12 + 2.*Box2Abb[4942]*m_z12_2 + 6.*Box2Abb[4939]*m_z12_3 - 2.*Box2Abb[230]*Box2Abb[4898]*m_z12_4*m_z2k + 2.*m_z12_5*m_z2k_2;

  Box2Abb[4944]=1. + m_z2k + m_z2k_2 + 2.*m_z2k_3;

  Box2Abb[4945]=10. - 18.*m_z2k + 11.*m_z2k_3;

  Box2Abb[4946]=9. + 2.*Box2Abb[3583]*m_z2k;

  Box2Abb[4947]=4. + Box2Abb[4946]*m_z2k;

  Box2Abb[4948]=-5. + Box2Abb[4947]*m_z2k;

  Box2Abb[4949]=30.*Box2Abb[4944] + 5.*Box2Abb[4948]*m_z12 + Box2Abb[4945]*m_z12_2 + 6.*Box2Abb[411]*m_z12_3*m_z2k;

  Box2Abb[4950]=Box2Abb[4924]*m_x_3 - Box2Abb[4911]*m_x_4 + Box2Abb[4949]*m_x_5 - Box2Abb[4935]*m_x_6 + Box2Abb[4932]*m_x_7 - Box2Abb[4926]*m_x_8 + m_x_9*m_z12 - Box2Abb[4943]*m_x_2*m_z2k + Box2Abb[4929]*m_x*m_z12*m_z2k_2 - Box2Abb[4904]*pow(Box2Abb[61],2.)*m_z12_2*m_z2k_3;

  Box2Abb[4951]=-67. - 11.*Box2Abb[201]*m_z12;

  Box2Abb[4952]=86. + 15.*m_z12;

  Box2Abb[4953]=-127. + Box2Abb[4952]*m_z12;

  Box2Abb[4954]=-16. + 17.*m_z12;

  Box2Abb[4955]=9. + Box2Abb[4954]*m_z12;

  Box2Abb[4956]=34. + Box2Abb[4951]*m_z12 + 66.*m_z2k + Box2Abb[4953]*m_z12*m_z2k + 7.*Box2Abb[4955]*m_z2k_2 + 84.*m_z12*m_z2k_3;

  Box2Abb[4957]=61. - 30.*m_z12;

  Box2Abb[4958]=-51. + Box2Abb[4957]*m_z12;

  Box2Abb[4959]=81. - 5.*m_z12;

  Box2Abb[4960]=-131. + Box2Abb[4959]*m_z12;

  Box2Abb[4961]=18. + Box2Abb[4960]*m_z12;

  Box2Abb[4962]=36. + 5.*m_z12;

  Box2Abb[4963]=-171. + 5.*Box2Abb[4962]*m_z12;

  Box2Abb[4964]=-20. + Box2Abb[4963]*m_z12;

  Box2Abb[4965]=30. + Box2Abb[4964]*m_z12;

  Box2Abb[4966]=-20. + 43.*m_z12;

  Box2Abb[4967]=9. + Box2Abb[4966]*m_z12;

  Box2Abb[4968]=4. + Box2Abb[4967]*m_z12;

  Box2Abb[4969]=19. + Box2Abb[4958]*m_z12 + 30.*m_z2k + Box2Abb[4961]*m_z12*m_z2k + Box2Abb[4965]*m_z2k_2 + 5.*Box2Abb[4968]*m_z2k_3 + 35.*Box2Abb[1479]*Box2Abb[4]*m_z2k_4 + 126.*m_z12*m_z2k_5;

  Box2Abb[4970]=-3. + 5.*m_z12 + 9.*m_z2k;

  Box2Abb[4971]=3. + Box2Abb[4970]*m_z12;

  Box2Abb[4972]=-7. + 9.*m_z2k;

  Box2Abb[4973]=22. + 37.*m_z2k;

  Box2Abb[4974]=-30. + Box2Abb[4973]*m_z12 - 2.*m_z12_2 + 4.*Box2Abb[4972]*m_z2k;

  Box2Abb[4975]=16. + Box2Abb[4974]*m_z12 + 21.*m_z2k;

  Box2Abb[4976]=3. + 2.*Box2Abb[841]*m_z2k;

  Box2Abb[4977]=-1. + Box2Abb[830]*m_z2k;

  Box2Abb[4978]=-3. + Box2Abb[2099]*m_z2k;

  Box2Abb[4979]=1. + Box2Abb[4978]*m_z2k;

  Box2Abb[4980]=13. + 5.*m_z2k;

  Box2Abb[4981]=-3. + Box2Abb[4980]*m_z2k;

  Box2Abb[4982]=1. + Box2Abb[4981]*m_z2k;

  Box2Abb[4983]=-Box2Abb[4977]*pow(Box2Abb[61],4.) - Box2Abb[4976]*pow(Box2Abb[61],4.)*m_z12 + 3.*Box2Abb[4979]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[4982]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[187]*m_z12_4*m_z2k_2;

  Box2Abb[4984]=-5. + 17.*m_z2k;

  Box2Abb[4985]=14. + 3.*Box2Abb[1567]*m_z2k;

  Box2Abb[4986]=7. + Box2Abb[4985]*m_z2k;

  Box2Abb[4987]=12. + 31.*m_z2k;

  Box2Abb[4988]=9. + 7.*Box2Abb[4987]*m_z2k;

  Box2Abb[4989]=59. + Box2Abb[4988]*m_z2k;

  Box2Abb[4990]=-9. + 7.*Box2Abb[179]*m_z2k;

  Box2Abb[4991]=-121. + 18.*Box2Abb[4990]*m_z2k;

  Box2Abb[4992]=-72. + Box2Abb[4991]*m_z2k;

  Box2Abb[4993]=5.*Box2Abb[4986] + Box2Abb[4992]*m_z12 + Box2Abb[4989]*m_z12_2 + Box2Abb[1916]*Box2Abb[4984]*m_z12_3 - m_z12_4*m_z2k;

  Box2Abb[4994]=-1. + m_z2k + 3.*m_z2k_2;

  Box2Abb[4995]=-10. + 13.*m_z2k;

  Box2Abb[4996]=-24. + Box2Abb[4995]*m_z2k;

  Box2Abb[4997]=22. + Box2Abb[4996]*m_z2k;

  Box2Abb[4998]=-5. + Box2Abb[4997]*m_z2k;

  Box2Abb[4999]=-25. + 9.*m_z2k;

  Box2Abb[5000]=36. + Box2Abb[4999]*m_z2k;

  Box2Abb[5001]=-5. + Box2Abb[5000]*m_z2k;

  Box2Abb[5002]=-8. + Box2Abb[5001]*m_z2k;

  Box2Abb[5003]=3. + Box2Abb[5002]*m_z2k;

  Box2Abb[5004]=-59. + 10.*m_z2k;

  Box2Abb[5005]=28. + Box2Abb[5004]*m_z2k;

  Box2Abb[5006]=46. + Box2Abb[5005]*m_z2k;

  Box2Abb[5007]=-14. + Box2Abb[5006]*m_z2k;

  Box2Abb[5008]=1. + Box2Abb[5007]*m_z2k;

  Box2Abb[5009]=-66. + 19.*m_z2k;

  Box2Abb[5010]=74. + Box2Abb[5009]*m_z2k;

  Box2Abb[5011]=21. + Box2Abb[5010]*m_z2k;

  Box2Abb[5012]=-15. + Box2Abb[5011]*m_z2k;

  Box2Abb[5013]=3. + Box2Abb[5012]*m_z2k;

  Box2Abb[5014]=Box2Abb[4994]*pow(Box2Abb[61],5.) + Box2Abb[5003]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[5013]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[5008]*Box2Abb[61]*m_z12_3 + Box2Abb[4998]*m_z12_4*m_z2k + 2.*Box2Abb[1941]*m_z12_5*m_z2k_3;

  Box2Abb[5015]=10. - 60.*m_z2k + 63.*m_z2k_2;

  Box2Abb[5016]=31. + 38.*m_z2k;

  Box2Abb[5017]=-5. + Box2Abb[5016]*m_z2k;

  Box2Abb[5018]=47. + 105.*m_z2k;

  Box2Abb[5019]=9. + Box2Abb[5018]*m_z2k;

  Box2Abb[5020]=33. + Box2Abb[5019]*m_z2k;

  Box2Abb[5021]=-10. + Box2Abb[5020]*m_z2k;

  Box2Abb[5022]=5. + 14.*Box2Abb[179]*m_z2k;

  Box2Abb[5023]=26. + Box2Abb[2011]*Box2Abb[5022]*m_z2k;

  Box2Abb[5024]=47. + Box2Abb[5023]*m_z2k;

  Box2Abb[5025]=-32. + 2.*Box2Abb[5024]*m_z2k;

  Box2Abb[5026]=-58. + 35.*m_z2k;

  Box2Abb[5027]=-74. + 5.*Box2Abb[5026]*m_z2k;

  Box2Abb[5028]=-130. + Box2Abb[5027]*m_z2k;

  Box2Abb[5029]=-145. + Box2Abb[5028]*m_z2k;

  Box2Abb[5030]=44. + Box2Abb[5029]*m_z2k;

  Box2Abb[5031]=8. + Box2Abb[5025]*m_z12 + Box2Abb[5030]*m_z12_2 + 2.*Box2Abb[5021]*m_z12_3 - 9.*m_z2k + 2.*Box2Abb[5017]*m_z12_4*m_z2k + Box2Abb[5015]*m_z2k_3;

  Box2Abb[5032]=27. + 30.*m_z2k + 34.*m_z2k_2;

  Box2Abb[5033]=-5. + Box2Abb[5032]*m_z2k;

  Box2Abb[5034]=4. + 3.*Box2Abb[1580]*m_z2k;

  Box2Abb[5035]=-4. + Box2Abb[5034]*m_z2k;

  Box2Abb[5036]=-83. + 93.*m_z2k;

  Box2Abb[5037]=-46. + Box2Abb[5036]*m_z2k;

  Box2Abb[5038]=-102. + Box2Abb[5037]*m_z2k;

  Box2Abb[5039]=37. + Box2Abb[5038]*m_z2k;

  Box2Abb[5040]=-7. + Box2Abb[5039]*m_z2k;

  Box2Abb[5041]=-11. + 3.*m_z2k;

  Box2Abb[5042]=61. + 4.*Box2Abb[5041]*m_z2k;

  Box2Abb[5043]=-52. + 3.*Box2Abb[5042]*m_z2k;

  Box2Abb[5044]=-10. + Box2Abb[5043]*m_z2k;

  Box2Abb[5045]=-46. + Box2Abb[5044]*m_z2k;

  Box2Abb[5046]=15. + Box2Abb[5045]*m_z2k;

  Box2Abb[5047]=-258. + 77.*m_z2k;

  Box2Abb[5048]=234. + Box2Abb[5047]*m_z2k;

  Box2Abb[5049]=-38. + Box2Abb[5048]*m_z2k;

  Box2Abb[5050]=63. + Box2Abb[5049]*m_z2k;

  Box2Abb[5051]=-72. + Box2Abb[5050]*m_z2k;

  Box2Abb[5052]=18. + Box2Abb[5051]*m_z2k;

  Box2Abb[5053]=Box2Abb[5035]*pow(Box2Abb[61],3.) + Box2Abb[5046]*Box2Abb[61]*m_z12 + Box2Abb[5052]*m_z12_2 + Box2Abb[5040]*m_z12_3 + 2.*Box2Abb[5033]*m_z12_4*m_z2k + 10.*m_z12_5*m_z2k_3;

  Box2Abb[5054]=Box2Abb[5014]*m_x - Box2Abb[5053]*m_x_2 + Box2Abb[5031]*m_x_3 - Box2Abb[4969]*m_x_4 + Box2Abb[4993]*m_x_5 - Box2Abb[4956]*m_x_6 + Box2Abb[4975]*m_x_7 - Box2Abb[4971]*m_x_8 + m_x_9*m_z12 + Box2Abb[4983]*Box2Abb[61]*m_z12*m_z2k;

  Box2Abb[5055]=-Box2Abb[4950]*Box2Abb[8]*m_cL + Box2Abb[5054]*m_cR*m_x;

  Box2Abb[5056]=-2. + 3.*m_z12 + 4.*m_z2k;

  Box2Abb[5057]=3. + 2.*Box2Abb[5056]*m_z12;

  Box2Abb[5058]=1. + 8.*m_z2k;

  Box2Abb[5059]=1. - 2.*m_z2k + 8.*m_z2k_3 - 7.*m_z2k_4;

  Box2Abb[5060]=1. + Box2Abb[797]*m_z2k_2;

  Box2Abb[5061]=-Box2Abb[514]*pow(Box2Abb[61],5.) - 4.*Box2Abb[230]*pow(Box2Abb[61],6.)*m_z12 - 3.*Box2Abb[5058]*pow(Box2Abb[61],5.)*m_z12_2 - 10.*Box2Abb[5060]*pow(Box2Abb[61],2.)*m_z12_3 + 5.*Box2Abb[5059]*m_z12_4 - 12.*m_z12_5*m_z2k_3;

  Box2Abb[5062]=6. + 13.*m_z2k;

  Box2Abb[5063]=-26. + 3.*Box2Abb[5062]*m_z12 - 5.*m_z12_2 + 4.*Box2Abb[2136]*m_z2k;

  Box2Abb[5064]=13. + Box2Abb[5063]*m_z12 + 18.*m_z2k;

  Box2Abb[5065]=22. - 8.*m_z2k;

  Box2Abb[5066]=7. + 9.*m_z2k;

  Box2Abb[5067]=4. + Box2Abb[5066]*m_z2k;

  Box2Abb[5068]=13. + 36.*m_z2k;

  Box2Abb[5069]=11. + Box2Abb[5068]*m_z2k;

  Box2Abb[5070]=27. - 14.*m_z2k;

  Box2Abb[5071]=18. + Box2Abb[5070]*m_z2k;

  Box2Abb[5072]=9. + Box2Abb[5071]*m_z2k;

  Box2Abb[5073]=-5.*Box2Abb[5067] + 4.*Box2Abb[5072]*m_z12 - 3.*Box2Abb[5069]*m_z12_2 + Box2Abb[5065]*m_z12_3 + m_z12_4;

  Box2Abb[5074]=1. + m_z2k + m_z2k_2 + 6.*m_z2k_3;

  Box2Abb[5075]=-20. + 7.*m_z2k;

  Box2Abb[5076]=-3. + Box2Abb[5075]*m_z2k_3;

  Box2Abb[5077]=22. + 85.*m_z2k;

  Box2Abb[5078]=-35. + Box2Abb[5077]*m_z2k;

  Box2Abb[5079]=-3. + 11.*m_z2k;

  Box2Abb[5080]=-3. + Box2Abb[5079]*m_z2k;

  Box2Abb[5081]=19. + 5.*Box2Abb[5080]*m_z2k;

  Box2Abb[5082]=10.*Box2Abb[5074] + 10.*Box2Abb[5076]*m_z12 + 3.*Box2Abb[5081]*m_z12_2 + Box2Abb[5078]*m_z12_3 - 5.*Box2Abb[12]*m_z12_4;

  Box2Abb[5083]=10. - 9.*m_z2k;

  Box2Abb[5084]=5. - 13.*m_z2k;

  Box2Abb[5085]=5. + Box2Abb[5084]*m_z2k;

  Box2Abb[5086]=-5. + 2.*m_z2k;

  Box2Abb[5087]=-10. + 7.*m_z2k;

  Box2Abb[5088]=-9. + 15.*m_z2k + Box2Abb[5086]*Box2Abb[5087]*m_z2k_3;

  Box2Abb[5089]=-1. + Box2Abb[680]*m_z2k;

  Box2Abb[5090]=9. + 40.*m_z2k;

  Box2Abb[5091]=9. + Box2Abb[5090]*m_z2k;

  Box2Abb[5092]=-5. + Box2Abb[5091]*m_z2k;

  Box2Abb[5093]=-5. - 4.*Box2Abb[5088]*m_z12 - 30.*Box2Abb[1069]*Box2Abb[5089]*m_z12_2 - 4.*Box2Abb[5092]*m_z12_3 + 2.*Box2Abb[5085]*m_z12_4 + 5.*Box2Abb[5083]*m_z2k_3;

  Box2Abb[5094]=-13. + 2.*Box2Abb[2136]*m_z2k;

  Box2Abb[5095]=-5. + 2.*Box2Abb[3911]*m_z2k_2;

  Box2Abb[5096]=-1. + 18.*m_z2k;

  Box2Abb[5097]=-7. + Box2Abb[5096]*m_z2k;

  Box2Abb[5098]=-18. + 5.*Box2Abb[258]*m_z2k;

  Box2Abb[5099]=4. + Box2Abb[5098]*m_z2k;

  Box2Abb[5100]=1. + Box2Abb[5099]*m_z2k;

  Box2Abb[5101]=-53. + 27.*m_z2k;

  Box2Abb[5102]=-3. + Box2Abb[5101]*m_z2k;

  Box2Abb[5103]=17. + Box2Abb[5102]*m_z2k;

  Box2Abb[5104]=-8. + Box2Abb[5103]*m_z2k;

  Box2Abb[5105]=Box2Abb[5097]*pow(Box2Abb[61],3.) + 2.*Box2Abb[5094]*pow(Box2Abb[61],4.)*m_z12 + 3.*Box2Abb[5104]*Box2Abb[61]*m_z12_2 + 5.*Box2Abb[5100]*m_z12_3 + 2.*Box2Abb[5095]*m_z12_4;

  Box2Abb[5106]=Box2Abb[5061]*m_x + Box2Abb[5105]*m_x_2 + Box2Abb[5093]*m_x_3 + Box2Abb[5082]*m_x_4 + Box2Abb[5073]*m_x_5 + Box2Abb[5064]*m_x_6 - Box2Abb[5057]*m_x_7 + pow(Box2Abb[5],3.)*pow(Box2Abb[61],5.)*m_z12 + m_x_8*m_z12;

  Box2Abb[5107]=3. + m_z12 + 2.*m_z12_2 + 9.*m_z12*m_z2k;

  Box2Abb[5108]=-1. + m_z2k + 9.*m_z2k_2;

  Box2Abb[5109]=3.*Box2Abb[719] + 4.*Box2Abb[5108]*m_z12 + 2.*Box2Abb[1580]*m_z12_2 - m_z12_3;

  Box2Abb[5110]=-18. + 13.*Box2Abb[12]*m_z2k;

  Box2Abb[5111]=8. + m_z2k;

  Box2Abb[5112]=-3. + Box2Abb[5111]*m_z2k;

  Box2Abb[5113]=-15. + Box2Abb[1567]*m_z2k;

  Box2Abb[5114]=-Box2Abb[222]*pow(Box2Abb[61],5.) - Box2Abb[5113]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[5110]*pow(Box2Abb[61],2.)*m_z12_2 - 2.*Box2Abb[5112]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[460]*m_z12_4*m_z2k;

  Box2Abb[5115]=m_z2k + 28.*m_z2k_3;

  Box2Abb[5116]=20. + 43.*m_z2k;

  Box2Abb[5117]=6. + Box2Abb[5116]*m_z2k;

  Box2Abb[5118]=15. + 3.*Box2Abb[5115]*m_z12 + Box2Abb[5117]*m_z12_2 + Box2Abb[4847]*m_z12_3 + 21.*Box2Abb[150]*m_z2k;

  Box2Abb[5119]=1. + Box2Abb[3116]*m_z2k;

  Box2Abb[5120]=7. + 19.*m_z2k;

  Box2Abb[5121]=-5. + Box2Abb[5120]*m_z2k;

  Box2Abb[5122]=18. + 77.*m_z2k;

  Box2Abb[5123]=-18. + Box2Abb[5122]*m_z2k;

  Box2Abb[5124]=4. + Box2Abb[5123]*m_z2k;

  Box2Abb[5125]=-2. + 9.*m_z2k;

  Box2Abb[5126]=36. + 7.*Box2Abb[5125]*m_z2k;

  Box2Abb[5127]=45. + 2.*Box2Abb[5126]*m_z2k;

  Box2Abb[5128]=19. + Box2Abb[5127]*m_z2k;

  Box2Abb[5129]=Box2Abb[5128]*m_z12 + Box2Abb[5124]*m_z12_2 + 2.*Box2Abb[5121]*m_z12_3 + 15.*Box2Abb[5119]*m_z2k;

  Box2Abb[5130]=16. - 9.*m_z2k;

  Box2Abb[5131]=11. + Box2Abb[5130]*m_z2k;

  Box2Abb[5132]=-12. + Box2Abb[5131]*m_z2k;

  Box2Abb[5133]=19. + m_z2k;

  Box2Abb[5134]=16. + Box2Abb[5133]*m_z2k;

  Box2Abb[5135]=-33. + Box2Abb[5134]*m_z2k;

  Box2Abb[5136]=3. + Box2Abb[5135]*m_z2k;

  Box2Abb[5137]=2. + Box2Abb[5066]*m_z2k;

  Box2Abb[5138]=-21. + Box2Abb[5137]*m_z2k;

  Box2Abb[5139]=15. + Box2Abb[5138]*m_z2k;

  Box2Abb[5140]=16. + 23.*m_z2k;

  Box2Abb[5141]=9. + Box2Abb[5140]*m_z2k;

  Box2Abb[5142]=-54. + Box2Abb[5141]*m_z2k;

  Box2Abb[5143]=18. + Box2Abb[5142]*m_z2k;

  Box2Abb[5144]=3.*pow(Box2Abb[61],6.) + Box2Abb[5139]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[5143]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[5136]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[5132]*m_z12_4*m_z2k - 2.*m_z12_5*m_z2k_3;

  Box2Abb[5145]=36. + 55.*m_z2k;

  Box2Abb[5146]=6. - Box2Abb[5145]*m_z2k;

  Box2Abb[5147]=10. + Box2Abb[5146]*m_z2k;

  Box2Abb[5148]=6. - 91.*m_z2k;

  Box2Abb[5149]=9. + Box2Abb[5148]*m_z2k;

  Box2Abb[5150]=16. + Box2Abb[5149]*m_z2k;

  Box2Abb[5151]=-6. + Box2Abb[5150]*m_z2k;

  Box2Abb[5152]=-5. + 9.*m_z2k;

  Box2Abb[5153]=25. + 2.*Box2Abb[5152]*m_z2k;

  Box2Abb[5154]=45. + 7.*Box2Abb[5153]*m_z2k;

  Box2Abb[5155]=55. + Box2Abb[5154]*m_z2k;

  Box2Abb[5156]=23. + Box2Abb[5155]*m_z2k;

  Box2Abb[5157]=13. - Box2Abb[5156]*m_z12 + Box2Abb[5151]*m_z12_2 + Box2Abb[5147]*m_z12_3 - 12.*m_z12_4*m_z2k_2 + 15.*Box2Abb[2001]*m_z2k_3;

  Box2Abb[5158]=-40. + 21.*m_z2k;

  Box2Abb[5159]=20. + Box2Abb[5158]*m_z2k;

  Box2Abb[5160]=15. + 7.*Box2Abb[61]*m_z2k;

  Box2Abb[5161]=-55. + 6.*Box2Abb[5160]*m_z2k;

  Box2Abb[5162]=4. + Box2Abb[5161]*m_z2k_3;

  Box2Abb[5163]=23. + 20.*m_z2k;

  Box2Abb[5164]=-18. + Box2Abb[5163]*m_z2k;

  Box2Abb[5165]=-32. + Box2Abb[5164]*m_z2k;

  Box2Abb[5166]=-5. + Box2Abb[5165]*m_z2k;

  Box2Abb[5167]=-20. + 77.*m_z2k;

  Box2Abb[5168]=99. + Box2Abb[5167]*m_z2k;

  Box2Abb[5169]=16. + Box2Abb[5168]*m_z2k;

  Box2Abb[5170]=50. + Box2Abb[5169]*m_z2k;

  Box2Abb[5171]=6. + Box2Abb[5170]*m_z2k;

  Box2Abb[5172]=-8. + 2.*Box2Abb[5162]*m_z12 + Box2Abb[5171]*m_z12_2 + Box2Abb[5166]*m_z12_3 + m_z2k + 4.*Box2Abb[559]*m_z12_4*m_z2k_2 + 3.*Box2Abb[5159]*m_z2k_3;

  Box2Abb[5173]=9. + Box2Abb[2800]*m_z2k;

  Box2Abb[5174]=-12. + 7.*m_z2k;

  Box2Abb[5175]=4. + Box2Abb[5174]*m_z2k;

  Box2Abb[5176]=-1. + 3.*Box2Abb[5175]*m_z2k_2;

  Box2Abb[5177]=61. + 4.*Box2Abb[5152]*m_z2k;

  Box2Abb[5178]=-89. + Box2Abb[5177]*m_z2k;

  Box2Abb[5179]=-3. + Box2Abb[5178]*m_z2k;

  Box2Abb[5180]=33. + Box2Abb[5179]*m_z2k;

  Box2Abb[5181]=-43. + 11.*m_z2k;

  Box2Abb[5182]=-40. + Box2Abb[5181]*m_z2k;

  Box2Abb[5183]=-38. + Box2Abb[5182]*m_z2k;

  Box2Abb[5184]=25. + Box2Abb[5183]*m_z2k;

  Box2Abb[5185]=1. + Box2Abb[5184]*m_z2k;

  Box2Abb[5186]=24. - 49.*m_z2k;

  Box2Abb[5187]=-96. + Box2Abb[5186]*m_z2k;

  Box2Abb[5188]=108. + Box2Abb[5187]*m_z2k;

  Box2Abb[5189]=51. + Box2Abb[5188]*m_z2k;

  Box2Abb[5190]=-60. + Box2Abb[5189]*m_z2k;

  Box2Abb[5191]=-2. + Box2Abb[5190]*m_z2k;

  Box2Abb[5192]=-Box2Abb[5176]*pow(Box2Abb[61],2.) + Box2Abb[5191]*m_z12_2 + Box2Abb[5185]*m_z12_3 - Box2Abb[5180]*Box2Abb[61]*m_z12*m_z2k + 2.*Box2Abb[5173]*m_z12_4*m_z2k_2 - 4.*m_z12_5*m_z2k_3;

  Box2Abb[5193]=Box2Abb[5192]*m_x_3 + Box2Abb[5172]*m_x_4 + Box2Abb[5157]*m_x_5 + Box2Abb[5129]*m_x_6 - Box2Abb[5118]*m_x_7 + Box2Abb[5109]*m_x_8 - Box2Abb[5107]*m_x_9 + m_x_10*m_z12 + Box2Abb[5144]*m_x_2*m_z2k + Box2Abb[5114]*Box2Abb[61]*m_x*m_z12*m_z2k_2 + Box2Abb[5]*Box2Abb[519]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_3;

  Box2Abb[5194]=Box2Abb[5193]*m_cL - Box2Abb[5106]*Box2Abb[8]*m_cR*m_x;

  Box2Abb[5195]=2. + 15.*Box2Abb[62]*m_z12;

  Box2Abb[5196]=16. + Box2Abb[5195]*m_z12;

  Box2Abb[5197]=8. + Box2Abb[533]*m_z12;

  Box2Abb[5198]=-119. + 5.*Box2Abb[5197]*m_z12;

  Box2Abb[5199]=17. + Box2Abb[5198]*m_z12;

  Box2Abb[5200]=-94. + Box2Abb[2695]*Box2Abb[734]*m_z12;

  Box2Abb[5201]=48. - Box2Abb[5200]*m_z12;

  Box2Abb[5202]=-261. + 46.*m_z12;

  Box2Abb[5203]=140. + Box2Abb[5202]*m_z12;

  Box2Abb[5204]=-95. + Box2Abb[5203]*m_z12;

  Box2Abb[5205]=20. + Box2Abb[5204]*m_z12;

  Box2Abb[5206]=-197. + 78.*m_z12;

  Box2Abb[5207]=162. + Box2Abb[5206]*m_z12;

  Box2Abb[5208]=-29. + Box2Abb[5207]*m_z12;

  Box2Abb[5209]=-56. + 39.*m_z12;

  Box2Abb[5210]=12. + Box2Abb[5209]*m_z12;

  Box2Abb[5211]=3. - Box2Abb[5196]*m_z12 - 12.*m_z2k - Box2Abb[5199]*m_z12*m_z2k + Box2Abb[5201]*m_z12*m_z2k_2 + Box2Abb[5205]*m_z2k_3 + 5.*Box2Abb[5208]*m_z2k_4 + 14.*Box2Abb[5210]*m_z2k_5 + 210.*m_z12*m_z2k_6;

  Box2Abb[5212]=-13. + 20.*m_z12;

  Box2Abb[5213]=-103. + 42.*m_z12 + 10.*m_z12_3;

  Box2Abb[5214]=50. + Box2Abb[5213]*m_z12;

  Box2Abb[5215]=-8. + Box2Abb[273]*m_z12;

  Box2Abb[5216]=-14. + Box2Abb[5215]*m_z12;

  Box2Abb[5217]=-3. + Box2Abb[5216]*m_z12;

  Box2Abb[5218]=-2. + Box2Abb[5217]*m_z12;

  Box2Abb[5219]=-32. + m_z12;

  Box2Abb[5220]=6. + Box2Abb[5219]*m_z12;

  Box2Abb[5221]=-66. + Box2Abb[5220]*m_z12;

  Box2Abb[5222]=71. + Box2Abb[5221]*m_z12;

  Box2Abb[5223]=-10. + Box2Abb[5222]*m_z12;

  Box2Abb[5224]=509. - 114.*m_z12;

  Box2Abb[5225]=-850. + Box2Abb[5224]*m_z12;

  Box2Abb[5226]=655. + Box2Abb[5225]*m_z12;

  Box2Abb[5227]=-110. + Box2Abb[5226]*m_z12;

  Box2Abb[5228]=1031. - 338.*m_z12;

  Box2Abb[5229]=-996. + Box2Abb[5228]*m_z12;

  Box2Abb[5230]=177. + Box2Abb[5229]*m_z12;

  Box2Abb[5231]=pow(Box2Abb[4],2.)*Box2Abb[5212]*m_z12 + m_z2k + Box2Abb[5214]*m_z12*m_z2k + 2.*Box2Abb[5218]*m_z2k_2 - 2.*Box2Abb[5223]*m_z2k_3 + Box2Abb[5227]*m_z2k_4 + Box2Abb[5230]*m_z2k_5 - 28.*Box2Abb[274]*Box2Abb[62]*m_z2k_6 - 120.*m_z12*m_z2k_7;

  Box2Abb[5232]=-2. + 3.*m_z12 + 5.*m_z2k;

  Box2Abb[5233]=3. + 2.*Box2Abb[5232]*m_z12;

  Box2Abb[5234]=14. + 51.*m_z2k;

  Box2Abb[5235]=-24. + Box2Abb[5234]*m_z12 - 7.*m_z12_2 + 5.*Box2Abb[4547]*m_z2k;

  Box2Abb[5236]=13. + Box2Abb[5235]*m_z12 + 24.*m_z2k;

  Box2Abb[5237]=31. + 18.*m_z2k;

  Box2Abb[5238]=61. + 84.*m_z2k;

  Box2Abb[5239]=43. + 192.*m_z2k;

  Box2Abb[5240]=22. + Box2Abb[5239]*m_z2k;

  Box2Abb[5241]=27. + 44.*m_z2k - 30.*m_z2k_2;

  Box2Abb[5242]=21. + 4.*Box2Abb[5241]*m_z2k;

  Box2Abb[5243]=-18. + Box2Abb[5242]*m_z12 - Box2Abb[5240]*m_z12_2 + Box2Abb[5237]*m_z12_3 + 2.*m_z12_4 - Box2Abb[5238]*m_z2k;

  Box2Abb[5244]=11. + m_z2k;

  Box2Abb[5245]=-9. + Box2Abb[5244]*m_z2k;

  Box2Abb[5246]=3. + Box2Abb[5245]*m_z2k;

  Box2Abb[5247]=-Box2Abb[4977]*pow(Box2Abb[61],3.) - Box2Abb[4976]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[5246]*Box2Abb[61]*m_z12_2 + Box2Abb[4979]*m_z12_3;

  Box2Abb[5248]=11. + 16.*m_z2k;

  Box2Abb[5249]=-81. + 40.*m_z2k;

  Box2Abb[5250]=-46. + Box2Abb[5249]*m_z2k;

  Box2Abb[5251]=31. + 56.*m_z2k;

  Box2Abb[5252]=38. + 3.*Box2Abb[5251]*m_z2k;

  Box2Abb[5253]=-61. + 420.*m_z2k;

  Box2Abb[5254]=-16. + Box2Abb[5253]*m_z2k;

  Box2Abb[5255]=50. + Box2Abb[5254]*m_z2k;

  Box2Abb[5256]=-32. + 15.*m_z2k;

  Box2Abb[5257]=-57. + 7.*Box2Abb[5256]*m_z2k;

  Box2Abb[5258]=5. + 2.*Box2Abb[5257]*m_z2k;

  Box2Abb[5259]=-1. + Box2Abb[5258]*m_z2k;

  Box2Abb[5260]=1. + Box2Abb[5259]*m_z12 + Box2Abb[5255]*m_z12_2 + Box2Abb[5250]*m_z12_3 - Box2Abb[5248]*m_z12_4 + Box2Abb[5252]*m_z2k;

  Box2Abb[5261]=25. + 3.*Box2Abb[2664]*m_z2k;

  Box2Abb[5262]=107. - 226.*m_z2k;

  Box2Abb[5263]=98. + Box2Abb[5262]*m_z2k;

  Box2Abb[5264]=10. + Box2Abb[5263]*m_z2k;

  Box2Abb[5265]=2. + m_z2k + 42.*m_z2k_2;

  Box2Abb[5266]=-1. + Box2Abb[5265]*m_z2k;

  Box2Abb[5267]=485. - 588.*m_z2k;

  Box2Abb[5268]=138. + Box2Abb[5267]*m_z2k;

  Box2Abb[5269]=-68. + Box2Abb[5268]*m_z2k;

  Box2Abb[5270]=-56. + Box2Abb[5269]*m_z2k;

  Box2Abb[5271]=-26. + 9.*m_z2k;

  Box2Abb[5272]=60. + 7.*Box2Abb[5271]*m_z2k;

  Box2Abb[5273]=111. + 4.*Box2Abb[5272]*m_z2k;

  Box2Abb[5274]=53. + Box2Abb[5273]*m_z2k;

  Box2Abb[5275]=-16. + Box2Abb[5274]*m_z2k;

  Box2Abb[5276]=6. - Box2Abb[5275]*m_z12 + Box2Abb[5270]*m_z12_2 + Box2Abb[5264]*m_z12_3 + Box2Abb[5261]*m_z12_4 - 5.*Box2Abb[5266]*m_z2k + m_z12_5*m_z2k;

  Box2Abb[5277]=6. - 7.*m_z2k;

  Box2Abb[5278]=10. + Box2Abb[5277]*m_z2k;

  Box2Abb[5279]=-10. + Box2Abb[5278]*m_z2k;

  Box2Abb[5280]=5. + Box2Abb[5279]*m_z2k;

  Box2Abb[5281]=22. + Box2Abb[4665]*m_z2k;

  Box2Abb[5282]=-11. + 2.*Box2Abb[5281]*m_z2k;

  Box2Abb[5283]=-9. + Box2Abb[5282]*m_z2k;

  Box2Abb[5284]=4. + Box2Abb[5283]*m_z2k;

  Box2Abb[5285]=-25. + 6.*m_z2k;

  Box2Abb[5286]=22. + Box2Abb[5285]*m_z2k;

  Box2Abb[5287]=19. + 5.*Box2Abb[5286]*m_z2k;

  Box2Abb[5288]=-26. + Box2Abb[5287]*m_z2k;

  Box2Abb[5289]=4. + Box2Abb[5288]*m_z2k;

  Box2Abb[5290]=-59. + 18.*m_z2k;

  Box2Abb[5291]=20. + Box2Abb[5290]*m_z2k;

  Box2Abb[5292]=28. + Box2Abb[5291]*m_z2k;

  Box2Abb[5293]=-18. + Box2Abb[5292]*m_z2k;

  Box2Abb[5294]=1. + Box2Abb[5293]*m_z2k;

  Box2Abb[5295]=-101. + 30.*m_z2k;

  Box2Abb[5296]=121. + Box2Abb[5295]*m_z2k;

  Box2Abb[5297]=-8. + Box2Abb[5296]*m_z2k;

  Box2Abb[5298]=-20. + Box2Abb[5297]*m_z2k;

  Box2Abb[5299]=6. + Box2Abb[5298]*m_z2k;

  Box2Abb[5300]=-Box2Abb[4994]*pow(Box2Abb[61],6.) - Box2Abb[5284]*pow(Box2Abb[61],4.)*m_z12 - Box2Abb[5299]*pow(Box2Abb[61],3.)*m_z12_2 - Box2Abb[5289]*pow(Box2Abb[61],2.)*m_z12_3 - Box2Abb[5294]*Box2Abb[61]*m_z12_4 + Box2Abb[5280]*m_z12_5*m_z2k;

  Box2Abb[5301]=1. + Box2Abb[150]*m_z2k;

  Box2Abb[5302]=-5. + Box2Abb[3910]*m_z2k_2;

  Box2Abb[5303]=-113. + 84.*m_z2k;

  Box2Abb[5304]=-4. + Box2Abb[5303]*m_z2k;

  Box2Abb[5305]=-36. + Box2Abb[5304]*m_z2k;

  Box2Abb[5306]=30. + Box2Abb[5305]*m_z2k;

  Box2Abb[5307]=-7. + Box2Abb[5306]*m_z2k;

  Box2Abb[5308]=-315. + 152.*m_z2k;

  Box2Abb[5309]=79. + Box2Abb[5308]*m_z2k;

  Box2Abb[5310]=-5. + Box2Abb[5309]*m_z2k;

  Box2Abb[5311]=65. + Box2Abb[5310]*m_z2k;

  Box2Abb[5312]=24. + Box2Abb[5311]*Box2Abb[61]*m_z2k;

  Box2Abb[5313]=-166. + 45.*m_z2k;

  Box2Abb[5314]=265. + Box2Abb[5313]*m_z2k;

  Box2Abb[5315]=-105. + Box2Abb[5314]*m_z2k;

  Box2Abb[5316]=-12. + Box2Abb[5315]*m_z2k;

  Box2Abb[5317]=-31. + Box2Abb[5316]*m_z2k;

  Box2Abb[5318]=16. + Box2Abb[5317]*m_z2k;

  Box2Abb[5319]=-475. + 132.*m_z2k;

  Box2Abb[5320]=557. + Box2Abb[5319]*m_z2k;

  Box2Abb[5321]=-161. + Box2Abb[5320]*m_z2k;

  Box2Abb[5322]=17. + Box2Abb[5321]*m_z2k;

  Box2Abb[5323]=-64. + Box2Abb[5322]*m_z2k;

  Box2Abb[5324]=30. + Box2Abb[5323]*m_z2k;

  Box2Abb[5325]=Box2Abb[1163]*Box2Abb[5301]*pow(Box2Abb[61],4.) + Box2Abb[5318]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[5324]*Box2Abb[61]*m_z12_2 + Box2Abb[5312]*m_z12_3 + Box2Abb[5307]*m_z12_4 + 2.*Box2Abb[5302]*m_z12_5*m_z2k;

  Box2Abb[5326]=Box2Abb[5300]*m_x + Box2Abb[5325]*m_x_2 + Box2Abb[5231]*m_x_3 + Box2Abb[5211]*m_x_4 + Box2Abb[5276]*m_x_5 + Box2Abb[5260]*m_x_6 + Box2Abb[5243]*m_x_7 + Box2Abb[5236]*m_x_8 - Box2Abb[5233]*m_x_9 + m_x_10*m_z12 - Box2Abb[5]*Box2Abb[5247]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[5327]=-Box2Abb[5193]*Box2Abb[8]*m_cL + Box2Abb[5326]*m_cR*m_x;

  Box2Abb[5328]=-12. + Box2Abb[734]*m_z12;

  Box2Abb[5329]=-9. + 5.*Box2Abb[5328]*m_z12;

  Box2Abb[5330]=32. + Box2Abb[5329]*m_z12;

  Box2Abb[5331]=-6. + 5.*Box2Abb[735]*m_z12;

  Box2Abb[5332]=-84. + Box2Abb[5331]*m_z12;

  Box2Abb[5333]=32. + Box2Abb[5332]*m_z12;

  Box2Abb[5334]=72. + 19.*m_z12;

  Box2Abb[5335]=-63. + Box2Abb[5334]*m_z12;

  Box2Abb[5336]=-5. + Box2Abb[5335]*m_z12;

  Box2Abb[5337]=62. - 17.*m_z12;

  Box2Abb[5338]=-46. + Box2Abb[5337]*m_z12;

  Box2Abb[5339]=10. + Box2Abb[5338]*m_z12;

  Box2Abb[5340]=-7. + Box2Abb[5330]*m_z12 + Box2Abb[5333]*m_z12*m_z2k + 2.*Box2Abb[5336]*m_z12*m_z2k_2 + 10.*Box2Abb[5339]*m_z2k_3 - 35.*Box2Abb[4]*Box2Abb[819]*m_z2k_4 - 126.*m_z12*m_z2k_5;

  Box2Abb[5341]=16. + Box2Abb[2678]*m_z12;

  Box2Abb[5342]=-60. + Box2Abb[5341]*m_z12;

  Box2Abb[5343]=10. + Box2Abb[5342]*m_z12;

  Box2Abb[5344]=77. + 34.*m_z12;

  Box2Abb[5345]=21. - Box2Abb[5344]*m_z12;

  Box2Abb[5346]=71. + Box2Abb[5345]*m_z12;

  Box2Abb[5347]=-74. + 7.*m_z12;

  Box2Abb[5348]=27. + Box2Abb[5347]*m_z12;

  Box2Abb[5349]=-52. + 47.*m_z12;

  Box2Abb[5350]=15. + Box2Abb[5349]*m_z12;

  Box2Abb[5351]=-5.*Box2Abb[187] - Box2Abb[5343]*m_z12 + Box2Abb[5346]*m_z12*m_z2k + 3.*Box2Abb[5348]*m_z12*m_z2k_2 + 7.*Box2Abb[5350]*m_z2k_3 + 126.*m_z12*m_z2k_4;

  Box2Abb[5352]=-5. + 7.*m_z12 + 9.*m_z2k;

  Box2Abb[5353]=3. + Box2Abb[5352]*m_z12;

  Box2Abb[5354]=8. + 53.*m_z2k;

  Box2Abb[5355]=-19. + Box2Abb[5354]*m_z12 - 11.*m_z12_2 + 4.*Box2Abb[1582]*m_z2k;

  Box2Abb[5356]=10. + Box2Abb[5355]*m_z12 + 21.*m_z2k;

  Box2Abb[5357]=2. - 25.*m_z2k;

  Box2Abb[5358]=-25. + 7.*Box2Abb[5357]*m_z2k;

  Box2Abb[5359]=4. - 7.*Box2Abb[179]*m_z2k;

  Box2Abb[5360]=-3. + 12.*Box2Abb[5359]*m_z2k;

  Box2Abb[5361]=-7. + Box2Abb[5360]*m_z12 + Box2Abb[5358]*m_z12_2 + 5.*Box2Abb[4695]*m_z12_3 + 6.*m_z12_4 - 3.*Box2Abb[888]*m_z2k;

  Box2Abb[5362]=1. + 26.*m_z2k;

  Box2Abb[5363]=-5. + Box2Abb[844]*m_z2k;

  Box2Abb[5364]=5. + Box2Abb[5363]*m_z2k;

  Box2Abb[5365]=-33. + 29.*m_z2k;

  Box2Abb[5366]=3. + Box2Abb[5365]*m_z2k;

  Box2Abb[5367]=9. + Box2Abb[5366]*m_z2k;

  Box2Abb[5368]=Box2Abb[514]*pow(Box2Abb[61],4.) + Box2Abb[5152]*pow(Box2Abb[61],5.)*m_z12 + Box2Abb[5362]*pow(Box2Abb[61],4.)*m_z12_2 + Box2Abb[5367]*Box2Abb[61]*m_z12_3 + Box2Abb[5364]*m_z12_4;

  Box2Abb[5369]=10. - 4.*m_z2k_3;

  Box2Abb[5370]=-5. + 3.*Box2Abb[1071]*m_z2k;

  Box2Abb[5371]=-9. + 8.*Box2Abb[3303]*m_z2k;

  Box2Abb[5372]=10. + Box2Abb[5371]*m_z2k;

  Box2Abb[5373]=5. + Box2Abb[5372]*m_z2k;

  Box2Abb[5374]=-29. + 9.*m_z2k;

  Box2Abb[5375]=30. + Box2Abb[5374]*m_z2k;

  Box2Abb[5376]=-1. + 4.*Box2Abb[5375]*m_z2k;

  Box2Abb[5377]=-29. + Box2Abb[5376]*m_z2k;

  Box2Abb[5378]=-316. + 133.*m_z2k;

  Box2Abb[5379]=156. + Box2Abb[5378]*m_z2k;

  Box2Abb[5380]=56. + Box2Abb[5379]*m_z2k;

  Box2Abb[5381]=-53. + Box2Abb[5380]*m_z2k;

  Box2Abb[5382]=-293. + 174.*m_z2k;

  Box2Abb[5383]=39. + Box2Abb[5382]*m_z2k;

  Box2Abb[5384]=81. + Box2Abb[5383]*m_z2k;

  Box2Abb[5385]=-29. + Box2Abb[5384]*m_z2k;

  Box2Abb[5386]=-Box2Abb[5370]*pow(Box2Abb[61],4.) - Box2Abb[5377]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[5381]*pow(Box2Abb[61],2.)*m_z12_2 - Box2Abb[5385]*Box2Abb[61]*m_z12_3 - 2.*Box2Abb[5373]*m_z12_4 + Box2Abb[5369]*m_z12_5;

  Box2Abb[5387]=-20. - 36.*m_z2k_2 + 34.*m_z2k_3;

  Box2Abb[5388]=-8. + 21.*m_z2k;

  Box2Abb[5389]=-11. + 3.*Box2Abb[5388]*m_z2k;

  Box2Abb[5390]=2. + Box2Abb[5389]*m_z2k;

  Box2Abb[5391]=-362. + 255.*m_z2k;

  Box2Abb[5392]=36. + Box2Abb[5391]*m_z2k;

  Box2Abb[5393]=-18. + Box2Abb[5392]*m_z2k;

  Box2Abb[5394]=85. + Box2Abb[5393]*m_z2k;

  Box2Abb[5395]=614. + 41.*Box2Abb[5075]*m_z2k;

  Box2Abb[5396]=-12. + Box2Abb[5395]*m_z2k;

  Box2Abb[5397]=27. + Box2Abb[5396]*m_z2k;

  Box2Abb[5398]=-72. + Box2Abb[5397]*m_z2k;

  Box2Abb[5399]=35. + 4.*Box2Abb[2534]*m_z2k;

  Box2Abb[5400]=-150. + 7.*Box2Abb[5399]*m_z2k;

  Box2Abb[5401]=8. + Box2Abb[5400]*m_z2k;

  Box2Abb[5402]=4. + Box2Abb[5401]*m_z2k;

  Box2Abb[5403]=5. + Box2Abb[5402]*m_z2k;

  Box2Abb[5404]=Box2Abb[5390]*pow(Box2Abb[61],2.) + 3.*Box2Abb[5403]*m_z12 + Box2Abb[5398]*m_z12_2 + Box2Abb[5394]*m_z12_3 + Box2Abb[5387]*m_z12_4 - 10.*Box2Abb[890]*m_z12_5;

  Box2Abb[5405]=Box2Abb[5]*Box2Abb[5368]*Box2Abb[61]*m_x + Box2Abb[5386]*m_x_2 + Box2Abb[5404]*m_x_3 + Box2Abb[5340]*m_x_4 + Box2Abb[5351]*m_x_5 + Box2Abb[5361]*m_x_6 + Box2Abb[5356]*m_x_7 - Box2Abb[5353]*m_x_8 - pow(Box2Abb[5],4.)*pow(Box2Abb[61],5.)*m_z12 + m_x_9*m_z12;

  Box2Abb[5406]=3. + 3.*m_z12_2 + 10.*m_z12*m_z2k;

  Box2Abb[5407]=4. - 45.*m_z2k;

  Box2Abb[5408]=2. - 5.*Box2Abb[510]*m_z12 + 3.*m_z12_2 + Box2Abb[5407]*m_z2k;

  Box2Abb[5409]=-3.*Box2Abb[2153] + Box2Abb[5408]*m_z12;

  Box2Abb[5410]=2. - 5.*m_z2k + 3.*m_z2k_3;

  Box2Abb[5411]=-9. + Box2Abb[3480]*m_z2k;

  Box2Abb[5412]=Box2Abb[222]*pow(Box2Abb[61],4.) + Box2Abb[5411]*pow(Box2Abb[61],2.)*m_z12 + 3.*Box2Abb[5410]*m_z12_2 + 8.*m_z12_3*m_z2k;

  Box2Abb[5413]=16. + 93.*m_z2k;

  Box2Abb[5414]=8. + Box2Abb[5413]*m_z2k;

  Box2Abb[5415]=-4. + 15.*m_z2k;

  Box2Abb[5416]=2. + Box2Abb[5415]*m_z2k;

  Box2Abb[5417]=2. + Box2Abb[5416]*m_z2k;

  Box2Abb[5418]=3. + 8.*Box2Abb[5417]*m_z12 + Box2Abb[5414]*m_z12_2 - 2.*Box2Abb[2099]*m_z12_3 - m_z12_4 + 33.*m_z2k + 84.*m_z2k_2;

  Box2Abb[5419]=9. - 22.*m_z2k;

  Box2Abb[5420]=5. + 2.*Box2Abb[5419]*m_z2k;

  Box2Abb[5421]=3. + 56.*m_z2k;

  Box2Abb[5422]=4. - Box2Abb[5421]*m_z2k;

  Box2Abb[5423]=5. + Box2Abb[5422]*m_z2k;

  Box2Abb[5424]=-5. + 204.*m_z2k;

  Box2Abb[5425]=13. + Box2Abb[5424]*m_z2k;

  Box2Abb[5426]=22. + Box2Abb[5425]*m_z2k;

  Box2Abb[5427]=69. + 7.*Box2Abb[1192]*m_z2k;

  Box2Abb[5428]=42. + Box2Abb[5427]*m_z2k;

  Box2Abb[5429]=11. + Box2Abb[5428]*m_z2k;

  Box2Abb[5430]=3.*Box2Abb[5423] - 2.*Box2Abb[5429]*m_z12 - Box2Abb[5426]*m_z12_2 + Box2Abb[5420]*m_z12_3 + 5.*Box2Abb[12]*m_z12_4;

  Box2Abb[5431]=6. - 12.*m_z2k + 4.*m_z2k_3 + 5.*m_z2k_4;

  Box2Abb[5432]=21. + 2.*m_z2k;

  Box2Abb[5433]=7. + Box2Abb[5432]*m_z2k;

  Box2Abb[5434]=-15. + Box2Abb[5433]*m_z2k;

  Box2Abb[5435]=23. + 17.*m_z2k;

  Box2Abb[5436]=8. + Box2Abb[5435]*m_z2k;

  Box2Abb[5437]=-33. + Box2Abb[5436]*m_z2k;

  Box2Abb[5438]=3. + Box2Abb[5437]*m_z2k;

  Box2Abb[5439]=22. + 39.*m_z2k;

  Box2Abb[5440]=2. + Box2Abb[5439]*m_z2k;

  Box2Abb[5441]=-54. + Box2Abb[5440]*m_z2k;

  Box2Abb[5442]=15. + Box2Abb[5441]*m_z2k;

  Box2Abb[5443]=3.*pow(Box2Abb[61],7.) + 2.*Box2Abb[5431]*pow(Box2Abb[61],4.)*m_z12 + Box2Abb[5442]*pow(Box2Abb[61],3.)*m_z12_2 + 2.*Box2Abb[5438]*pow(Box2Abb[61],2.)*m_z12_3 + 2.*Box2Abb[5434]*Box2Abb[61]*m_z12_4*m_z2k + 4.*Box2Abb[444]*m_z12_5*m_z2k_2;

  Box2Abb[5444]=-5. + Box2Abb[725]*m_z2k;

  Box2Abb[5445]=-1. + Box2Abb[3704]*m_z2k;

  Box2Abb[5446]=-1. + Box2Abb[5445]*m_z2k;

  Box2Abb[5447]=-9. + 71.*m_z2k;

  Box2Abb[5448]=-9. + Box2Abb[5447]*m_z2k;

  Box2Abb[5449]=10. + Box2Abb[5448]*m_z2k;

  Box2Abb[5450]=-31. + 98.*m_z2k;

  Box2Abb[5451]=21. + Box2Abb[5450]*m_z2k;

  Box2Abb[5452]=31. + Box2Abb[5451]*m_z2k;

  Box2Abb[5453]=5. + Box2Abb[5452]*m_z2k;

  Box2Abb[5454]=95. + 7.*Box2Abb[4547]*m_z2k;

  Box2Abb[5455]=12. + Box2Abb[5454]*m_z2k;

  Box2Abb[5456]=13. + 2.*Box2Abb[5455]*m_z2k;

  Box2Abb[5457]=-8. + 2.*Box2Abb[5456]*m_z2k;

  Box2Abb[5458]=-9. + Box2Abb[5457]*m_z12 + 3.*Box2Abb[5453]*m_z12_2 + 2.*Box2Abb[5449]*m_z12_3 + 2.*Box2Abb[5444]*m_z12_4 + 15.*Box2Abb[5446]*m_z2k;

  Box2Abb[5459]=10. + 6.*m_z2k_2 - 51.*m_z2k_3;

  Box2Abb[5460]=-95. + 56.*m_z2k;

  Box2Abb[5461]=40. + Box2Abb[5460]*m_z2k;

  Box2Abb[5462]=25. - 97.*m_z2k;

  Box2Abb[5463]=15. + 2.*Box2Abb[5462]*m_z2k;

  Box2Abb[5464]=-4. + Box2Abb[5463]*m_z2k;

  Box2Abb[5465]=-35. + Box2Abb[5464]*m_z2k;

  Box2Abb[5466]=-187. + 294.*m_z2k;

  Box2Abb[5467]=241. + Box2Abb[5466]*m_z2k;

  Box2Abb[5468]=96. + Box2Abb[5467]*m_z2k;

  Box2Abb[5469]=91. + Box2Abb[5468]*m_z2k;

  Box2Abb[5470]=-19. + Box2Abb[5469]*m_z2k;

  Box2Abb[5471]=53. + 7.*Box2Abb[2800]*m_z2k;

  Box2Abb[5472]=-30. + Box2Abb[5471]*m_z2k;

  Box2Abb[5473]=-11. - 5.*Box2Abb[5472]*m_z2k;

  Box2Abb[5474]=22. + Box2Abb[5473]*m_z2k;

  Box2Abb[5475]=10. + 2.*Box2Abb[5474]*m_z2k;

  Box2Abb[5476]=-5. + Box2Abb[5475]*m_z12 - Box2Abb[5470]*m_z12_2 + Box2Abb[5465]*m_z12_3 + Box2Abb[5459]*m_z12_4 + 14.*m_z2k - 3.*Box2Abb[5461]*m_z2k_3;

  Box2Abb[5477]=-19. + 64.*m_z2k;

  Box2Abb[5478]=-12. + Box2Abb[5477]*m_z2k;

  Box2Abb[5479]=10. + Box2Abb[5478]*m_z2k;

  Box2Abb[5480]=-5. + Box2Abb[5479]*m_z2k;

  Box2Abb[5481]=-28. + 15.*m_z2k;

  Box2Abb[5482]=51. + Box2Abb[5481]*m_z2k;

  Box2Abb[5483]=-70. + Box2Abb[5482]*m_z2k;

  Box2Abb[5484]=35. + Box2Abb[5483]*m_z2k;

  Box2Abb[5485]=2. - 9.*m_z2k + 2.*Box2Abb[5484]*m_z2k_3;

  Box2Abb[5486]=-87. + 28.*m_z2k;

  Box2Abb[5487]=95. + Box2Abb[5486]*m_z2k;

  Box2Abb[5488]=-40. + Box2Abb[5487]*m_z2k;

  Box2Abb[5489]=8. + 3.*Box2Abb[5488]*m_z2k;

  Box2Abb[5490]=3. + Box2Abb[5489]*m_z2k;

  Box2Abb[5491]=-27. + 71.*m_z2k;

  Box2Abb[5492]=40. + Box2Abb[5491]*m_z2k;

  Box2Abb[5493]=40. + Box2Abb[5492]*m_z2k;

  Box2Abb[5494]=3. + Box2Abb[5493]*m_z2k;

  Box2Abb[5495]=11. + Box2Abb[5494]*m_z2k;

  Box2Abb[5496]=-41. + 42.*m_z2k;

  Box2Abb[5497]=67. + Box2Abb[5496]*m_z2k;

  Box2Abb[5498]=-217. + 5.*Box2Abb[5497]*m_z2k;

  Box2Abb[5499]=-11. + Box2Abb[5498]*m_z2k;

  Box2Abb[5500]=22. + Box2Abb[5499]*m_z2k;

  Box2Abb[5501]=-26. + Box2Abb[5500]*m_z2k;

  Box2Abb[5502]=1. + 4.*Box2Abb[5485]*m_z12 + Box2Abb[5501]*m_z12_2 + 2.*Box2Abb[5495]*m_z12_3 + Box2Abb[5480]*m_z12_4 + Box2Abb[5490]*m_z2k + 6.*m_z12_5*m_z2k_3;

  Box2Abb[5503]=-1. + Box2Abb[3244]*Box2Abb[61]*m_z2k;

  Box2Abb[5504]=1. + 3.*Box2Abb[5503]*m_z2k;

  Box2Abb[5505]=6. + m_z2k + 5.*m_z2k_2;

  Box2Abb[5506]=26. + 5.*Box2Abb[5505]*m_z2k;

  Box2Abb[5507]=5. + Box2Abb[5506]*m_z2k;

  Box2Abb[5508]=-1. + Box2Abb[5507]*m_z2k;

  Box2Abb[5509]=-69. + 6.*m_z2k - 68.*m_z2k_2;

  Box2Abb[5510]=112. + Box2Abb[5509]*m_z2k;

  Box2Abb[5511]=6. + Box2Abb[5510]*m_z2k;

  Box2Abb[5512]=-6. + Box2Abb[5511]*m_z2k;

  Box2Abb[5513]=-5. + Box2Abb[5512]*m_z2k;

  Box2Abb[5514]=-17. + 36.*m_z2k;

  Box2Abb[5515]=44. + Box2Abb[5514]*m_z2k;

  Box2Abb[5516]=-218. + 3.*Box2Abb[5515]*m_z2k;

  Box2Abb[5517]=36. + Box2Abb[5516]*m_z2k;

  Box2Abb[5518]=21. + Box2Abb[5517]*m_z2k;

  Box2Abb[5519]=8. + Box2Abb[5518]*m_z2k;

  Box2Abb[5520]=-22. + 45.*m_z2k;

  Box2Abb[5521]=77. + Box2Abb[5520]*m_z2k;

  Box2Abb[5522]=-148. + Box2Abb[5521]*m_z2k;

  Box2Abb[5523]=45. + Box2Abb[5522]*m_z2k;

  Box2Abb[5524]=6. + Box2Abb[5523]*m_z2k;

  Box2Abb[5525]=5. + Box2Abb[5524]*m_z2k;

  Box2Abb[5526]=-Box2Abb[5504]*pow(Box2Abb[61],3.) - Box2Abb[5525]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[5519]*Box2Abb[61]*m_z12_2 + Box2Abb[5513]*m_z12_3 - Box2Abb[5508]*m_z12_4 - 8.*Box2Abb[61]*m_z12_5*m_z2k_3;

  Box2Abb[5527]=Box2Abb[5526]*m_x_3 + Box2Abb[5502]*m_x_4 + Box2Abb[5476]*m_x_5 + Box2Abb[5458]*m_x_6 + Box2Abb[5430]*m_x_7 + Box2Abb[5418]*m_x_8 + Box2Abb[5409]*m_x_9 + Box2Abb[5406]*m_x_10 - m_x_11*m_z12 + Box2Abb[5443]*m_x_2*m_z2k - Box2Abb[5]*Box2Abb[5412]*pow(Box2Abb[61],2.)*m_x*m_z12*m_z2k_2 + pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*Box2Abb[65]*m_z12_2*m_z2k_3;

  Box2Abb[5528]=Box2Abb[5527]*m_cL + Box2Abb[5405]*Box2Abb[8]*m_cR*m_x;

  Box2Abb[5529]=39. - 19.*m_z12;

  Box2Abb[5530]=15. + Box2Abb[5529]*m_z12;

  Box2Abb[5531]=50. + Box2Abb[5530]*m_z12;

  Box2Abb[5532]=30. - 37.*m_z12;

  Box2Abb[5533]=18. + Box2Abb[5532]*m_z12;

  Box2Abb[5534]=50. + Box2Abb[5533]*m_z12;

  Box2Abb[5535]=40. + 2.*Box2Abb[2331]*Box2Abb[9]*m_z12 + Box2Abb[5531]*m_z2k + Box2Abb[5534]*m_z2k_2 - 14.*Box2Abb[1629]*m_z12*m_z2k_3;

  Box2Abb[5536]=12. - 5.*m_z12;

  Box2Abb[5537]=-45. + 4.*Box2Abb[5536]*m_z12;

  Box2Abb[5538]=2. + Box2Abb[5537]*m_z12;

  Box2Abb[5539]=25. + Box2Abb[3731]*m_z12;

  Box2Abb[5540]=7. + Box2Abb[5539]*m_z12;

  Box2Abb[5541]=8. + m_z12 + m_z12_2;

  Box2Abb[5542]=-27. + Box2Abb[5541]*m_z12;

  Box2Abb[5543]=16. + Box2Abb[5542]*m_z12;

  Box2Abb[5544]=-10. + Box2Abb[5543]*m_z12;

  Box2Abb[5545]=43. - 12.*m_z12;

  Box2Abb[5546]=66. + Box2Abb[5545]*m_z12;

  Box2Abb[5547]=-60. + Box2Abb[5546]*m_z12;

  Box2Abb[5548]=40. + Box2Abb[5547]*m_z12;

  Box2Abb[5549]=-9. + 14.*m_z12;

  Box2Abb[5550]=9. + 4.*Box2Abb[5549]*m_z12;

  Box2Abb[5551]=14. + Box2Abb[5550]*m_z12;

  Box2Abb[5552]=-2. + 2.*m_z12_2 - m_z12_3 + Box2Abb[5538]*m_z2k - 2.*Box2Abb[4]*Box2Abb[5540]*m_z2k_2 + 4.*Box2Abb[5544]*m_z2k_3 + Box2Abb[5548]*m_z2k_4 - Box2Abb[5551]*m_z2k_5 - 14.*Box2Abb[4]*m_z12*m_z2k_6;

  Box2Abb[5553]=-6. + m_z12 - 11.*m_z2k;

  Box2Abb[5554]=14. + Box2Abb[5553]*m_z12 + 26.*m_z2k;

  Box2Abb[5555]=4. + Box2Abb[5554]*m_z12;

  Box2Abb[5556]=-35. + 11.*m_z12;

  Box2Abb[5557]=-10. + Box2Abb[1087]*m_z12 - 43.*m_z2k + 8.*Box2Abb[9]*m_z12*m_z2k + 2.*Box2Abb[5556]*m_z2k_2;

  Box2Abb[5558]=-20. + Box2Abb[5557]*m_z12 - 22.*m_z2k;

  Box2Abb[5559]=-3. + Box2Abb[833]*m_z2k;

  Box2Abb[5560]=Box2Abb[460]*pow(Box2Abb[61],4.) + 2.*Box2Abb[222]*Box2Abb[460]*pow(Box2Abb[61],3.)*m_z12 + 6.*Box2Abb[5559]*pow(Box2Abb[61],2.)*m_z12_2 + 2.*Box2Abb[1163]*Box2Abb[61]*m_z12_3 + 2.*Box2Abb[231]*m_z12_4*m_z2k;

  Box2Abb[5561]=1. + Box2Abb[5066]*m_z2k;

  Box2Abb[5562]=-5. + 2.*Box2Abb[5561]*m_z2k;

  Box2Abb[5563]=33. + 7.*m_z2k;

  Box2Abb[5564]=27. + Box2Abb[5563]*m_z2k;

  Box2Abb[5565]=20. + Box2Abb[5564]*m_z2k;

  Box2Abb[5566]=-3. + Box2Abb[5565]*m_z2k;

  Box2Abb[5567]=-13. + 14.*m_z2k;

  Box2Abb[5568]=9. + Box2Abb[5567]*m_z2k;

  Box2Abb[5569]=14. + 5.*Box2Abb[5568]*m_z2k;

  Box2Abb[5570]=22. - Box2Abb[5569]*m_z2k;

  Box2Abb[5571]=-38. + Box2Abb[5570]*m_z12 - 2.*Box2Abb[5566]*m_z12_2 + 2.*Box2Abb[5562]*m_z12_3 - 20.*Box2Abb[474]*m_z2k + 12.*m_z12_4*m_z2k_2;

  Box2Abb[5572]=-5. - 8.*m_z2k + 6.*m_z2k_2;

  Box2Abb[5573]=6. + Box2Abb[5572]*m_z2k;

  Box2Abb[5574]=-9. + 2.*Box2Abb[258]*m_z2k;

  Box2Abb[5575]=15. + Box2Abb[5574]*m_z2k;

  Box2Abb[5576]=-3. + 2.*Box2Abb[5244]*m_z2k;

  Box2Abb[5577]=-60. + Box2Abb[5576]*m_z2k;

  Box2Abb[5578]=15. + Box2Abb[5577]*m_z2k;

  Box2Abb[5579]=-50. + 17.*m_z2k;

  Box2Abb[5580]=-36. + Box2Abb[5579]*m_z2k;

  Box2Abb[5581]=62. + Box2Abb[5580]*m_z2k;

  Box2Abb[5582]=-5. + Box2Abb[5581]*m_z2k;

  Box2Abb[5583]=2.*pow(Box2Abb[61],5.) - Box2Abb[5575]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[5578]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[5582]*Box2Abb[61]*m_z12_3 + 4.*Box2Abb[5573]*m_z12_4*m_z2k + 2.*m_z12_5*m_z2k_3;

  Box2Abb[5584]=8. + 25.*m_z2k;

  Box2Abb[5585]=36. + Box2Abb[5584]*m_z2k;

  Box2Abb[5586]=22. + Box2Abb[5585]*m_z2k;

  Box2Abb[5587]=5. + Box2Abb[5586]*m_z2k;

  Box2Abb[5588]=55. + 7.*Box2Abb[2534]*m_z2k;

  Box2Abb[5589]=11. + Box2Abb[5588]*m_z2k;

  Box2Abb[5590]=23. + Box2Abb[5589]*m_z2k;

  Box2Abb[5591]=-4. + Box2Abb[5590]*m_z2k;

  Box2Abb[5592]=-69. + 4.*Box2Abb[3239]*m_z2k;

  Box2Abb[5593]=8. + Box2Abb[5592]*m_z2k;

  Box2Abb[5594]=-29. + Box2Abb[5593]*m_z2k;

  Box2Abb[5595]=-6. + Box2Abb[5594]*m_z2k;

  Box2Abb[5596]=-4.*Box2Abb[1359] + 2.*Box2Abb[5591]*m_z12 + Box2Abb[5595]*m_z12_2 + Box2Abb[5587]*m_z12_3 - 6.*Box2Abb[150]*m_z12_4*m_z2k_2 + 40.*Box2Abb[61]*m_z2k_3;

  Box2Abb[5597]=Box2Abb[5552]*m_x_3 + Box2Abb[5596]*m_x_4 + Box2Abb[5571]*m_x_5 + Box2Abb[5535]*m_x_6 + Box2Abb[5558]*m_x_7 + Box2Abb[5555]*m_x_8 + 2.*Box2Abb[72]*m_x_9*m_z12 + Box2Abb[5583]*m_x_2*m_z2k + Box2Abb[5560]*Box2Abb[61]*m_x*m_z12*m_z2k_2 - Box2Abb[4904]*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_3;

  Box2Abb[5598]=78. - 5.*Box2Abb[703]*m_z12;

  Box2Abb[5599]=12. + Box2Abb[5598]*m_z12;

  Box2Abb[5600]=-6. + Box2Abb[872]*m_z12;

  Box2Abb[5601]=-42. - 5.*Box2Abb[5600]*m_z12;

  Box2Abb[5602]=-3. + 77.*m_z12;

  Box2Abb[5603]=30. + Box2Abb[5602]*m_z12;

  Box2Abb[5604]=36. + Box2Abb[5603]*m_z12;

  Box2Abb[5605]=30. + Box2Abb[5604]*m_z12;

  Box2Abb[5606]=-2. + 15.*m_z12;

  Box2Abb[5607]=192. + 13.*Box2Abb[5606]*m_z12;

  Box2Abb[5608]=-424. + Box2Abb[5607]*m_z12;

  Box2Abb[5609]=76. + Box2Abb[5608]*m_z12;

  Box2Abb[5610]=32. + 5.*Box2Abb[1682]*m_z12;

  Box2Abb[5611]=-8. + Box2Abb[5610]*m_z12;

  Box2Abb[5612]=-1. + 9.*m_z12;

  Box2Abb[5613]=-46. + Box2Abb[5599]*m_z12 + 20.*m_z2k + 2.*Box2Abb[5601]*m_z12*m_z2k + Box2Abb[5605]*m_z2k_2 + Box2Abb[5609]*m_z2k_3 + 10.*Box2Abb[5611]*m_z2k_4 + 28.*Box2Abb[5612]*m_z12*m_z2k_5;

  Box2Abb[5614]=16. + Box2Abb[2459]*m_z12 + 52.*m_z2k;

  Box2Abb[5615]=8. + Box2Abb[5614]*m_z12;

  Box2Abb[5616]=-5. + 9.*m_z12;

  Box2Abb[5617]=44. + Box2Abb[2480]*m_z12 - 14.*m_z2k + 10.*Box2Abb[735]*m_z12*m_z2k + 28.*Box2Abb[5616]*m_z2k_2;

  Box2Abb[5618]=-44.*Box2Abb[12] + Box2Abb[5617]*m_z12;

  Box2Abb[5619]=17. + 10.*m_z2k;

  Box2Abb[5620]=1. - Box2Abb[5619]*m_z2k;

  Box2Abb[5621]=2.*Box2Abb[230]*pow(Box2Abb[61],3.) + 3.*Box2Abb[183]*pow(Box2Abb[61],2.)*m_z12 - 6.*Box2Abb[4994]*Box2Abb[61]*m_z12_2 + Box2Abb[5620]*m_z12_3;

  Box2Abb[5622]=57. + 50.*m_z2k;

  Box2Abb[5623]=50. + Box2Abb[5622]*m_z2k;

  Box2Abb[5624]=15. + 106.*m_z2k;

  Box2Abb[5625]=85. - 2.*Box2Abb[5624]*m_z2k;

  Box2Abb[5626]=6. + 35.*m_z2k;

  Box2Abb[5627]=5. + Box2Abb[5626]*m_z2k;

  Box2Abb[5628]=43. + 4.*Box2Abb[5627]*m_z2k;

  Box2Abb[5629]=-36. + 49.*m_z2k;

  Box2Abb[5630]=-37. + Box2Abb[5629]*m_z2k;

  Box2Abb[5631]=-57. + 2.*Box2Abb[5630]*m_z2k;

  Box2Abb[5632]=2.*Box2Abb[5623] + 2.*Box2Abb[5631]*m_z12 - 3.*Box2Abb[5628]*m_z12_2 + Box2Abb[5625]*m_z12_3 + Box2Abb[2486]*m_z12_4;

  Box2Abb[5633]=-4. + 5.*Box2Abb[230]*m_z2k;

  Box2Abb[5634]=25. + 4.*Box2Abb[556]*m_z2k;

  Box2Abb[5635]=49. + Box2Abb[5634]*m_z2k;

  Box2Abb[5636]=107. + 261.*m_z2k;

  Box2Abb[5637]=41. + Box2Abb[5636]*m_z2k;

  Box2Abb[5638]=-60. + 2.*Box2Abb[5637]*m_z2k;

  Box2Abb[5639]=37. - 14.*m_z2k;

  Box2Abb[5640]=-8. + 5.*Box2Abb[5639]*m_z2k;

  Box2Abb[5641]=29. + Box2Abb[5640]*m_z2k;

  Box2Abb[5642]=36. + Box2Abb[5641]*m_z2k;

  Box2Abb[5643]=-9. + 20.*Box2Abb[560]*m_z2k;

  Box2Abb[5644]=27. + Box2Abb[5643]*m_z2k;

  Box2Abb[5645]=50. + Box2Abb[5644]*m_z2k;

  Box2Abb[5646]=-2.*Box2Abb[5635] + 2.*Box2Abb[5642]*m_z12 + 3.*Box2Abb[5645]*m_z12_2 + Box2Abb[5638]*m_z12_3 + 5.*Box2Abb[5633]*m_z12_4;

  Box2Abb[5647]=1. + 2.*Box2Abb[291]*Box2Abb[61]*m_z2k;

  Box2Abb[5648]=69. - 46.*m_z2k + 20.*m_z2k_2;

  Box2Abb[5649]=8. + Box2Abb[5648]*m_z2k;

  Box2Abb[5650]=5. + Box2Abb[5649]*m_z2k;

  Box2Abb[5651]=9. + Box2Abb[222]*m_z2k;

  Box2Abb[5652]=-9. + 4.*Box2Abb[5651]*m_z2k;

  Box2Abb[5653]=1. + Box2Abb[5652]*m_z2k;

  Box2Abb[5654]=-90. + 59.*m_z2k;

  Box2Abb[5655]=8. + Box2Abb[5654]*m_z2k;

  Box2Abb[5656]=20. + Box2Abb[5655]*m_z2k;

  Box2Abb[5657]=5. + Box2Abb[5656]*m_z2k;

  Box2Abb[5658]=-2. + Box2Abb[5657]*m_z2k;

  Box2Abb[5659]=-2.*Box2Abb[230]*pow(Box2Abb[61],5.) + 2.*Box2Abb[5647]*pow(Box2Abb[61],4.)*m_z12 + 3.*Box2Abb[5653]*pow(Box2Abb[61],3.)*m_z12_2 + Box2Abb[5650]*pow(Box2Abb[61],2.)*m_z12_3 + Box2Abb[5658]*m_z12_4 + 24.*m_z12_5*m_z2k_4;

  Box2Abb[5660]=55. + 14.*m_z2k;

  Box2Abb[5661]=-82. + Box2Abb[5660]*m_z2k;

  Box2Abb[5662]=-37. + Box2Abb[5661]*m_z2k;

  Box2Abb[5663]=2. + Box2Abb[5662]*m_z2k;

  Box2Abb[5664]=2. + 101.*m_z2k;

  Box2Abb[5665]=15. + Box2Abb[5664]*m_z2k;

  Box2Abb[5666]=5. + Box2Abb[5665]*m_z2k;

  Box2Abb[5667]=-5. + Box2Abb[5666]*m_z2k;

  Box2Abb[5668]=81. + 4.*Box2Abb[4847]*m_z2k;

  Box2Abb[5669]=-13. + Box2Abb[5668]*m_z2k;

  Box2Abb[5670]=-17. + Box2Abb[5669]*m_z2k;

  Box2Abb[5671]=-3. + Box2Abb[5670]*m_z2k;

  Box2Abb[5672]=-147. + 127.*m_z2k;

  Box2Abb[5673]=22. + Box2Abb[5672]*m_z2k;

  Box2Abb[5674]=16. + Box2Abb[5673]*m_z2k;

  Box2Abb[5675]=15. + Box2Abb[5674]*m_z2k;

  Box2Abb[5676]=7. + Box2Abb[5675]*m_z2k;

  Box2Abb[5677]=-14.*Box2Abb[2307]*pow(Box2Abb[61],3.) + 2.*Box2Abb[5663]*pow(Box2Abb[61],2.)*m_z12 + 3.*Box2Abb[5671]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[5676]*m_z12_3 + 2.*Box2Abb[5667]*m_z12_4 + 24.*m_z12_5*m_z2k_3;

  Box2Abb[5678]=-Box2Abb[5659]*m_x_2 + Box2Abb[5677]*m_x_3 - Box2Abb[5613]*m_x_4 + Box2Abb[5646]*m_x_5 + Box2Abb[5632]*m_x_6 + Box2Abb[5618]*m_x_7 + Box2Abb[5615]*m_x_8 + 4.*Box2Abb[155]*m_x_9*m_z12 + Box2Abb[5621]*pow(Box2Abb[61],3.)*m_x*m_z12*m_z2k + 3.*Box2Abb[5]*pow(Box2Abb[61],5.)*m_z12_3*m_z2k_2;

  Box2Abb[5679]=-2.*Box2Abb[5597]*m_cL + Box2Abb[5678]*m_cR;

  Box2Abb[5680]=11. + 5.*Box2Abb[1058]*m_z12;

  Box2Abb[5681]=-58. + Box2Abb[5680]*m_z12;

  Box2Abb[5682]=58. + 7.*m_z12;

  Box2Abb[5683]=156. + Box2Abb[5682]*m_z12;

  Box2Abb[5684]=-100. + Box2Abb[5683]*m_z12;

  Box2Abb[5685]=-25. + 9.*Box2Abb[3403]*m_z12;

  Box2Abb[5686]=60. + Box2Abb[5681]*m_z12 + 82.*m_z2k + Box2Abb[5684]*m_z12*m_z2k - 4.*Box2Abb[5685]*m_z2k_2 + 28.*Box2Abb[500]*m_z12*m_z2k_3;

  Box2Abb[5687]=4. + 7.*m_z12 - 20.*m_z2k;

  Box2Abb[5688]=16. + Box2Abb[5687]*m_z12 + 52.*m_z2k;

  Box2Abb[5689]=8. + Box2Abb[5688]*m_z12;

  Box2Abb[5690]=23. + 2.*Box2Abb[3437]*m_z12;

  Box2Abb[5691]=40. + 3.*m_z12;

  Box2Abb[5692]=14.*Box2Abb[179] + Box2Abb[5690]*m_z12 + 2.*Box2Abb[5691]*m_z12*m_z2k - 28.*Box2Abb[815]*m_z2k_2;

  Box2Abb[5693]=36. + Box2Abb[5692]*m_z12 + 44.*m_z2k;

  Box2Abb[5694]=9. + 4.*Box2Abb[460]*m_z2k;

  Box2Abb[5695]=-9. + 17.*m_z2k;

  Box2Abb[5696]=5. + Box2Abb[5695]*m_z2k;

  Box2Abb[5697]=-11. + 42.*m_z2k;

  Box2Abb[5698]=3. + Box2Abb[5697]*m_z2k;

  Box2Abb[5699]=-2.*Box2Abb[230]*pow(Box2Abb[61],4.) + Box2Abb[5694]*pow(Box2Abb[61],3.)*m_z12 + 2.*Box2Abb[5696]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[5698]*Box2Abb[61]*m_z12_3 + 16.*m_z12_4*m_z2k_2;

  Box2Abb[5700]=-20. + Box2Abb[1125]*m_z2k;

  Box2Abb[5701]=5. + 4.*m_z2k + 60.*m_z2k_2;

  Box2Abb[5702]=13. + Box2Abb[5701]*m_z2k;

  Box2Abb[5703]=-97. + 13.*m_z2k;

  Box2Abb[5704]=-63. + Box2Abb[5703]*m_z2k;

  Box2Abb[5705]=10. + Box2Abb[5704]*m_z2k;

  Box2Abb[5706]=-24. + 5.*Box2Abb[5639]*m_z2k;

  Box2Abb[5707]=5. + Box2Abb[5706]*m_z2k;

  Box2Abb[5708]=4. + Box2Abb[5707]*m_z2k;

  Box2Abb[5709]=19. + 7.*m_z2k;

  Box2Abb[5710]=63. + 20.*Box2Abb[5709]*m_z2k;

  Box2Abb[5711]=-29. + Box2Abb[5710]*m_z2k;

  Box2Abb[5712]=10. + Box2Abb[5711]*m_z2k;

  Box2Abb[5713]=-2.*Box2Abb[5702] + 2.*Box2Abb[5708]*m_z12 - Box2Abb[5712]*m_z12_2 + 2.*Box2Abb[5705]*m_z12_3 + Box2Abb[5700]*m_z12_4;

  Box2Abb[5714]=-5. + Box2Abb[1151]*m_z2k;

  Box2Abb[5715]=5. + 2.*Box2Abb[5714]*m_z2k;

  Box2Abb[5716]=-11. + 13.*m_z2k;

  Box2Abb[5717]=-27. + Box2Abb[5716]*m_z2k;

  Box2Abb[5718]=63. + 4.*Box2Abb[5717]*m_z2k;

  Box2Abb[5719]=-15. + Box2Abb[5718]*m_z2k;

  Box2Abb[5720]=17. + 35.*m_z2k;

  Box2Abb[5721]=-85. + 3.*Box2Abb[5720]*m_z2k;

  Box2Abb[5722]=23. + Box2Abb[5721]*m_z2k;

  Box2Abb[5723]=-2. + Box2Abb[5722]*m_z2k;

  Box2Abb[5724]=-5. + 78.*m_z2k;

  Box2Abb[5725]=-169. + 2.*Box2Abb[5724]*m_z2k;

  Box2Abb[5726]=72. + Box2Abb[5725]*m_z2k;

  Box2Abb[5727]=-9. + Box2Abb[5726]*m_z2k;

  Box2Abb[5728]=2.*Box2Abb[230]*pow(Box2Abb[61],5.) - 2.*Box2Abb[5715]*pow(Box2Abb[61],4.)*m_z12 + Box2Abb[5719]*pow(Box2Abb[61],3.)*m_z12_2 + Box2Abb[5727]*pow(Box2Abb[61],2.)*m_z12_3 + Box2Abb[5723]*Box2Abb[61]*m_z12_4 + 8.*Box2Abb[411]*m_z12_5*m_z2k_3;

  Box2Abb[5729]=1. + 2.*Box2Abb[1071]*m_z2k;

  Box2Abb[5730]=-58. + Box2Abb[5660]*m_z2k;

  Box2Abb[5731]=-13. + Box2Abb[5730]*m_z2k;

  Box2Abb[5732]=10. + Box2Abb[5731]*m_z2k;

  Box2Abb[5733]=70. + 27.*m_z2k;

  Box2Abb[5734]=77. + Box2Abb[5733]*m_z2k;

  Box2Abb[5735]=-25. + Box2Abb[5734]*m_z2k;

  Box2Abb[5736]=5. + Box2Abb[5735]*m_z2k;

  Box2Abb[5737]=102. + m_z2k - 121.*m_z2k_2;

  Box2Abb[5738]=128. + Box2Abb[5737]*m_z2k;

  Box2Abb[5739]=-89. + Box2Abb[5738]*m_z2k;

  Box2Abb[5740]=19. + Box2Abb[5739]*m_z2k;

  Box2Abb[5741]=-19. + 35.*m_z2k;

  Box2Abb[5742]=-343. + 4.*Box2Abb[5741]*m_z2k;

  Box2Abb[5743]=27. + Box2Abb[5742]*m_z2k;

  Box2Abb[5744]=135. + Box2Abb[5743]*m_z2k;

  Box2Abb[5745]=-51. + Box2Abb[5744]*m_z2k;

  Box2Abb[5746]=-2.*Box2Abb[5729]*pow(Box2Abb[61],3.) + 2.*Box2Abb[5732]*pow(Box2Abb[61],2.)*m_z12 - Box2Abb[5745]*Box2Abb[61]*m_z12_2 + 2.*Box2Abb[5740]*m_z12_3 - 2.*Box2Abb[5736]*m_z12_4 + 8.*m_z12_5*m_z2k_3;

  Box2Abb[5747]=47. - 31.*m_z2k;

  Box2Abb[5748]=-30. + Box2Abb[5747]*m_z2k;

  Box2Abb[5749]=20. + Box2Abb[5748]*m_z2k;

  Box2Abb[5750]=-7. + 20.*m_z2k;

  Box2Abb[5751]=-9. + 2.*Box2Abb[5750]*m_z2k;

  Box2Abb[5752]=5. + Box2Abb[5751]*m_z2k;

  Box2Abb[5753]=161. + 65.*m_z2k;

  Box2Abb[5754]=191. + 2.*Box2Abb[5753]*m_z2k;

  Box2Abb[5755]=172. + Box2Abb[5754]*m_z2k;

  Box2Abb[5756]=-55. + Box2Abb[5755]*m_z2k;

  Box2Abb[5757]=-80. + 7.*m_z2k;

  Box2Abb[5758]=98. + Box2Abb[5757]*m_z2k;

  Box2Abb[5759]=-13. + Box2Abb[5758]*m_z2k;

  Box2Abb[5760]=13. + Box2Abb[5759]*m_z2k;

  Box2Abb[5761]=1. + Box2Abb[5760]*m_z2k;

  Box2Abb[5762]=25. + 49.*m_z2k;

  Box2Abb[5763]=-92. + Box2Abb[5762]*m_z2k;

  Box2Abb[5764]=-59. + 2.*Box2Abb[5763]*m_z2k;

  Box2Abb[5765]=-86. + Box2Abb[5764]*m_z2k;

  Box2Abb[5766]=25. + Box2Abb[5765]*m_z2k;

  Box2Abb[5767]=2.*Box2Abb[5752]*Box2Abb[61] + 4.*Box2Abb[5761]*m_z12 + 2.*Box2Abb[5766]*m_z12_2 + Box2Abb[5756]*m_z12_3 + Box2Abb[5749]*m_z12_4;

  Box2Abb[5768]=Box2Abb[5728]*m_x_2 + Box2Abb[5746]*m_x_3 + Box2Abb[5767]*m_x_4 + Box2Abb[5713]*m_x_5 + Box2Abb[5686]*m_x_6 - Box2Abb[5693]*m_x_7 + Box2Abb[5689]*m_x_8 + 4.*Box2Abb[72]*m_x_9*m_z12 - Box2Abb[5699]*pow(Box2Abb[61],2.)*m_x*m_z12*m_z2k - Box2Abb[5]*pow(Box2Abb[61],5.)*m_z12_3*m_z2k_2;

  Box2Abb[5769]=32. + m_z12;

  Box2Abb[5770]=53. + Box2Abb[5769]*m_z12;

  Box2Abb[5771]=-82. + Box2Abb[5770]*m_z12;

  Box2Abb[5772]=-36. + Box2Abb[5771]*m_z12;

  Box2Abb[5773]=42. + 11.*m_z12;

  Box2Abb[5774]=-55. + 4.*Box2Abb[5773]*m_z12;

  Box2Abb[5775]=-183. + Box2Abb[5774]*m_z12;

  Box2Abb[5776]=274. + 55.*m_z12;

  Box2Abb[5777]=-191. + Box2Abb[5776]*m_z12;

  Box2Abb[5778]=72. + Box2Abb[5777]*m_z12;

  Box2Abb[5779]=14. - 25.*m_z12;

  Box2Abb[5780]=40. + Box2Abb[5772]*m_z12 + 68.*m_z2k + Box2Abb[5775]*m_z12*m_z2k + Box2Abb[5778]*m_z2k_2 + 12.*Box2Abb[5779]*m_z12*m_z2k_3;

  Box2Abb[5781]=-3. + Box2Abb[2741]*m_z12;

  Box2Abb[5782]=-6. + 5.*Box2Abb[5781]*m_z12;

  Box2Abb[5783]=8. + 3.*m_z12;

  Box2Abb[5784]=-22. + 5.*Box2Abb[5783]*m_z12;

  Box2Abb[5785]=-53. + Box2Abb[5784]*m_z12;

  Box2Abb[5786]=7. + Box2Abb[5785]*m_z12;

  Box2Abb[5787]=16. + Box2Abb[5786]*m_z12;

  Box2Abb[5788]=37. + 7.*m_z12;

  Box2Abb[5789]=-125. + 2.*Box2Abb[5788]*m_z12;

  Box2Abb[5790]=31. + Box2Abb[5789]*m_z12;

  Box2Abb[5791]=-707. + 474.*m_z12 + 4.*m_z12_2;

  Box2Abb[5792]=315. + Box2Abb[5791]*m_z12;

  Box2Abb[5793]=-40. + Box2Abb[5792]*m_z12;

  Box2Abb[5794]=842. - 277.*m_z12;

  Box2Abb[5795]=-555. + Box2Abb[5794]*m_z12;

  Box2Abb[5796]=100. + Box2Abb[5795]*m_z12;

  Box2Abb[5797]=2.*Box2Abb[4]*Box2Abb[5782]*m_z12 + Box2Abb[5787]*m_z2k + Box2Abb[5790]*Box2Abb[9]*m_z12*m_z2k_2 + Box2Abb[5793]*m_z2k_3 + Box2Abb[5796]*m_z2k_4 + 42.*Box2Abb[879]*m_z12*m_z2k_5;

  Box2Abb[5798]=-5. + Box2Abb[3432]*m_z12;

  Box2Abb[5799]=69. - 5.*Box2Abb[5798]*m_z12;

  Box2Abb[5800]=-24. + Box2Abb[5799]*m_z12;

  Box2Abb[5801]=65. + 6.*Box2Abb[2481]*m_z12;

  Box2Abb[5802]=-212. + Box2Abb[5801]*m_z12;

  Box2Abb[5803]=9. + Box2Abb[5802]*m_z12;

  Box2Abb[5804]=279. + 76.*m_z12;

  Box2Abb[5805]=227. - Box2Abb[5804]*m_z12;

  Box2Abb[5806]=140. + Box2Abb[5805]*m_z12;

  Box2Abb[5807]=-60. + Box2Abb[5806]*m_z12;

  Box2Abb[5808]=-684. + 101.*m_z12;

  Box2Abb[5809]=467. + Box2Abb[5808]*m_z12;

  Box2Abb[5810]=-110. + Box2Abb[5809]*m_z12;

  Box2Abb[5811]=6. + 25.*m_z2k;

  Box2Abb[5812]=-2.*Box2Abb[5811] + Box2Abb[5800]*m_z12 - Box2Abb[5803]*m_z12*m_z2k + Box2Abb[5807]*m_z2k_2 + Box2Abb[5810]*m_z2k_3 + 168.*Box2Abb[693]*m_z12*m_z2k_4;

  Box2Abb[5813]=-2. + 14.*m_z12 - 47.*m_z2k;

  Box2Abb[5814]=2. + Box2Abb[5813]*m_z12 + 30.*m_z2k;

  Box2Abb[5815]=4. + Box2Abb[5814]*m_z12;

  Box2Abb[5816]=-9. + 32.*m_z2k;

  Box2Abb[5817]=52. + 58.*m_z2k;

  Box2Abb[5818]=40. - 159.*m_z2k;

  Box2Abb[5819]=-19. + Box2Abb[5818]*m_z2k;

  Box2Abb[5820]=-46. + Box2Abb[5819]*m_z12 + Box2Abb[5817]*m_z12_2 + 7.*m_z12_3 + 3.*Box2Abb[5816]*m_z2k;

  Box2Abb[5821]=20. + Box2Abb[5820]*m_z12 + 26.*m_z2k;

  Box2Abb[5822]=-1. + Box2Abb[2152]*m_z2k;

  Box2Abb[5823]=-4. + Box2Abb[5125]*m_z2k;

  Box2Abb[5824]=-5. + m_z2k + 3.*m_z2k_2;

  Box2Abb[5825]=1. + Box2Abb[5824]*m_z2k;

  Box2Abb[5826]=-Box2Abb[444]*pow(Box2Abb[61],4.) + Box2Abb[5823]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[4894]*pow(Box2Abb[61],2.)*m_z12_2 + 4.*Box2Abb[5825]*m_z12_3 + Box2Abb[5822]*m_z12_4;

  Box2Abb[5827]=-7. + 3.*Box2Abb[1941]*m_z2k;

  Box2Abb[5828]=1. + m_z2k + Box2Abb[1133]*Box2Abb[444]*m_z2k_2;

  Box2Abb[5829]=-69. + 4.*m_z2k + 30.*m_z2k_2;

  Box2Abb[5830]=8. + Box2Abb[5829]*m_z2k;

  Box2Abb[5831]=-1. + Box2Abb[5830]*m_z2k;

  Box2Abb[5832]=-16. + 3.*m_z2k;

  Box2Abb[5833]=15. + 4.*Box2Abb[5832]*m_z2k;

  Box2Abb[5834]=2. + Box2Abb[5833]*m_z2k;

  Box2Abb[5835]=3. + Box2Abb[5834]*m_z2k;

  Box2Abb[5836]=68. + 25.*m_z2k;

  Box2Abb[5837]=-72. + Box2Abb[5836]*m_z2k;

  Box2Abb[5838]=2. + Box2Abb[5837]*m_z2k;

  Box2Abb[5839]=-3. + Box2Abb[5838]*m_z2k;

  Box2Abb[5840]=-Box2Abb[5831]*pow(Box2Abb[61],3.)*m_z12_2 - Box2Abb[5839]*pow(Box2Abb[61],2.)*m_z12_3 + Box2Abb[5835]*Box2Abb[61]*m_z12_4 + Box2Abb[5828]*m_z12_5 - 2.*pow(Box2Abb[61],5.)*m_z2k + Box2Abb[5827]*pow(Box2Abb[61],4.)*m_z12*m_z2k;

  Box2Abb[5841]=-2. + m_z2k_2 - 25.*m_z2k_3 + 48.*m_z2k_4 - 23.*m_z2k_5;

  Box2Abb[5842]=12. - 17.*m_z2k;

  Box2Abb[5843]=-10. + Box2Abb[5842]*m_z2k;

  Box2Abb[5844]=5. + Box2Abb[5843]*m_z2k_2;

  Box2Abb[5845]=109. + 24.*m_z2k;

  Box2Abb[5846]=-64. + Box2Abb[5845]*m_z2k;

  Box2Abb[5847]=-3. + Box2Abb[5846]*m_z2k;

  Box2Abb[5848]=-4. + Box2Abb[5847]*m_z2k;

  Box2Abb[5849]=-301. + 6.*Box2Abb[2664]*m_z2k;

  Box2Abb[5850]=90. + Box2Abb[5849]*m_z2k;

  Box2Abb[5851]=-13. + Box2Abb[5850]*m_z2k;

  Box2Abb[5852]=10. + Box2Abb[5851]*m_z2k;

  Box2Abb[5853]=-305. + 67.*m_z2k;

  Box2Abb[5854]=284. + Box2Abb[5853]*m_z2k;

  Box2Abb[5855]=-40. + Box2Abb[5854]*m_z2k;

  Box2Abb[5856]=17. + Box2Abb[5855]*m_z2k;

  Box2Abb[5857]=-3. + Box2Abb[5856]*m_z2k;

  Box2Abb[5858]=-Box2Abb[5848]*pow(Box2Abb[61],3.)*m_z12 + Box2Abb[5852]*pow(Box2Abb[61],2.)*m_z12_2 - Box2Abb[5857]*Box2Abb[61]*m_z12_3 + 4.*Box2Abb[5841]*m_z12_4 + Box2Abb[5844]*m_z12_5 + 4.*Box2Abb[183]*pow(Box2Abb[61],4.)*m_z2k;

  Box2Abb[5859]=4. + 27.*m_z2k;

  Box2Abb[5860]=1. + Box2Abb[5859]*m_z2k;

  Box2Abb[5861]=10. + Box2Abb[488]*m_z2k;

  Box2Abb[5862]=10. + Box2Abb[5861]*m_z2k;

  Box2Abb[5863]=-57. + 53.*m_z2k;

  Box2Abb[5864]=20. + m_z2k + 2.*Box2Abb[5863]*m_z2k_2;

  Box2Abb[5865]=-5. + Box2Abb[5864]*m_z2k;

  Box2Abb[5866]=965. - 488.*m_z2k + 42.*m_z2k_2;

  Box2Abb[5867]=-512. + Box2Abb[5866]*m_z2k;

  Box2Abb[5868]=-39. + 72.*m_z2k + Box2Abb[5867]*m_z2k_3;

  Box2Abb[5869]=-297. + 313.*m_z2k;

  Box2Abb[5870]=35. + Box2Abb[5869]*m_z2k;

  Box2Abb[5871]=11. + Box2Abb[5870]*m_z2k;

  Box2Abb[5872]=-8. + Box2Abb[5871]*m_z2k;

  Box2Abb[5873]=-625. + 241.*m_z2k;

  Box2Abb[5874]=375. + Box2Abb[5873]*m_z2k;

  Box2Abb[5875]=35. + Box2Abb[5874]*m_z2k;

  Box2Abb[5876]=-60. + Box2Abb[5875]*m_z2k;

  Box2Abb[5877]=46. + Box2Abb[5876]*m_z2k;

  Box2Abb[5878]=Box2Abb[5872]*Box2Abb[61]*m_z12 + Box2Abb[5868]*m_z12_2 + Box2Abb[5877]*m_z12_3 + Box2Abb[5865]*m_z12_4 - Box2Abb[5862]*m_z12_5 - 2.*Box2Abb[5860]*pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[5879]=Box2Abb[5840]*Box2Abb[61]*m_x_2 + Box2Abb[5858]*m_x_3 + Box2Abb[5878]*m_x_4 + Box2Abb[5797]*m_x_5 + Box2Abb[5812]*m_x_6 + Box2Abb[5780]*m_x_7 - Box2Abb[5821]*m_x_8 + Box2Abb[5815]*m_x_9 + 2.*Box2Abb[155]*m_x_10*m_z12 + Box2Abb[5826]*pow(Box2Abb[61],3.)*m_x*m_z12*m_z2k - pow(Box2Abb[5],3.)*pow(Box2Abb[61],5.)*m_z12_2*m_z2k_2;

  Box2Abb[5880]=Box2Abb[5768]*Box2Abb[8]*m_cL - 2.*Box2Abb[5879]*m_cR;

  Box2Abb[5945]=(-4.*m_m)/(m_s_2*m_z12);

  Box2Abb[5946]=(4.*m_cR*m_m*m_x)/m_s_2;

  Box2Abb[5947]=(4.*m_m*m_x)/m_s_2;

  Box2Abb[5948]=(4.*m_m*m_x)/(m_s_2*m_z2k);

  Box2Abb[5949]=(4.*m_m)/(m_s_2*m_z12);

  Box2Abb[5950]=-m_m;

  Box2Abb[5951]=(-4.*m_m)/m_s;

  Box2Abb[5952]=(4.*m_m)/m_s;

  Box2Abb[5953]=-m_m*m_s;

  Box2Abb[5954]=-4.*m_m;

  Box2Abb[5955]=4.*m_m;

  Box2Abb[5956]=(-4.*m_m)/(m_s_2*m_z2k);

  Box2Abb[5957]=2.*Box2Abb[2]*m_cL + Box2Abb[72]*m_cR;

  Box2Abb[5958]=-Box2Abb[609]*Box2Abb[72]*m_x + Box2Abb[5957]*m_z12*m_z2k + 2.*Box2Abb[50]*m_z12*m_z2k_2;

  Box2Abb[5959]=-10. + 11.*m_z12;

  Box2Abb[5960]=5. + m_z12 - 4.*m_z12_2 + 20.*m_z2k - 23.*m_z12*m_z2k;

  Box2Abb[5961]=-2. + Box2Abb[5960]*m_z12;

  Box2Abb[5962]=Box2Abb[5961]*m_x_2 + Box2Abb[3125]*m_x*m_z12 + Box2Abb[5959]*m_x_3*m_z12 - Box2Abb[3123]*m_z12_2*m_z2k;

  Box2Abb[5963]=-Box2Abb[3120]*m_x_3 + Box2Abb[3115]*m_x_2*m_z12 - Box2Abb[872]*m_x_4*m_z12 - Box2Abb[3118]*m_x*m_z12_2*m_z2k + Box2Abb[397]*m_z12_3*m_z2k_2;

  Box2Abb[5964]=Box2Abb[5963]*m_cR + Box2Abb[5962]*m_cL*m_x;

  Box2Abb[5965]=10. - 3.*m_z12;

  Box2Abb[5966]=3. - 7.*m_z2k;

  Box2Abb[5967]=4. - Box2Abb[3134]*m_z12 + Box2Abb[5966]*m_z2k;

  Box2Abb[5968]=9. + m_z2k;

  Box2Abb[5969]=3. - Box2Abb[5968]*m_z12 + 2.*m_z12_2 - 20.*m_z2k;

  Box2Abb[5970]=2. + Box2Abb[5969]*m_z12;

  Box2Abb[5971]=Box2Abb[5970]*m_x_3 + Box2Abb[3133]*m_x_2*m_z12 + Box2Abb[5965]*m_x_4*m_z12 + Box2Abb[5967]*m_x*m_z12_2*m_z2k + 2.*Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[5972]=-Box2Abb[3146]*m_x_2 + Box2Abb[3143]*m_x*m_z12 + Box2Abb[3124]*m_x_3*m_z12 - 3.*pow(Box2Abb[5],2.)*m_z12_2*m_z2k;

  Box2Abb[5973]=Box2Abb[5971]*m_cR + Box2Abb[5972]*m_cL*m_x;

  Box2Abb[5974]=-Box2Abb[3156]*m_x_3 + Box2Abb[3151]*m_x_2*m_z12 + Box2Abb[3124]*m_x_4*m_z12 - Box2Abb[3153]*Box2Abb[5]*m_x*m_z12_2 + 2.*pow(Box2Abb[5],2.)*m_z12_3*m_z2k;

  Box2Abb[5975]=Box2Abb[5974]*m_cL + Box2Abb[5971]*m_cR;

  Box2Abb[5976]=-pow(Box2Abb[4],2.)*Box2Abb[62] - 3.*Box2Abb[4]*Box2Abb[815]*m_z2k + 5.*Box2Abb[70]*m_z2k_2;

  Box2Abb[5977]=Box2Abb[3162]*m_x_2 + Box2Abb[5976]*m_x*m_z12 - 5.*Box2Abb[72]*m_x_3*m_z12 - Box2Abb[3159]*Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[5978]=15. + 13.*Box2Abb[72]*m_z12;

  Box2Abb[5979]=-26. + 21.*m_z12;

  Box2Abb[5980]=-2. + Box2Abb[5978]*m_z12 + 13.*m_z2k + 2.*Box2Abb[5979]*m_z12*m_z2k + 5.*Box2Abb[550]*m_z2k_2;

  Box2Abb[5981]=-Box2Abb[3170]*m_x_3 + Box2Abb[5980]*m_x_2*m_z12 + 5.*Box2Abb[155]*m_x_4*m_z12 - Box2Abb[3164]*Box2Abb[5]*m_x*m_z12_2 + pow(Box2Abb[5],2.)*m_z12_3*m_z2k;

  Box2Abb[5982]=Box2Abb[5981]*m_cL + Box2Abb[5977]*m_cR*m_x;

  Box2Abb[5983]=-5. - 4.*m_x + 2.*m_z12;

  Box2Abb[5984]=4. + Box2Abb[5983]*m_z12;

  Box2Abb[5985]=-1. + Box2Abb[5984]*m_x + m_z12;

  Box2Abb[5986]=4. + Box2Abb[166]*m_z12_2;

  Box2Abb[5987]=Box2Abb[5986]*m_x - 2.*Box2Abb[4]*Box2Abb[72]*m_z12 - 12.*m_x_2*m_z12;

  Box2Abb[5988]=3. - 4.*m_z12;

  Box2Abb[5989]=-2. + Box2Abb[5988]*m_z12;

  Box2Abb[5990]=Box2Abb[5989]*m_x + Box2Abb[4]*Box2Abb[72]*m_z12 + 12.*m_x_2*m_z12;

  Box2Abb[5991]=Box2Abb[1025]*m_x + Box2Abb[155]*m_z12;

  Box2Abb[5992]=Box2Abb[5985]*Box2Abb[72]*m_x_2 + Box2Abb[5987]*m_x*m_z2k + Box2Abb[5990]*m_z12*m_z2k_2 + Box2Abb[5991]*m_z12*m_z2k_3;

  Box2Abb[5993]=11. - 12.*m_z12 - 48.*m_z2k;

  Box2Abb[5994]=2.*Box2Abb[1173] + Box2Abb[5993]*m_z12;

  Box2Abb[5995]=8. + Box2Abb[5994]*m_z12;

  Box2Abb[5996]=-7. + Box2Abb[3186]*m_z12 - 8.*m_z2k;

  Box2Abb[5997]=2. + Box2Abb[5996]*m_z12 - 4.*m_z2k;

  Box2Abb[5998]=Box2Abb[5997]*m_x_2 + Box2Abb[5995]*m_x_3 + 4.*Box2Abb[550]*m_x_4*m_z12 - Box2Abb[3182]*m_x*m_z12*m_z2k + 3.*Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[5999]=Box2Abb[5998]*m_cL + Box2Abb[5992]*m_cR;

  Box2Abb[6000]=-27. + 8.*m_z12;

  Box2Abb[6001]=26. + Box2Abb[6000]*m_z12 - 12.*m_z2k;

  Box2Abb[6002]=-8. + Box2Abb[6001]*m_z12;

  Box2Abb[6003]=2. - 2.*Box2Abb[1627]*m_z12 - 6.*m_z2k + Box2Abb[1048]*m_z12*m_z2k + 4.*Box2Abb[1424]*m_z2k_2;

  Box2Abb[6004]=-11. - 2.*Box2Abb[201]*m_z12 + 2.*m_z2k + 5.*m_z12*m_z2k + 12.*m_z2k_2;

  Box2Abb[6005]=7. + Box2Abb[6004]*m_z12 - 8.*m_z2k;

  Box2Abb[6006]=-2. + Box2Abb[6005]*m_z12 + 4.*m_z2k;

  Box2Abb[6007]=Box2Abb[6006]*m_x_2 + Box2Abb[6002]*m_x_3 - 4.*Box2Abb[72]*m_x_4*m_z12 + Box2Abb[6003]*m_x*m_z12*m_z2k + Box2Abb[5]*m_z12_3*m_z2k_2;

  Box2Abb[6008]=45. - 34.*m_z12 - 48.*m_z2k;

  Box2Abb[6009]=-22. + Box2Abb[6008]*m_z12 + 12.*m_z2k;

  Box2Abb[6010]=8. + Box2Abb[6009]*m_z12;

  Box2Abb[6011]=-17. + Box2Abb[3206]*m_z12;

  Box2Abb[6012]=6. + Box2Abb[6011]*m_z12 - 4.*m_z2k;

  Box2Abb[6013]=Box2Abb[6012]*m_x_2 + Box2Abb[6010]*m_x_3 - Box2Abb[3201]*Box2Abb[5]*m_x*m_z12 + 4.*Box2Abb[550]*m_x_4*m_z12 + pow(Box2Abb[5],2.)*Box2Abb[70]*m_z12_2*m_z2k;

  Box2Abb[6014]=Box2Abb[6013]*m_cL + Box2Abb[6007]*m_cR;

  Box2Abb[6015]=Box2Abb[429]*m_cL + Box2Abb[435]*m_cR;

  Box2Abb[6016]=Box2Abb[1340]*m_cL - Box2Abb[429]*pow(Box2Abb[8],2.)*m_cR*m_z2k;

  Box2Abb[6017]=Box2Abb[3215]*Box2Abb[3226]*m_cL + Box2Abb[3217]*m_cR;

  Box2Abb[6018]=m_x - 3.*Box2Abb[4]*m_z12 + 3.*m_x*m_z12;

  Box2Abb[6019]=Box2Abb[6018]*m_cR - Box2Abb[9]*m_cL*m_x + Box2Abb[72]*m_cL*m_z12;

  Box2Abb[6020]=m_cL - 2.*m_cR;

  Box2Abb[6021]=pow(Box2Abb[1],2.)*Box2Abb[3214]*Box2Abb[4]*m_x - Box2Abb[6017]*m_z2k + Box2Abb[6019]*m_z2k_2 + Box2Abb[6020]*m_z12*m_z2k_3;

  Box2Abb[6022]=1. - m_x - 2.*m_z12 + 3.*m_x*m_z12;

  Box2Abb[6023]=Box2Abb[1]*Box2Abb[6022] + Box2Abb[3234]*m_z2k + Box2Abb[1485]*m_z2k_2;

  Box2Abb[6024]=-pow(Box2Abb[1],2.)*Box2Abb[4]*m_x + Box2Abb[3215]*Box2Abb[3216]*m_z2k + Box2Abb[3220]*m_z2k_2 - m_z12*m_z2k_3;

  Box2Abb[6025]=Box2Abb[6024]*m_cR + Box2Abb[6023]*m_cL*m_x;

  Box2Abb[6026]=-pow(Box2Abb[1],2.)*Box2Abb[4]*m_x + Box2Abb[3215]*Box2Abb[3226]*m_z2k + Box2Abb[3220]*m_z2k_2 - m_z12*m_z2k_3;

  Box2Abb[6027]=Box2Abb[6026]*m_cR + Box2Abb[3235]*m_cL*m_x;

  Box2Abb[6028]=-1. + 4.*m_z12;

  Box2Abb[6029]=2. - Box2Abb[3237]*m_z12 + 4.*m_z2k;

  Box2Abb[6030]=2. + Box2Abb[854]*m_z2k;

  Box2Abb[6031]=Box2Abb[3245]*Box2Abb[61]*m_z12 + Box2Abb[6030]*m_z12_2 - pow(Box2Abb[61],2.)*m_z2k;

  Box2Abb[6032]=-5. + 28.*m_z2k;

  Box2Abb[6033]=-4. + Box2Abb[6032]*m_z2k;

  Box2Abb[6034]=-1. + Box2Abb[6033]*m_z12 + Box2Abb[3239]*m_z12_2 - 2.*Box2Abb[256]*m_z2k;

  Box2Abb[6035]=4. - 9.*m_z2k;

  Box2Abb[6036]=-3. + Box2Abb[6035]*m_z2k;

  Box2Abb[6037]=19. - 22.*m_z2k;

  Box2Abb[6038]=-11. + Box2Abb[6037]*m_z2k;

  Box2Abb[6039]=4. + Box2Abb[6038]*m_z2k;

  Box2Abb[6040]=Box2Abb[6039]*m_z12 + Box2Abb[6036]*m_z12_2 + 2.*Box2Abb[4189]*m_z2k;

  Box2Abb[6041]=Box2Abb[6040]*m_x_2 + Box2Abb[6034]*m_x_3 + Box2Abb[6029]*m_x_4 + Box2Abb[6028]*m_x_5 + Box2Abb[6031]*m_x*m_z2k - Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[6042]=Box2Abb[6041]*m_cL - Box2Abb[3235]*pow(Box2Abb[8],2.)*m_cR*m_x;

  Box2Abb[6043]=-Box2Abb[1341]*m_cR + Box2Abb[1342]*m_cL*m_x;

  Box2Abb[6044]=Box2Abb[3262]*m_cL + Box2Abb[1342]*pow(Box2Abb[8],2.)*m_cR;

  Box2Abb[6045]=Box2Abb[3274]*m_cL + Box2Abb[3266]*Box2Abb[8]*m_cR;

  Box2Abb[6046]=-Box2Abb[3280]*pow(Box2Abb[8],2.)*m_cR + Box2Abb[3287]*m_cL*m_x;

  Box2Abb[6047]=pow(Box2Abb[4],2.)*Box2Abb[550] + 3.*Box2Abb[1485]*Box2Abb[4]*m_z2k + 3.*m_z12*m_z2k_2;

  Box2Abb[6048]=Box2Abb[6047]*m_x - 3.*Box2Abb[3290]*m_x_2 - pow(Box2Abb[5],3.)*m_z12 + m_x_3*m_z12;

  Box2Abb[6049]=-Box2Abb[148]*Box2Abb[3280]*m_cR + Box2Abb[6048]*m_cL*m_x;

  Box2Abb[6050]=-15. + 2.*m_z12;

  Box2Abb[6051]=15. + Box2Abb[6050]*m_z12;

  Box2Abb[6052]=-2. + Box2Abb[6051]*m_z12 + 6.*m_z2k - 6.*Box2Abb[70]*m_z12*m_z2k - 4.*m_z2k_2;

  Box2Abb[6053]=4. - Box2Abb[3300]*m_z12 + 12.*m_z2k;

  Box2Abb[6054]=-24. + 11.*m_z12 + 36.*m_z2k;

  Box2Abb[6055]=24. + Box2Abb[6054]*m_z12 - 20.*m_z2k;

  Box2Abb[6056]=-8. + Box2Abb[6055]*m_z12;

  Box2Abb[6057]=Box2Abb[6053]*m_x_3 + Box2Abb[6056]*m_x_4 + 4.*Box2Abb[699]*m_x_5*m_z12 + Box2Abb[3308]*m_x_2*m_z2k + Box2Abb[6052]*m_x*m_z12*m_z2k_2 + 3.*Box2Abb[5]*m_z12_3*m_z2k_3;

  Box2Abb[6058]=Box2Abb[6057]*m_cL + 2.*Box2Abb[3296]*Box2Abb[382]*Box2Abb[8]*m_cR;

  Box2Abb[6059]=-2.*Box2Abb[3332]*Box2Abb[388]*m_cL + Box2Abb[3325]*m_cR;

  Box2Abb[6060]=-2. + 7.*m_z12;

  Box2Abb[6061]=-2.*Box2Abb[12] + Box2Abb[3851]*m_z12 + 2.*Box2Abb[1116]*m_z12*m_z2k + 2.*Box2Abb[70]*m_z2k_2;

  Box2Abb[6062]=4. - Box2Abb[3340]*m_z12;

  Box2Abb[6063]=Box2Abb[6062]*m_x_2 + Box2Abb[6061]*m_x*m_z12 + 2.*Box2Abb[6060]*m_x_3*m_z12 - Box2Abb[3336]*Box2Abb[5]*m_z12_2*m_z2k;

  Box2Abb[6064]=Box2Abb[6063]*m_cL - Box2Abb[3334]*Box2Abb[382]*m_cR;

  Box2Abb[6065]=-Box2Abb[3346]*Box2Abb[382]*Box2Abb[8]*m_cR + Box2Abb[3356]*m_cL*m_x;

  Box2Abb[6066]=Box2Abb[148]*Box2Abb[3346]*Box2Abb[382]*m_cR - Box2Abb[3362]*Box2Abb[8]*m_cL*m_x;

  Box2Abb[6067]=Box2Abb[3384]*m_cL - 2.*Box2Abb[3371]*m_cR;

  Box2Abb[6068]=pow(Box2Abb[4],3.)*m_z12_2 + Box2Abb[3401]*pow(Box2Abb[4],2.)*m_z2k + Box2Abb[3403]*Box2Abb[4]*m_z2k_2 + 6.*Box2Abb[699]*m_z12*m_z2k_3 - 6.*m_z12*m_z2k_4;

  Box2Abb[6069]=89. - 40.*m_z12;

  Box2Abb[6070]=-63. + Box2Abb[6069]*m_z12;

  Box2Abb[6071]=16. + Box2Abb[6070]*m_z12;

  Box2Abb[6072]=-3.*Box2Abb[155]*pow(Box2Abb[4],2.)*m_z12 + Box2Abb[6071]*m_z2k - 4.*Box2Abb[3410]*m_z2k_2 + 4.*m_z12*m_z2k_3;

  Box2Abb[6073]=Box2Abb[6072]*m_z12 - 2.*m_z2k;

  Box2Abb[6074]=Box2Abb[6073]*m_x_3 + Box2Abb[3417]*m_x_4 + Box2Abb[6068]*m_x_2*m_z12 - 6.*Box2Abb[3406]*m_x_5*m_z12 + 2.*m_x_6*m_z12_2 - Box2Abb[3400]*Box2Abb[5]*m_x*m_z12_2*m_z2k - pow(Box2Abb[5],3.)*m_z12_3*m_z2k_2;

  Box2Abb[6075]=2.*Box2Abb[6074]*m_cL - Box2Abb[3399]*Box2Abb[8]*m_cR;

  Box2Abb[6076]=1. + Box2Abb[240]*m_x + 2.*m_x_2 - m_z12 + 3.*Box2Abb[4]*m_z2k + 2.*m_z2k_2;

  Box2Abb[6077]=-Box2Abb[3423]*m_cR + Box2Abb[6076]*m_cL*m_x;

  Box2Abb[6078]=Box2Abb[3495]*m_cL + Box2Abb[3463]*m_cR;

  Box2Abb[6079]=-11. - 2.*m_z12 + 3.*Box2Abb[3402]*m_z2k;

  Box2Abb[6080]=-1. + 15.*m_z12;

  Box2Abb[6081]=13. - Box2Abb[3638]*m_z12;

  Box2Abb[6082]=16. + Box2Abb[6081]*m_z12;

  Box2Abb[6083]=-29. + 22.*m_z12 + 6.*m_z12_2;

  Box2Abb[6084]=-9. + 11.*m_z12 + m_z2k + 2.*Box2Abb[6080]*m_z12*m_z2k + Box2Abb[6082]*m_z2k_2 + 2.*Box2Abb[6083]*m_z2k_3 + 15.*Box2Abb[3452]*m_z2k_4;

  Box2Abb[6085]=-1. + m_z2k - 40.*m_z2k_2;

  Box2Abb[6086]=56. + 25.*m_z2k;

  Box2Abb[6087]=41. + Box2Abb[6086]*m_z2k;

  Box2Abb[6088]=20. + Box2Abb[6087]*m_z2k;

  Box2Abb[6089]=9. - Box2Abb[6088]*m_z12 + 2.*Box2Abb[6085]*m_z2k + 2.*Box2Abb[2534]*m_z12_2*m_z2k;

  Box2Abb[6090]=-Box2Abb[3545]*m_x_2 + Box2Abb[6084]*m_x_3 + Box2Abb[6089]*m_x_4 + Box2Abb[3534]*m_x_5 + Box2Abb[6079]*m_x_6 - Box2Abb[201]*m_x_7 + Box2Abb[3552]*m_x*m_z2k - Box2Abb[3546]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[6091]=11. - Box2Abb[3500]*m_z12 + 24.*m_z2k;

  Box2Abb[6092]=3. + Box2Abb[4068]*m_z2k;

  Box2Abb[6093]=31. - 115.*m_z2k;

  Box2Abb[6094]=9. + Box2Abb[6093]*m_z2k;

  Box2Abb[6095]=5. + Box2Abb[6094]*m_z2k;

  Box2Abb[6096]=-1. + Box2Abb[6095]*m_z12 - 2.*Box2Abb[3515]*m_z12_2 + 2.*Box2Abb[6092]*m_z2k;

  Box2Abb[6097]=-15. - 22.*m_z2k + 68.*m_z2k_2;

  Box2Abb[6098]=6. + Box2Abb[6097]*m_z2k;

  Box2Abb[6099]=-8. + Box2Abb[3511]*m_z2k;

  Box2Abb[6100]=-2.*Box2Abb[3509] + Box2Abb[6099]*m_z12 + Box2Abb[6098]*m_z12_2 + 5.*m_z12_3*m_z2k_2;

  Box2Abb[6101]=Box2Abb[3522]*pow(Box2Abb[61],2.) - Box2Abb[3525]*Box2Abb[61]*m_z12 - 2.*Box2Abb[3529]*m_z12_2 + Box2Abb[1612]*m_z12_3*m_z2k_2;

  Box2Abb[6102]=Box2Abb[6101]*m_x_2 + Box2Abb[6100]*m_x_3 + Box2Abb[6096]*m_x_4 + Box2Abb[3504]*m_x_5 + Box2Abb[6091]*m_x_6 + Box2Abb[1324]*m_x_7 + Box2Abb[3499]*Box2Abb[5]*Box2Abb[61]*m_x*m_z2k - pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[6103]=Box2Abb[6102]*m_cL + Box2Abb[6090]*m_cR;

  Box2Abb[6104]=Box2Abb[3600]*m_cL + Box2Abb[3562]*m_cR;

  Box2Abb[6105]=Box2Abb[3635]*m_cL + Box2Abb[3666]*m_cR;

  Box2Abb[6106]=5. + Box2Abb[3720]*m_z12;

  Box2Abb[6107]=60. + Box2Abb[3567]*m_z12;

  Box2Abb[6108]=-14. + Box2Abb[6107]*m_z12;

  Box2Abb[6109]=11. + Box2Abb[6108]*m_z12;

  Box2Abb[6110]=-13. + 20.*m_z12 - 10.*m_z12_2 + 6.*m_z12_3 - 2.*Box2Abb[6106]*m_z2k + Box2Abb[6109]*m_z2k_2 + 4.*Box2Abb[3726]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[3727]*m_z12*m_z2k_4;

  Box2Abb[6111]=-7. + Box2Abb[1046]*m_z12;

  Box2Abb[6112]=21. + Box2Abb[4905]*m_z12;

  Box2Abb[6113]=-20. + Box2Abb[6112]*m_z12;

  Box2Abb[6114]=56. - 11.*m_z12;

  Box2Abb[6115]=-86. + Box2Abb[6114]*m_z12;

  Box2Abb[6116]=34. + Box2Abb[6115]*m_z12;

  Box2Abb[6117]=-29. + Box2Abb[6116]*m_z12;

  Box2Abb[6118]=51. - 26.*m_z12;

  Box2Abb[6119]=-2. + Box2Abb[6118]*m_z12;

  Box2Abb[6120]=7. + Box2Abb[6119]*m_z12;

  Box2Abb[6121]=4. + Box2Abb[6111]*m_z12 - 9.*m_z2k + 9.*Box2Abb[201]*Box2Abb[4]*m_z12*m_z2k + Box2Abb[4]*Box2Abb[6113]*m_z2k_2 + Box2Abb[6117]*m_z2k_3 + 2.*Box2Abb[6120]*m_z2k_4 - 9.*Box2Abb[1637]*m_z12*m_z2k_5;

  Box2Abb[6122]=6. - Box2Abb[3734]*m_z12;

  Box2Abb[6123]=32. - 115.*m_z2k;

  Box2Abb[6124]=-20. + Box2Abb[6123]*m_z2k;

  Box2Abb[6125]=5. + Box2Abb[6124]*m_z2k;

  Box2Abb[6126]=-27. + Box2Abb[1667]*m_z2k;

  Box2Abb[6127]=-37. + 2.*Box2Abb[6126]*m_z2k;

  Box2Abb[6128]=26. + Box2Abb[6127]*m_z12 + Box2Abb[6125]*m_z12_2 - Box2Abb[1590]*m_z12_3 + 11.*Box2Abb[234]*m_z2k;

  Box2Abb[6129]=Box2Abb[6121]*m_x_2 + Box2Abb[6110]*m_x_3 + Box2Abb[6128]*m_x_4 + Box2Abb[3733]*m_x_5 + Box2Abb[6122]*m_x_6 + Box2Abb[550]*m_x_7*m_z12 + Box2Abb[3740]*Box2Abb[61]*m_x*m_z2k - Box2Abb[3729]*Box2Abb[5]*pow(Box2Abb[61],3.)*m_z12*m_z2k_2;

  Box2Abb[6130]=17. - 27.*m_z12 + 11.*m_z12_2 + 2.*Box2Abb[3669]*m_z2k - Box2Abb[3673]*m_z2k_2 + 4.*Box2Abb[3675]*Box2Abb[4]*m_z2k_3 + 5.*Box2Abb[3731]*m_z12*m_z2k_4;

  Box2Abb[6131]=-31. + 14.*m_z12 - 2.*Box2Abb[3680]*m_z2k + 3.*Box2Abb[669]*m_z2k_2;

  Box2Abb[6132]=23. + Box2Abb[6131]*m_z12 + 26.*m_z2k;

  Box2Abb[6133]=4. + Box2Abb[514]*m_z12 - 10.*m_z2k;

  Box2Abb[6134]=-6. + Box2Abb[6133]*m_z12;

  Box2Abb[6135]=32. - 5.*m_z2k;

  Box2Abb[6136]=24. + Box2Abb[6135]*m_z2k;

  Box2Abb[6137]=47. + 2.*Box2Abb[6136]*m_z2k;

  Box2Abb[6138]=-31. + Box2Abb[6137]*m_z12 - Box2Abb[3689]*m_z12_2 - Box2Abb[3687]*m_z2k + 2.*Box2Abb[2534]*m_z12_3*m_z2k;

  Box2Abb[6139]=-Box2Abb[3716]*m_x_2 + Box2Abb[6130]*m_x_3 + Box2Abb[6138]*m_x_4 + Box2Abb[6132]*m_x_5 + Box2Abb[6134]*m_x_6 - Box2Abb[72]*m_x_7*m_z12 + Box2Abb[3703]*m_x*m_z2k - Box2Abb[3686]*pow(Box2Abb[61],2.)*m_z12*m_z2k_2;

  Box2Abb[6140]=Box2Abb[6129]*m_cL + Box2Abb[6139]*m_cR;

  Box2Abb[6141]=6. - Box2Abb[3783]*m_z12;

  Box2Abb[6142]=19. + 44.*m_z2k;

  Box2Abb[6143]=-62. + 5.*m_z2k;

  Box2Abb[6144]=9. + Box2Abb[6143]*m_z2k;

  Box2Abb[6145]=29. + 2.*Box2Abb[1558]*m_z2k;

  Box2Abb[6146]=62. - 115.*m_z2k;

  Box2Abb[6147]=-63. + Box2Abb[6146]*m_z2k;

  Box2Abb[6148]=15. + Box2Abb[6147]*m_z2k;

  Box2Abb[6149]=3. + m_z12 + Box2Abb[6148]*m_z12_2 - Box2Abb[6145]*m_z12_3 - 2.*m_z12_4 + Box2Abb[6142]*m_z2k + 2.*Box2Abb[6144]*m_z12*m_z2k;

  Box2Abb[6150]=Box2Abb[3792]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x - Box2Abb[3803]*m_x_2 + Box2Abb[3781]*m_x_3 + Box2Abb[6149]*m_x_4 + Box2Abb[3788]*m_x_5 + Box2Abb[6141]*m_x_6 + Box2Abb[550]*m_x_7*m_z12 - pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*Box2Abb[70]*m_z12*m_z2k;

  Box2Abb[6151]=Box2Abb[6150]*m_cL - Box2Abb[3847]*m_cR;

  Box2Abb[6152]=Box2Abb[431]*m_cL*m_x - Box2Abb[429]*Box2Abb[8]*m_cR*m_z2k;

  Box2Abb[6153]=Box2Abb[3874]*m_cL - Box2Abb[3888]*Box2Abb[8]*m_cR;

  Box2Abb[6154]=Box2Abb[3927]*Box2Abb[8]*m_cL + Box2Abb[3908]*m_cR*m_z2k;

  Box2Abb[6155]=Box2Abb[3961]*m_cL + Box2Abb[3908]*m_cR*m_z2k;

  Box2Abb[6156]=Box2Abb[4033]*m_cL + Box2Abb[3991]*m_cR*m_z2k;

  Box2Abb[6157]=Box2Abb[4081]*m_cL + Box2Abb[4053]*m_cR*m_z2k;

  Box2Abb[6158]=Box2Abb[4160]*m_cL + Box2Abb[4107]*Box2Abb[8]*m_cR*m_z2k;

  Box2Abb[6159]=4. - 7.*m_z12 - 6.*m_z12*m_z2k - 4.*Box2Abb[4]*m_z2k_2;

  Box2Abb[6160]=-2.*pow(Box2Abb[61],2.) + Box2Abb[4163]*Box2Abb[61]*m_z12 + 2.*Box2Abb[187]*m_z12_2*m_z2k;

  Box2Abb[6161]=Box2Abb[6160]*m_x + Box2Abb[6159]*m_x_2 + Box2Abb[783]*m_x_3 - 2.*Box2Abb[70]*m_x_4 - 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[6162]=4. + 6.*m_z12;

  Box2Abb[6163]=6. - Box2Abb[4166]*m_z12;

  Box2Abb[6164]=-2.*Box2Abb[4168]*m_z12 + 3.*m_z12_2 - 4.*Box2Abb[841]*m_z2k;

  Box2Abb[6165]=-3. + Box2Abb[2153]*m_z2k;

  Box2Abb[6166]=-2.*pow(Box2Abb[61],2.) + Box2Abb[61]*Box2Abb[6165]*m_z12 + Box2Abb[4056]*m_z12_2;

  Box2Abb[6167]=Box2Abb[6166]*m_x + Box2Abb[6164]*m_x_2 + Box2Abb[6163]*m_x_3 + Box2Abb[6162]*m_x_4 - 2.*Box2Abb[5]*pow(Box2Abb[61],2.)*m_z12*m_z2k;

  Box2Abb[6168]=Box2Abb[6167]*m_cL + Box2Abb[6161]*m_cR;

  Box2Abb[6169]=-6. + Box2Abb[1629]*Box2Abb[4]*m_z12;

  Box2Abb[6170]=-22. + Box2Abb[4230]*m_z12;

  Box2Abb[6171]=-4. + Box2Abb[6170]*m_z12;

  Box2Abb[6172]=-8. + Box2Abb[6171]*m_z12;

  Box2Abb[6173]=8. - 9.*m_z12;

  Box2Abb[6174]=2.*Box2Abb[6169] + Box2Abb[6172]*m_z2k - 2.*Box2Abb[4235]*m_z2k_2 + 5.*Box2Abb[6173]*m_z12*m_z2k_3;

  Box2Abb[6175]=-1. + Box2Abb[487]*m_z2k;

  Box2Abb[6176]=Box2Abb[4241]*pow(Box2Abb[61],4.) + Box2Abb[4243]*pow(Box2Abb[61],2.)*m_z12 + 2.*Box2Abb[4245]*Box2Abb[61]*m_z12_2 + Box2Abb[6175]*m_z12_3*m_z2k;

  Box2Abb[6177]=9. + Box2Abb[4248]*m_z12 - Box2Abb[285]*m_z12_2 + Box2Abb[4247]*m_z2k;

  Box2Abb[6178]=8.*Box2Abb[12] + Box2Abb[6177]*m_z12;

  Box2Abb[6179]=-2. + 3.*Box2Abb[1071]*m_z2k;

  Box2Abb[6180]=2. + Box2Abb[6179]*m_z2k;

  Box2Abb[6181]=46. - 33.*m_z2k;

  Box2Abb[6182]=46. + Box2Abb[6181]*m_z2k;

  Box2Abb[6183]=-108. + Box2Abb[6182]*m_z2k;

  Box2Abb[6184]=51. + Box2Abb[6183]*m_z2k;

  Box2Abb[6185]=-2. + Box2Abb[6184]*m_z2k_2;

  Box2Abb[6186]=-Box2Abb[230]*Box2Abb[510]*pow(Box2Abb[61],4.) + Box2Abb[6185]*m_z12 - Box2Abb[12]*Box2Abb[4268]*m_z12_2 + 2.*Box2Abb[6180]*m_z12_3*m_z2k + Box2Abb[424]*m_z12_4*m_z2k_2;

  Box2Abb[6187]=Box2Abb[4289]*m_x_3 + Box2Abb[4260]*m_x_4 + Box2Abb[6174]*m_x_5 + Box2Abb[6178]*m_x_6 + Box2Abb[4239]*m_x_7 + Box2Abb[6186]*m_x_2*m_z12 - 3.*Box2Abb[72]*m_x_8*m_z12 + Box2Abb[6176]*m_x*m_z12_2*m_z2k - Box2Abb[470]*pow(Box2Abb[5],2.)*pow(Box2Abb[61],2.)*m_z12_3*m_z2k_2;

  Box2Abb[6188]=Box2Abb[6187]*m_cR + Box2Abb[4229]*m_cL*m_x;

  Box2Abb[6189]=-2. + Box2Abb[4346]*m_z12;

  Box2Abb[6190]=28. - 33.*m_z12;

  Box2Abb[6191]=-14. + 6.*Box2Abb[203]*m_z12 - 11.*m_z2k + Box2Abb[6190]*m_z12*m_z2k + Box2Abb[2287]*m_z2k_2;

  Box2Abb[6192]=8.*Box2Abb[12] + Box2Abb[6191]*m_z12;

  Box2Abb[6193]=-20. + Box2Abb[4363]*m_z2k;

  Box2Abb[6194]=38. + 2.*Box2Abb[4365]*m_z2k;

  Box2Abb[6195]=-4.*Box2Abb[892] + Box2Abb[6194]*m_z12 - Box2Abb[4368]*m_z12_2 + Box2Abb[6193]*m_z12_3 + 8.*m_z12_4*m_z2k;

  Box2Abb[6196]=-72. + 25.*m_z2k;

  Box2Abb[6197]=-46. + Box2Abb[6196]*m_z2k;

  Box2Abb[6198]=-11. + Box2Abb[6197]*m_z2k_2;

  Box2Abb[6199]=8.*Box2Abb[12]*pow(Box2Abb[61],2.) - 2.*Box2Abb[4356]*Box2Abb[61]*m_z12 + Box2Abb[6198]*m_z12_2 + 2.*Box2Abb[4361]*m_z12_3 - 2.*Box2Abb[4354]*m_z12_4*m_z2k;

  Box2Abb[6200]=-2.*pow(Box2Abb[61],4.) + Box2Abb[4386]*pow(Box2Abb[61],2.)*m_z12 + Box2Abb[4395]*m_z12_2 - 2.*Box2Abb[4390]*m_z12_3 - 2.*Box2Abb[4383]*m_z12_4*m_z2k + 4.*m_z12_5*m_z2k_2;

  Box2Abb[6201]=Box2Abb[6200]*m_x_3 + Box2Abb[6199]*m_x_4 + Box2Abb[6195]*m_x_5 + Box2Abb[6192]*m_x_6 + Box2Abb[6189]*m_x_7 + Box2Abb[4382]*m_x_2*m_z12 + Box2Abb[1087]*m_x_8*m_z12 - Box2Abb[4353]*Box2Abb[5]*Box2Abb[61]*m_x*m_z12_2*m_z2k - 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12_3*m_z2k_2;

  Box2Abb[6202]=Box2Abb[6201]*m_cR + Box2Abb[4345]*m_cL*m_x;

  Box2Abb[6203]=124. - 195.*m_z2k;

  Box2Abb[6204]=144. + Box2Abb[6203]*m_z2k;

  Box2Abb[6205]=108. + Box2Abb[6204]*m_z2k;

  Box2Abb[6206]=153. + Box2Abb[6205]*m_z2k;

  Box2Abb[6207]=-9. + 5.*Box2Abb[4451]*m_z2k;

  Box2Abb[6208]=-12. + Box2Abb[6207]*m_z2k;

  Box2Abb[6209]=-8.*Box2Abb[12]*pow(Box2Abb[61],2.) + 2.*Box2Abb[6208]*m_z12 + Box2Abb[6206]*m_z12_2 - 2.*Box2Abb[4450]*m_z12_3 - 2.*Box2Abb[4447]*m_z12_4 - 4.*m_z12_5*m_z2k;

  Box2Abb[6210]=Box2Abb[4445]*m_x_3 + Box2Abb[6209]*m_x_4 + Box2Abb[4431]*m_x_5 - Box2Abb[4407]*m_x_6 + Box2Abb[4402]*m_x_7 - Box2Abb[4422]*m_x_2*m_z12 + Box2Abb[3347]*m_x_8*m_z12 + Box2Abb[4410]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12_2 - 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_3*m_z2k;

  Box2Abb[6211]=Box2Abb[6210]*m_cL + Box2Abb[6201]*m_cR;

  Box2Abb[6212]=4. - 3.*m_z2k;

  Box2Abb[6213]=-26. + 7.*Box2Abb[6212]*m_z2k;

  Box2Abb[6214]=10. + Box2Abb[6213]*m_z2k;

  Box2Abb[6215]=-7. + Box2Abb[6214]*m_z2k;

  Box2Abb[6216]=-Box2Abb[4591]*pow(Box2Abb[61],3.)*m_z12 - Box2Abb[4584]*pow(Box2Abb[61],3.)*m_z12_2 - Box2Abb[4595]*Box2Abb[61]*m_z12_3 + Box2Abb[6215]*m_z12_4 + Box2Abb[4581]*pow(Box2Abb[61],4.)*m_z2k;

  Box2Abb[6217]=Box2Abb[4644]*m_x_3 - Box2Abb[4621]*m_x_4 + Box2Abb[4541]*m_x_5 - Box2Abb[4579]*m_x_6 + Box2Abb[4563]*m_x_7 - Box2Abb[4554]*m_x_8 + Box2Abb[4544]*m_x_9 + Box2Abb[710]*m_x_10*m_z12 + Box2Abb[6216]*m_x_2*m_z12*m_z2k + Box2Abb[4549]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12_2*m_z2k_2 - pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_3*m_z2k_3;

  Box2Abb[6218]=Box2Abb[6217]*m_cL + Box2Abb[4520]*pow(Box2Abb[8],2.)*m_cR*m_x;

  Box2Abb[6219]=-2.*Box2Abb[9] + 3.*Box2Abb[70]*m_z2k - 4.*m_z2k_2;

  Box2Abb[6220]=-2.*Box2Abb[230]*pow(Box2Abb[61],4.) + Box2Abb[4660]*pow(Box2Abb[61],3.)*m_z12 - 2.*Box2Abb[4663]*pow(Box2Abb[61],2.)*m_z12_2 + Box2Abb[1880]*Box2Abb[61]*m_z12_3 + 12.*m_z12_4*m_z2k_3;

  Box2Abb[6221]=-1. + m_z2k - 54.*m_z2k_2 + 136.*m_z2k_3 - 80.*m_z2k_4;

  Box2Abb[6222]=-7. + 31.*m_z2k;

  Box2Abb[6223]=-31. + 2.*Box2Abb[6222]*m_z2k;

  Box2Abb[6224]=-20. + Box2Abb[6223]*m_z2k;

  Box2Abb[6225]=32. + Box2Abb[6221]*m_z12 + 2.*Box2Abb[6224]*m_z12_2 + 3.*Box2Abb[1871]*m_z12_3 - 8.*Box2Abb[2193]*m_z2k_2;

  Box2Abb[6226]=28. - 127.*m_z2k;

  Box2Abb[6227]=84. + Box2Abb[6226]*m_z2k;

  Box2Abb[6228]=19. + Box2Abb[2153]*m_z2k;

  Box2Abb[6229]=-4.*Box2Abb[6228] + Box2Abb[4670]*m_z12 + Box2Abb[6227]*m_z12_2 + Box2Abb[4509]*m_z12_3;

  Box2Abb[6230]=1. - Box2Abb[4681]*m_z2k;

  Box2Abb[6231]=4.*Box2Abb[4672]*pow(Box2Abb[61],2.) + Box2Abb[4678]*m_z12 + Box2Abb[6230]*m_z12_2 + Box2Abb[1886]*m_z12_3 - 12.*m_z12_4*m_z2k_2;

  Box2Abb[6232]=Box2Abb[6220]*m_x + Box2Abb[6231]*m_x_2 + Box2Abb[6225]*m_x_3 + Box2Abb[6229]*m_x_4 - Box2Abb[4652]*m_x_5 + Box2Abb[4649]*m_x_6 + Box2Abb[5]*pow(Box2Abb[61],3.)*Box2Abb[6219]*m_z12*m_z2k;

  Box2Abb[6233]=-4. + Box2Abb[4687]*m_z12;

  Box2Abb[6234]=28. - 10.*Box2Abb[9]*m_z12 + Box2Abb[6233]*m_z2k + Box2Abb[4688]*m_z2k_2 - 10.*m_z12*m_z2k_3;

  Box2Abb[6235]=6. + Box2Abb[1065]*m_z12;

  Box2Abb[6236]=3. - 17.*m_z2k;

  Box2Abb[6237]=13. + Box2Abb[6236]*m_z12 - 2.*Box2Abb[841]*m_z2k;

  Box2Abb[6238]=-23. + Box2Abb[6237]*m_z12 - 14.*m_z2k;

  Box2Abb[6239]=-5. + Box2Abb[4696]*m_z2k;

  Box2Abb[6240]=-8. + 5.*Box2Abb[179]*m_z2k;

  Box2Abb[6241]=-4. + Box2Abb[6240]*m_z2k;

  Box2Abb[6242]=-7. + 2.*Box2Abb[6241]*m_z2k;

  Box2Abb[6243]=2.*Box2Abb[6239] + Box2Abb[6242]*m_z12 + Box2Abb[4698]*m_z12_2 - 2.*Box2Abb[3116]*m_z12_3*m_z2k;

  Box2Abb[6244]=3. + 4.*m_z2k_2;

  Box2Abb[6245]=23. + 4.*Box2Abb[5083]*m_z2k;

  Box2Abb[6246]=-10. + Box2Abb[6245]*m_z2k;

  Box2Abb[6247]=-6. + Box2Abb[6246]*m_z2k;

  Box2Abb[6248]=-2.*Box2Abb[4714]*pow(Box2Abb[61],2.) + Box2Abb[4722]*m_z12 + Box2Abb[6247]*m_z12_2 - 2.*Box2Abb[6244]*m_z12_3*m_z2k + 2.*m_z12_4*m_z2k_2;

  Box2Abb[6249]=Box2Abb[4713]*m_x_2 + Box2Abb[6248]*m_x_3 + Box2Abb[6243]*m_x_4 + Box2Abb[6234]*m_x_5 + Box2Abb[6238]*m_x_6 + Box2Abb[6235]*m_x_7 - 2.*m_x_8*m_z12 + Box2Abb[4691]*Box2Abb[61]*m_x*m_z12*m_z2k - pow(Box2Abb[5],2.)*pow(Box2Abb[61],3.)*m_z12_2*m_z2k_2;

  Box2Abb[6250]=2.*Box2Abb[6249]*m_cR + Box2Abb[6232]*m_cL*m_x;

  Box2Abb[6251]=-Box2Abb[4864]*m_x_3 - Box2Abb[4845]*m_x_4 + Box2Abb[4809]*m_x_5 - Box2Abb[4778]*m_x_6 + Box2Abb[4794]*m_x_7 - Box2Abb[4787]*m_x_8 + Box2Abb[4780]*m_x_9 + Box2Abb[4823]*m_x_2*m_z2k - Box2Abb[4783]*Box2Abb[5]*pow(Box2Abb[61],2.)*m_x*m_z12*m_z2k_2 + 2.*pow(Box2Abb[5],2.)*pow(Box2Abb[61],4.)*m_z12_2*m_z2k_3;

  Box2Abb[6252]=Box2Abb[6251]*m_cL - Box2Abb[4767]*pow(Box2Abb[8],2.)*m_cR*m_x;

  Box2Abb[6253]=Box2Abb[4883]*m_cL + Box2Abb[4902]*m_cR;

  Box2Abb[6254]=-Box2Abb[4950]*Box2Abb[8]*m_cR + Box2Abb[5054]*m_cL*m_x;

  Box2Abb[6255]=-5. + Box2Abb[1377]*m_z2k;

  Box2Abb[6256]=5. + 4.*Box2Abb[5088]*m_z12 + 30.*Box2Abb[1069]*Box2Abb[5089]*m_z12_2 + 4.*Box2Abb[5092]*m_z12_3 + 2.*Box2Abb[6255]*m_z12_4 + 5.*Box2Abb[2303]*m_z2k_3;

  Box2Abb[6257]=Box2Abb[5061]*m_x + Box2Abb[5105]*m_x_2 - Box2Abb[6256]*m_x_3 + Box2Abb[5082]*m_x_4 + Box2Abb[5073]*m_x_5 + Box2Abb[5064]*m_x_6 - Box2Abb[5057]*m_x_7 + pow(Box2Abb[5],3.)*pow(Box2Abb[61],5.)*m_z12 + m_x_8*m_z12;

  Box2Abb[6258]=-Box2Abb[5193]*m_cR + Box2Abb[6257]*Box2Abb[8]*m_cL*m_x;

  Box2Abb[6259]=-Box2Abb[5193]*Box2Abb[8]*m_cR + Box2Abb[5326]*m_cL*m_x;

  Box2Abb[6260]=Box2Abb[5527]*m_cR + Box2Abb[5405]*Box2Abb[8]*m_cL*m_x;

  Box2Abb[6261]=Box2Abb[5678]*m_cL - 2.*Box2Abb[5597]*m_cR;

  Box2Abb[6262]=-2.*Box2Abb[5879]*m_cL + Box2Abb[5768]*Box2Abb[8]*m_cR;

  Box2Abb[6326]=(4.*m_cL*m_m*m_x)/m_s_2;

  return;
}

DivArrC Z_Decay_RV_Box_2::RV_Box_2(const int& ME, const int& LR,
				   const Vec4C& epsV, const Vec4C& epsP)
{
  // Box corrections for emission off leg 2
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One
	*((m_p1*epsP)*((Box2Abb[3223]*Box2Abb[5949]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3231]*Box2Abb[5949]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[3231]*Box2Abb[5949]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3253]*Box2Abb[5949]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[3255]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3263]*Box2Abb[5947]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))
	
	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[3112]*Box2Abb[86]*(epsP*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (m_p1*epsP)*((Box2Abb[3130]*Box2Abb[5945]*(m_p1*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (Box2Abb[3148]*Box2Abb[5945]*(m_p2*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405])) + (Box2Abb[3158]*Box2Abb[5945]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (Box2Abb[3172]*Box2Abb[5945]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (m_pP*epsV)*((Box2Abb[319]*Box2Abb[3190]*(m_p1*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[405]) + (Box2Abb[319]*Box2Abb[3210]*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[405])))

	+B_0(0.,0.,m_m2,m_mu2)*((m_p1*epsP)*((Box2Abb[5946]*Box2Abb[707]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3211]*Box2Abb[600]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[3212]*Box2Abb[5947]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3213]*Box2Abb[5948]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[3211]*Box2Abb[600]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3213]*Box2Abb[5948]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[3424]*Box2Abb[4]*Box2Abb[86]*(epsP*epsV))/(Box2Abb[0]*Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[322]*Box2Abb[3496]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[3563]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))) + (Box2Abb[322]*Box2Abb[3601]*Box2Abb[4]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[3667]*Box2Abb[4]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[3761]*(m_p1*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[3848]*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[3849]*Box2Abb[86]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[321]*Box2Abb[3889]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[3928]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))) + (Box2Abb[322]*Box2Abb[3962]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[4034]*Box2Abb[5956]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[4082]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[4161]*Box2Abb[5956]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8])))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[4176]*Box2Abb[86]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)*Box2Abb[51]) + (m_p1*epsP)*((Box2Abb[4291]*Box2Abb[5949]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[4398]*Box2Abb[5949]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51])) + (Box2Abb[4460]*Box2Abb[5949]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[4646]*Box2Abb[5949]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[4725]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[322]*Box2Abb[4866]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*pow(Box2Abb[8],2.))))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[4903]*Box2Abb[5950]*Box2Abb[8]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],3.)) + (m_p1*epsP)*((Box2Abb[5055]*Box2Abb[5952]*Box2Abb[8]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[340]*Box2Abb[5194]*Box2Abb[5952]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))) + (Box2Abb[5327]*Box2Abb[5952]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[5528]*Box2Abb[5952]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[5679]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[5880]*Box2Abb[86]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[1365]*Box2Abb[3275]*Box2Abb[5950]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[0]*Box2Abb[3288]*Box2Abb[5951]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[3294]*Box2Abb[340]*Box2Abb[5951]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[0]*Box2Abb[148]*Box2Abb[3288]*Box2Abb[5951]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[148]*Box2Abb[3294]*Box2Abb[340]*Box2Abb[5951]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[0]*Box2Abb[3310]*Box2Abb[86]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[3333]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[3275]*Box2Abb[5950]*Box2Abb[8]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[3288]*Box2Abb[5951]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[3294]*Box2Abb[340]*Box2Abb[5951]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[3288]*Box2Abb[5951]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[3294]*Box2Abb[340]*Box2Abb[5951]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[3310]*Box2Abb[8]*Box2Abb[86]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[3333]*Box2Abb[8]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[3343]*Box2Abb[5950]*Box2Abb[8]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[3357]*Box2Abb[5951]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[3363]*Box2Abb[5952]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[3357]*Box2Abb[5951]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[3363]*Box2Abb[5952]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[3385]*Box2Abb[8]*Box2Abb[86]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[3419]*Box2Abb[86]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[3343]*Box2Abb[5953]*pow(Box2Abb[8],2.)*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[3357]*Box2Abb[5954]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[3363]*Box2Abb[5955]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[3357]*Box2Abb[5954]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[3363]*Box2Abb[5955]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[320]*Box2Abb[3385]*pow(Box2Abb[8],2.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[320]*Box2Abb[3419]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],3.)));
    }
    else if (LR == 1) {
      return One
	*((m_p1*epsP)*((Box2Abb[5945]*Box2Abb[6021]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5945]*Box2Abb[6025]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[5945]*Box2Abb[6027]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5945]*Box2Abb[6042]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[319]*Box2Abb[6043]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5947]*Box2Abb[6044]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[5958]*Box2Abb[86]*(epsP*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (m_p1*epsP)*((Box2Abb[5949]*Box2Abb[5964]*(m_p1*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (Box2Abb[5949]*Box2Abb[5973]*(m_p2*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[405])) + (Box2Abb[5949]*Box2Abb[5975]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (Box2Abb[5949]*Box2Abb[5982]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[405]) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[5999]*(m_p1*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[405]) + (Box2Abb[322]*Box2Abb[6014]*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[405])))

	+B_0(0.,0.,m_m2,m_mu2)*((m_p1*epsP)*((Box2Abb[6326]*Box2Abb[707]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5947]*Box2Abb[6015]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[5947]*Box2Abb[6015]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5948]*Box2Abb[6016]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[5947]*Box2Abb[6015]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[5948]*Box2Abb[6016]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[4]*Box2Abb[6077]*Box2Abb[85]*(epsP*epsV))/(Box2Abb[0]*Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[322]*Box2Abb[4]*Box2Abb[6078]*(m_p1*epsV))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[319]*Box2Abb[4]*Box2Abb[6103]*(m_p2*epsV))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[6104]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[6105]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (m_pP*epsV)*((Box2Abb[319]*Box2Abb[6140]*(m_p1*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[319]*Box2Abb[6151]*(m_p2*epsP))/(Box2Abb[0]*pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[6152]*Box2Abb[86]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[321]*Box2Abb[6153]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[322]*Box2Abb[6154]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.))) + (Box2Abb[322]*Box2Abb[6155]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[5956]*Box2Abb[6156]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (m_pP*epsV)*((Box2Abb[322]*Box2Abb[6157]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[5956]*Box2Abb[6158]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8])))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[6168]*Box2Abb[85]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)*Box2Abb[51]) + (m_p1*epsP)*((Box2Abb[5945]*Box2Abb[6188]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[5945]*Box2Abb[6202]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51])) + (Box2Abb[5945]*Box2Abb[6211]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[5945]*Box2Abb[6218]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[319]*Box2Abb[6250]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]) + (Box2Abb[319]*Box2Abb[6252]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*pow(Box2Abb[8],2.))))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[5950]*Box2Abb[6253]*Box2Abb[8]*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],3.)) + (m_p1*epsP)*((Box2Abb[5952]*Box2Abb[6254]*Box2Abb[8]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[5952]*Box2Abb[6258]*Box2Abb[8]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))) + (Box2Abb[5952]*Box2Abb[6259]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[5952]*Box2Abb[6260]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[6261]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[6262]*Box2Abb[86]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[6064]*Box2Abb[8]*Box2Abb[87]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[5951]*Box2Abb[6065]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[5952]*Box2Abb[6066]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[5951]*Box2Abb[6065]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[5952]*Box2Abb[6066]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[6067]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[6075]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[5950]*Box2Abb[6045]*Box2Abb[8]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[5951]*Box2Abb[6046]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[5951]*Box2Abb[6049]*pow(Box2Abb[8],3.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[5951]*Box2Abb[6046]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[5951]*Box2Abb[6049]*pow(Box2Abb[8],2.)*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[6058]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[6059]*Box2Abb[8]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[1365]*Box2Abb[5950]*Box2Abb[6045]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[0]*Box2Abb[5951]*Box2Abb[6046]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[5951]*Box2Abb[6049]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[0]*Box2Abb[148]*Box2Abb[5951]*Box2Abb[6046]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[148]*Box2Abb[5951]*Box2Abb[6049]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[0]*Box2Abb[6058]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[0]*Box2Abb[6059]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[6064]*pow(Box2Abb[8],2.)*Box2Abb[88]*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[5954]*Box2Abb[6065]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[5955]*Box2Abb[6066]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[148]*Box2Abb[5954]*Box2Abb[6065]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[148]*Box2Abb[5955]*Box2Abb[6066]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[598]*Box2Abb[6067]*pow(Box2Abb[8],2.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[598]*Box2Abb[6075]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],3.)));
    }
  }
  // ubar1 \slashed{k} P_i v2
  if (ME == 2) {
    if (LR == 0) {
      return One
	*((m_p1*epsP)*((Box2Abb[1341]*Box2Abb[2983]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1342]*Box2Abb[2985]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3])) + (Box2Abb[1342]*Box2Abb[2985]*Box2Abb[4]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1352]*Box2Abb[2985]*Box2Abb[4]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[1353]*Box2Abb[2983]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1363]*Box2Abb[2983]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[2981]*Box2Abb[382]*(epsP*epsV))/(Box2Abb[3]*Box2Abb[4]) + (m_p1*epsP)*((Box2Abb[1300]*Box2Abb[2982]*(m_p1*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1308]*Box2Abb[2983]*(m_p2*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[51])) + (Box2Abb[1317]*Box2Abb[2982]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[51]) + (Box2Abb[1323]*Box2Abb[2982]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]) + (m_pP*epsV)*((Box2Abb[1330]*Box2Abb[2984]*(m_p1*epsP))/(pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.)) + (Box2Abb[1337]*Box2Abb[2984]*(m_p2*epsP))/(pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.))))

	+B_0(0.,0.,m_m2,m_mu2)*((m_p1*epsP)*((Box2Abb[2985]*Box2Abb[707]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[2985]*Box2Abb[429]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[2985]*Box2Abb[429]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1340]*Box2Abb[2986]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[2985]*Box2Abb[429]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1340]*Box2Abb[2986]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1251]*Box2Abb[1526]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[4]) + (m_p1*epsP)*((Box2Abb[1578]*Box2Abb[2983]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1626]*Box2Abb[2983]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4])) + (Box2Abb[1677]*Box2Abb[2983]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1721]*Box2Abb[2983]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (m_pP*epsV)*((Box2Abb[1773]*Box2Abb[2983]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.)) + (Box2Abb[1831]*Box2Abb[2983]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1246]*Box2Abb[1833]*Box2Abb[218]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (m_p1*epsP)*((Box2Abb[1856]*Box2Abb[3004]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[1896]*Box2Abb[3005]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51])) + (Box2Abb[1974]*Box2Abb[3005]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[2037]*Box2Abb[3005]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[2059]*Box2Abb[3004]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (Box2Abb[2115]*Box2Abb[3005]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.))))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1245]*Box2Abb[2118]*Box2Abb[340]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_p1*epsP)*((Box2Abb[2143]*Box2Abb[3006]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[2174]*Box2Abb[2982]*Box2Abb[340]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8])) + (Box2Abb[2224]*Box2Abb[2982]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (Box2Abb[2279]*Box2Abb[3007]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (m_pP*epsV)*((Box2Abb[2315]*Box2Abb[2982]*Box2Abb[340]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.)) + (Box2Abb[2369]*Box2Abb[3007]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.))))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1253]*Box2Abb[2373]*Box2Abb[429]*pow(Box2Abb[8],2.)*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[2458]*Box2Abb[3008]*Box2Abb[8]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[1245]*Box2Abb[2545]*Box2Abb[8]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))) + (Box2Abb[1245]*Box2Abb[2677]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[2781]*Box2Abb[3009]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[2818]*Box2Abb[3008]*pow(Box2Abb[8],2.)*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[2909]*Box2Abb[3009]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[2987]*pow(Box2Abb[405],2.)*pow(Box2Abb[8],2.)*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[1371]*Box2Abb[2988]*Box2Abb[51]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1245]*Box2Abb[1388]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1245]*Box2Abb[1409]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1412]*Box2Abb[148]*Box2Abb[2989]*Box2Abb[51]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[2990]*Box2Abb[382]*Box2Abb[405]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1414]*Box2Abb[2991]*Box2Abb[372]*Box2Abb[4]*Box2Abb[51]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[1249]*Box2Abb[1422]*Box2Abb[4]*(epsP*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[1432]*Box2Abb[2998]*pow(Box2Abb[4],2.)*(m_p1*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1448]*Box2Abb[2999]*pow(Box2Abb[4],2.)*(m_p2*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],3.))) + (Box2Abb[1464]*Box2Abb[2999]*pow(Box2Abb[4],2.)*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1469]*Box2Abb[148]*Box2Abb[3000]*pow(Box2Abb[4],2.)*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[1478]*Box2Abb[3001]*(m_p1*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1491]*Box2Abb[3000]*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.))))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1253]*Box2Abb[1496]*Box2Abb[382]*(epsP*epsV))/pow(Box2Abb[3],2.) + (m_p1*epsP)*((Box2Abb[1432]*Box2Abb[3002]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1251]*Box2Abb[1448]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1251]*Box2Abb[1464]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1469]*Box2Abb[148]*Box2Abb[3003]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[1503]*Box2Abb[3002]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1521]*Box2Abb[3003]*(m_p2*epsP))/(pow(Box2Abb[3],3.)*Box2Abb[8])))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[2992]*pow(Box2Abb[405],2.)*pow(Box2Abb[8],3.)*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[1371]*Box2Abb[2993]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1388]*Box2Abb[2994]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1409]*Box2Abb[2994]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1412]*Box2Abb[148]*Box2Abb[2995]*Box2Abb[51]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[2996]*Box2Abb[382]*Box2Abb[405]*Box2Abb[51]*pow(Box2Abb[8],3.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1414]*Box2Abb[2997]*Box2Abb[372]*Box2Abb[4]*Box2Abb[51]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],3.)));
    }
    else if (LR == 1) {
      return One
	*((m_p1*epsP)*((Box2Abb[1341]*Box2Abb[3083]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1342]*Box2Abb[3085]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3])) + (Box2Abb[1342]*Box2Abb[3085]*Box2Abb[4]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1352]*Box2Abb[3085]*Box2Abb[4]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[1353]*Box2Abb[3083]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1363]*Box2Abb[3083]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[3081]*Box2Abb[382]*(epsP*epsV))/(Box2Abb[3]*Box2Abb[4]) + (m_p1*epsP)*((Box2Abb[1300]*Box2Abb[3082]*(m_p1*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1308]*Box2Abb[3083]*(m_p2*epsV))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[51])) + (Box2Abb[1317]*Box2Abb[3082]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]*Box2Abb[51]) + (Box2Abb[1323]*Box2Abb[3082]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[3],2.)*Box2Abb[4]) + (m_pP*epsV)*((Box2Abb[1330]*Box2Abb[3084]*(m_p1*epsP))/(pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.)) + (Box2Abb[1337]*Box2Abb[3084]*(m_p2*epsP))/(pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.))))

	+B_0(0.,0.,m_m2,m_mu2)*((m_p1*epsP)*((Box2Abb[3085]*Box2Abb[707]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[3085]*Box2Abb[429]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3])) + (Box2Abb[3085]*Box2Abb[429]*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1340]*Box2Abb[3086]*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[3085]*Box2Abb[429]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1340]*Box2Abb[3086]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.))))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1293]*Box2Abb[1526]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[4]) + (m_p1*epsP)*((Box2Abb[1578]*Box2Abb[3083]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1626]*Box2Abb[3083]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4])) + (Box2Abb[1677]*Box2Abb[3083]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (Box2Abb[1721]*Box2Abb[3083]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[4]) + (m_pP*epsV)*((Box2Abb[1773]*Box2Abb[3083]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.)) + (Box2Abb[1831]*Box2Abb[3083]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[4],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1288]*Box2Abb[1833]*Box2Abb[218]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (m_p1*epsP)*((Box2Abb[1856]*Box2Abb[3104]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[1896]*Box2Abb[3105]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51])) + (Box2Abb[1974]*Box2Abb[3105]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[2037]*Box2Abb[3105]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.)) + (m_pP*epsV)*((Box2Abb[2059]*Box2Abb[3104]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (Box2Abb[2115]*Box2Abb[3105]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.))))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1287]*Box2Abb[2118]*Box2Abb[340]*(epsP*epsV))/(Box2Abb[142]*Box2Abb[3]*pow(Box2Abb[8],2.)) + (m_p1*epsP)*((Box2Abb[2143]*Box2Abb[3106]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)) + (Box2Abb[2174]*Box2Abb[3082]*Box2Abb[340]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8])) + (Box2Abb[2224]*Box2Abb[3082]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (Box2Abb[2279]*Box2Abb[3107]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*Box2Abb[8]) + (m_pP*epsV)*((Box2Abb[2315]*Box2Abb[3082]*Box2Abb[340]*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.)) + (Box2Abb[2369]*Box2Abb[3107]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],2.)*pow(Box2Abb[8],2.))))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1295]*Box2Abb[2373]*Box2Abb[429]*pow(Box2Abb[8],2.)*(epsP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[2458]*Box2Abb[3108]*Box2Abb[8]*(m_p1*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[1287]*Box2Abb[2545]*Box2Abb[8]*(m_p2*epsV))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))) + (Box2Abb[1287]*Box2Abb[2677]*(m_p1*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[2781]*Box2Abb[3109]*(m_p2*epsV)*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[2818]*Box2Abb[3108]*pow(Box2Abb[8],2.)*(m_p1*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.)) + (Box2Abb[2909]*Box2Abb[3109]*(m_p2*epsP))/(pow(Box2Abb[142],2.)*pow(Box2Abb[3],3.))))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1365]*Box2Abb[3087]*pow(Box2Abb[405],2.)*pow(Box2Abb[8],2.)*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[1371]*Box2Abb[3088]*Box2Abb[51]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1287]*Box2Abb[1388]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1287]*Box2Abb[1409]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1412]*Box2Abb[148]*Box2Abb[3089]*Box2Abb[51]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[3090]*Box2Abb[382]*Box2Abb[405]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1414]*Box2Abb[3091]*Box2Abb[372]*Box2Abb[4]*Box2Abb[51]*(m_p2*epsP))/pow(Box2Abb[3],3.)))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[1291]*Box2Abb[1422]*Box2Abb[4]*(epsP*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],2.)) + (m_p1*epsP)*((Box2Abb[1432]*Box2Abb[3098]*pow(Box2Abb[4],2.)*(m_p1*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1448]*Box2Abb[3099]*pow(Box2Abb[4],2.)*(m_p2*epsV))/(Box2Abb[1415]*pow(Box2Abb[3],3.))) + (Box2Abb[1464]*Box2Abb[3099]*pow(Box2Abb[4],2.)*(m_p1*epsV)*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1469]*Box2Abb[148]*Box2Abb[3100]*pow(Box2Abb[4],2.)*(m_p2*epsV)*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (m_pP*epsV)*((Box2Abb[1478]*Box2Abb[3101]*(m_p1*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.)) + (Box2Abb[1491]*Box2Abb[3100]*(m_p2*epsP))/(Box2Abb[1415]*pow(Box2Abb[3],3.))))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1295]*Box2Abb[1496]*Box2Abb[382]*(epsP*epsV))/pow(Box2Abb[3],2.) + (m_p1*epsP)*((Box2Abb[1432]*Box2Abb[3102]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1293]*Box2Abb[1448]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1293]*Box2Abb[1464]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1469]*Box2Abb[148]*Box2Abb[3103]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[1503]*Box2Abb[3102]*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1521]*Box2Abb[3103]*(m_p2*epsP))/(pow(Box2Abb[3],3.)*Box2Abb[8])))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1365]*Box2Abb[3092]*pow(Box2Abb[405],2.)*pow(Box2Abb[8],3.)*(epsP*epsV))/pow(Box2Abb[3],3.) + (m_p1*epsP)*((Box2Abb[1371]*Box2Abb[3093]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],3.) + (Box2Abb[1388]*Box2Abb[3094]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],3.)) + (Box2Abb[1409]*Box2Abb[3094]*Box2Abb[8]*(m_p1*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1412]*Box2Abb[148]*Box2Abb[3095]*Box2Abb[51]*Box2Abb[8]*(m_p2*epsV)*(m_p2*epsP))/pow(Box2Abb[3],3.) + (m_pP*epsV)*((Box2Abb[3096]*Box2Abb[382]*Box2Abb[405]*Box2Abb[51]*pow(Box2Abb[8],3.)*(m_p1*epsP))/pow(Box2Abb[3],3.) + (Box2Abb[1414]*Box2Abb[3097]*Box2Abb[372]*Box2Abb[4]*Box2Abb[51]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],3.)));
    }
  }
  // ubar1 \slashed{epsP} P_i v2
  if (ME == 3) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[1009]*Box2Abb[1244]*Box2Abb[405]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[1012]*Box2Abb[1245]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[1014]*Box2Abb[1246]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1056]*Box2Abb[1251]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1063]*Box2Abb[1251]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1067]*Box2Abb[1251]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1078]*Box2Abb[1079]*Box2Abb[1245]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1085]*Box2Abb[1252]*Box2Abb[340]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1089]*Box2Abb[1245]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1096]*Box2Abb[1246]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1106]*Box2Abb[1245]*Box2Abb[340]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[1115]*Box2Abb[1246]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1016]*Box2Abb[1247]*Box2Abb[51]*Box2Abb[8]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1022]*Box2Abb[1248]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1027]*Box2Abb[1248]*Box2Abb[148]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[1033]*Box2Abb[1248]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1041]*Box2Abb[1248]*Box2Abb[4]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1044]*Box2Abb[1248]*Box2Abb[148]*Box2Abb[4]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1077]*Box2Abb[1248]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1041]*Box2Abb[1248]*Box2Abb[340]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1044]*Box2Abb[1248]*Box2Abb[148]*Box2Abb[340]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1138]*Box2Abb[1253]*Box2Abb[340]*(m_pP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[1182]*Box2Abb[1253]*(m_p1*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[1211]*Box2Abb[1253]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1016]*Box2Abb[1249]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1022]*Box2Abb[1250]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1027]*Box2Abb[1250]*Box2Abb[148]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[1009]*Box2Abb[1286]*Box2Abb[405]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[1012]*Box2Abb[1287]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[1014]*Box2Abb[1288]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[1056]*Box2Abb[1293]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1063]*Box2Abb[1293]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[1067]*Box2Abb[1293]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[1078]*Box2Abb[1079]*Box2Abb[1287]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1085]*Box2Abb[1294]*Box2Abb[340]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1089]*Box2Abb[1287]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[1096]*Box2Abb[1288]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[1106]*Box2Abb[1287]*Box2Abb[340]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[1115]*Box2Abb[1288]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[1016]*Box2Abb[1289]*Box2Abb[51]*Box2Abb[8]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1022]*Box2Abb[1290]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1027]*Box2Abb[1290]*Box2Abb[148]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[1033]*Box2Abb[1290]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1041]*Box2Abb[1290]*Box2Abb[4]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1044]*Box2Abb[1290]*Box2Abb[148]*Box2Abb[4]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1077]*Box2Abb[1290]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1041]*Box2Abb[1290]*Box2Abb[340]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1044]*Box2Abb[1290]*Box2Abb[148]*Box2Abb[340]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1138]*Box2Abb[1295]*Box2Abb[340]*(m_pP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[1182]*Box2Abb[1295]*(m_p1*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[1211]*Box2Abb[1295]*Box2Abb[340]*(m_p2*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[1016]*Box2Abb[1291]*Box2Abb[51]*pow(Box2Abb[8],2.)*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1022]*Box2Abb[1292]*pow(Box2Abb[8],2.)*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[1027]*Box2Abb[1292]*Box2Abb[148]*pow(Box2Abb[8],2.)*(m_p2*epsV))/pow(Box2Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsV} P_i v2
  if (ME == 4) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[698]*Box2Abb[938]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[706]*Box2Abb[938]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[713]*Box2Abb[939]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[722]*Box2Abb[939]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[752]*Box2Abb[940]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[767]*Box2Abb[939]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[782]*Box2Abb[941]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[805]*Box2Abb[938]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((-Box2Abb[340]*Box2Abb[4]*Box2Abb[673]*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[4]*Box2Abb[673]*Box2Abb[674]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((-Box2Abb[673]*pow(Box2Abb[8],2.)*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[148]*Box2Abb[673]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((-Box2Abb[729]*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[748]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((-Box2Abb[849]*(m_p1*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)) - (Box2Abb[914]*(m_p2*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[37]*Box2Abb[691]*Box2Abb[8]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[148]*Box2Abb[37]*Box2Abb[691]*(m_p2*epsP))/pow(Box2Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[938]*Box2Abb[944]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[938]*Box2Abb[945]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[939]*Box2Abb[946]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[939]*Box2Abb[947]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[940]*Box2Abb[950]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[939]*Box2Abb[951]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[941]*Box2Abb[952]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[939]*Box2Abb[957]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((-Box2Abb[340]*Box2Abb[4]*Box2Abb[942]*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[4]*Box2Abb[674]*Box2Abb[942]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((-pow(Box2Abb[8],2.)*Box2Abb[942]*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[148]*Box2Abb[8]*Box2Abb[942]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((-Box2Abb[948]*(m_p1*epsP))/pow(Box2Abb[3],2.) - (Box2Abb[949]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((-Box2Abb[958]*(m_p1*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)) - (Box2Abb[983]*(m_p2*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[37]*Box2Abb[8]*Box2Abb[943]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[148]*Box2Abb[37]*Box2Abb[943]*(m_p2*epsP))/pow(Box2Abb[3],2.));
    }
  }
  // ubar1 \slashed{k} \slashed{epsP} P_i v2
  if (ME == 5) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[369]*Box2Abb[595]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[371]*Box2Abb[4]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[322]*Box2Abb[373]*Box2Abb[4]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[322]*Box2Abb[396]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[402]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[404]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[319]*Box2Abb[432]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[340]*Box2Abb[436]*Box2Abb[599]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[319]*Box2Abb[340]*Box2Abb[439]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[443]*Box2Abb[600]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[322]*Box2Abb[453]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[322]*Box2Abb[467]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[381]*Box2Abb[86]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[387]*Box2Abb[86]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[369]*Box2Abb[388]*Box2Abb[8]*Box2Abb[86]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[381]*Box2Abb[8]*Box2Abb[86]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[387]*Box2Abb[8]*Box2Abb[86]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[369]*Box2Abb[388]*pow(Box2Abb[8],2.)*Box2Abb[86]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[405]*Box2Abb[406]*Box2Abb[596]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[415]*Box2Abb[85]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[428]*Box2Abb[85]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[498]*Box2Abb[86]*(m_pP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[530]*Box2Abb[86]*(m_p1*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[565]*Box2Abb[86]*(m_p2*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[405]*Box2Abb[406]*Box2Abb[597]*Box2Abb[8]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[415]*Box2Abb[598]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[428]*Box2Abb[598]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[595]*Box2Abb[601]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[602]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[603]*(m_p2*epsV))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[322]*Box2Abb[616]*(m_pP*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[4]*Box2Abb[617]*(m_p1*epsV))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[618]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[319]*Box2Abb[622]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[340]*Box2Abb[599]*Box2Abb[623]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[319]*Box2Abb[340]*Box2Abb[624]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[600]*Box2Abb[625]*(m_pP*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]) + (Box2Abb[322]*Box2Abb[340]*Box2Abb[628]*(m_p1*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[322]*Box2Abb[340]*Box2Abb[631]*(m_p2*epsV))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[612]*Box2Abb[85]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[615]*Box2Abb[85]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[340]*Box2Abb[388]*Box2Abb[601]*Box2Abb[85]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[612]*Box2Abb[8]*Box2Abb[85]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[615]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[340]*Box2Abb[388]*Box2Abb[601]*Box2Abb[8]*Box2Abb[85]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[405]*Box2Abb[619]*Box2Abb[664]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[620]*Box2Abb[85]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[621]*Box2Abb[85]*(m_p2*epsV))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[632]*Box2Abb[86]*(m_pP*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[633]*Box2Abb[86]*(m_p1*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[634]*Box2Abb[86]*(m_p2*epsV))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[405]*Box2Abb[619]*Box2Abb[665]*Box2Abb[8]*(m_pP*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[598]*Box2Abb[620]*Box2Abb[8]*(m_p1*epsV))/pow(Box2Abb[3],2.) + (Box2Abb[598]*Box2Abb[621]*Box2Abb[8]*(m_p2*epsV))/pow(Box2Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsV} \slashed{k} P_i v2
  if (ME == 6) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[135]*Box2Abb[319]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[141]*Box2Abb[319]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[147]*Box2Abb[319]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[153]*Box2Abb[319]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[209]*Box2Abb[321]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[217]*Box2Abb[322]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[225]*Box2Abb[322]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[239]*Box2Abb[322]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[165]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[173]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[165]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[173]*Box2Abb[8]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[191]*Box2Abb[86]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[207]*Box2Abb[86]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[265]*Box2Abb[86]*(m_p1*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[299]*Box2Abb[86]*(m_p2*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[191]*Box2Abb[320]*Box2Abb[8]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[207]*Box2Abb[320]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box2Abb[319]*Box2Abb[325]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]) + (Box2Abb[319]*Box2Abb[327]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[3]*Box2Abb[51]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box2Abb[322]*Box2Abb[331]*(m_p1*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[333]*(m_p2*epsP))/(Box2Abb[0]*Box2Abb[142]*Box2Abb[3]))

	+B_0(m_s2k,0.,m_m2,m_mu2)*((Box2Abb[321]*Box2Abb[338]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]) + (Box2Abb[322]*Box2Abb[339]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[8]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box2Abb[322]*Box2Abb[340]*Box2Abb[345]*(m_p1*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]) + (Box2Abb[322]*Box2Abb[346]*(m_p2*epsP))/(Box2Abb[142]*Box2Abb[3]*Box2Abb[51]*Box2Abb[8]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box2Abb[0]*Box2Abb[334]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[0]*Box2Abb[335]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[334]*Box2Abb[8]*Box2Abb[85]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[335]*Box2Abb[8]*Box2Abb[85]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box2Abb[336]*Box2Abb[86]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[337]*Box2Abb[86]*(m_p2*epsP))/pow(Box2Abb[3],2.))

	+C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)*((Box2Abb[347]*Box2Abb[86]*(m_p1*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)) + (Box2Abb[348]*Box2Abb[86]*(m_p2*epsP))/(Box2Abb[142]*pow(Box2Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box2Abb[320]*Box2Abb[336]*Box2Abb[8]*(m_p1*epsP))/pow(Box2Abb[3],2.) + (Box2Abb[320]*Box2Abb[337]*Box2Abb[8]*(m_p2*epsP))/pow(Box2Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsV} \slashed{epsP} P_i v2
  if (ME == 7) {
    if (LR == 0) {
      return (Box2Abb[50]*Box2Abb[85]*B_0(m_m2,0.,m_m2,m_mu2))/Box2Abb[51]

	+(Box2Abb[50]*Box2Abb[86]*B_0(m_s12,m_m2,m_m2,m_mu2))/Box2Abb[51]

	+(Box2Abb[0]*Box2Abb[58]*Box2Abb[87]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[58]*Box2Abb[8]*Box2Abb[87]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[69]*Box2Abb[87]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[76]*Box2Abb[87]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[76]*Box2Abb[8]*Box2Abb[88]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[3];
    }
    else if (LR == 1) {
      return (Box2Abb[50]*Box2Abb[85]*B_0(m_m2,0.,m_m2,m_mu2))/Box2Abb[51]

	+(Box2Abb[50]*Box2Abb[86]*B_0(m_s12,m_m2,m_m2,m_mu2))/Box2Abb[51]

	+(Box2Abb[0]*Box2Abb[87]*Box2Abb[89]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[8]*Box2Abb[87]*Box2Abb[89]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[87]*Box2Abb[90]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[87]*Box2Abb[91]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[8]*Box2Abb[88]*Box2Abb[91]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[3];
    }
  }
  // ubar1  \slashed{epsV}  \slashed{k} \slashed{epsP} P_i v2
  if (ME == 8) {
    if (LR == 0) {
      return (Box2Abb[0]*Box2Abb[7]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[7]*Box2Abb[8]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(-Box2Abb[16]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[20]*Box2Abb[36]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[25]*Box2Abb[37]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[3];
    }
    else if (LR == 1) {
      return (Box2Abb[0]*Box2Abb[39]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[39]*Box2Abb[8]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(-Box2Abb[40]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[3]

	+(Box2Abb[36]*Box2Abb[41]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[3]

	+(Box2Abb[37]*Box2Abb[42]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[3];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
		<< "Values range from 1 to 8.";
  }
  return Zero;
}
