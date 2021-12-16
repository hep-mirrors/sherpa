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

Z_Decay_RV_Box_1::Z_Decay_RV_Box_1
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Z_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Box1Abb = new Complex[6241];
  Init_Coefficients();
}

Z_Decay_RV_Box_1::~Z_Decay_RV_Box_1()
{
  delete [] Box1Abb;
}

void Z_Decay_RV_Box_1::Init_Coefficients() 
{
  Init_Box_1_Coefficients();
  return;
}

// Set up coefficients multiplying master integrals for emission from leg 2
// Calculated using FeynCalc 9.0.1
void Z_Decay_RV_Box_1::Init_Box_1_Coefficients() 
{
  Box1Abb[0]=-1. + m_z12;

  Box1Abb[1]=-1. + m_x;

  Box1Abb[2]=-1. - 2.*m_x + m_z12;

  Box1Abb[3]=m_x + Box1Abb[1]*m_x*m_z12 + Box1Abb[2]*m_z12*m_z1k + m_z12*m_z1k_2;

  Box1Abb[4]=-1. + m_z12 + m_z1k;

  Box1Abb[5]=m_x + Box1Abb[4]*m_z12 - 2.*m_x*m_z12;

  Box1Abb[6]=Box1Abb[5]*m_cL - Box1Abb[0]*m_cR*m_x;

  Box1Abb[7]=m_x - m_z1k;

  Box1Abb[8]=-1. + 2.*m_z12;

  Box1Abb[9]=m_z12 + m_z1k;

  Box1Abb[10]=Box1Abb[8]*m_x + m_z12 - Box1Abb[9]*m_z12;

  Box1Abb[11]=Box1Abb[10]*m_cL + Box1Abb[0]*m_cR*m_x;

  Box1Abb[12]=1. + m_z12;

  Box1Abb[13]=1. + m_x - m_z12 + m_x*m_z12 - Box1Abb[12]*m_z1k;

  Box1Abb[14]=1. + 2.*m_z12;

  Box1Abb[15]=1. + m_z1k;

  Box1Abb[16]=3. + m_z12 + 3.*m_z1k;

  Box1Abb[17]=-3. + Box1Abb[16]*m_z12 + m_z1k;

  Box1Abb[18]=-Box1Abb[17]*m_x + Box1Abb[14]*m_x_2 + Box1Abb[15]*Box1Abb[4]*m_z12;

  Box1Abb[19]=Box1Abb[18]*m_cL + Box1Abb[13]*m_cR*m_x;

  Box1Abb[20]=1. + 2.*m_x - m_z12 - 2.*m_z1k;

  Box1Abb[21]=3. - 4.*m_z12 - 2.*m_z1k;

  Box1Abb[22]=Box1Abb[21]*m_x + 2.*m_x_2 + Box1Abb[4]*m_z12;

  Box1Abb[23]=Box1Abb[22]*m_cL + Box1Abb[20]*m_cR*m_x;

  Box1Abb[24]=2. - 3.*Box1Abb[15]*m_z12 + m_z12_2 + m_z1k - 2.*m_z1k_2;

  Box1Abb[25]=7. - 2.*m_z12 + 4.*m_z1k;

  Box1Abb[26]=-4. + Box1Abb[25]*m_z12;

  Box1Abb[27]=Box1Abb[26]*m_x_2 + Box1Abb[24]*m_x*m_z12 - 2.*m_x_3*m_z12 + Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[28]=Box1Abb[27]*m_cL + Box1Abb[20]*Box1Abb[7]*m_cR*m_x*m_z12;

  Box1Abb[37]=m_s*m_z12;

  Box1Abb[38]=m_s*m_x;

  Box1Abb[39]=m_s*m_z1k;

  Box1Abb[40]=-m_z12;

  Box1Abb[41]=-m_s;

  Box1Abb[42]=m_x - m_x*m_z12;

  Box1Abb[43]=Box1Abb[42]*m_cL + Box1Abb[5]*m_cR;

  Box1Abb[44]=Box1Abb[10]*m_cR + Box1Abb[0]*m_cL*m_x;

  Box1Abb[45]=Box1Abb[18]*m_cR + Box1Abb[13]*m_cL*m_x;

  Box1Abb[46]=Box1Abb[22]*m_cR + Box1Abb[20]*m_cL*m_x;

  Box1Abb[47]=Box1Abb[27]*m_cR + Box1Abb[20]*Box1Abb[7]*m_cL*m_x*m_z12;

  Box1Abb[56]=m_cL + m_cR;

  Box1Abb[57]=4.*m_x - m_z12;

  Box1Abb[58]=1. - m_z12;

  Box1Abb[59]=1. - 2.*m_z12 - 4.*m_z1k;

  Box1Abb[60]=-3. + 4.*m_z12;

  Box1Abb[61]=pow(Box1Abb[0],2.) + Box1Abb[59]*m_x + 2.*m_x_2 + Box1Abb[60]*m_z1k + 2.*m_z1k_2;

  Box1Abb[62]=-2. + m_z12 - 4.*m_z1k;

  Box1Abb[63]=m_z12 + 2.*m_z1k;

  Box1Abb[64]=Box1Abb[62]*m_x + 2.*m_x_2 + Box1Abb[63]*m_z1k;

  Box1Abb[65]=Box1Abb[61]*m_cL - Box1Abb[64]*m_cR;

  Box1Abb[66]=1. + 6.*m_z1k;

  Box1Abb[67]=-2. + m_z12_2 + 4.*m_z1k - 2.*m_z12*m_z1k - 6.*m_z1k_2;

  Box1Abb[68]=-1. + m_z1k;

  Box1Abb[69]=-3. + 2.*m_z12;

  Box1Abb[70]=pow(Box1Abb[0],2.) + Box1Abb[69]*m_z1k + 2.*m_z1k_2;

  Box1Abb[71]=Box1Abb[68]*Box1Abb[70] + Box1Abb[67]*m_x + Box1Abb[66]*m_x_2 - 2.*m_x_3;

  Box1Abb[72]=-1. + 2.*m_z12 + m_z1k;

  Box1Abb[73]=4. + 3.*m_z12 + 6.*m_z1k;

  Box1Abb[74]=-2. + 3.*m_z12 + 2.*m_z1k + 8.*m_z12*m_z1k + 6.*m_z1k_2;

  Box1Abb[75]=Box1Abb[74]*m_x - Box1Abb[73]*m_x_2 + 2.*m_x_3 - Box1Abb[63]*Box1Abb[72]*m_z1k;

  Box1Abb[76]=Box1Abb[71]*m_cL + Box1Abb[75]*m_cR;

  Box1Abb[77]=2. + m_z12;

  Box1Abb[78]=2.*Box1Abb[77]*m_x - Box1Abb[63]*m_z12;

  Box1Abb[79]=-2. + m_z12;

  Box1Abb[80]=7. - 4.*m_z12 - 4.*m_z1k;

  Box1Abb[81]=-2. + Box1Abb[80]*m_z12 + 4.*m_z1k;

  Box1Abb[82]=Box1Abb[81]*m_x + 2.*Box1Abb[79]*m_x_2 + Box1Abb[70]*m_z12;

  Box1Abb[83]=Box1Abb[82]*m_cL + Box1Abb[7]*Box1Abb[78]*m_cR;

  Box1Abb[92]=(-2.*m_m)/m_s;

  Box1Abb[93]=(2.*m_m)/m_s;

  Box1Abb[94]=m_m;

  Box1Abb[95]=m_m*m_s;

  Box1Abb[96]=-Box1Abb[64]*m_cL + Box1Abb[61]*m_cR;

  Box1Abb[97]=Box1Abb[75]*m_cL + Box1Abb[71]*m_cR;

  Box1Abb[98]=Box1Abb[7]*Box1Abb[78]*m_cL + Box1Abb[82]*m_cR;

  Box1Abb[107]=m_x + 2.*Box1Abb[12]*m_x_2 - 2.*Box1Abb[72]*m_x*m_z12 + Box1Abb[4]*m_z12_2 - 2.*m_x*m_z1k;

  Box1Abb[108]=3. + m_z12 + 2.*m_z1k;

  Box1Abb[109]=-3. + Box1Abb[108]*m_z12 + 2.*m_z1k;

  Box1Abb[110]=-Box1Abb[109]*m_x + 2.*Box1Abb[12]*m_x_2 + Box1Abb[4]*m_z12;

  Box1Abb[111]=Box1Abb[110]*m_cL - Box1Abb[107]*m_cR;

  Box1Abb[112]=1. - 2.*Box1Abb[12]*m_x + Box1Abb[0]*m_z12;

  Box1Abb[113]=2.*Box1Abb[12]*m_x - m_z12;

  Box1Abb[114]=m_x + 2.*Box1Abb[12]*m_x_2 - 2.*m_x*m_z12 - 2.*Box1Abb[12]*m_x*m_z1k + m_z12_2*m_z1k;

  Box1Abb[115]=Box1Abb[114]*m_cL + Box1Abb[112]*m_cR*m_x + Box1Abb[113]*m_cR*m_z1k;

  Box1Abb[116]=pow(Box1Abb[68],2.) - 2.*Box1Abb[15]*m_x + m_x_2;

  Box1Abb[117]=1. + m_x - m_z12 - m_z1k;

  Box1Abb[118]=m_z12 - m_z1k;

  Box1Abb[119]=2. + 3.*m_z1k;

  Box1Abb[120]=-1. + m_z12 + m_z12_2 + 3.*m_z1k_2;

  Box1Abb[121]=-Box1Abb[118]*Box1Abb[4]*Box1Abb[68] - Box1Abb[120]*m_x + Box1Abb[119]*m_x_2 - m_x_3;

  Box1Abb[122]=Box1Abb[116]*Box1Abb[117]*m_cL + Box1Abb[121]*m_cR;

  Box1Abb[123]=-2. + m_x + m_z12;

  Box1Abb[124]=3. - 2.*m_z12;

  Box1Abb[125]=2.*Box1Abb[0]*m_x + 3.*m_x_2 + Box1Abb[124]*m_z12;

  Box1Abb[126]=1. + 3.*m_x + m_z12;

  Box1Abb[127]=-Box1Abb[1]*Box1Abb[123]*m_x + Box1Abb[125]*m_z1k - Box1Abb[126]*m_z1k_2 + m_z1k_3;

  Box1Abb[128]=Box1Abb[116]*Box1Abb[7]*m_cL + Box1Abb[127]*m_cR;

  Box1Abb[129]=-4. + 3.*m_z12;

  Box1Abb[130]=6. - 4.*m_z12 - 3.*m_z1k;

  Box1Abb[131]=2.*Box1Abb[68] + Box1Abb[130]*m_z12;

  Box1Abb[132]=-Box1Abb[79]*m_x_3 + Box1Abb[131]*m_x*m_z1k + Box1Abb[129]*m_x_2*m_z1k + pow(Box1Abb[4],2.)*m_z12*m_z1k;

  Box1Abb[133]=-2. + 3.*m_z12;

  Box1Abb[134]=-5. + 3.*m_z12 - 3.*m_z1k;

  Box1Abb[135]=2. + Box1Abb[134]*m_z12 + 4.*m_z1k;

  Box1Abb[136]=Box1Abb[135]*m_x_2 + Box1Abb[79]*m_x_3 - pow(Box1Abb[0],2.)*m_x*m_z12 + Box1Abb[133]*m_x*m_z1k_2 - Box1Abb[4]*m_z12*m_z1k_2;

  Box1Abb[137]=Box1Abb[132]*m_cL + Box1Abb[136]*m_cR;

  Box1Abb[138]=m_z12_2 + 4.*m_z1k + 3.*m_z12*m_z1k;

  Box1Abb[139]=pow(Box1Abb[0],2.) + Box1Abb[133]*m_z1k + m_z1k_2;

  Box1Abb[140]=-1. + m_z12 + m_z12_2;

  Box1Abb[141]=2. + 3.*m_z12;

  Box1Abb[142]=pow(Box1Abb[0],2.) + 2.*Box1Abb[140]*m_z1k + Box1Abb[141]*m_z1k_2;

  Box1Abb[143]=Box1Abb[142]*m_x - Box1Abb[138]*m_x_2 + Box1Abb[77]*m_x_3 - Box1Abb[139]*m_z12*m_z1k;

  Box1Abb[144]=2.*m_z12 + 2.*m_z1k + 3.*m_z12*m_z1k;

  Box1Abb[145]=-2. + m_z12 + 3.*m_z1k;

  Box1Abb[146]=2. + Box1Abb[145]*m_z12 + 4.*m_z1k;

  Box1Abb[147]=Box1Abb[146]*m_x_2 - Box1Abb[77]*m_x_3 - Box1Abb[144]*m_x*m_z1k + Box1Abb[9]*m_z12*m_z1k_2;

  Box1Abb[148]=Box1Abb[143]*m_cL + Box1Abb[147]*m_cR;

  Box1Abb[149]=-5. + 3.*m_z12;

  Box1Abb[150]=4. + 2.*Box1Abb[149]*m_z12 - 6.*m_z1k + 9.*m_z12*m_z1k + 2.*m_z1k_2;

  Box1Abb[151]=-2. + m_z12 + 6.*m_z1k;

  Box1Abb[152]=2. + Box1Abb[151]*m_z12;

  Box1Abb[153]=-8. + 9.*m_z12 + 6.*m_z1k;

  Box1Abb[154]=2. + Box1Abb[153]*m_z12;

  Box1Abb[155]=Box1Abb[152]*m_x_3 - 2.*m_x_4*m_z12 - Box1Abb[154]*m_x_2*m_z1k + Box1Abb[150]*m_x*m_z12*m_z1k - pow(Box1Abb[4],2.)*m_z12_2*m_z1k;

  Box1Abb[156]=m_z12 - 2.*m_z1k;

  Box1Abb[157]=-2. + 3.*Box1Abb[156]*m_z12;

  Box1Abb[158]=-4. + m_z12;

  Box1Abb[159]=pow(Box1Abb[0],2.)*m_z12 + 4.*Box1Abb[0]*m_z1k - Box1Abb[158]*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[160]=4. + m_z12;

  Box1Abb[161]=6. - 9.*m_z12 + 5.*m_z12_2 + Box1Abb[160]*m_z1k - 6.*m_z1k_2;

  Box1Abb[162]=2.*Box1Abb[15] - Box1Abb[161]*m_z12;

  Box1Abb[163]=Box1Abb[162]*m_x_2 + Box1Abb[157]*m_x_3 + Box1Abb[159]*m_x*m_z12 + 2.*m_x_4*m_z12 - Box1Abb[4]*m_z12_2*m_z1k_2;

  Box1Abb[164]=Box1Abb[155]*m_cL + Box1Abb[163]*m_cR;

  Box1Abb[165]=-6. + 3.*m_z12 - 4.*m_z1k;

  Box1Abb[166]=-2. + Box1Abb[165]*m_z12;

  Box1Abb[167]=2. + 4.*m_z1k;

  Box1Abb[168]=Box1Abb[167]*m_z12 - m_z12_2 + 2.*Box1Abb[15]*m_z1k;

  Box1Abb[169]=Box1Abb[166]*m_x_2 + Box1Abb[168]*m_x*m_z12 + 2.*m_x_3*m_z12 - Box1Abb[9]*m_z12_2*m_z1k;

  Box1Abb[170]=-2. + m_z1k;

  Box1Abb[171]=pow(Box1Abb[68],2.) + Box1Abb[170]*m_z12 + m_z12_2;

  Box1Abb[172]=4. - 3.*m_z12 + 6.*m_z1k;

  Box1Abb[173]=2. + Box1Abb[172]*m_z12;

  Box1Abb[174]=1. + 4.*m_z1k;

  Box1Abb[175]=6. + m_z1k;

  Box1Abb[176]=2. + Box1Abb[175]*m_z1k;

  Box1Abb[177]=-1. + Box1Abb[176]*m_z12 - Box1Abb[174]*m_z12_2 + 2.*m_z1k_3;

  Box1Abb[178]=2. + m_z1k;

  Box1Abb[179]=-6. + Box1Abb[178]*m_z12 + m_z12_2 - 2.*Box1Abb[119]*m_z1k;

  Box1Abb[180]=2. + Box1Abb[179]*m_z12 - 2.*m_z1k;

  Box1Abb[181]=Box1Abb[180]*m_x_2 + Box1Abb[173]*m_x_3 + Box1Abb[177]*m_x*m_z12 - 2.*m_x_4*m_z12 + Box1Abb[171]*m_z12_2*m_z1k;

  Box1Abb[182]=Box1Abb[181]*m_cL + Box1Abb[169]*Box1Abb[7]*m_cR;

  Box1Abb[183]=1. + m_x;

  Box1Abb[184]=pow(Box1Abb[1],2.) - 2.*Box1Abb[183]*m_z1k + m_z1k_2;

  Box1Abb[185]=1. + Box1Abb[1]*m_z12;

  Box1Abb[186]=-2. + Box1Abb[123]*m_x + m_z12;

  Box1Abb[187]=-3.*m_x + m_z12;

  Box1Abb[188]=-1. + Box1Abb[187]*m_x + m_z12;

  Box1Abb[189]=2. + 3.*m_x - m_z12;

  Box1Abb[190]=-Box1Abb[1]*Box1Abb[185]*m_x + Box1Abb[186]*m_x*m_z1k + Box1Abb[188]*m_z1k_2 + Box1Abb[189]*m_z1k_3 - m_z1k_4;

  Box1Abb[191]=Box1Abb[190]*m_cR - Box1Abb[117]*Box1Abb[184]*m_cL*m_z1k;

  Box1Abb[192]=-m_x + m_z1k;

  Box1Abb[193]=-Box1Abb[108]*m_x + m_x_2 + m_z12 + Box1Abb[0]*m_z1k + m_z1k_2;

  Box1Abb[194]=Box1Abb[116]*m_cL - Box1Abb[193]*m_cR;

  Box1Abb[195]=-Box1Abb[0]*m_x_2 + 2.*m_x_3 + Box1Abb[21]*m_x*m_z1k + Box1Abb[4]*m_z12*m_z1k;

  Box1Abb[196]=5. + 8.*m_z12 + 4.*m_z1k;

  Box1Abb[197]=1. + Box1Abb[0]*m_z12;

  Box1Abb[198]=-1. + 8.*m_z12;

  Box1Abb[199]=-3. + m_z12 + 3.*m_z12_2 + 4.*Box1Abb[197]*m_z1k + Box1Abb[198]*m_z1k_2 - 4.*m_z1k_3;

  Box1Abb[200]=3. + 2.*m_z12;

  Box1Abb[201]=m_z12 - Box1Abb[12]*m_z1k + Box1Abb[200]*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[202]=6. + m_z12 + 8.*m_z1k;

  Box1Abb[203]=8. - 2.*Box1Abb[202]*m_z12 + m_z1k;

  Box1Abb[204]=-Box1Abb[201]*Box1Abb[4] + Box1Abb[199]*m_x + Box1Abb[203]*m_x_2 + Box1Abb[196]*m_x_3 - 2.*m_x_4;

  Box1Abb[205]=Box1Abb[116]*Box1Abb[195]*m_cL + Box1Abb[204]*m_cR*m_x;

  Box1Abb[206]=-1. + 2.*m_x + 2.*m_z1k;

  Box1Abb[207]=7. + m_z12 + 2.*m_z1k;

  Box1Abb[208]=3. - 5.*m_z12;

  Box1Abb[209]=-1. + 2.*m_z12 - Box1Abb[0]*m_z12*m_z1k + Box1Abb[208]*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[210]=3. + m_z1k;

  Box1Abb[211]=-4. + 3.*Box1Abb[68]*m_z12 + 2.*Box1Abb[210]*m_z1k;

  Box1Abb[212]=Box1Abb[209]*m_x + Box1Abb[211]*m_x_2 + Box1Abb[207]*m_x_3 - 2.*m_x_4 + Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[213]=Box1Abb[212]*m_cR + Box1Abb[116]*Box1Abb[206]*m_cL*m_x;

  Box1Abb[214]=2. - 5.*m_z12 + 6.*Box1Abb[0]*m_z1k;

  Box1Abb[215]=-2. + m_z12 - 6.*m_z1k + 4.*m_z12*m_z1k;

  Box1Abb[216]=9. - 4.*m_z1k;

  Box1Abb[217]=-5. + Box1Abb[216]*m_z1k;

  Box1Abb[218]=2.*pow(Box1Abb[68],2.) + Box1Abb[217]*m_z12 - 2.*Box1Abb[68]*m_z12_2 + m_z12_3;

  Box1Abb[219]=-Box1Abb[215]*m_x_3 + Box1Abb[79]*m_x_4 + Box1Abb[218]*m_x*m_z1k + Box1Abb[214]*m_x_2*m_z1k + pow(Box1Abb[4],2.)*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[220]=-4. + 7.*m_z12;

  Box1Abb[221]=27. + 5.*m_z12;

  Box1Abb[222]=-1. + 6.*m_z12;

  Box1Abb[223]=5. + 2.*Box1Abb[222]*m_z12;

  Box1Abb[224]=3.*Box1Abb[0]*m_z12_2 + Box1Abb[220]*m_z12*m_z1k + Box1Abb[0]*Box1Abb[221]*m_z12*m_z1k_2 + 2.*Box1Abb[223]*m_z1k_3 + 5.*Box1Abb[133]*m_z1k_4;

  Box1Abb[225]=8. + 3.*m_z12 + 6.*m_z1k;

  Box1Abb[226]=8. - Box1Abb[225]*m_z12 + 10.*m_z1k;

  Box1Abb[227]=7. + 13.*m_z1k;

  Box1Abb[228]=8. + 5.*m_z1k;

  Box1Abb[229]=5. + 3.*Box1Abb[228]*m_z1k;

  Box1Abb[230]=9. + 10.*m_z1k;

  Box1Abb[231]=6. + Box1Abb[230]*m_z1k;

  Box1Abb[232]=-2.*Box1Abb[231] + Box1Abb[229]*m_z12 + Box1Abb[227]*m_z12_2 + m_z12_3;

  Box1Abb[233]=3. + 4.*m_z1k;

  Box1Abb[234]=3. + 10.*m_z1k;

  Box1Abb[235]=3. + Box1Abb[234]*m_z1k;

  Box1Abb[236]=17. - 20.*Box1Abb[15]*m_z1k;

  Box1Abb[237]=3. + Box1Abb[236]*m_z1k;

  Box1Abb[238]=2. + Box1Abb[237]*m_z12 - Box1Abb[233]*Box1Abb[66]*m_z12_2 - 3.*Box1Abb[15]*m_z12_3 + 2.*Box1Abb[235]*m_z1k;

  Box1Abb[239]=1. - 5.*m_z1k + 6.*m_z1k_3;

  Box1Abb[240]=9. + 5.*m_z1k;

  Box1Abb[241]=-3. + Box1Abb[240]*m_z1k;

  Box1Abb[242]=1. + Box1Abb[241]*m_z1k;

  Box1Abb[243]=11. + 13.*m_z1k;

  Box1Abb[244]=-8. + Box1Abb[243]*m_z1k;

  Box1Abb[245]=2. + Box1Abb[244]*m_z1k;

  Box1Abb[246]=Box1Abb[239]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[245]*Box1Abb[68]*m_z12_2 + Box1Abb[242]*m_z12_3 - 2.*pow(Box1Abb[68],3.)*m_z1k_2;

  Box1Abb[247]=-Box1Abb[246]*m_x + Box1Abb[224]*m_x_2 + Box1Abb[238]*m_x_3 + Box1Abb[232]*m_x_4 + Box1Abb[226]*m_x_5 + Box1Abb[79]*m_x_6 + Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[72]*m_z12*m_z1k_2;

  Box1Abb[248]=Box1Abb[116]*Box1Abb[219]*m_cL - Box1Abb[247]*m_cR;

  Box1Abb[249]=-1. + m_z12 - 4.*m_z1k;

  Box1Abb[250]=-2. + Box1Abb[249]*m_z12 - 6.*m_z1k;

  Box1Abb[251]=3. + m_z12;

  Box1Abb[252]=1. - 2.*Box1Abb[0]*m_z12 - 2.*m_z1k + Box1Abb[251]*m_z12*m_z1k + 6.*Box1Abb[12]*m_z1k_2;

  Box1Abb[253]=-3. + m_z12 - m_z12_2 + m_z12_3;

  Box1Abb[254]=4. + m_z12 - 3.*m_z12_2;

  Box1Abb[255]=pow(Box1Abb[0],2.) + Box1Abb[253]*m_z1k + Box1Abb[254]*m_z1k_2 - 2.*Box1Abb[14]*m_z1k_3;

  Box1Abb[256]=Box1Abb[255]*m_x + Box1Abb[252]*m_x_2 + Box1Abb[250]*m_x_3 + Box1Abb[77]*m_x_4 + Box1Abb[171]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[257]=5. + Box1Abb[77]*m_z12;

  Box1Abb[258]=6. + 5.*m_z12;

  Box1Abb[259]=12. + m_z12 - 3.*m_z12_2 + 2.*Box1Abb[257]*m_z1k + 2.*Box1Abb[258]*m_z1k_2;

  Box1Abb[260]=-2. + m_z12 + 5.*m_z12_2;

  Box1Abb[261]=4. + 5.*m_z12;

  Box1Abb[262]=-3.*Box1Abb[12]*Box1Abb[79] + Box1Abb[133]*Box1Abb[251]*m_z1k + 2.*Box1Abb[260]*m_z1k_2 + 2.*Box1Abb[261]*m_z1k_3;

  Box1Abb[263]=-1. + m_z12 - 5.*m_z1k;

  Box1Abb[264]=-8.*Box1Abb[15] + Box1Abb[263]*m_z12;

  Box1Abb[265]=1. + 3.*m_z1k;

  Box1Abb[266]=Box1Abb[265]*Box1Abb[68]*m_z12 + pow(Box1Abb[68],2.)*m_z1k + 2.*m_z12_2*m_z1k;

  Box1Abb[267]=-4. + 5.*m_z1k;

  Box1Abb[268]=-2. + Box1Abb[267]*m_z1k;

  Box1Abb[269]=-7. + 10.*m_z1k;

  Box1Abb[270]=2. + Box1Abb[269]*m_z1k;

  Box1Abb[271]=-1. + Box1Abb[270]*m_z1k;

  Box1Abb[272]=Box1Abb[15]*Box1Abb[268]*Box1Abb[68]*m_z12 + Box1Abb[271]*m_z12_2 + 2.*pow(Box1Abb[68],3.)*m_z1k + 2.*m_z12_3*m_z1k_2;

  Box1Abb[273]=Box1Abb[272]*m_x - Box1Abb[262]*m_x_2 + Box1Abb[259]*m_x_3 + Box1Abb[264]*m_x_4 + Box1Abb[77]*m_x_5 - Box1Abb[266]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[274]=-Box1Abb[116]*Box1Abb[256]*m_cL + Box1Abb[273]*Box1Abb[7]*m_cR;

  Box1Abb[295]=(-4.*m_m)/m_s_2;

  Box1Abb[296]=(4.*m_m)/m_s_2;

  Box1Abb[297]=2.*m_m;

  Box1Abb[298]=(-4.*m_m*m_z1k)/m_s_2;

  Box1Abb[299]=Box1Abb[109]*m_x - 2.*Box1Abb[12]*m_x_2 + m_z12 - Box1Abb[9]*m_z12;

  Box1Abb[300]=Box1Abb[107]*m_cL + Box1Abb[299]*m_cR;

  Box1Abb[301]=1. + 2.*Box1Abb[12]*m_x - 2.*m_z12;

  Box1Abb[302]=-1. + 2.*Box1Abb[12]*m_x + m_z12 - m_z12_2;

  Box1Abb[303]=-2.*Box1Abb[12]*m_x + m_z12;

  Box1Abb[304]=2.*Box1Abb[12]*m_x - m_z12_2;

  Box1Abb[305]=Box1Abb[302]*m_cL*m_x - Box1Abb[301]*m_cR*m_x + Box1Abb[303]*m_cL*m_z1k + Box1Abb[304]*m_cR*m_z1k;

  Box1Abb[306]=Box1Abb[121]*m_cL + Box1Abb[116]*Box1Abb[117]*m_cR;

  Box1Abb[307]=-2.*Box1Abb[0]*m_x - 3.*m_x_2 + Box1Abb[69]*m_z12;

  Box1Abb[308]=Box1Abb[1]*Box1Abb[123]*m_x + Box1Abb[307]*m_z1k + Box1Abb[126]*m_z1k_2 - m_z1k_3;

  Box1Abb[309]=-Box1Abb[308]*m_cL + Box1Abb[116]*Box1Abb[7]*m_cR;

  Box1Abb[310]=Box1Abb[136]*m_cL + Box1Abb[132]*m_cR;

  Box1Abb[311]=Box1Abb[147]*m_cL + Box1Abb[143]*m_cR;

  Box1Abb[312]=Box1Abb[163]*m_cL + Box1Abb[155]*m_cR;

  Box1Abb[313]=Box1Abb[169]*Box1Abb[7]*m_cL + Box1Abb[181]*m_cR;

  Box1Abb[314]=Box1Abb[190]*m_cL - Box1Abb[117]*Box1Abb[184]*m_cR*m_z1k;

  Box1Abb[315]=-Box1Abb[193]*m_cL + Box1Abb[116]*m_cR;

  Box1Abb[316]=-8. + 2.*Box1Abb[202]*m_z12 - m_z1k;

  Box1Abb[317]=-1. + 2.*m_z1k;

  Box1Abb[318]=3. - 2.*m_z1k;

  Box1Abb[319]=-1. + Box1Abb[318]*m_z1k;

  Box1Abb[320]=m_z12 + Box1Abb[319]*m_z1k + Box1Abb[317]*m_z12*m_z1k;

  Box1Abb[321]=-1. + 4.*m_z1k - 8.*m_z1k_2;

  Box1Abb[322]=-4. + m_z1k + 4.*m_z1k_2;

  Box1Abb[323]=3. + Box1Abb[321]*m_z12 - Box1Abb[233]*m_z12_2 + Box1Abb[322]*m_z1k;

  Box1Abb[324]=Box1Abb[320]*Box1Abb[4] + Box1Abb[323]*m_x + Box1Abb[316]*m_x_2 - Box1Abb[196]*m_x_3 + 2.*m_x_4;

  Box1Abb[325]=-Box1Abb[116]*Box1Abb[195]*m_cR + Box1Abb[324]*m_cL*m_x;

  Box1Abb[326]=-3. + 5.*m_z12;

  Box1Abb[327]=1. - 2.*m_z12 + Box1Abb[0]*m_z12*m_z1k + Box1Abb[326]*m_z1k_2 + 2.*m_z1k_3;

  Box1Abb[328]=4. - 3.*Box1Abb[68]*m_z12 - 2.*Box1Abb[210]*m_z1k;

  Box1Abb[329]=Box1Abb[327]*m_x + Box1Abb[328]*m_x_2 - Box1Abb[207]*m_x_3 + 2.*m_x_4 - Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[330]=Box1Abb[329]*m_cL - Box1Abb[116]*Box1Abb[206]*m_cR*m_x;

  Box1Abb[331]=Box1Abb[247]*m_cL - Box1Abb[116]*Box1Abb[219]*m_cR;

  Box1Abb[332]=-Box1Abb[273]*Box1Abb[7]*m_cL + Box1Abb[116]*Box1Abb[256]*m_cR;

  Box1Abb[353]=Box1Abb[117]*m_cL + Box1Abb[192]*m_cR;

  Box1Abb[354]=-1. - 2.*m_x + m_z12 + 2.*m_z1k;

  Box1Abb[355]=Box1Abb[22]*m_cL + Box1Abb[354]*m_cR*m_x;

  Box1Abb[356]=m_x - 2.*m_x_2 + 2.*m_x*m_z1k - m_z12*m_z1k;

  Box1Abb[357]=Box1Abb[356]*m_cR + Box1Abb[20]*m_cL*m_x;

  Box1Abb[358]=-2. + m_z12 + m_z1k;

  Box1Abb[359]=pow(Box1Abb[0],2.) + 3.*m_z12*m_z1k;

  Box1Abb[360]=1. + m_z12 + 3.*m_z12*m_z1k;

  Box1Abb[361]=m_x - Box1Abb[360]*m_x_2 + m_x_3*m_z12 + Box1Abb[359]*m_x*m_z1k - Box1Abb[358]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[362]=1. + m_z12 + Box1Abb[79]*m_z1k + 3.*m_z1k_2;

  Box1Abb[363]=-1. + Box1Abb[362]*m_z12 + m_z1k;

  Box1Abb[364]=-Box1Abb[363]*m_x + Box1Abb[360]*m_x_2 + Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12 - m_x_3*m_z12;

  Box1Abb[365]=Box1Abb[364]*m_cL + Box1Abb[361]*m_cR;

  Box1Abb[366]=-Box1Abb[4]*pow(Box1Abb[68],2.) + m_x_3 + Box1Abb[15]*m_x*m_z12 + 3.*Box1Abb[68]*m_x*m_z1k - 3.*m_x_2*m_z1k;

  Box1Abb[367]=-1. + m_z12 + 3.*m_z1k;

  Box1Abb[368]=Box1Abb[265]*m_x_2 - m_x_3 + Box1Abb[4]*Box1Abb[68]*m_z1k - Box1Abb[367]*m_x*m_z1k;

  Box1Abb[369]=Box1Abb[366]*m_cL + Box1Abb[368]*m_cR;

  Box1Abb[370]=1. + m_z12 + 3.*m_z1k;

  Box1Abb[371]=m_z12 + Box1Abb[0]*m_z1k + m_z1k_2;

  Box1Abb[372]=m_x - Box1Abb[119]*m_x_2 + m_x_3 - Box1Abb[371]*m_z1k + Box1Abb[370]*m_x*m_z1k;

  Box1Abb[373]=Box1Abb[368]*m_cL + Box1Abb[372]*m_cR;

  Box1Abb[374]=2. - 2.*Box1Abb[183]*m_z12 + m_z12_2;

  Box1Abb[375]=3. + m_x;

  Box1Abb[376]=1. + 2.*m_x;

  Box1Abb[377]=-2.*m_x + m_z12 + 2.*Box1Abb[375]*m_x*m_z12 - 2.*Box1Abb[376]*m_z12_2 + m_z12_3;

  Box1Abb[378]=Box1Abb[377]*m_cL + Box1Abb[374]*m_cR*m_x;

  Box1Abb[379]=m_cL + m_cR + 2.*m_cL*m_x - 2.*m_cR*m_x - m_cL*m_z12;

  Box1Abb[380]=-2.*m_x + m_z12;

  Box1Abb[381]=Box1Abb[380]*m_cR + 2.*m_cL*m_x;

  Box1Abb[382]=-Box1Abb[378]*m_x + 2.*Box1Abb[379]*m_x*m_z12*m_z1k - Box1Abb[381]*m_z12*m_z1k_2;

  Box1Abb[383]=Box1Abb[133]*m_x + m_z12 - Box1Abb[9]*m_z12;

  Box1Abb[384]=Box1Abb[79]*m_x + m_z12*m_z1k;

  Box1Abb[385]=2. + Box1Abb[170]*m_z12 - 4.*m_z1k;

  Box1Abb[386]=-1. + m_z12_2;

  Box1Abb[387]=pow(Box1Abb[0],2.) + 2.*Box1Abb[386]*m_z1k + Box1Abb[77]*m_z1k_2;

  Box1Abb[388]=Box1Abb[387]*m_x + Box1Abb[385]*m_x_2 - Box1Abb[79]*m_x_3 - pow(Box1Abb[4],2.)*m_z12*m_z1k;

  Box1Abb[389]=Box1Abb[388]*m_cL + Box1Abb[384]*pow(Box1Abb[7],2.)*m_cR;

  Box1Abb[390]=-4.*m_x + m_z12;

  Box1Abb[391]=-Box1Abb[384]*Box1Abb[7]*m_cR + Box1Abb[0]*Box1Abb[20]*m_cL*m_x;

  Box1Abb[392]=-2. + 3.*m_z12 + m_z1k;

  Box1Abb[393]=-10. + 9.*m_z12 + 4.*m_z1k;

  Box1Abb[394]=2. + Box1Abb[393]*m_z12;

  Box1Abb[395]=-Box1Abb[394]*m_x_2 + 2.*Box1Abb[392]*Box1Abb[4]*m_x*m_z12 + 2.*m_x_3*m_z12 - pow(Box1Abb[4],2.)*m_z12_2;

  Box1Abb[396]=-6. + m_z1k;

  Box1Abb[397]=4. + Box1Abb[396]*m_z12 + 2.*m_z12_2 + 2.*Box1Abb[170]*m_z1k;

  Box1Abb[398]=-8. + 5.*m_z12 + 6.*m_z1k;

  Box1Abb[399]=2. + Box1Abb[398]*m_z12;

  Box1Abb[400]=-1. + 7.*m_z1k;

  Box1Abb[401]=-2. + Box1Abb[400]*m_z12 + m_z12_2 + 6.*Box1Abb[170]*m_z1k;

  Box1Abb[402]=2.*Box1Abb[15] + Box1Abb[401]*m_z12;

  Box1Abb[403]=-Box1Abb[402]*m_x_2 + Box1Abb[399]*m_x_3 - 2.*m_x_4*m_z12 + Box1Abb[397]*m_x*m_z12*m_z1k + Box1Abb[4]*m_z12_2*m_z1k_2;

  Box1Abb[404]=Box1Abb[395]*Box1Abb[7]*m_cL + Box1Abb[403]*m_cR;

  Box1Abb[405]=-2. + m_z12 + 4.*m_z1k;

  Box1Abb[406]=2. + Box1Abb[405]*m_z12;

  Box1Abb[407]=-Box1Abb[406]*m_x_2 + 2.*m_x_3*m_z12 + 2.*Box1Abb[15]*m_x*m_z12*m_z1k - m_z12_2*m_z1k_2;

  Box1Abb[408]=pow(Box1Abb[0],2.) + 4.*Box1Abb[0]*m_z12*m_z1k + 3.*m_z12*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[409]=-4. + m_z12 + 6.*m_z1k;

  Box1Abb[410]=2. + Box1Abb[409]*m_z12;

  Box1Abb[411]=4. + m_z1k;

  Box1Abb[412]=-6. + Box1Abb[411]*m_z12 + 4.*m_z1k - 6.*m_z1k_2;

  Box1Abb[413]=2. + Box1Abb[412]*m_z12 - 2.*m_z1k;

  Box1Abb[414]=-Box1Abb[413]*m_x_2 - Box1Abb[410]*m_x_3 + Box1Abb[408]*m_x*m_z12 + 2.*m_x_4*m_z12 - pow(Box1Abb[4],2.)*m_z12_2*m_z1k;

  Box1Abb[415]=Box1Abb[414]*m_cL - Box1Abb[407]*Box1Abb[7]*m_cR;

  Box1Abb[416]=Box1Abb[77]*m_x + Box1Abb[68]*m_z12;

  Box1Abb[417]=3. + 2.*m_x - 2.*m_z12;

  Box1Abb[418]=-1. + m_x + pow(Box1Abb[1],2.)*m_z12 + Box1Abb[417]*m_z1k - Box1Abb[77]*m_z1k_2;

  Box1Abb[419]=-Box1Abb[418]*m_cL*m_x + Box1Abb[416]*Box1Abb[7]*m_cR*m_z1k;

  Box1Abb[420]=Box1Abb[4]*Box1Abb[68] + Box1Abb[156]*m_x + m_x_2;

  Box1Abb[421]=1. - 2.*m_z12 + m_z1k - 2.*m_z1k_2;

  Box1Abb[422]=Box1Abb[4]*pow(Box1Abb[68],2.) + Box1Abb[421]*m_x + Box1Abb[9]*m_x_2;

  Box1Abb[423]=Box1Abb[422]*m_cL - Box1Abb[420]*m_cR*m_z1k;

  Box1Abb[424]=2.*m_z12 + m_z1k;

  Box1Abb[425]=m_x - m_x_2 + m_z1k - Box1Abb[424]*m_z1k + 2.*m_x*m_z1k;

  Box1Abb[426]=Box1Abb[420]*m_cL + Box1Abb[425]*m_cR;

  Box1Abb[427]=-1. + m_x + m_z1k + 2.*m_z12*m_z1k;

  Box1Abb[428]=2. - 3.*m_z12 - 2.*m_z12*m_z1k;

  Box1Abb[429]=-Box1Abb[4]*Box1Abb[68] + Box1Abb[428]*m_x + Box1Abb[14]*m_x_2;

  Box1Abb[430]=Box1Abb[429]*m_cL - Box1Abb[427]*Box1Abb[7]*m_cR;

  Box1Abb[431]=1. + 2.*m_z1k;

  Box1Abb[432]=3. + 2.*m_z12 - 2.*m_z1k;

  Box1Abb[433]=-1. + Box1Abb[432]*m_z1k;

  Box1Abb[434]=m_z12 + 2.*Box1Abb[15]*m_z1k - 8.*m_z12*m_z1k;

  Box1Abb[435]=Box1Abb[4]*Box1Abb[433] + Box1Abb[434]*m_x + Box1Abb[431]*m_x_2 - 2.*m_x_3;

  Box1Abb[436]=1. + 7.*m_z12 - 2.*m_z1k;

  Box1Abb[437]=-3. + 2.*m_z1k;

  Box1Abb[438]=-1. + 5.*m_z1k;

  Box1Abb[439]=Box1Abb[437]*pow(Box1Abb[68],2.) + Box1Abb[438]*Box1Abb[68]*m_z12 + Box1Abb[210]*m_z12_2;

  Box1Abb[440]=11. + 3.*m_z1k;

  Box1Abb[441]=-4. + m_z1k + m_z1k_2;

  Box1Abb[442]=2.*Box1Abb[441] + Box1Abb[440]*m_z12 + 2.*m_z12_2;

  Box1Abb[443]=Box1Abb[439]*m_x - Box1Abb[442]*m_x_2 + Box1Abb[436]*m_x_3 + 2.*m_x_4 - Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12;

  Box1Abb[444]=Box1Abb[443]*m_cL + Box1Abb[435]*m_cR*m_x;

  Box1Abb[445]=1. + m_x + 2.*m_x_2 - m_z12;

  Box1Abb[446]=m_x + m_x_2 - 4.*m_x*m_z12 + m_z12_2;

  Box1Abb[447]=-4. + 2.*Box1Abb[446] + m_z12;

  Box1Abb[448]=5. + 2.*m_x;

  Box1Abb[449]=Box1Abb[1]*Box1Abb[445] - Box1Abb[447]*m_z1k - Box1Abb[448]*m_z1k_2 + 2.*m_z1k_3;

  Box1Abb[450]=-1. + m_z12 + Box1Abb[77]*m_z12*m_z1k + Box1Abb[208]*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[451]=2. + 5.*m_z1k;

  Box1Abb[452]=-Box1Abb[451]*m_z12 + 2.*Box1Abb[15]*m_z1k;

  Box1Abb[453]=Box1Abb[450]*m_x + Box1Abb[452]*m_x_2 + Box1Abb[108]*m_x_3 - 2.*m_x_4 + Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[454]=Box1Abb[453]*m_cR + Box1Abb[449]*m_cL*m_x;

  Box1Abb[455]=-4. + m_z12 - 8.*m_z1k;

  Box1Abb[456]=-2. + Box1Abb[455]*m_z12;

  Box1Abb[457]=-1. + 3.*m_z1k;

  Box1Abb[458]=1. + m_z1k_2;

  Box1Abb[459]=2.*Box1Abb[458]*pow(Box1Abb[68],2.) + Box1Abb[317]*Box1Abb[457]*Box1Abb[68]*m_z12 + 4.*m_z12_2*m_z1k_2;

  Box1Abb[460]=-3. + 6.*m_z1k;

  Box1Abb[461]=1. + m_z1k + 3.*m_z1k_2;

  Box1Abb[462]=6. + 4.*Box1Abb[461]*m_z12 + Box1Abb[460]*m_z12_2 + 4.*m_z1k;

  Box1Abb[463]=-3. + 14.*m_z1k;

  Box1Abb[464]=4. + 4.*Box1Abb[431]*Box1Abb[68]*m_z1k;

  Box1Abb[465]=2.*Box1Abb[170]*Box1Abb[68] + Box1Abb[464]*m_z12 + Box1Abb[15]*Box1Abb[463]*m_z12_2;

  Box1Abb[466]=-Box1Abb[465]*m_x_2 + Box1Abb[462]*m_x_3 + Box1Abb[456]*m_x_4 + Box1Abb[459]*m_x*m_z12 + 2.*m_x_5*m_z12 + pow(Box1Abb[68],3.)*m_z12_2*m_z1k;

  Box1Abb[467]=-1. + 2.*m_z12_2 - 5.*m_z12*m_z1k;

  Box1Abb[468]=-11. + 3.*Box1Abb[77]*m_z12;

  Box1Abb[469]=11. + 3.*m_z12;

  Box1Abb[470]=1. + Box1Abb[158]*m_z12;

  Box1Abb[471]=2. + Box1Abb[468]*m_z12 + Box1Abb[0]*Box1Abb[469]*m_z12*m_z1k - 6.*Box1Abb[470]*m_z1k_2 - 20.*m_z12*m_z1k_3;

  Box1Abb[472]=-7. + 3.*m_z12;

  Box1Abb[473]=pow(Box1Abb[0],2.) - 2.*Box1Abb[0]*Box1Abb[79]*m_z1k + Box1Abb[0]*Box1Abb[472]*m_z1k_2 + 6.*Box1Abb[0]*m_z1k_3 + 2.*m_z1k_4;

  Box1Abb[474]=2. - 5.*m_z1k;

  Box1Abb[475]=5. + 3.*m_z1k;

  Box1Abb[476]=-5. + 2.*Box1Abb[475]*m_z12 + m_z12_2 + 4.*Box1Abb[474]*m_z1k;

  Box1Abb[477]=-6.*Box1Abb[15] + Box1Abb[476]*m_z12;

  Box1Abb[478]=-3. + m_z1k_2;

  Box1Abb[479]=1. + m_z1k + 7.*Box1Abb[68]*m_z1k_2;

  Box1Abb[480]=-7. + 5.*m_z1k;

  Box1Abb[481]=5. + 2.*Box1Abb[480]*m_z1k;

  Box1Abb[482]=-3. + Box1Abb[481]*m_z1k;

  Box1Abb[483]=2.*pow(Box1Abb[68],3.) + Box1Abb[482]*Box1Abb[68]*m_z12 + 2.*Box1Abb[479]*m_z12_2 + Box1Abb[478]*m_z12_3;

  Box1Abb[484]=Box1Abb[483]*m_x + Box1Abb[471]*m_x_2 - Box1Abb[477]*m_x_3 + 2.*Box1Abb[467]*m_x_4 - Box1Abb[473]*Box1Abb[68]*m_z12 + 2.*m_x_5*m_z12;

  Box1Abb[485]=Box1Abb[466]*Box1Abb[7]*m_cR - Box1Abb[484]*m_cL*m_x;

  Box1Abb[486]=-3. + 4.*m_z12 - 13.*m_z1k;

  Box1Abb[487]=4. + Box1Abb[486]*m_z12 + 8.*m_z1k;

  Box1Abb[488]=10. + m_z12;

  Box1Abb[489]=-8. + Box1Abb[488]*m_z12;

  Box1Abb[490]=6. - 11.*m_z12;

  Box1Abb[491]=Box1Abb[489]*m_z12 + 2.*Box1Abb[160]*m_z12*m_z1k + 2.*Box1Abb[490]*m_z1k_2;

  Box1Abb[492]=-4. + Box1Abb[77]*m_z12;

  Box1Abb[493]=-3. + m_z12;

  Box1Abb[494]=6. + 5.*Box1Abb[493]*m_z12;

  Box1Abb[495]=4. - 9.*m_z12;

  Box1Abb[496]=4. + 3.*Box1Abb[492]*m_z12 + 3.*m_z12_3*m_z1k - 2.*Box1Abb[494]*m_z1k_2 + 2.*Box1Abb[495]*m_z1k_3;

  Box1Abb[497]=-2.*pow(Box1Abb[68],2.) + Box1Abb[400]*Box1Abb[68]*m_z12 + 3.*Box1Abb[15]*m_z12_2;

  Box1Abb[498]=Box1Abb[4]*Box1Abb[497]*Box1Abb[68]*m_x + Box1Abb[496]*m_x_2 - Box1Abb[491]*m_x_3 + Box1Abb[487]*m_x_4 + Box1Abb[133]*m_x_5 - pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12;

  Box1Abb[499]=-8. + m_z12 - 16.*m_z1k;

  Box1Abb[500]=6. + Box1Abb[499]*m_z12 + 10.*m_z1k;

  Box1Abb[501]=7. + 3.*m_z12;

  Box1Abb[502]=3. + Box1Abb[501]*m_z12;

  Box1Abb[503]=-14. + 19.*m_z12;

  Box1Abb[504]=6. + Box1Abb[503]*m_z12;

  Box1Abb[505]=-2. + 5.*m_z12;

  Box1Abb[506]=-2. - Box1Abb[493]*m_z12 + 6.*m_z1k + 3.*Box1Abb[129]*m_z12*m_z1k + 2.*Box1Abb[0]*Box1Abb[502]*m_z1k_2 + 2.*Box1Abb[504]*m_z1k_3 + 5.*Box1Abb[505]*m_z1k_4;

  Box1Abb[507]=-1. - 2.*m_z1k + 4.*m_z1k_2;

  Box1Abb[508]=-2. + m_z1k + 19.*m_z1k_2;

  Box1Abb[509]=2.*Box1Abb[507]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[508]*Box1Abb[68]*m_z12_2 - 2.*pow(Box1Abb[68],3.)*m_z1k + 4.*Box1Abb[431]*m_z12_3*m_z1k;

  Box1Abb[510]=-8. + 6.*m_z1k + 4.*m_z1k_2 - 40.*m_z1k_3;

  Box1Abb[511]=7. + 15.*m_z1k;

  Box1Abb[512]=3. - 2.*Box1Abb[511]*m_z1k;

  Box1Abb[513]=4. + Box1Abb[510]*m_z12 + Box1Abb[512]*m_z12_2 + 6.*m_z1k + 20.*m_z1k_3;

  Box1Abb[514]=3. + 5.*m_z1k;

  Box1Abb[515]=-3. + 7.*m_z1k;

  Box1Abb[516]=16. + 35.*m_z1k;

  Box1Abb[517]=10. + Box1Abb[516]*m_z1k;

  Box1Abb[518]=-6. + Box1Abb[517]*m_z12 + Box1Abb[515]*m_z12_2 - 4.*Box1Abb[514]*m_z1k;

  Box1Abb[519]=Box1Abb[506]*m_x_2 + Box1Abb[513]*m_x_3 + Box1Abb[518]*m_x_4 + Box1Abb[500]*m_x_5 + Box1Abb[133]*m_x_6 - Box1Abb[509]*m_x*m_z1k + Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[72]*m_z12*m_z1k_2;

  Box1Abb[520]=-Box1Abb[498]*Box1Abb[7]*m_cL + Box1Abb[519]*m_cR;

  Box1Abb[521]=-4. + m_z12 - 10.*m_z1k + 4.*m_z12*m_z1k;

  Box1Abb[522]=7. + m_z12;

  Box1Abb[523]=-11. + Box1Abb[522]*m_z12;

  Box1Abb[524]=-5. + 8.*m_z12;

  Box1Abb[525]=6. + Box1Abb[524]*m_z12;

  Box1Abb[526]=1. + 3.*Box1Abb[0]*m_z12 - 3.*m_z1k + Box1Abb[523]*m_z12*m_z1k + 2.*Box1Abb[525]*m_z1k_2 - 20.*m_z1k_3;

  Box1Abb[527]=4. - 3.*m_z12;

  Box1Abb[528]=15. + m_z12 + m_z12_2 - 5.*m_z12_3;

  Box1Abb[529]=-10. + m_z12 - 6.*m_z12_2;

  Box1Abb[530]=-1. + Box1Abb[527]*m_z12 - 4.*m_z1k + Box1Abb[527]*m_z12_2*m_z1k + Box1Abb[528]*m_z1k_2 + 2.*Box1Abb[529]*m_z1k_3 + 5.*Box1Abb[77]*m_z1k_4;

  Box1Abb[531]=3. - 5.*m_z1k;

  Box1Abb[532]=1. + 5.*m_z1k;

  Box1Abb[533]=1. + m_z12 - Box1Abb[66]*m_z12_2 + 4.*Box1Abb[532]*m_z1k + Box1Abb[531]*m_z12*m_z1k;

  Box1Abb[534]=1. + 2.*Box1Abb[68]*m_z1k;

  Box1Abb[535]=-1. + Box1Abb[175]*m_z1k;

  Box1Abb[536]=-2. + 3.*m_z1k;

  Box1Abb[537]=3. + Box1Abb[536]*m_z1k;

  Box1Abb[538]=-2. + Box1Abb[233]*m_z1k;

  Box1Abb[539]=-Box1Abb[534]*pow(Box1Abb[68],3.) - Box1Abb[538]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[535]*Box1Abb[68]*m_z12_2 + Box1Abb[537]*m_z12_3*m_z1k;

  Box1Abb[540]=Box1Abb[539]*m_x + Box1Abb[530]*m_x_2 + Box1Abb[526]*m_x_3 + Box1Abb[533]*m_x_4 + Box1Abb[521]*m_x_5 - Box1Abb[79]*m_x_6 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[541]=pow(Box1Abb[68],2.) + 4.*Box1Abb[68]*m_z12 + 2.*m_z12_2;

  Box1Abb[542]=-6. + m_z12;

  Box1Abb[543]=-6.*Box1Abb[15] + 3.*m_z12 + 4.*m_z12_2*m_z1k + 2.*Box1Abb[542]*m_z1k_2;

  Box1Abb[544]=6. - 3.*Box1Abb[15]*m_z12 + 8.*m_z1k;

  Box1Abb[545]=-1. + m_z1k + m_z1k_2;

  Box1Abb[546]=-8. + 11.*m_z1k - 3.*m_z1k_3;

  Box1Abb[547]=-2.*pow(Box1Abb[68],3.) + Box1Abb[546]*m_z12 - 4.*Box1Abb[545]*m_z12_2 + 2.*m_z12_3*m_z1k;

  Box1Abb[548]=11. + 2.*Box1Abb[210]*m_z1k;

  Box1Abb[549]=-1. + Box1Abb[548]*m_z1k;

  Box1Abb[550]=2. + Box1Abb[549]*m_z12 - 4.*Box1Abb[178]*m_z12_2*m_z1k - 6.*m_z1k_2 + 8.*m_z1k_3;

  Box1Abb[551]=Box1Abb[550]*m_x_2 + Box1Abb[543]*m_x_3 + Box1Abb[544]*m_x_4 + Box1Abb[79]*m_x_5 + Box1Abb[547]*m_x*m_z1k + Box1Abb[541]*Box1Abb[68]*m_z12*m_z1k_2;

  Box1Abb[552]=Box1Abb[540]*m_cL + Box1Abb[551]*Box1Abb[7]*m_cR;

  Box1Abb[582]=(4.*m_m*m_z12)/m_s_2;

  Box1Abb[583]=(-2.*m_m*m_z12)/m_s;

  Box1Abb[584]=-2.*m_m*m_z12;

  Box1Abb[585]=-2.*m_m;

  Box1Abb[586]=(-4.*m_m*m_x)/m_s_2;

  Box1Abb[587]=Box1Abb[192]*m_cL + Box1Abb[117]*m_cR;

  Box1Abb[588]=Box1Abb[22]*m_cR + Box1Abb[354]*m_cL*m_x;

  Box1Abb[589]=Box1Abb[356]*m_cL + Box1Abb[20]*m_cR*m_x;

  Box1Abb[590]=pow(Box1Abb[0],2.) + 2.*Box1Abb[0]*m_z1k + 2.*m_z1k_2;

  Box1Abb[591]=-3. + 2.*m_z12 + 2.*m_z1k;

  Box1Abb[592]=1. + Box1Abb[591]*m_z12;

  Box1Abb[593]=-2.*Box1Abb[592]*m_x + Box1Abb[590]*m_z12 + 2.*m_x_2*m_z12;

  Box1Abb[594]=Box1Abb[406]*m_x_2 - 2.*m_x_3*m_z12 - 2.*Box1Abb[15]*m_x*m_z12*m_z1k + m_z12_2*m_z1k_2;

  Box1Abb[595]=Box1Abb[594]*m_cL + Box1Abb[593]*m_cR*m_x;

  Box1Abb[596]=-2. - Box1Abb[170]*m_z12 + 4.*m_z1k;

  Box1Abb[597]=-Box1Abb[387]*m_x + Box1Abb[596]*m_x_2 + Box1Abb[79]*m_x_3 + pow(Box1Abb[4],2.)*m_z12*m_z1k;

  Box1Abb[598]=-Box1Abb[384]*pow(Box1Abb[7],2.)*m_cL + Box1Abb[597]*m_cR;

  Box1Abb[599]=Box1Abb[361]*m_cL + Box1Abb[364]*m_cR;

  Box1Abb[600]=Box1Abb[368]*m_cL + Box1Abb[366]*m_cR;

  Box1Abb[601]=Box1Abb[372]*m_cL + Box1Abb[368]*m_cR;

  Box1Abb[602]=Box1Abb[384]*Box1Abb[7]*m_cL + Box1Abb[0]*Box1Abb[354]*m_cR*m_x;

  Box1Abb[603]=Box1Abb[403]*m_cL + Box1Abb[395]*Box1Abb[7]*m_cR;

  Box1Abb[604]=-Box1Abb[407]*Box1Abb[7]*m_cL + Box1Abb[414]*m_cR;

  Box1Abb[605]=-Box1Abb[418]*m_cR*m_x + Box1Abb[416]*Box1Abb[7]*m_cL*m_z1k;

  Box1Abb[606]=Box1Abb[422]*m_cR - Box1Abb[420]*m_cL*m_z1k;

  Box1Abb[607]=Box1Abb[425]*m_cL + Box1Abb[420]*m_cR;

  Box1Abb[608]=-Box1Abb[427]*Box1Abb[7]*m_cL + Box1Abb[429]*m_cR;

  Box1Abb[609]=-1. - 7.*m_z12 + 2.*m_z1k;

  Box1Abb[610]=-Box1Abb[439]*m_x + Box1Abb[442]*m_x_2 + Box1Abb[609]*m_x_3 - 2.*m_x_4 + Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12;

  Box1Abb[611]=Box1Abb[610]*m_cR + Box1Abb[449]*m_cL*m_x;

  Box1Abb[612]=Box1Abb[451]*m_z12 - 2.*Box1Abb[15]*m_z1k;

  Box1Abb[613]=-Box1Abb[450]*m_x + Box1Abb[612]*m_x_2 - Box1Abb[108]*m_x_3 + 2.*m_x_4 - Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[614]=Box1Abb[613]*m_cL + Box1Abb[435]*m_cR*m_x;

  Box1Abb[615]=Box1Abb[466]*Box1Abb[7]*m_cL - Box1Abb[484]*m_cR*m_x;

  Box1Abb[616]=Box1Abb[519]*m_cL - Box1Abb[498]*Box1Abb[7]*m_cR;

  Box1Abb[617]=Box1Abb[551]*Box1Abb[7]*m_cL + Box1Abb[540]*m_cR;

  Box1Abb[647]=(2.*m_m*m_z12)/m_s;

  Box1Abb[648]=2.*m_m*m_z12;

  Box1Abb[649]=-1. - m_x + m_z12 + m_z1k;

  Box1Abb[650]=pow(Box1Abb[0],2.) + 2.*m_z12*m_z1k;

  Box1Abb[651]=-Box1Abb[650]*m_x + m_x_2*m_z12 + m_z12*m_z1k_2;

  Box1Abb[652]=-2. + 2.*m_z12 + m_z1k;

  Box1Abb[653]=6. + m_z12;

  Box1Abb[654]=-6. + Box1Abb[653]*m_z12 + 8.*m_z1k;

  Box1Abb[655]=Box1Abb[15]*m_z12 - Box1Abb[265]*m_z12_2 + m_z1k - 2.*m_z1k_2;

  Box1Abb[656]=2.*Box1Abb[655]*m_x + Box1Abb[654]*m_x_2 - 4.*m_x_3 + Box1Abb[652]*m_z12_2*m_z1k;

  Box1Abb[657]=Box1Abb[656]*m_cL + 2.*Box1Abb[651]*m_cR*m_x;

  Box1Abb[658]=pow(Box1Abb[0],2.) + Box1Abb[527]*m_x;

  Box1Abb[659]=Box1Abb[658]*m_x + 2.*Box1Abb[79]*m_x*m_z1k + m_z12*m_z1k_2;

  Box1Abb[660]=-2. + m_z12 + 8.*m_z1k;

  Box1Abb[661]=2. + Box1Abb[660]*m_z12;

  Box1Abb[662]=8. + 7.*m_z1k;

  Box1Abb[663]=1. + m_z1k - m_z1k_2;

  Box1Abb[664]=6.*Box1Abb[663]*m_z12 - Box1Abb[662]*m_z12_2 + 2.*m_z12_3 - 2.*Box1Abb[534]*m_z1k;

  Box1Abb[665]=-5. + 4.*m_z1k;

  Box1Abb[666]=-2. + Box1Abb[431]*Box1Abb[665]*m_z1k;

  Box1Abb[667]=3. + Box1Abb[662]*m_z1k;

  Box1Abb[668]=2.*Box1Abb[666]*m_z12 + 2.*Box1Abb[667]*m_z12_2 - Box1Abb[119]*m_z12_3 + 4.*m_z1k_2;

  Box1Abb[669]=2. - 10.*m_z1k;

  Box1Abb[670]=-1. + m_z1k - 2.*m_z1k_2;

  Box1Abb[671]=12.*Box1Abb[670] + Box1Abb[669]*m_z12 + m_z12_2;

  Box1Abb[672]=8. + Box1Abb[671]*m_z12 - 8.*m_z1k;

  Box1Abb[673]=Box1Abb[668]*m_x_2 + Box1Abb[672]*m_x_3 + 2.*Box1Abb[661]*m_x_4 - 4.*m_x_5*m_z12 + Box1Abb[664]*m_x*m_z12*m_z1k + Box1Abb[652]*m_z12_3*m_z1k_2;

  Box1Abb[674]=Box1Abb[673]*m_cL - 2.*Box1Abb[659]*Box1Abb[7]*m_cR*m_x*m_z12;

  Box1Abb[675]=7. - 9.*m_z12 + m_z12_2 - 4.*Box1Abb[12]*m_z1k;

  Box1Abb[676]=7. + Box1Abb[542]*m_z12;

  Box1Abb[677]=1. + 4.*m_z12;

  Box1Abb[678]=2. - Box1Abb[676]*m_z12 - 6.*m_z1k + 2.*Box1Abb[677]*m_z12*m_z1k + 4.*Box1Abb[12]*m_z1k_2;

  Box1Abb[679]=Box1Abb[678]*m_x + 2.*Box1Abb[675]*m_x_2 + 4.*Box1Abb[12]*m_x_3 - 2.*Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[680]=Box1Abb[679]*m_cL + 2.*Box1Abb[0]*Box1Abb[383]*m_cR*m_x;

  Box1Abb[681]=3. - 3.*m_z12 - 2.*m_z1k;

  Box1Abb[682]=3.*Box1Abb[0]*m_z12 - 2.*m_z1k + 6.*Box1Abb[8]*m_z12*m_z1k + 4.*Box1Abb[12]*m_z1k_2;

  Box1Abb[683]=5. + m_z12 + 4.*m_z1k;

  Box1Abb[684]=-5. + Box1Abb[683]*m_z12 + 4.*m_z1k;

  Box1Abb[685]=Box1Abb[682]*m_x - 2.*Box1Abb[684]*m_x_2 + 4.*Box1Abb[12]*m_x_3 + Box1Abb[681]*m_z12_2*m_z1k;

  Box1Abb[686]=Box1Abb[685]*m_cL + 2.*Box1Abb[0]*Box1Abb[384]*m_cR*m_x;

  Box1Abb[687]=2. - 6.*m_z12 - 8.*m_z1k;

  Box1Abb[688]=12. - 7.*m_z12;

  Box1Abb[689]=7. - 5.*m_z12;

  Box1Abb[690]=3. + Box1Abb[158]*m_z12 - 9.*m_z1k + Box1Abb[688]*m_z12*m_z1k + 2.*Box1Abb[689]*m_z1k_2 - 8.*m_z1k_3;

  Box1Abb[691]=5. - 6.*m_z1k;

  Box1Abb[692]=4. + 7.*m_z1k;

  Box1Abb[693]=5. - 2.*Box1Abb[692]*m_z12 + m_z12_2 + 2.*Box1Abb[691]*m_z1k;

  Box1Abb[694]=Box1Abb[690]*m_x - Box1Abb[693]*m_x_2 + Box1Abb[687]*m_x_3 + 2.*m_x_4 + 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[695]=Box1Abb[694]*m_cL - 2.*Box1Abb[0]*Box1Abb[420]*m_cR*m_x;

  Box1Abb[696]=-Box1Abb[431]*m_x + m_x_2 + Box1Abb[72]*m_z1k;

  Box1Abb[697]=-3.*Box1Abb[0]*m_z12 + 2.*Box1Abb[68]*m_z1k;

  Box1Abb[698]=-5. + 7.*m_z12 - 4.*m_z1k + 8.*m_z12*m_z1k + 12.*m_z1k_2;

  Box1Abb[699]=-6. + 5.*m_z12;

  Box1Abb[700]=-3. + 3.*m_z12 + 5.*m_z1k + Box1Abb[699]*m_z12*m_z1k + 4.*Box1Abb[79]*m_z1k_2 + 8.*m_z1k_3;

  Box1Abb[701]=-Box1Abb[700]*m_x + Box1Abb[698]*m_x_2 - 4.*Box1Abb[63]*m_x_3 + 2.*m_x_4 + Box1Abb[68]*Box1Abb[697]*m_z1k;

  Box1Abb[702]=Box1Abb[701]*m_cL - 2.*Box1Abb[0]*Box1Abb[696]*m_cR*m_x;

  Box1Abb[703]=-pow(Box1Abb[0],2.)*m_x + Box1Abb[129]*m_x_2 - 2.*Box1Abb[79]*m_x*m_z1k - m_z12*m_z1k_2;

  Box1Abb[704]=4. + m_z12 - 8.*m_z1k;

  Box1Abb[705]=-2. + Box1Abb[704]*m_z12;

  Box1Abb[706]=-2. + 2.*m_z12 + 3.*m_z1k;

  Box1Abb[707]=8. + m_z12;

  Box1Abb[708]=12. + m_z12;

  Box1Abb[709]=-36. + 3.*Box1Abb[707]*m_z12 + 2.*Box1Abb[708]*m_z1k - 24.*m_z1k_2;

  Box1Abb[710]=12. + Box1Abb[709]*m_z12 - 8.*m_z1k;

  Box1Abb[711]=4. - 5.*m_z1k;

  Box1Abb[712]=4. + 3.*Box1Abb[711]*m_z1k;

  Box1Abb[713]=-1. + m_z1k + m_z1k_3;

  Box1Abb[714]=2.*Box1Abb[713]*m_z12 + Box1Abb[712]*m_z12_2 - 2.*Box1Abb[66]*m_z12_3 + 2.*Box1Abb[534]*Box1Abb[68]*m_z1k;

  Box1Abb[715]=11. + 23.*m_z1k;

  Box1Abb[716]=10. + m_z1k_2;

  Box1Abb[717]=-9. + 4.*Box1Abb[318]*m_z1k;

  Box1Abb[718]=4. + Box1Abb[717]*m_z1k;

  Box1Abb[719]=2.*Box1Abb[718]*m_z12 - 2.*Box1Abb[716]*m_z12_2 + Box1Abb[715]*m_z12_3 + m_z12_4 - 4.*Box1Abb[68]*m_z1k;

  Box1Abb[720]=Box1Abb[719]*m_x_2 - Box1Abb[710]*m_x_3 + 2.*Box1Abb[705]*m_x_4 + Box1Abb[714]*m_x*m_z12 + 4.*m_x_5*m_z12 + Box1Abb[4]*Box1Abb[706]*m_z12_3*m_z1k;

  Box1Abb[721]=-Box1Abb[720]*m_cL + 2.*Box1Abb[117]*Box1Abb[703]*m_cR*m_x*m_z12;

  Box1Abb[722]=1. - Box1Abb[178]*m_z12 + m_z1k;

  Box1Abb[723]=-Box1Abb[4]*Box1Abb[68] + Box1Abb[722]*m_x + m_x_2*m_z12;

  Box1Abb[724]=3.*m_z12 + 2.*m_z1k;

  Box1Abb[725]=8. + 9.*m_z1k;

  Box1Abb[726]=3. - Box1Abb[725]*m_z12 + 2.*m_z1k - 8.*m_z1k_2;

  Box1Abb[727]=-7. + 4.*m_z1k;

  Box1Abb[728]=7. + 2.*Box1Abb[727]*m_z1k;

  Box1Abb[729]=1. + 11.*m_z1k;

  Box1Abb[730]=2. + Box1Abb[68]*Box1Abb[729]*m_z1k;

  Box1Abb[731]=-2. + Box1Abb[730]*m_z12 + m_z1k + 4.*Box1Abb[15]*m_z12_2*m_z1k + Box1Abb[728]*m_z1k_2;

  Box1Abb[732]=8. + 13.*m_z1k;

  Box1Abb[733]=7. + Box1Abb[732]*m_z1k;

  Box1Abb[734]=-5. + Box1Abb[733]*m_z12 + 2.*Box1Abb[66]*Box1Abb[68]*m_z1k + 2.*m_z12_2*m_z1k;

  Box1Abb[735]=-Box1Abb[731]*m_x + Box1Abb[734]*m_x_2 + Box1Abb[726]*m_x_3 + Box1Abb[724]*m_x_4 + 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[9]*m_z1k;

  Box1Abb[736]=Box1Abb[735]*m_cL + 2.*Box1Abb[7]*Box1Abb[723]*m_cR*m_x;

  Box1Abb[737]=-1. - Box1Abb[12]*m_x + m_z12 + m_z1k + m_z12*m_z1k;

  Box1Abb[738]=m_z12 - 6.*m_z1k;

  Box1Abb[739]=-2. + 3.*m_z12 + 2.*m_z1k;

  Box1Abb[740]=m_z12 + 4.*m_z1k - 2.*m_z12*m_z1k - 6.*m_z1k_2;

  Box1Abb[741]=Box1Abb[740]*m_x - Box1Abb[738]*m_x_2 - 2.*m_x_3 + Box1Abb[68]*Box1Abb[739]*m_z1k;

  Box1Abb[742]=Box1Abb[741]*m_cL + 2.*Box1Abb[737]*m_cR*m_x;

  Box1Abb[743]=2. + 7.*m_z12;

  Box1Abb[744]=-2.*pow(Box1Abb[68],2.) + 3.*Box1Abb[15]*m_z12_2 + 5.*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[745]=10. + 2.*m_z12 + 11.*m_z1k;

  Box1Abb[746]=-8. + Box1Abb[745]*m_z12;

  Box1Abb[747]=Box1Abb[744]*m_x - Box1Abb[746]*m_x_2 + Box1Abb[743]*m_x_3 - Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12;

  Box1Abb[748]=1. + 7.*m_z12 - 6.*m_z1k;

  Box1Abb[749]=20. + 9.*m_z12;

  Box1Abb[750]=43. + 18.*m_z12;

  Box1Abb[751]=-1. + 10.*m_z12;

  Box1Abb[752]=-18.*Box1Abb[15] + Box1Abb[749]*m_z12 + Box1Abb[750]*m_z12*m_z1k + 4.*Box1Abb[751]*m_z1k_2 + 8.*m_z1k_3;

  Box1Abb[753]=33. + 36.*m_z1k;

  Box1Abb[754]=-14. + Box1Abb[753]*m_z12 + 4.*m_z12_2 - 8.*Box1Abb[68]*m_z1k;

  Box1Abb[755]=3. + m_z1k + 6.*Box1Abb[68]*m_z1k_2;

  Box1Abb[756]=7. + 26.*Box1Abb[15]*m_z1k;

  Box1Abb[757]=-3. + 4.*m_z1k;

  Box1Abb[758]=-2. + Box1Abb[757]*m_z1k;

  Box1Abb[759]=2.*Box1Abb[68]*Box1Abb[755] + m_z12 + Box1Abb[756]*m_z12_2 + 7.*Box1Abb[758]*m_z12*m_z1k + 2.*m_z12_3*m_z1k;

  Box1Abb[760]=-2. + 7.*Box1Abb[317]*m_z1k;

  Box1Abb[761]=-11. + 10.*m_z1k;

  Box1Abb[762]=-2. + Box1Abb[761]*m_z1k;

  Box1Abb[763]=pow(Box1Abb[68],2.)*Box1Abb[762]*m_z12 + Box1Abb[68]*Box1Abb[760]*m_z12_2 + 2.*Box1Abb[317]*pow(Box1Abb[68],3.)*m_z1k + 4.*Box1Abb[15]*m_z12_3*m_z1k;

  Box1Abb[764]=Box1Abb[763]*m_x - Box1Abb[759]*m_x_2 + Box1Abb[752]*m_x_3 - Box1Abb[754]*m_x_4 + 2.*Box1Abb[748]*m_x_5 + 4.*m_x_6 - 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k;

  Box1Abb[765]=Box1Abb[764]*m_cL + 2.*Box1Abb[7]*Box1Abb[747]*m_cR*m_x;

  Box1Abb[766]=-1. + 2.*m_x;

  Box1Abb[767]=2. + 2.*m_x - m_z12;

  Box1Abb[768]=-2. + m_x + 4.*m_x_2;

  Box1Abb[769]=-1. + 9.*m_x;

  Box1Abb[770]=1. - 2.*m_x;

  Box1Abb[771]=2. - 2.*Box1Abb[768]*m_x - 5.*m_z12 + 2.*Box1Abb[769]*m_x*m_z12 + 2.*Box1Abb[770]*m_z12_2;

  Box1Abb[772]=7. + 4.*m_z12;

  Box1Abb[773]=3. + 7.*m_z12;

  Box1Abb[774]=-8. - 2.*Box1Abb[773]*m_x + Box1Abb[772]*m_z12;

  Box1Abb[775]=5. + 4.*m_x - m_z12;

  Box1Abb[776]=Box1Abb[1]*Box1Abb[766]*Box1Abb[767]*m_x + Box1Abb[771]*m_z1k + Box1Abb[774]*m_z1k_2 + 2.*Box1Abb[775]*m_z1k_3 - 4.*m_z1k_4;

  Box1Abb[777]=2. + m_z1k - 5.*m_z1k_2;

  Box1Abb[778]=-2.*pow(Box1Abb[68],2.) + Box1Abb[777]*m_z12 - m_z12_2*m_z1k;

  Box1Abb[779]=Box1Abb[778]*m_x + Box1Abb[77]*m_x_3 + 3.*Box1Abb[68]*m_x_2*m_z12 + Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[780]=Box1Abb[776]*m_cL + 2.*Box1Abb[779]*m_cR;

  Box1Abb[781]=1. - 2.*m_z12_2 + 6.*m_z12*m_z1k;

  Box1Abb[782]=4. - 15.*m_z1k;

  Box1Abb[783]=-5. + Box1Abb[514]*m_z12 + m_z12_2 + Box1Abb[782]*m_z1k;

  Box1Abb[784]=-2. + Box1Abb[783]*m_z12 - 4.*m_z1k;

  Box1Abb[785]=-3. + 5.*m_z1k;

  Box1Abb[786]=1. + 2.*m_z1k_2;

  Box1Abb[787]=3. - 4.*m_z1k + 6.*m_z1k_2;

  Box1Abb[788]=-3. + Box1Abb[787]*m_z1k;

  Box1Abb[789]=pow(Box1Abb[68],4.) + pow(Box1Abb[68],2.)*Box1Abb[788]*m_z12 + Box1Abb[68]*Box1Abb[785]*Box1Abb[786]*m_z12_2 + Box1Abb[457]*Box1Abb[458]*m_z12_3;

  Box1Abb[790]=3. - 4.*m_z1k;

  Box1Abb[791]=15. + 4.*Box1Abb[267]*m_z1k;

  Box1Abb[792]=5. + Box1Abb[791]*m_z1k;

  Box1Abb[793]=-4. + Box1Abb[792]*m_z12 + Box1Abb[790]*m_z12_2 - 3.*Box1Abb[15]*m_z12_3 + 6.*m_z1k_2;

  Box1Abb[794]=3. + m_z1k_2;

  Box1Abb[795]=5. + 2.*Box1Abb[267]*m_z1k;

  Box1Abb[796]=7. + Box1Abb[795]*m_z1k;

  Box1Abb[797]=-8. + 5.*m_z1k;

  Box1Abb[798]=7. + Box1Abb[797]*m_z1k;

  Box1Abb[799]=4. - 3.*Box1Abb[798]*m_z1k;

  Box1Abb[800]=4. + Box1Abb[799]*m_z1k;

  Box1Abb[801]=Box1Abb[800]*m_z12 - Box1Abb[796]*m_z12_2 + Box1Abb[794]*m_z12_3 - 2.*Box1Abb[317]*Box1Abb[68]*m_z1k;

  Box1Abb[802]=Box1Abb[789]*m_x + Box1Abb[801]*m_x_2 + Box1Abb[793]*m_x_3 + Box1Abb[784]*m_x_4 + Box1Abb[781]*m_x_5 - m_x_6*m_z12 - Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[72]*m_z12*m_z1k_2;

  Box1Abb[803]=-10. + m_z12;

  Box1Abb[804]=-2. + Box1Abb[803]*m_z12 - 24.*m_z1k;

  Box1Abb[805]=22. + 9.*m_z12;

  Box1Abb[806]=40. - Box1Abb[805]*m_z12;

  Box1Abb[807]=17. + m_z12;

  Box1Abb[808]=-14. + Box1Abb[807]*m_z12;

  Box1Abb[809]=7. + Box1Abb[808]*m_z12;

  Box1Abb[810]=2. - 5.*m_z12;

  Box1Abb[811]=-8. + 3.*Box1Abb[810]*m_z12;

  Box1Abb[812]=30. + Box1Abb[811]*m_z12;

  Box1Abb[813]=-13. + 11.*m_z12;

  Box1Abb[814]=25. + Box1Abb[813]*m_z12;

  Box1Abb[815]=-12.*Box1Abb[15] + Box1Abb[806]*m_z12 - 2.*Box1Abb[809]*m_z12*m_z1k + 2.*Box1Abb[812]*m_z1k_2 - 4.*Box1Abb[814]*m_z1k_3 + 60.*m_z1k_4;

  Box1Abb[816]=-2. + 2.*Box1Abb[15]*m_z12 + m_z1k + m_z1k_2;

  Box1Abb[817]=1. - 6.*m_z1k;

  Box1Abb[818]=7. + 5.*m_z1k;

  Box1Abb[819]=10. + 21.*m_z1k;

  Box1Abb[820]=12. - 2.*Box1Abb[819]*m_z12 + 2.*Box1Abb[818]*m_z12_2 + m_z12_3 + 10.*Box1Abb[817]*m_z1k;

  Box1Abb[821]=1. + m_z1k + m_z1k_2;

  Box1Abb[822]=2. + m_z1k - 17.*m_z1k_2 + 18.*m_z1k_3 + m_z1k_4 - 5.*m_z1k_5;

  Box1Abb[823]=3. + Box1Abb[119]*m_z1k;

  Box1Abb[824]=2. + 17.*m_z1k;

  Box1Abb[825]=-2. + Box1Abb[824]*m_z1k;

  Box1Abb[826]=2.*pow(Box1Abb[68],3.)*Box1Abb[821]*m_z12 + 2.*Box1Abb[822]*m_z12_2 - Box1Abb[15]*Box1Abb[68]*Box1Abb[825]*m_z12_3 + 2.*Box1Abb[317]*pow(Box1Abb[68],4.)*m_z1k - 2.*Box1Abb[823]*m_z12_4*m_z1k;

  Box1Abb[827]=5. + 11.*m_z1k;

  Box1Abb[828]=1. + 5.*Box1Abb[790]*m_z1k;

  Box1Abb[829]=6. + 17.*m_z1k;

  Box1Abb[830]=17. + 2.*Box1Abb[829]*m_z1k;

  Box1Abb[831]=14. + 31.*m_z1k;

  Box1Abb[832]=32. + Box1Abb[831]*m_z1k;

  Box1Abb[833]=10. - 2.*Box1Abb[830]*m_z12 + Box1Abb[832]*m_z12_2 + Box1Abb[827]*m_z12_3 + 4.*Box1Abb[828]*m_z1k;

  Box1Abb[834]=-1. + 4.*m_z1k;

  Box1Abb[835]=2. + Box1Abb[834]*m_z1k;

  Box1Abb[836]=22. - 3.*m_z1k + 9.*m_z1k_2;

  Box1Abb[837]=-7. + Box1Abb[836]*m_z1k;

  Box1Abb[838]=25. + 4.*m_z1k + 34.*m_z1k_2;

  Box1Abb[839]=7. + Box1Abb[838]*m_z1k;

  Box1Abb[840]=-12. + 31.*m_z1k;

  Box1Abb[841]=2. + Box1Abb[840]*m_z1k;

  Box1Abb[842]=-48. + Box1Abb[841]*m_z1k;

  Box1Abb[843]=-1. + Box1Abb[842]*m_z1k;

  Box1Abb[844]=-2.*Box1Abb[536]*pow(Box1Abb[68],2.)*Box1Abb[835] - 2.*Box1Abb[68]*Box1Abb[837]*m_z12 + Box1Abb[843]*m_z12_2 + Box1Abb[839]*m_z12_3 + 6.*Box1Abb[15]*m_z12_4*m_z1k;

  Box1Abb[845]=Box1Abb[826]*m_x + Box1Abb[844]*m_x_2 + Box1Abb[815]*m_x_3 + Box1Abb[833]*m_x_4 - Box1Abb[820]*m_x_5 + Box1Abb[804]*m_x_6 + 4.*m_x_7 + Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[816]*m_z12_2*m_z1k;

  Box1Abb[846]=Box1Abb[845]*m_cL + 2.*Box1Abb[802]*m_cR*m_x;

  Box1Abb[847]=-6. - 6.*m_z12 + m_z12_2 - 24.*m_z1k;

  Box1Abb[848]=-31. + 6.*m_z12;

  Box1Abb[849]=36. + Box1Abb[848]*m_z12;

  Box1Abb[850]=16. + Box1Abb[129]*m_z12;

  Box1Abb[851]=-8. + Box1Abb[77]*Box1Abb[8]*m_z12;

  Box1Abb[852]=2. + 31.*m_z12;

  Box1Abb[853]=50. + Box1Abb[852]*m_z12;

  Box1Abb[854]=4. + 2.*Box1Abb[493]*m_z12 - 14.*m_z1k + Box1Abb[849]*m_z12*m_z1k - 2.*Box1Abb[0]*Box1Abb[850]*m_z1k_2 + 6.*Box1Abb[851]*m_z1k_3 + Box1Abb[853]*m_z1k_4 - 24.*m_z1k_5;

  Box1Abb[855]=-20. + 7.*m_z12;

  Box1Abb[856]=-5. + m_z12;

  Box1Abb[857]=2. + Box1Abb[856]*m_z12;

  Box1Abb[858]=-5. + 4.*m_z12;

  Box1Abb[859]=-18. + Box1Abb[77]*Box1Abb[858]*m_z12;

  Box1Abb[860]=-3. + 11.*m_z12;

  Box1Abb[861]=15. + Box1Abb[860]*m_z12;

  Box1Abb[862]=6.*Box1Abb[178] + Box1Abb[855]*m_z12 + 6.*Box1Abb[857]*m_z12*m_z1k + 2.*Box1Abb[859]*m_z1k_2 + 4.*Box1Abb[861]*m_z1k_3 - 60.*m_z1k_4;

  Box1Abb[863]=20. + 22.*m_z1k;

  Box1Abb[864]=-6. + Box1Abb[863]*m_z12 - 5.*Box1Abb[431]*m_z12_2 + 10.*Box1Abb[66]*m_z1k;

  Box1Abb[865]=1. - 4.*m_z1k;

  Box1Abb[866]=9. + 7.*m_z1k;

  Box1Abb[867]=7. + Box1Abb[866]*m_z1k;

  Box1Abb[868]=3. + 31.*m_z1k;

  Box1Abb[869]=9. + Box1Abb[868]*m_z1k;

  Box1Abb[870]=16. - 4.*Box1Abb[867]*m_z12 + Box1Abb[869]*m_z12_2 + 2.*m_z12_3*m_z1k + 20.*Box1Abb[865]*m_z1k_2;

  Box1Abb[871]=3. + Box1Abb[170]*m_z1k;

  Box1Abb[872]=-1. + Box1Abb[178]*m_z1k;

  Box1Abb[873]=3. + 2.*m_z1k;

  Box1Abb[874]=-19. + 5.*Box1Abb[873]*m_z1k;

  Box1Abb[875]=8. + Box1Abb[874]*m_z1k;

  Box1Abb[876]=2.*pow(Box1Abb[68],2.)*Box1Abb[871]*m_z12 + Box1Abb[68]*Box1Abb[875]*m_z12_2 + 2.*Box1Abb[834]*Box1Abb[872]*m_z12_3 - 2.*Box1Abb[317]*pow(Box1Abb[68],3.)*m_z1k;

  Box1Abb[877]=Box1Abb[854]*m_x_2 - Box1Abb[862]*m_x_3 + Box1Abb[870]*m_x_4 + Box1Abb[864]*m_x_5 + Box1Abb[847]*m_x_6 + 4.*m_x_7 - Box1Abb[876]*m_x*m_z1k + Box1Abb[652]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_2;

  Box1Abb[878]=-1. + Box1Abb[263]*m_z12;

  Box1Abb[879]=1. + Box1Abb[133]*m_z12;

  Box1Abb[880]=-1. + 3.*Box1Abb[0]*m_z12 + m_z12*m_z1k - 3.*Box1Abb[879]*m_z1k_2 - 10.*m_z12*m_z1k_3;

  Box1Abb[881]=-3. + m_z1k;

  Box1Abb[882]=3.*Box1Abb[15] + m_z12 + Box1Abb[881]*m_z12_2 + 10.*m_z12*m_z1k_2;

  Box1Abb[883]=2. + m_z1k_2 - 8.*m_z1k_3 + 5.*m_z1k_4;

  Box1Abb[884]=-5. + 11.*m_z1k;

  Box1Abb[885]=-1. + Box1Abb[884]*m_z1k;

  Box1Abb[886]=-1. + Box1Abb[885]*m_z1k;

  Box1Abb[887]=pow(Box1Abb[68],3.) + Box1Abb[883]*m_z12 + Box1Abb[886]*m_z12_2 + 2.*m_z12_3*m_z1k_2;

  Box1Abb[888]=Box1Abb[887]*m_x + Box1Abb[880]*m_x_2 + Box1Abb[882]*m_x_3 + Box1Abb[878]*m_x_4 + m_x_5*m_z12 - Box1Abb[541]*Box1Abb[68]*m_z12*m_z1k_2;

  Box1Abb[889]=Box1Abb[877]*m_cL - 2.*Box1Abb[7]*Box1Abb[888]*m_cR*m_x;

  Box1Abb[912]=m_s;

  Box1Abb[913]=-2./m_s;

  Box1Abb[914]=2./m_s;

  Box1Abb[915]=(2.*m_z1k)/m_s;

  Box1Abb[916]=(-2.*m_x)/m_s;

  Box1Abb[917]=Box1Abb[656]*m_cR + 2.*Box1Abb[651]*m_cL*m_x;

  Box1Abb[918]=Box1Abb[673]*m_cR - 2.*Box1Abb[659]*Box1Abb[7]*m_cL*m_x*m_z12;

  Box1Abb[919]=Box1Abb[679]*m_cR + 2.*Box1Abb[0]*Box1Abb[383]*m_cL*m_x;

  Box1Abb[920]=Box1Abb[685]*m_cR + 2.*Box1Abb[0]*Box1Abb[384]*m_cL*m_x;

  Box1Abb[921]=Box1Abb[694]*m_cR - 2.*Box1Abb[0]*Box1Abb[420]*m_cL*m_x;

  Box1Abb[922]=Box1Abb[701]*m_cR - 2.*Box1Abb[0]*Box1Abb[696]*m_cL*m_x;

  Box1Abb[923]=-Box1Abb[720]*m_cR + 2.*Box1Abb[117]*Box1Abb[703]*m_cL*m_x*m_z12;

  Box1Abb[924]=Box1Abb[735]*m_cR + 2.*Box1Abb[7]*Box1Abb[723]*m_cL*m_x;

  Box1Abb[925]=Box1Abb[741]*m_cR + 2.*Box1Abb[737]*m_cL*m_x;

  Box1Abb[926]=Box1Abb[764]*m_cR + 2.*Box1Abb[7]*Box1Abb[747]*m_cL*m_x;

  Box1Abb[927]=2.*Box1Abb[779]*m_cL + Box1Abb[776]*m_cR;

  Box1Abb[928]=2. - Box1Abb[803]*m_z12 + 24.*m_z1k;

  Box1Abb[929]=-40. + Box1Abb[805]*m_z12;

  Box1Abb[930]=8. + 3.*Box1Abb[505]*m_z12;

  Box1Abb[931]=-30. + Box1Abb[930]*m_z12;

  Box1Abb[932]=12.*Box1Abb[15] + Box1Abb[929]*m_z12 + 2.*Box1Abb[809]*m_z12*m_z1k + 2.*Box1Abb[931]*m_z1k_2 + 4.*Box1Abb[814]*m_z1k_3 - 60.*m_z1k_4;

  Box1Abb[933]=-5. + Box1Abb[240]*m_z1k;

  Box1Abb[934]=-2. + Box1Abb[933]*m_z1k;

  Box1Abb[935]=-2.*pow(Box1Abb[68],3.)*Box1Abb[821]*m_z12 + 2.*pow(Box1Abb[68],2.)*Box1Abb[934]*m_z12_2 + Box1Abb[15]*Box1Abb[68]*Box1Abb[825]*m_z12_3 - 2.*Box1Abb[317]*pow(Box1Abb[68],4.)*m_z1k + 2.*Box1Abb[823]*m_z12_4*m_z1k;

  Box1Abb[936]=12. - 31.*m_z1k;

  Box1Abb[937]=-2. + Box1Abb[936]*m_z1k;

  Box1Abb[938]=48. + Box1Abb[937]*m_z1k;

  Box1Abb[939]=1. + Box1Abb[938]*m_z1k;

  Box1Abb[940]=2.*Box1Abb[536]*pow(Box1Abb[68],2.)*Box1Abb[835] + 2.*Box1Abb[68]*Box1Abb[837]*m_z12 + Box1Abb[939]*m_z12_2 - Box1Abb[839]*m_z12_3 - 6.*Box1Abb[15]*m_z12_4*m_z1k;

  Box1Abb[941]=Box1Abb[935]*m_x + Box1Abb[940]*m_x_2 + Box1Abb[932]*m_x_3 - Box1Abb[833]*m_x_4 + Box1Abb[820]*m_x_5 + Box1Abb[928]*m_x_6 - 4.*m_x_7 - Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[816]*m_z12_2*m_z1k;

  Box1Abb[942]=m_z12 - 3.*m_z1k;

  Box1Abb[943]=-1. + 2.*Box1Abb[942]*m_z12;

  Box1Abb[944]=-5. + 3.*Box1Abb[0]*m_z12;

  Box1Abb[945]=-3. + 8.*m_z12;

  Box1Abb[946]=4. + Box1Abb[944]*m_z12 + Box1Abb[149]*Box1Abb[251]*m_z12*m_z1k + 2.*Box1Abb[945]*m_z1k_2 - 20.*m_z12*m_z1k_3;

  Box1Abb[947]=-pow(Box1Abb[68],4.) - pow(Box1Abb[68],2.)*Box1Abb[788]*m_z12 - Box1Abb[68]*Box1Abb[785]*Box1Abb[786]*m_z12_2 - Box1Abb[457]*Box1Abb[458]*m_z12_3;

  Box1Abb[948]=-4. + 3.*Box1Abb[798]*m_z1k;

  Box1Abb[949]=-4. + Box1Abb[948]*m_z1k;

  Box1Abb[950]=Box1Abb[949]*m_z12 + Box1Abb[796]*m_z12_2 - Box1Abb[794]*m_z12_3 + 2.*Box1Abb[317]*Box1Abb[68]*m_z1k;

  Box1Abb[951]=Box1Abb[947]*m_x + Box1Abb[950]*m_x_2 + Box1Abb[946]*m_x_3 - Box1Abb[784]*m_x_4 + Box1Abb[943]*m_x_5 + m_x_6*m_z12 + Box1Abb[4]*pow(Box1Abb[68],2.)*Box1Abb[72]*m_z12*m_z1k_2;

  Box1Abb[952]=-Box1Abb[941]*m_cR - 2.*Box1Abb[951]*m_cL*m_x;

  Box1Abb[953]=Box1Abb[877]*m_cR - 2.*Box1Abb[7]*Box1Abb[888]*m_cL*m_x;

  Box1Abb[976]=2. - 3.*m_z12 - 2.*m_z1k;

  Box1Abb[977]=m_z12 + 4.*m_z1k;

  Box1Abb[978]=Box1Abb[977]*m_x - 2.*m_x_2 + Box1Abb[976]*m_z1k;

  Box1Abb[979]=-1. + m_z12 + 4.*m_z1k;

  Box1Abb[980]=2. - 2.*Box1Abb[979]*m_x + 4.*m_x_2 - 3.*m_z12 + m_z12_2 + 6.*Box1Abb[0]*m_z1k + 4.*m_z1k_2;

  Box1Abb[981]=1. + m_z12 + 4.*m_z1k;

  Box1Abb[982]=m_z12 - 2.*m_z1k + 6.*m_z12*m_z1k + 4.*m_z1k_2;

  Box1Abb[983]=Box1Abb[982]*m_x - 2.*Box1Abb[981]*m_x_2 + 4.*m_x_3 - m_z12_2*m_z1k;

  Box1Abb[984]=2. + Box1Abb[158]*m_z12 - 4.*m_z1k;

  Box1Abb[985]=Box1Abb[984]*m_x + 4.*m_x_2 + m_z12_2*m_z1k;

  Box1Abb[986]=2.*m_x - m_z12;

  Box1Abb[987]=-2. + Box1Abb[767]*m_z12;

  Box1Abb[988]=4. - 8.*m_z12;

  Box1Abb[989]=Box1Abb[988]*m_x + 8.*m_x_2 + m_z12_2;

  Box1Abb[990]=Box1Abb[986]*Box1Abb[987]*m_x - Box1Abb[989]*m_z12*m_z1k + 4.*m_x*m_z12*m_z1k_2;

  Box1Abb[991]=-1. + m_z12 + m_z12_2 - 3.*m_z12*m_z1k;

  Box1Abb[992]=-1. + Box1Abb[326]*m_z12;

  Box1Abb[993]=Box1Abb[0]*m_z12 + Box1Abb[992]*m_z1k + 2.*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[994]=-10. + 10.*m_z12 + m_z12_2 + 4.*Box1Abb[77]*m_z1k - 12.*m_z1k_2;

  Box1Abb[995]=Box1Abb[994]*m_z12 - 4.*m_z1k;

  Box1Abb[996]=-Box1Abb[995]*m_x_2 + 4.*Box1Abb[991]*m_x_3 + 2.*Box1Abb[993]*m_x*m_z12 + 4.*m_x_4*m_z12 - Box1Abb[652]*m_z12_3*m_z1k;

  Box1Abb[997]=6. + Box1Abb[856]*m_z12 - 6.*m_z1k;

  Box1Abb[998]=-2. + Box1Abb[997]*m_z12;

  Box1Abb[999]=-1. + m_z12 + Box1Abb[0]*Box1Abb[79]*m_z1k - 2.*m_z1k_2;

  Box1Abb[1000]=Box1Abb[0]*pow(Box1Abb[79],2.) + 4.*pow(Box1Abb[79],2.)*m_z1k - 12.*m_z1k_2;

  Box1Abb[1001]=Box1Abb[1000]*m_z12 - 4.*m_z1k;

  Box1Abb[1002]=Box1Abb[1001]*m_x_2 - 2.*Box1Abb[998]*m_x_3 - 4.*m_x_4*m_z12 - 2.*Box1Abb[999]*m_x*m_z12*m_z1k + Box1Abb[0]*m_z12_3*m_z1k_2;

  Box1Abb[1003]=-3. + m_z12 + m_z1k;

  Box1Abb[1004]=1. + Box1Abb[1003]*m_z12 - 2.*m_z1k;

  Box1Abb[1005]=-Box1Abb[79]*Box1Abb[977]*m_x_2 + 2.*Box1Abb[79]*m_x_3 + 2.*Box1Abb[1004]*m_x*m_z1k + m_z12_2*m_z1k_2;

  Box1Abb[1006]=2.*m_x + m_z12;

  Box1Abb[1007]=2. + Box1Abb[1006]*m_x - 2.*m_z12;

  Box1Abb[1008]=2. + m_z12_2;

  Box1Abb[1009]=Box1Abb[1008]*m_x - 6.*Box1Abb[79]*m_x_2 + 2.*Box1Abb[0]*Box1Abb[493]*m_z12;

  Box1Abb[1010]=2. - 3.*m_z12;

  Box1Abb[1011]=-2. + Box1Abb[1010]*m_z12;

  Box1Abb[1012]=Box1Abb[1011]*m_x + 6.*Box1Abb[79]*m_x_2 + 2.*Box1Abb[0]*m_z12_2;

  Box1Abb[1013]=-2.*Box1Abb[79]*m_x + m_z12_2;

  Box1Abb[1014]=Box1Abb[1007]*Box1Abb[79]*m_x_2 + Box1Abb[1009]*m_x*m_z1k + Box1Abb[1012]*m_z1k_2 + Box1Abb[1013]*m_z1k_3;

  Box1Abb[1015]=2. + m_z12 + 2.*m_z12_2 + 8.*m_z12*m_z1k;

  Box1Abb[1016]=5. - 2.*m_z12;

  Box1Abb[1017]=-8. + Box1Abb[1016]*m_z12;

  Box1Abb[1018]=13. - 10.*m_z12;

  Box1Abb[1019]=-2. + Box1Abb[1018]*m_z12;

  Box1Abb[1020]=1. + Box1Abb[493]*m_z12 + m_z1k + Box1Abb[1017]*m_z12*m_z1k + Box1Abb[1019]*m_z1k_2 - 8.*m_z12*m_z1k_3;

  Box1Abb[1021]=Box1Abb[493]*m_z12 + 2.*m_z12*m_z1k + m_z1k_2;

  Box1Abb[1022]=3. + 2.*Box1Abb[1021] - 5.*m_z1k;

  Box1Abb[1023]=-5. + 12.*m_z1k;

  Box1Abb[1024]=m_z12 + Box1Abb[1023]*m_z1k + 8.*m_z12*m_z1k;

  Box1Abb[1025]=3. + Box1Abb[1024]*m_z12 + 4.*m_z1k;

  Box1Abb[1026]=Box1Abb[1020]*m_x + Box1Abb[1025]*m_x_2 - Box1Abb[1015]*m_x_3 + 2.*m_x_4*m_z12 + Box1Abb[1022]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[1027]=1. - m_z12 + m_z1k + Box1Abb[69]*m_z12*m_z1k + 10.*Box1Abb[0]*m_z1k_2 + 8.*m_z1k_3;

  Box1Abb[1028]=-1. + 6.*m_z1k;

  Box1Abb[1029]=-1. + m_z12 + 2.*Box1Abb[1028]*m_z1k + 8.*m_z12*m_z1k;

  Box1Abb[1030]=-Box1Abb[1027]*m_x + Box1Abb[1029]*m_x_2 - 2.*Box1Abb[981]*m_x_3 + 2.*m_x_4 + 2.*pow(Box1Abb[4],2.)*Box1Abb[68]*m_z1k;

  Box1Abb[1031]=2. + m_z12 + 4.*m_z1k;

  Box1Abb[1032]=5. + 2.*m_z12;

  Box1Abb[1033]=1. - 2.*m_z12 + m_z1k - Box1Abb[1032]*m_z12*m_z1k + 2.*Box1Abb[810]*m_z1k_2 - 8.*m_z1k_3;

  Box1Abb[1034]=Box1Abb[174]*Box1Abb[68]*m_z12 + 2.*Box1Abb[15]*m_z12_2 + 2.*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[1035]=4. + 8.*m_z1k;

  Box1Abb[1036]=1. + Box1Abb[1035]*m_z12 + 4.*Box1Abb[265]*m_z1k;

  Box1Abb[1037]=Box1Abb[1033]*m_x + Box1Abb[1036]*m_x_2 - 2.*Box1Abb[1031]*m_x_3 + 2.*m_x_4 + Box1Abb[1034]*m_z1k;

  Box1Abb[1038]=2. + m_z12_2 + 8.*m_z12*m_z1k;

  Box1Abb[1039]=-2. + 5.*m_z1k;

  Box1Abb[1040]=2. + 2.*Box1Abb[1039]*m_z12 + 5.*m_z12_2 + 4.*Box1Abb[68]*m_z1k;

  Box1Abb[1041]=-2. + 7.*m_z1k;

  Box1Abb[1042]=8. + 2.*Box1Abb[1041]*m_z12 + m_z12_2 + 4.*Box1Abb[1028]*m_z1k;

  Box1Abb[1043]=-4. + Box1Abb[1042]*m_z12 + 8.*m_z1k;

  Box1Abb[1044]=4. + 22.*m_z1k;

  Box1Abb[1045]=-2. + Box1Abb[1044]*m_z12 + m_z12_2 + 8.*Box1Abb[317]*m_z1k;

  Box1Abb[1046]=Box1Abb[1045]*m_z12 + 4.*m_z1k;

  Box1Abb[1047]=-Box1Abb[1043]*m_x_3 + 2.*Box1Abb[1038]*m_x_4 - 4.*m_x_5*m_z12 + Box1Abb[1046]*m_x_2*m_z1k - Box1Abb[1040]*m_x*m_z12*m_z1k_2 + m_z12_3*m_z1k_3;

  Box1Abb[1048]=Box1Abb[431]*m_x - m_x_2 - Box1Abb[72]*m_z1k;

  Box1Abb[1049]=-Box1Abb[979]*m_x + m_x_2*m_z12 - Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[1050]=1. - 2.*m_z12 + 6.*m_z1k_2;

  Box1Abb[1051]=6. - 5.*m_z12;

  Box1Abb[1052]=-1. + m_z12 + m_z1k + 2.*Box1Abb[0]*m_z12*m_z1k + Box1Abb[1051]*m_z1k_2 - 6.*m_z1k_3;

  Box1Abb[1053]=Box1Abb[1052]*m_x + Box1Abb[1050]*m_x_2 + Box1Abb[156]*m_x_3 + 2.*pow(Box1Abb[4],2.)*Box1Abb[68]*m_z1k;

  Box1Abb[1054]=2. - 2.*m_x + 3.*m_z12;

  Box1Abb[1055]=2. + Box1Abb[1054]*m_x - 3.*m_z12;

  Box1Abb[1056]=2. - 8.*m_z12;

  Box1Abb[1057]=2. + Box1Abb[1056]*m_x + 6.*m_x_2 + Box1Abb[858]*m_z12;

  Box1Abb[1058]=-4. - 6.*m_x + 5.*m_z12;

  Box1Abb[1059]=Box1Abb[1055]*m_x + Box1Abb[1057]*m_z1k + Box1Abb[1058]*m_z1k_2 + 2.*m_z1k_3;

  Box1Abb[1060]=-2. + m_x;

  Box1Abb[1061]=-2.*Box1Abb[1060]*m_x + Box1Abb[1]*Box1Abb[376]*m_z12;

  Box1Abb[1062]=m_x - 3.*Box1Abb[183]*m_z12 + 2.*m_z12_2;

  Box1Abb[1063]=-2. + 2.*Box1Abb[1062]*m_x + 3.*m_z12;

  Box1Abb[1064]=2. + 6.*m_z12;

  Box1Abb[1065]=4. + Box1Abb[1064]*m_x - Box1Abb[677]*m_z12;

  Box1Abb[1066]=Box1Abb[1061]*m_x + Box1Abb[1063]*m_z1k + Box1Abb[1065]*m_z1k_2 - 2.*Box1Abb[12]*m_z1k_3;

  Box1Abb[1067]=1. + 6.*m_z12 + 4.*m_z1k;

  Box1Abb[1068]=-2.*Box1Abb[317]*pow(Box1Abb[68],2.) + Box1Abb[68]*m_z12 + 4.*m_z12_2*m_z1k;

  Box1Abb[1069]=13. + 2.*m_z12 + 28.*m_z1k;

  Box1Abb[1070]=-2.*Box1Abb[475] + Box1Abb[1069]*m_z12;

  Box1Abb[1071]=1. + 22.*m_z1k;

  Box1Abb[1072]=1. + 6.*Box1Abb[170]*m_z1k;

  Box1Abb[1073]=-1. + m_z1k + 4.*m_z1k_2;

  Box1Abb[1074]=-2.*Box1Abb[1073]*Box1Abb[68] + 2.*Box1Abb[1072]*m_z12 + Box1Abb[1071]*m_z12_2;

  Box1Abb[1075]=Box1Abb[1068]*Box1Abb[4] - Box1Abb[1074]*m_x + Box1Abb[1070]*m_x_2 - 2.*Box1Abb[1067]*m_x_3 + 4.*m_x_4;

  Box1Abb[1076]=3. + 4.*m_z12 + 4.*m_z1k;

  Box1Abb[1077]=15. + 2.*m_z12 + 16.*m_z1k;

  Box1Abb[1078]=-2.*Box1Abb[178] + Box1Abb[1077]*m_z12;

  Box1Abb[1079]=3. + Box1Abb[665]*m_z1k;

  Box1Abb[1080]=-2.*Box1Abb[1079]*Box1Abb[15] + 5.*m_z12 + 4.*Box1Abb[265]*m_z12_2;

  Box1Abb[1081]=2. - 8.*Box1Abb[68]*m_z1k;

  Box1Abb[1082]=-9. + 8.*m_z1k;

  Box1Abb[1083]=-2. + Box1Abb[1082]*m_z1k;

  Box1Abb[1084]=-Box1Abb[1083]*Box1Abb[68]*m_z12 + Box1Abb[1081]*m_z12_2 - 2.*Box1Abb[317]*pow(Box1Abb[68],2.)*m_z1k + 2.*m_z12_3*m_z1k;

  Box1Abb[1085]=Box1Abb[1084]*m_x - Box1Abb[1080]*m_x_2 + Box1Abb[1078]*m_x_3 - 2.*Box1Abb[1076]*m_x_4 + 4.*m_x_5 + 2.*Box1Abb[4]*Box1Abb[68]*m_z12_2*m_z1k;

  Box1Abb[1086]=-1. + 5.*m_z12;

  Box1Abb[1087]=1. + Box1Abb[1086]*m_z12;

  Box1Abb[1088]=4. + 12.*m_z12 + 14.*m_z12_2 - 3.*m_z12_3 + 2.*Box1Abb[140]*Box1Abb[743]*m_z1k + 12.*Box1Abb[1087]*m_z1k_2 + 40.*m_z12*m_z1k_3;

  Box1Abb[1089]=6. + 3.*m_z12 + 10.*m_z1k;

  Box1Abb[1090]=2. + Box1Abb[1089]*m_z12;

  Box1Abb[1091]=5. + 8.*m_z12;

  Box1Abb[1092]=12. + 14.*m_z12 - m_z12_2 + 4.*Box1Abb[1091]*m_z1k + 40.*m_z1k_2;

  Box1Abb[1093]=8. + Box1Abb[1092]*m_z12 + 12.*m_z1k;

  Box1Abb[1094]=2. + Box1Abb[400]*m_z1k;

  Box1Abb[1095]=-5. + 18.*m_z1k;

  Box1Abb[1096]=1. + Box1Abb[1095]*m_z1k;

  Box1Abb[1097]=-2.*Box1Abb[534]*pow(Box1Abb[68],3.) - 2.*Box1Abb[1094]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[1096]*Box1Abb[68]*m_z12_2 - 8.*m_z12_3*m_z1k_2;

  Box1Abb[1098]=19. + 34.*m_z1k;

  Box1Abb[1099]=-3. + Box1Abb[1098]*m_z1k;

  Box1Abb[1100]=-13. + 24.*m_z1k;

  Box1Abb[1101]=-14. + Box1Abb[1100]*m_z1k;

  Box1Abb[1102]=5. + Box1Abb[1101]*m_z1k;

  Box1Abb[1103]=-9. + 5.*m_z1k;

  Box1Abb[1104]=1. + Box1Abb[1103]*m_z1k;

  Box1Abb[1105]=3. + 2.*Box1Abb[1104]*m_z1k;

  Box1Abb[1106]=6. + 2.*Box1Abb[1105]*m_z1k;

  Box1Abb[1107]=4.*pow(Box1Abb[68],3.) + Box1Abb[1106]*m_z12 + 2.*Box1Abb[1102]*m_z12_2 + Box1Abb[1099]*m_z12_3;

  Box1Abb[1108]=Box1Abb[1107]*m_x_2 - Box1Abb[1088]*m_x_3 + Box1Abb[1093]*m_x_4 - 2.*Box1Abb[1090]*m_x_5 + Box1Abb[1097]*m_x*m_z12 + 4.*m_x_6*m_z12 - pow(Box1Abb[68],3.)*m_z12_3*m_z1k;

  Box1Abb[1109]=2. - 5.*m_z12 - 10.*m_z1k;

  Box1Abb[1110]=8. + Box1Abb[1109]*m_z12 + 20.*m_z1k;

  Box1Abb[1111]=20. + Box1Abb[803]*m_z12;

  Box1Abb[1112]=-8. + 9.*m_z12;

  Box1Abb[1113]=-2. + Box1Abb[1112]*m_z12;

  Box1Abb[1114]=-Box1Abb[1111]*m_z12 + 3.*Box1Abb[1113]*m_z1k + 20.*Box1Abb[79]*m_z1k_2;

  Box1Abb[1115]=2. + Box1Abb[396]*m_z1k;

  Box1Abb[1116]=19. - 20.*m_z1k;

  Box1Abb[1117]=1. + Box1Abb[1116]*m_z1k_2;

  Box1Abb[1118]=-8. + 15.*m_z1k;

  Box1Abb[1119]=3. + Box1Abb[1118]*m_z1k;

  Box1Abb[1120]=2.*Box1Abb[317]*pow(Box1Abb[68],4.) - 2.*Box1Abb[1115]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[1119]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[1117]*m_z12_3 - 8.*m_z12_4*m_z1k_2;

  Box1Abb[1121]=3. - 12.*m_z1k;

  Box1Abb[1122]=4. + m_z1k + 26.*m_z1k_2;

  Box1Abb[1123]=-2. + 5.*Box1Abb[757]*m_z1k;

  Box1Abb[1124]=19. - 10.*Box1Abb[881]*m_z1k;

  Box1Abb[1125]=20. + 2.*Box1Abb[1124]*m_z1k;

  Box1Abb[1126]=-10. + Box1Abb[1125]*m_z12 - 2.*Box1Abb[1122]*m_z12_2 + Box1Abb[1121]*m_z12_3 + 2.*Box1Abb[1123]*m_z1k;

  Box1Abb[1127]=-3. + 10.*m_z1k;

  Box1Abb[1128]=-2. + Box1Abb[1127]*m_z1k;

  Box1Abb[1129]=13. + 32.*m_z1k;

  Box1Abb[1130]=-3. + Box1Abb[1129]*m_z1k;

  Box1Abb[1131]=-23. + 5.*m_z1k;

  Box1Abb[1132]=-9. + Box1Abb[1131]*m_z1k;

  Box1Abb[1133]=4. + Box1Abb[1132]*m_z1k;

  Box1Abb[1134]=-21. + 22.*m_z1k;

  Box1Abb[1135]=-22. + Box1Abb[1134]*m_z1k;

  Box1Abb[1136]=6. + 2.*Box1Abb[1135]*m_z1k;

  Box1Abb[1137]=-2.*Box1Abb[1128]*pow(Box1Abb[68],2.) + 2.*Box1Abb[1133]*Box1Abb[68]*m_z12 + Box1Abb[1136]*m_z12_2 + Box1Abb[1130]*m_z12_3;

  Box1Abb[1138]=Box1Abb[1120]*m_x + Box1Abb[1137]*m_x_2 + Box1Abb[1126]*m_x_3 + Box1Abb[1114]*m_x_4 + Box1Abb[1110]*m_x_5 + 2.*Box1Abb[79]*m_x_6 + Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k;

  Box1Abb[1139]=-2. + 2.*m_z12 + m_z1k + m_z1k_2;

  Box1Abb[1140]=9. + 23.*m_z12;

  Box1Abb[1141]=-21. + Box1Abb[1140]*m_z12;

  Box1Abb[1142]=3. + Box1Abb[1141]*m_z12;

  Box1Abb[1143]=5. + 14.*m_z12;

  Box1Abb[1144]=-3. + 2.*Box1Abb[12]*m_z12;

  Box1Abb[1145]=-4. + 29.*m_z12;

  Box1Abb[1146]=9. + Box1Abb[1145]*m_z12;

  Box1Abb[1147]=6. - 32.*m_z12 + 27.*m_z12_2 + 2.*Box1Abb[1142]*m_z1k + 2.*Box1Abb[1143]*Box1Abb[1144]*m_z1k_2 + 4.*Box1Abb[1146]*m_z1k_3 + 30.*Box1Abb[79]*m_z1k_4;

  Box1Abb[1148]=8. + 9.*m_z12 + 12.*m_z1k;

  Box1Abb[1149]=12. - Box1Abb[1148]*m_z12 + 24.*m_z1k;

  Box1Abb[1150]=33. + 50.*m_z1k;

  Box1Abb[1151]=4. + 5.*m_z1k;

  Box1Abb[1152]=-2. + 6.*Box1Abb[1151]*m_z1k;

  Box1Abb[1153]=17. + 30.*m_z1k;

  Box1Abb[1154]=8. + Box1Abb[1153]*m_z1k;

  Box1Abb[1155]=-2.*Box1Abb[1154] + Box1Abb[1152]*m_z12 + Box1Abb[1150]*m_z12_2;

  Box1Abb[1156]=3. + 8.*m_z1k_2;

  Box1Abb[1157]=95. + 109.*m_z1k;

  Box1Abb[1158]=45. + Box1Abb[1157]*m_z1k;

  Box1Abb[1159]=-5. + Box1Abb[451]*m_z1k;

  Box1Abb[1160]=30. - 8.*Box1Abb[1159]*m_z1k;

  Box1Abb[1161]=2.*Box1Abb[1156]*Box1Abb[532] + Box1Abb[1160]*m_z12 - Box1Abb[1158]*m_z12_2 - 18.*m_z12_3*m_z1k;

  Box1Abb[1162]=3. + Box1Abb[170]*m_z1k_2;

  Box1Abb[1163]=-3. + 8.*m_z1k;

  Box1Abb[1164]=-8. + Box1Abb[1163]*m_z1k;

  Box1Abb[1165]=5. + Box1Abb[1164]*m_z1k;

  Box1Abb[1166]=-5. + Box1Abb[1127]*m_z1k;

  Box1Abb[1167]=16. + Box1Abb[1166]*m_z1k;

  Box1Abb[1168]=2.*Box1Abb[1162]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[1167]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[1165]*m_z12_3 - 2.*Box1Abb[317]*pow(Box1Abb[68],3.)*m_z1k + 4.*Box1Abb[873]*m_z12_4*m_z1k;

  Box1Abb[1169]=19. + 7.*Box1Abb[233]*m_z1k;

  Box1Abb[1170]=43. + 32.*m_z1k - 59.*m_z1k_3;

  Box1Abb[1171]=-6. + Box1Abb[1170]*m_z1k;

  Box1Abb[1172]=-5. + 6.*m_z1k;

  Box1Abb[1173]=7. + 2.*Box1Abb[1172]*m_z1k;

  Box1Abb[1174]=-4. + Box1Abb[1173]*m_z1k;

  Box1Abb[1175]=2. + Box1Abb[1174]*m_z1k;

  Box1Abb[1176]=4. - 3.*Box1Abb[170]*m_z1k;

  Box1Abb[1177]=5. + 2.*Box1Abb[1176]*m_z1k;

  Box1Abb[1178]=-10. + Box1Abb[1177]*m_z1k;

  Box1Abb[1179]=5. + Box1Abb[1178]*m_z1k;

  Box1Abb[1180]=2.*Box1Abb[1175]*Box1Abb[68] + 2.*Box1Abb[1179]*m_z12 + Box1Abb[1171]*m_z12_2 - 2.*Box1Abb[1169]*m_z12_3*m_z1k - 8.*m_z12_4*m_z1k_2;

  Box1Abb[1181]=-Box1Abb[1180]*m_x_2 - Box1Abb[1147]*m_x_3 - Box1Abb[1161]*m_x_4 - Box1Abb[1155]*m_x_5 - Box1Abb[1149]*m_x_6 - 2.*Box1Abb[79]*m_x_7 - Box1Abb[1168]*m_x*m_z1k - Box1Abb[1139]*Box1Abb[68]*Box1Abb[72]*m_z12_2*m_z1k_2;

  Box1Abb[1214]=(-2.*m_cL*m_z12)/m_s;

  Box1Abb[1215]=(-2.*m_cL*m_x)/m_s;

  Box1Abb[1216]=(-2.*m_cL)/m_s;

  Box1Abb[1217]=m_cL*m_z12;

  Box1Abb[1218]=m_cL;

  Box1Abb[1219]=m_cL*m_s*m_z12;

  Box1Abb[1220]=m_cL*m_s;

  Box1Abb[1221]=(2.*m_cL)/m_s;

  Box1Abb[1222]=(-2.*m_cL*m_z1k)/m_s;

  Box1Abb[1223]=-m_cL;

  Box1Abb[1256]=(-2.*m_cR*m_z12)/m_s;

  Box1Abb[1257]=(-2.*m_cR*m_x)/m_s;

  Box1Abb[1258]=(-2.*m_cR)/m_s;

  Box1Abb[1259]=m_cR*m_z12;

  Box1Abb[1260]=m_cR;

  Box1Abb[1261]=m_cR*m_s*m_z12;

  Box1Abb[1262]=m_cR*m_s;

  Box1Abb[1263]=(2.*m_cR)/m_s;

  Box1Abb[1264]=(-2.*m_cR*m_z1k)/m_s;

  Box1Abb[1265]=-m_cR;

  Box1Abb[1266]=1. - 3.*m_x + m_z12 + m_x_2*m_z12 + Box1Abb[1]*m_z12_2;

  Box1Abb[1267]=2. + 3.*m_x - 3.*m_z12;

  Box1Abb[1268]=m_x + m_z12 + Box1Abb[1267]*m_x*m_z12;

  Box1Abb[1269]=2. + 3.*m_x;

  Box1Abb[1270]=Box1Abb[1266]*m_x - Box1Abb[1268]*m_z1k + Box1Abb[1269]*m_z12*m_z1k_2 - m_z12*m_z1k_3;

  Box1Abb[1271]=-1. + Box1Abb[62]*m_z12;

  Box1Abb[1272]=4. + Box1Abb[68]*m_z12 + m_z1k + 6.*m_z1k_2;

  Box1Abb[1273]=-2. + Box1Abb[1272]*m_z12;

  Box1Abb[1274]=-1. + m_z1k + 4.*m_z1k_2 - 4.*m_z1k_3;

  Box1Abb[1275]=pow(Box1Abb[68],2.) + Box1Abb[1274]*m_z12 - Box1Abb[265]*m_z12_2*m_z1k;

  Box1Abb[1276]=Box1Abb[1275]*m_x + Box1Abb[1273]*m_x_2 + Box1Abb[1271]*m_x_3 + m_x_4*m_z12 + Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[1277]=1. + Box1Abb[249]*m_x + 2.*m_x_2 - m_z12 + 3.*Box1Abb[0]*m_z1k + 2.*m_z1k_2;

  Box1Abb[1278]=-1. + 3.*m_z12 + 2.*m_z1k;

  Box1Abb[1279]=m_x - Box1Abb[233]*m_x_2 + 2.*m_x_3 + Box1Abb[1278]*m_x*m_z1k + Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[1280]=-4. + 5.*m_z12;

  Box1Abb[1281]=-3. + 3.*m_z12 + m_z1k;

  Box1Abb[1282]=4. + m_z12 + 7.*m_z1k;

  Box1Abb[1283]=-5. + Box1Abb[1282]*m_z12 - 4.*m_z1k;

  Box1Abb[1284]=5. - 2.*m_z12 - 11.*m_z1k;

  Box1Abb[1285]=-7. + Box1Abb[1284]*m_z12 + 8.*m_z1k;

  Box1Abb[1286]=4. + Box1Abb[1285]*m_z12;

  Box1Abb[1287]=Box1Abb[1286]*m_x_2 + Box1Abb[1280]*m_x_3*m_z12 + Box1Abb[1283]*m_x*m_z12*m_z1k - Box1Abb[1281]*m_z12_2*m_z1k_2;

  Box1Abb[1288]=11. - 2.*m_z12;

  Box1Abb[1289]=2. - 2.*pow(Box1Abb[79],2.)*m_z12 - 6.*m_z1k + Box1Abb[1288]*m_z12*m_z1k + 4.*m_z1k_2;

  Box1Abb[1290]=-9. + 4.*m_z12 + 8.*m_z1k;

  Box1Abb[1291]=-14. - 3.*Box1Abb[1290]*m_z12 + 36.*m_z1k;

  Box1Abb[1292]=-1. + 2.*Box1Abb[79]*m_z12;

  Box1Abb[1293]=-26. + 11.*m_z12;

  Box1Abb[1294]=5. + Box1Abb[1292]*m_z12 + 8.*m_z1k + Box1Abb[1293]*m_z12*m_z1k + 12.*Box1Abb[79]*m_z1k_2;

  Box1Abb[1295]=-2. + Box1Abb[1294]*m_z12 + 4.*m_z1k;

  Box1Abb[1296]=Box1Abb[1295]*m_x_2 + Box1Abb[1291]*m_x_3*m_z12 + 4.*Box1Abb[129]*m_x_4*m_z12 + Box1Abb[1289]*m_x*m_z12*m_z1k - Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[1297]=-4. + 2.*Box1Abb[370]*m_z12 - 3.*m_z1k;

  Box1Abb[1298]=1. - 2.*m_z12 + 6.*m_z1k - 9.*m_z12*m_z1k;

  Box1Abb[1299]=1. + Box1Abb[1298]*m_z12;

  Box1Abb[1300]=Box1Abb[1299]*m_x_2 + Box1Abb[60]*m_x_3*m_z12 + Box1Abb[1297]*m_x*m_z12*m_z1k - Box1Abb[1281]*m_z12_2*m_z1k_2;

  Box1Abb[1301]=2. + Box1Abb[69]*m_z12;

  Box1Abb[1302]=-2. + m_z12 + m_z12_3;

  Box1Abb[1303]=pow(Box1Abb[0],2.)*Box1Abb[1301] + Box1Abb[1302]*m_z1k + Box1Abb[77]*m_z12*m_z1k_2;

  Box1Abb[1304]=18. - 7.*m_z12 - 5.*m_z1k;

  Box1Abb[1305]=-17. + Box1Abb[1304]*m_z12 + 2.*m_z1k;

  Box1Abb[1306]=6. + Box1Abb[1305]*m_z12;

  Box1Abb[1307]=Box1Abb[1303]*m_x + Box1Abb[1306]*m_x_2 + Box1Abb[133]*m_x_3*m_z12 + pow(Box1Abb[4],2.)*Box1Abb[79]*m_z12*m_z1k;

  Box1Abb[1308]=-5. + 2.*m_z12;

  Box1Abb[1309]=3. + Box1Abb[1308]*m_z12 - 3.*m_z1k;

  Box1Abb[1310]=13. - 5.*m_z12 - 3.*m_z1k;

  Box1Abb[1311]=-9. + Box1Abb[1310]*m_z12 + 6.*m_z1k;

  Box1Abb[1312]=1. + Box1Abb[1311]*m_z12;

  Box1Abb[1313]=Box1Abb[1312]*m_x_2 + Box1Abb[1309]*Box1Abb[4]*m_x*m_z12 + Box1Abb[69]*m_x_3*m_z12 + pow(Box1Abb[4],2.)*m_z12_2*m_z1k;

  Box1Abb[1314]=1. + m_z12 - 8.*m_z12_2 + 6.*m_z12_3;

  Box1Abb[1315]=2. + m_z12 - 4.*m_z12_2;

  Box1Abb[1316]=-2.*pow(Box1Abb[0],2.)*m_z12 - 2.*Box1Abb[1314]*m_z1k + 3.*Box1Abb[1315]*m_z1k_2 - 4.*m_z1k_3;

  Box1Abb[1317]=-7. + 4.*m_z12_2 + 6.*Box1Abb[437]*m_z1k + 25.*m_z12*m_z1k;

  Box1Abb[1318]=2. + m_z12 + Box1Abb[1317]*m_z12_2 - 4.*m_z1k;

  Box1Abb[1319]=-9. + 6.*m_z12 + 8.*m_z1k;

  Box1Abb[1320]=6. + Box1Abb[1319]*m_z12 - 4.*m_z1k;

  Box1Abb[1321]=8. - 3.*Box1Abb[1320]*m_z12;

  Box1Abb[1322]=Box1Abb[1318]*m_x_2 + Box1Abb[1321]*m_x_3 + Box1Abb[1316]*m_x*m_z12 + 4.*Box1Abb[133]*m_x_4*m_z12 + Box1Abb[4]*Box1Abb[652]*m_z12_3*m_z1k;

  Box1Abb[1323]=-1. + m_z12 + 2.*m_z12*m_z1k;

  Box1Abb[1324]=-Box1Abb[1323]*m_x + m_x_2*m_z12 + Box1Abb[4]*m_z12*m_z1k;

  Box1Abb[1325]=-24. + 15.*m_z12 + 20.*m_z1k;

  Box1Abb[1326]=12. + Box1Abb[1325]*m_z12;

  Box1Abb[1327]=-2. + 11.*m_z1k;

  Box1Abb[1328]=3. + 14.*Box1Abb[68]*m_z1k;

  Box1Abb[1329]=-2. + Box1Abb[1328]*m_z12 + Box1Abb[1327]*m_z12_2 + m_z12_3 + 6.*m_z1k + 4.*Box1Abb[170]*m_z1k_2;

  Box1Abb[1330]=-5. + 13.*m_z1k;

  Box1Abb[1331]=3. + 4.*Box1Abb[536]*m_z1k;

  Box1Abb[1332]=-8. + 17.*m_z1k;

  Box1Abb[1333]=3. + 2.*Box1Abb[1332]*m_z1k;

  Box1Abb[1334]=2. + Box1Abb[1333]*m_z12 + Box1Abb[1330]*m_z12_2 + 2.*m_z12_3 + 2.*Box1Abb[1331]*m_z1k;

  Box1Abb[1335]=-2. + Box1Abb[1334]*m_z12;

  Box1Abb[1336]=-7. + 6.*m_z1k;

  Box1Abb[1337]=-34. + 50.*m_z1k;

  Box1Abb[1338]=24. + Box1Abb[1337]*m_z12 + 15.*m_z12_2 + 8.*Box1Abb[1336]*m_z1k;

  Box1Abb[1339]=Box1Abb[1338]*m_z12 + 24.*m_z1k;

  Box1Abb[1340]=-4. + Box1Abb[1339]*m_z12;

  Box1Abb[1341]=Box1Abb[1340]*m_x_3 - Box1Abb[1335]*m_x_2*m_z12 - 2.*Box1Abb[1326]*m_x_4*m_z12 + 12.*m_x_5*m_z12_2 + Box1Abb[1329]*m_x*m_z12_2*m_z1k - Box1Abb[4]*m_z12_4*m_z1k_2;

  Box1Abb[1342]=3. - Box1Abb[411]*m_z12 + m_z12_2 + m_z1k + m_z1k_2;

  Box1Abb[1343]=-5. + 2.*m_z12 + 3.*m_z1k;

  Box1Abb[1344]=3. + Box1Abb[1343]*m_z12;

  Box1Abb[1345]=-3. + m_z12 + Box1Abb[158]*m_z1k + 3.*m_z1k_2;

  Box1Abb[1346]=2. + Box1Abb[1345]*m_z12 + 3.*m_z1k;

  Box1Abb[1347]=Box1Abb[1346]*m_x_2 - Box1Abb[1344]*m_x_3 + m_x_4*m_z12 - Box1Abb[1342]*m_x*m_z12*m_z1k - Box1Abb[0]*m_z12_2*m_z1k_2;

  Box1Abb[1348]=-1. + m_z12 + 2.*m_z1k;

  Box1Abb[1349]=1. + Box1Abb[1348]*m_z12;

  Box1Abb[1350]=-Box1Abb[1349]*m_x + 3.*m_x_2*m_z12 - Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[1351]=pow(Box1Abb[0],3.) + 3.*pow(Box1Abb[0],2.)*m_z1k + 3.*m_z12*m_z1k_2;

  Box1Abb[1352]=1. + Box1Abb[358]*m_z12;

  Box1Abb[1353]=Box1Abb[1351]*m_x - 3.*Box1Abb[1352]*m_x_2 + m_x_3*m_z12 - m_z12*m_z1k_3;

  Box1Abb[1354]=-2. + m_z12 + 4.*m_z12_2;

  Box1Abb[1355]=-3. + Box1Abb[200]*m_z12;

  Box1Abb[1356]=-2.*pow(Box1Abb[0],3.)*m_z12 - Box1Abb[0]*Box1Abb[1354]*m_z1k - 2.*Box1Abb[1355]*m_z1k_2 + 2.*Box1Abb[158]*m_z1k_3 + 4.*m_z1k_4;

  Box1Abb[1357]=-8. + 5.*m_z12 + 2.*m_z1k;

  Box1Abb[1358]=4. + Box1Abb[1357]*m_z12;

  Box1Abb[1359]=6. + Box1Abb[472]*m_z12;

  Box1Abb[1360]=-12. + 5.*Box1Abb[1359]*m_z12;

  Box1Abb[1361]=2.*pow(Box1Abb[0],3.) + Box1Abb[1360]*m_z1k + 6.*Box1Abb[60]*m_z12*m_z1k_2 + 24.*m_z1k_3 - 12.*m_z1k_4;

  Box1Abb[1362]=Box1Abb[1361]*m_z12 + 2.*m_z1k;

  Box1Abb[1363]=36. + 23.*Box1Abb[79]*m_z12 - 88.*m_z1k + 64.*m_z12*m_z1k + 8.*m_z1k_2;

  Box1Abb[1364]=-16. + Box1Abb[1363]*m_z12 + 48.*m_z1k;

  Box1Abb[1365]=4. + Box1Abb[1364]*m_z12;

  Box1Abb[1366]=1. - 44.*m_z1k;

  Box1Abb[1367]=2. - 3.*m_z1k;

  Box1Abb[1368]=15. + 70.*m_z1k - 36.*m_z1k_2;

  Box1Abb[1369]=9. + 2.*m_z1k;

  Box1Abb[1370]=3. + Box1Abb[1369]*m_z1k;

  Box1Abb[1371]=6. + 2.*Box1Abb[1370]*Box1Abb[437]*m_z12 + Box1Abb[1368]*m_z12_2 + Box1Abb[1366]*m_z12_3 - 4.*m_z12_4 + 8.*Box1Abb[1367]*m_z1k;

  Box1Abb[1372]=Box1Abb[1371]*m_z12 - 4.*m_z1k;

  Box1Abb[1373]=Box1Abb[1372]*m_x_3 + Box1Abb[1365]*m_x_4 + Box1Abb[1362]*m_x_2*m_z12 - 6.*Box1Abb[1358]*m_x_5*m_z12 + 4.*m_x_6*m_z12_2 + Box1Abb[1356]*m_x*m_z12_2*m_z1k + Box1Abb[4]*m_z12_4*m_z1k_3;

  Box1Abb[1374]=m_s - m_s*m_z12;

  Box1Abb[1375]=2. + Box1Abb[881]*m_z12 + m_z12_2 + Box1Abb[317]*m_z1k;

  Box1Abb[1376]=-5. + 3.*m_z12 + 6.*m_z1k;

  Box1Abb[1377]=2. + Box1Abb[1376]*m_z12;

  Box1Abb[1378]=-1. + 8.*m_z1k;

  Box1Abb[1379]=-2. + Box1Abb[1378]*m_z12 + m_z12_2 + 12.*Box1Abb[68]*m_z1k;

  Box1Abb[1380]=2. + Box1Abb[1379]*m_z12 + 4.*m_z1k;

  Box1Abb[1381]=Box1Abb[1380]*m_x_2 - 2.*Box1Abb[1377]*m_x_3 + 4.*m_x_4*m_z12 - 2.*Box1Abb[1375]*m_x*m_z12*m_z1k - Box1Abb[0]*m_z12_2*m_z1k_2;

  Box1Abb[1382]=-9. + 3.*m_z12 - 5.*m_z1k;

  Box1Abb[1383]=6. + Box1Abb[1382]*m_z12;

  Box1Abb[1384]=-3.*pow(Box1Abb[0],4.) - 6.*Box1Abb[0]*Box1Abb[493]*m_z12*m_z1k + 10.*m_z12_2*m_z1k_2;

  Box1Abb[1385]=pow(Box1Abb[0],2.)*Box1Abb[158] + 6.*Box1Abb[0]*m_z1k_2 + 2.*m_z1k_3;

  Box1Abb[1386]=6. + Box1Abb[158]*m_z12;

  Box1Abb[1387]=pow(Box1Abb[0],2.)*Box1Abb[1386] - 12.*Box1Abb[386]*m_z1k_2 - 10.*m_z12*m_z1k_3;

  Box1Abb[1388]=2. + Box1Abb[79]*m_z12_2;

  Box1Abb[1389]=pow(Box1Abb[0],2.)*Box1Abb[1388] + 6.*pow(Box1Abb[0],4.)*m_z1k - 36.*Box1Abb[0]*m_z12*m_z1k_2 - 20.*m_z12_2*m_z1k_3;

  Box1Abb[1390]=Box1Abb[1389]*m_x_3 + 2.*Box1Abb[1384]*m_x_4 + 2.*Box1Abb[1383]*m_x_5*m_z12 + 2.*m_x_6*m_z12_2 - Box1Abb[1387]*m_x_2*m_z12*m_z1k - Box1Abb[1385]*m_x*m_z12_2*m_z1k_2 - pow(Box1Abb[0],2.)*m_z12_3*m_z1k_3;

  Box1Abb[1391]=-3. + Box1Abb[129]*m_x + Box1Abb[1016]*m_z12;

  Box1Abb[1392]=9. + 26.*m_x;

  Box1Abb[1393]=-3. + Box1Abb[1392]*m_x;

  Box1Abb[1394]=4. + 3.*m_x;

  Box1Abb[1395]=7. + 22.*m_x;

  Box1Abb[1396]=2. + 4.*m_x + 2.*Box1Abb[1393]*m_z12 - 12.*Box1Abb[1394]*m_x*m_z12_2 + Box1Abb[1395]*m_z12_3 - 3.*m_z12_4;

  Box1Abb[1397]=1. + 5.*m_x;

  Box1Abb[1398]=1. + 6.*Box1Abb[1397]*m_x;

  Box1Abb[1399]=9. + 14.*m_x;

  Box1Abb[1400]=m_z12 + 2.*m_x*m_z12;

  Box1Abb[1401]=9.*pow(Box1Abb[1400],2.) - 4.*m_x - 2.*Box1Abb[1398]*m_z12 - Box1Abb[1399]*m_z12_3 + 2.*m_z12_4;

  Box1Abb[1402]=7. - 3.*m_z12;

  Box1Abb[1403]=3. + Box1Abb[542]*m_z12;

  Box1Abb[1404]=2.*Box1Abb[1403]*m_x + 4.*Box1Abb[1402]*m_x_2 + Box1Abb[0]*m_z12_2;

  Box1Abb[1405]=-4.*m_x + m_z12_2;

  Box1Abb[1406]=Box1Abb[1391]*Box1Abb[57]*m_x_3*m_z12 + Box1Abb[1396]*m_x_2*m_z1k + Box1Abb[1401]*m_x*m_z1k_2 + Box1Abb[1404]*m_z12*m_z1k_3 + Box1Abb[1405]*m_z12*m_z1k_4;

  Box1Abb[1407]=6. + Box1Abb[493]*m_z12;

  Box1Abb[1408]=1. - 2.*m_z12;

  Box1Abb[1409]=4. - Box1Abb[1407]*m_z12 + 3.*m_z1k + 3.*Box1Abb[1408]*m_z12*m_z1k - 2.*Box1Abb[493]*m_z1k_2;

  Box1Abb[1410]=-3. + 3.*Box1Abb[15]*m_z12 - Box1Abb[210]*m_z1k;

  Box1Abb[1411]=3. - 2.*m_z12 - 14.*m_z1k;

  Box1Abb[1412]=-3. + Box1Abb[1411]*m_z12 + 18.*m_z1k;

  Box1Abb[1413]=2. + Box1Abb[1412]*m_z12;

  Box1Abb[1414]=4. + Box1Abb[493]*m_z12 + Box1Abb[326]*m_z12*m_z1k + 6.*Box1Abb[69]*m_z1k_2;

  Box1Abb[1415]=-2.*Box1Abb[15] + Box1Abb[1414]*m_z12;

  Box1Abb[1416]=Box1Abb[1415]*m_x_3 + Box1Abb[1413]*m_x_4 + Box1Abb[699]*m_x_5*m_z12 + Box1Abb[1409]*m_x_2*m_z12*m_z1k + Box1Abb[1410]*m_x*m_z12_2*m_z1k_2 - Box1Abb[0]*m_z12_3*m_z1k_3;

  Box1Abb[1417]=-3. + 7.*m_z12;

  Box1Abb[1418]=-9. + 7.*m_z12;

  Box1Abb[1419]=6. + Box1Abb[1418]*m_z12;

  Box1Abb[1420]=-2. + Box1Abb[1419]*m_z12;

  Box1Abb[1421]=pow(Box1Abb[0],3.)*Box1Abb[1420]*m_x - pow(Box1Abb[0],4.)*Box1Abb[197]*m_z12 - 2.*pow(Box1Abb[0],3.)*Box1Abb[1417]*m_x_2*m_z12 + 2.*Box1Abb[0]*Box1Abb[149]*m_x_3*m_z12_2 + 2.*m_x_4*m_z12_3;

  Box1Abb[1422]=-1. + 3.*m_z12;

  Box1Abb[1423]=6.*pow(Box1Abb[0],3.)*Box1Abb[1422]*m_x_2 + pow(Box1Abb[0],4.)*m_z12 - 4.*Box1Abb[0]*Box1Abb[472]*m_x_3*m_z12 - 4.*pow(Box1Abb[0],3.)*m_x*m_z12_2 - 10.*m_x_4*m_z12_2;

  Box1Abb[1424]=pow(Box1Abb[0],3.)*Box1Abb[493] + 24.*Box1Abb[0]*m_x_2 - 20.*m_x_3*m_z12;

  Box1Abb[1425]=pow(Box1Abb[0],3.) + m_x + Box1Abb[1010]*m_x*m_z12 + 5.*m_x_2*m_z12;

  Box1Abb[1426]=4. + 5.*m_x - 3.*m_z12;

  Box1Abb[1427]=-1. + Box1Abb[1426]*m_z12;

  Box1Abb[1428]=Box1Abb[1421]*m_x + Box1Abb[1423]*m_z12*m_z1k - Box1Abb[1424]*m_z12_2*m_z1k_2 - 4.*Box1Abb[1425]*m_z12_2*m_z1k_3 + 2.*Box1Abb[1427]*m_z12_2*m_z1k_4 - 2.*m_z12_3*m_z1k_5;

  Box1Abb[1429]=pow(Box1Abb[0],3.) + 3.*pow(Box1Abb[0],2.)*m_z1k + 3.*Box1Abb[79]*m_z1k_2;

  Box1Abb[1430]=12. - 5.*m_z12 - 9.*m_z1k;

  Box1Abb[1431]=-9. + Box1Abb[1430]*m_z12 + 12.*m_z1k;

  Box1Abb[1432]=2. + Box1Abb[1431]*m_z12;

  Box1Abb[1433]=Box1Abb[1432]*m_x_2 + Box1Abb[1429]*m_x*m_z12 + Box1Abb[699]*m_x_3*m_z12 + m_z12_2*m_z1k_3;

  Box1Abb[1434]=-9. + 4.*m_z12;

  Box1Abb[1435]=3. + Box1Abb[1434]*m_z12;

  Box1Abb[1436]=3. + Box1Abb[79]*m_z12;

  Box1Abb[1437]=2.*pow(Box1Abb[0],3.)*m_z12 + 2.*m_z1k + Box1Abb[1435]*m_z12*m_z1k - 2.*Box1Abb[1436]*m_z1k_2 - 4.*Box1Abb[0]*m_z1k_3;

  Box1Abb[1438]=-12. + 11.*m_z12;

  Box1Abb[1439]=-2. + Box1Abb[1438]*m_z12_2;

  Box1Abb[1440]=2. + Box1Abb[1051]*m_z12_2;

  Box1Abb[1441]=-2.*pow(Box1Abb[0],3.)*m_z12 - Box1Abb[0]*Box1Abb[1439]*m_z1k + 2.*Box1Abb[1440]*m_z1k_2 + 4.*Box1Abb[0]*m_z12*m_z1k_3;

  Box1Abb[1442]=-24. + 32.*m_z12 - 15.*m_z12_2 - 20.*Box1Abb[0]*m_z1k;

  Box1Abb[1443]=8. + Box1Abb[1442]*m_z12;

  Box1Abb[1444]=5. - 2.*m_z1k;

  Box1Abb[1445]=-5. + 26.*m_z1k;

  Box1Abb[1446]=-4. + m_z1k;

  Box1Abb[1447]=-1. + 4.*Box1Abb[1446]*m_z1k;

  Box1Abb[1448]=4. + 3.*Box1Abb[1447]*m_z12 + Box1Abb[1445]*m_z12_2 + 4.*m_z12_3 + 6.*Box1Abb[1444]*m_z1k;

  Box1Abb[1449]=Box1Abb[1448]*m_z12 - 12.*m_z1k;

  Box1Abb[1450]=Box1Abb[1441]*m_x_2 + Box1Abb[1449]*m_x_3 + Box1Abb[1443]*m_x_4 + 8.*Box1Abb[0]*m_x_5*m_z12 + Box1Abb[1437]*m_x*m_z12*m_z1k + Box1Abb[4]*m_z12_3*m_z1k_3;

  Box1Abb[1451]=-2. + 6.*m_x - m_z12;

  Box1Abb[1452]=2. + Box1Abb[1451]*m_z12;

  Box1Abb[1453]=1. + 6.*m_x - 2.*m_z12;

  Box1Abb[1454]=6.*m_x - m_z12;

  Box1Abb[1455]=Box1Abb[1452]*m_x_2 - 2.*Box1Abb[1453]*m_x*m_z12*m_z1k + Box1Abb[1454]*m_z12*m_z1k_2;

  Box1Abb[1456]=m_z12 + 2.*Box1Abb[12]*m_z1k + m_z1k_2;

  Box1Abb[1457]=6. - 3.*m_z12 + m_z12_2 + 6.*Box1Abb[12]*m_z1k + 12.*m_z1k_2;

  Box1Abb[1458]=-2. + Box1Abb[1457]*m_z12 + 6.*m_z1k;

  Box1Abb[1459]=m_z12 + 6.*m_z1k;

  Box1Abb[1460]=4. + Box1Abb[1459]*Box1Abb[151]*m_z12 + 24.*m_z1k;

  Box1Abb[1461]=-2. + Box1Abb[1460]*m_z12;

  Box1Abb[1462]=Box1Abb[1461]*m_x_4 - 6.*Box1Abb[406]*m_x_5*m_z12 + 6.*m_x_6*m_z12_2 - 2.*Box1Abb[1458]*m_x_3*m_z12*m_z1k + 6.*Box1Abb[1456]*m_x_2*m_z12_2*m_z1k_2 - 2.*Box1Abb[370]*m_x*m_z12_3*m_z1k_3 + m_z12_4*m_z1k_4;

  Box1Abb[1463]=-9. + 10.*m_z12;

  Box1Abb[1464]=Box1Abb[0]*Box1Abb[1418]*m_z12 + 4.*Box1Abb[0]*Box1Abb[1463]*m_z1k + 30.*m_z12*m_z1k_2;

  Box1Abb[1465]=-13. + 7.*m_z12 + 9.*m_z1k;

  Box1Abb[1466]=6. + Box1Abb[1465]*m_z12;

  Box1Abb[1467]=Box1Abb[534]*Box1Abb[68] + 4.*pow(Box1Abb[68],2.)*m_z12 + 5.*Box1Abb[68]*m_z12_2 + 2.*m_z12_3;

  Box1Abb[1468]=-3. + Box1Abb[69]*m_z12;

  Box1Abb[1469]=4. + Box1Abb[1468]*m_z12 + 6.*m_z1k + 8.*Box1Abb[493]*m_z12*m_z1k - 4.*Box1Abb[251]*m_z1k_2;

  Box1Abb[1470]=Box1Abb[1469]*m_z12 - 2.*m_z1k;

  Box1Abb[1471]=7. - Box1Abb[158]*Box1Abb[493]*m_z12 + 8.*Box1Abb[0]*m_z12*m_z1k + 14.*Box1Abb[0]*m_z1k_2 + 6.*m_z1k_3;

  Box1Abb[1472]=-1. + Box1Abb[1471]*m_z12;

  Box1Abb[1473]=4. + Box1Abb[493]*m_z12;

  Box1Abb[1474]=-5. + Box1Abb[1473]*m_z12;

  Box1Abb[1475]=-8. + 3.*m_z12;

  Box1Abb[1476]=6. + Box1Abb[1475]*m_z12;

  Box1Abb[1477]=-9. + 8.*m_z12;

  Box1Abb[1478]=5. + Box1Abb[1474]*m_z12 - 8.*m_z1k + 6.*Box1Abb[1476]*m_z12*m_z1k + 4.*Box1Abb[0]*Box1Abb[1477]*m_z1k_2 + 20.*m_z12*m_z1k_3;

  Box1Abb[1479]=2. - Box1Abb[1478]*m_z12 - 2.*m_z1k;

  Box1Abb[1480]=Box1Abb[1479]*m_x_3 + Box1Abb[1464]*m_x_4*m_z12 - 2.*Box1Abb[1466]*m_x_5*m_z12 + 4.*m_x_6*m_z12_2 + Box1Abb[0]*Box1Abb[1470]*m_x_2*m_z1k + Box1Abb[1472]*m_x*m_z12*m_z1k_2 - Box1Abb[1467]*m_z12_2*m_z1k_3;

  Box1Abb[1481]=-1. + Box1Abb[129]*Box1Abb[493]*m_z12;

  Box1Abb[1482]=-26. + 7.*m_z12;

  Box1Abb[1483]=32. + Box1Abb[1482]*m_z12;

  Box1Abb[1484]=-10. + Box1Abb[1483]*m_z12;

  Box1Abb[1485]=2. + Box1Abb[1484]*m_z12;

  Box1Abb[1486]=-13. + 6.*m_z12;

  Box1Abb[1487]=4. + Box1Abb[1486]*m_z12;

  Box1Abb[1488]=-4. + 11.*m_z12;

  Box1Abb[1489]=2. + Box1Abb[1403]*m_z12 - 4.*m_z1k + Box1Abb[1481]*m_z12*m_z1k + Box1Abb[1485]*m_z1k_2 + 3.*Box1Abb[1487]*m_z12*m_z1k_3 + Box1Abb[1488]*m_z12*m_z1k_4;

  Box1Abb[1490]=2. + 2.*Box1Abb[170]*m_z12 + m_z12_2 + Box1Abb[170]*m_z1k;

  Box1Abb[1491]=1. + Box1Abb[1490]*m_z12 - m_z1k;

  Box1Abb[1492]=5. - 6.*m_z12;

  Box1Abb[1493]=-27. + Box1Abb[1492]*m_z12;

  Box1Abb[1494]=14. + 9.*m_z12;

  Box1Abb[1495]=59. - 42.*m_z12;

  Box1Abb[1496]=-109. + Box1Abb[1495]*m_z12;

  Box1Abb[1497]=72. + Box1Abb[1496]*m_z12;

  Box1Abb[1498]=16. - 23.*m_z12;

  Box1Abb[1499]=17. + Box1Abb[1493]*m_z12 + 55.*m_z1k - 3.*Box1Abb[1494]*m_z12*m_z1k + Box1Abb[1497]*m_z1k_2 + 5.*Box1Abb[1498]*m_z12*m_z1k_3;

  Box1Abb[1500]=7. + Box1Abb[1499]*m_z12;

  Box1Abb[1501]=8. + 3.*Box1Abb[79]*m_z12;

  Box1Abb[1502]=88. - 21.*m_z12;

  Box1Abb[1503]=-41. + Box1Abb[1502]*m_z12;

  Box1Abb[1504]=82. + 5.*Box1Abb[542]*m_z12;

  Box1Abb[1505]=-30. + Box1Abb[1504]*m_z12;

  Box1Abb[1506]=-15. + Box1Abb[1505]*m_z12;

  Box1Abb[1507]=-22. + 17.*m_z12;

  Box1Abb[1508]=12. + Box1Abb[1507]*m_z12;

  Box1Abb[1509]=-12. + 19.*m_z12;

  Box1Abb[1510]=-1. + 2.*Box1Abb[1501]*m_z12 - 34.*m_z1k + Box1Abb[1503]*m_z12*m_z1k + Box1Abb[1506]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[1508]*m_z1k_3 + 5.*Box1Abb[1509]*m_z12*m_z1k_4;

  Box1Abb[1511]=-8. + Box1Abb[1510]*m_z12 + 16.*m_z1k;

  Box1Abb[1512]=11. + 31.*m_z1k;

  Box1Abb[1513]=5. + Box1Abb[1512]*m_z12 - 24.*m_z1k;

  Box1Abb[1514]=12. - Box1Abb[1513]*m_z12;

  Box1Abb[1515]=14. + 81.*m_z1k;

  Box1Abb[1516]=7. + Box1Abb[1515]*m_z1k;

  Box1Abb[1517]=28. + Box1Abb[1516]*m_z12 + 2.*Box1Abb[532]*m_z12_2 + 20.*Box1Abb[1367]*m_z1k;

  Box1Abb[1518]=-31. + Box1Abb[1517]*m_z12 - 48.*m_z1k;

  Box1Abb[1519]=9. - 11.*m_z1k;

  Box1Abb[1520]=-3. + m_z1k_2 + 12.*m_z1k_3;

  Box1Abb[1521]=-31. + 26.*m_z1k;

  Box1Abb[1522]=20. + Box1Abb[1521]*m_z1k;

  Box1Abb[1523]=-6. + Box1Abb[1522]*m_z1k;

  Box1Abb[1524]=1. + Box1Abb[1523]*m_z1k;

  Box1Abb[1525]=139. - 45.*m_z1k;

  Box1Abb[1526]=-122. + Box1Abb[1525]*m_z1k;

  Box1Abb[1527]=58. + Box1Abb[1526]*m_z1k;

  Box1Abb[1528]=-54. + Box1Abb[1527]*m_z1k;

  Box1Abb[1529]=6. + Box1Abb[1528]*m_z1k;

  Box1Abb[1530]=-79. + 24.*m_z1k;

  Box1Abb[1531]=70. + Box1Abb[1530]*m_z1k;

  Box1Abb[1532]=-27. + Box1Abb[1531]*m_z1k;

  Box1Abb[1533]=44. + Box1Abb[1532]*m_z1k;

  Box1Abb[1534]=-8. + Box1Abb[1533]*m_z1k;

  Box1Abb[1535]=pow(Box1Abb[68],2.) + Box1Abb[1520]*Box1Abb[68]*m_z12 + Box1Abb[1534]*m_z12_2 + Box1Abb[1529]*m_z12_3 - 2.*Box1Abb[1524]*m_z12_4 + Box1Abb[1519]*m_z12_5*m_z1k_2;

  Box1Abb[1536]=Box1Abb[1535]*m_x_2 + Box1Abb[1511]*m_x_3 + Box1Abb[1500]*m_x_4 + Box1Abb[1518]*m_x_5*m_z12 + Box1Abb[1514]*m_x_6*m_z12 + Box1Abb[1280]*m_x_7*m_z12_2 + Box1Abb[1489]*Box1Abb[68]*m_x*m_z12*m_z1k - Box1Abb[1491]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_2;

  Box1Abb[1537]=10. - 17.*m_z12;

  Box1Abb[1538]=12. - 2.*Box1Abb[522]*m_z12 + Box1Abb[1537]*m_z12*m_z1k;

  Box1Abb[1539]=7. - 8.*m_z1k + 10.*m_z1k_3;

  Box1Abb[1540]=-23. + 2.*Box1Abb[411]*m_z1k;

  Box1Abb[1541]=7. + Box1Abb[1540]*m_z1k;

  Box1Abb[1542]=10. - 3.*m_z1k;

  Box1Abb[1543]=16. + Box1Abb[1542]*m_z1k;

  Box1Abb[1544]=-28. + Box1Abb[1543]*m_z1k;

  Box1Abb[1545]=10. + Box1Abb[1544]*m_z1k;

  Box1Abb[1546]=2.*pow(Box1Abb[68],2.) + Box1Abb[1539]*Box1Abb[68]*m_z12 + Box1Abb[1545]*m_z12_2 - Box1Abb[1541]*m_z12_3 + Box1Abb[1115]*m_z12_4;

  Box1Abb[1547]=-23. + 72.*m_z1k;

  Box1Abb[1548]=27. + 4.*m_z1k + 22.*m_z1k_2;

  Box1Abb[1549]=74. - 45.*m_z1k;

  Box1Abb[1550]=-45. + Box1Abb[1549]*m_z1k;

  Box1Abb[1551]=44. + Box1Abb[1550]*m_z1k;

  Box1Abb[1552]=-72. + 5.*m_z1k;

  Box1Abb[1553]=26. + Box1Abb[1552]*m_z1k;

  Box1Abb[1554]=-7. + 2.*Box1Abb[1553]*m_z1k;

  Box1Abb[1555]=-32. + Box1Abb[1554]*m_z12 + Box1Abb[1551]*m_z12_2 - Box1Abb[1548]*m_z12_3 - 2.*m_z12_4 + Box1Abb[1547]*m_z1k;

  Box1Abb[1556]=20. + Box1Abb[1555]*m_z12 + 6.*m_z1k;

  Box1Abb[1557]=2. + Box1Abb[178]*m_z1k;

  Box1Abb[1558]=-19. + 28.*Box1Abb[170]*m_z1k;

  Box1Abb[1559]=11. + Box1Abb[1558]*m_z1k;

  Box1Abb[1560]=-21. + 16.*m_z1k;

  Box1Abb[1561]=10. + Box1Abb[1560]*m_z1k;

  Box1Abb[1562]=3. + Box1Abb[1561]*m_z1k;

  Box1Abb[1563]=58. + 5.*m_z1k;

  Box1Abb[1564]=-83. + Box1Abb[1563]*m_z1k;

  Box1Abb[1565]=7. + 2.*Box1Abb[1564]*m_z1k;

  Box1Abb[1566]=43. + Box1Abb[1565]*m_z1k;

  Box1Abb[1567]=-96. + 25.*m_z1k;

  Box1Abb[1568]=174. + Box1Abb[1567]*m_z1k;

  Box1Abb[1569]=34. + Box1Abb[1568]*m_z1k;

  Box1Abb[1570]=-46. + Box1Abb[1569]*m_z1k;

  Box1Abb[1571]=-4. - 3.*Box1Abb[1562]*m_z12 + Box1Abb[1566]*m_z12_2 + Box1Abb[1570]*m_z12_3 + Box1Abb[1559]*m_z12_4 + 3.*Box1Abb[1557]*m_z12_5 + 10.*m_z1k - 6.*m_z1k_2;

  Box1Abb[1572]=37. - 9.*m_z1k;

  Box1Abb[1573]=-16. + 39.*m_z1k;

  Box1Abb[1574]=-7. + Box1Abb[1573]*m_z1k;

  Box1Abb[1575]=8. + Box1Abb[1574]*m_z12 + 6.*Box1Abb[178]*m_z12_2 + 2.*Box1Abb[1572]*m_z1k;

  Box1Abb[1576]=-5. + Box1Abb[1575]*m_z12 - 48.*m_z1k;

  Box1Abb[1577]=-2. + Box1Abb[1576]*m_z12;

  Box1Abb[1578]=-13. + 12.*m_z1k;

  Box1Abb[1579]=8. + Box1Abb[1578]*m_z1k;

  Box1Abb[1580]=16. - 5.*m_z1k;

  Box1Abb[1581]=3. + Box1Abb[1580]*m_z1k;

  Box1Abb[1582]=-6. + Box1Abb[1581]*m_z1k;

  Box1Abb[1583]=16. - 3.*m_z1k;

  Box1Abb[1584]=-91. + 4.*Box1Abb[1583]*m_z1k;

  Box1Abb[1585]=12. + Box1Abb[1584]*m_z1k;

  Box1Abb[1586]=11. + Box1Abb[1585]*m_z1k;

  Box1Abb[1587]=16. + 9.*m_z1k;

  Box1Abb[1588]=-41. + Box1Abb[1587]*m_z1k;

  Box1Abb[1589]=59. + 2.*Box1Abb[1588]*m_z1k;

  Box1Abb[1590]=-9. + Box1Abb[1589]*m_z1k;

  Box1Abb[1591]=-34. + 3.*m_z1k;

  Box1Abb[1592]=150. + Box1Abb[1591]*m_z1k;

  Box1Abb[1593]=-176. + Box1Abb[1592]*m_z1k;

  Box1Abb[1594]=63. + Box1Abb[1593]*m_z1k;

  Box1Abb[1595]=2. + Box1Abb[1594]*m_z1k;

  Box1Abb[1596]=2.*pow(Box1Abb[68],3.) + Box1Abb[1579]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[1590]*Box1Abb[68]*m_z12_2 - Box1Abb[1595]*m_z12_3 + Box1Abb[1586]*m_z12_4 + Box1Abb[1582]*m_z12_5;

  Box1Abb[1597]=Box1Abb[1596]*m_x_2 + Box1Abb[1571]*m_x_3 + Box1Abb[1556]*m_x_4 + Box1Abb[1577]*m_x_5 + Box1Abb[1546]*pow(Box1Abb[68],2.)*m_x*m_z12 + Box1Abb[1538]*m_x_6*m_z12 + Box1Abb[133]*m_x_7*m_z12_2 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*Box1Abb[79]*m_z12_2*m_z1k;

  Box1Abb[1598]=21. + 2.*m_z12;

  Box1Abb[1599]=-45. + Box1Abb[1598]*m_z12;

  Box1Abb[1600]=15. + Box1Abb[1599]*m_z12;

  Box1Abb[1601]=12. + Box1Abb[1308]*m_z12;

  Box1Abb[1602]=-6. + 11.*m_z12;

  Box1Abb[1603]=16. - 9.*m_z12;

  Box1Abb[1604]=-3. - Box1Abb[1600]*m_z12 + 5.*m_z1k + Box1Abb[1601]*m_z12*m_z1k - 2.*Box1Abb[1602]*Box1Abb[857]*m_z1k_2 + 5.*Box1Abb[1603]*m_z12*m_z1k_3;

  Box1Abb[1605]=5. + 6.*m_z12;

  Box1Abb[1606]=-41. + Box1Abb[1605]*m_z12;

  Box1Abb[1607]=36. + Box1Abb[1606]*m_z12;

  Box1Abb[1608]=3. + Box1Abb[493]*m_z12;

  Box1Abb[1609]=-7. + 3.*Box1Abb[1608]*m_z12;

  Box1Abb[1610]=-21. + m_z12;

  Box1Abb[1611]=37. + Box1Abb[1610]*m_z12;

  Box1Abb[1612]=-23. + Box1Abb[1611]*m_z12;

  Box1Abb[1613]=-30. + 7.*m_z12;

  Box1Abb[1614]=4. + Box1Abb[1613]*m_z12;

  Box1Abb[1615]=-12. + 5.*m_z12;

  Box1Abb[1616]=-7. + Box1Abb[1607]*m_z12 + 2.*Box1Abb[1609]*m_z12*m_z1k + Box1Abb[1422]*Box1Abb[1612]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[1614]*m_z1k_3 + 5.*Box1Abb[1615]*m_z12*m_z1k_4;

  Box1Abb[1617]=4. - 17.*m_z1k;

  Box1Abb[1618]=-4. + Box1Abb[1617]*m_z12 + 24.*m_z1k;

  Box1Abb[1619]=4. + Box1Abb[1618]*m_z12;

  Box1Abb[1620]=-9. + 5.*m_z12;

  Box1Abb[1621]=-7. + m_z12;

  Box1Abb[1622]=-20. + 13.*m_z12;

  Box1Abb[1623]=13. + 2.*Box1Abb[1620]*m_z12 + 36.*m_z1k + 6.*Box1Abb[1621]*m_z12*m_z1k + 3.*Box1Abb[1622]*m_z1k_2;

  Box1Abb[1624]=-11. + Box1Abb[1623]*m_z12 - 16.*m_z1k;

  Box1Abb[1625]=-7. - 4.*Box1Abb[170]*m_z1k;

  Box1Abb[1626]=10. + 3.*m_z1k;

  Box1Abb[1627]=5. + Box1Abb[1626]*Box1Abb[68]*m_z1k;

  Box1Abb[1628]=3. - Box1Abb[1627]*m_z12 + Box1Abb[1115]*m_z12_2 + Box1Abb[1625]*m_z1k;

  Box1Abb[1629]=7. - 9.*m_z1k + 6.*m_z1k_2;

  Box1Abb[1630]=-11. + 4.*Box1Abb[1629]*m_z1k;

  Box1Abb[1631]=-23. + m_z1k;

  Box1Abb[1632]=31. + Box1Abb[1631]*m_z1k;

  Box1Abb[1633]=-40. + 3.*Box1Abb[1632]*m_z1k;

  Box1Abb[1634]=-1. + Box1Abb[1633]*m_z1k;

  Box1Abb[1635]=-73. - 12.*Box1Abb[396]*m_z1k;

  Box1Abb[1636]=8. + Box1Abb[1635]*m_z1k;

  Box1Abb[1637]=13. + Box1Abb[1636]*m_z1k;

  Box1Abb[1638]=Box1Abb[665]*pow(Box1Abb[68],3.) + Box1Abb[1630]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[1634]*Box1Abb[68]*m_z12_2 + Box1Abb[1637]*m_z12_3 + Box1Abb[1582]*m_z12_4;

  Box1Abb[1639]=Box1Abb[1638]*m_x_2 + Box1Abb[1616]*m_x_3 + Box1Abb[1604]*m_x_4 + Box1Abb[1624]*m_x_5 + Box1Abb[1619]*m_x_6 + Box1Abb[1628]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_x*m_z12 + Box1Abb[129]*m_x_7*m_z12 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12_2*m_z1k;

  Box1Abb[1640]=2. + Box1Abb[79]*m_z12;

  Box1Abb[1641]=2.*Box1Abb[0]*m_z12 + Box1Abb[1640]*m_z1k + Box1Abb[79]*m_z1k_2;

  Box1Abb[1642]=32. + Box1Abb[1621]*m_z12;

  Box1Abb[1643]=-23. + Box1Abb[1642]*m_z12;

  Box1Abb[1644]=-22. + 3.*m_z12;

  Box1Abb[1645]=-62. + Box1Abb[158]*Box1Abb[1644]*m_z12;

  Box1Abb[1646]=17. + Box1Abb[1645]*m_z12;

  Box1Abb[1647]=-10. + 7.*m_z12;

  Box1Abb[1648]=9. + Box1Abb[1647]*m_z12;

  Box1Abb[1649]=2. + 5.*m_z12;

  Box1Abb[1650]=1. + 18.*Box1Abb[0]*m_z12_2 + 2.*m_z1k + 2.*Box1Abb[1643]*m_z12*m_z1k + Box1Abb[1646]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[1648]*m_z1k_3 + 5.*Box1Abb[1649]*m_z12*m_z1k_4;

  Box1Abb[1651]=-3. + 5.*Box1Abb[124]*m_z12;

  Box1Abb[1652]=-23. + 9.*m_z12 - 6.*m_z12_2;

  Box1Abb[1653]=16. + Box1Abb[1652]*m_z12;

  Box1Abb[1654]=-28. + 9.*m_z12;

  Box1Abb[1655]=43. + Box1Abb[1654]*m_z12;

  Box1Abb[1656]=-18. + Box1Abb[1655]*m_z12;

  Box1Abb[1657]=36. - 5.*m_z12;

  Box1Abb[1658]=-82. + Box1Abb[1657]*m_z12;

  Box1Abb[1659]=70. + Box1Abb[1658]*m_z12;

  Box1Abb[1660]=-31. + Box1Abb[1659]*m_z12;

  Box1Abb[1661]=11. - 6.*m_z12;

  Box1Abb[1662]=-2. + Box1Abb[1661]*m_z12;

  Box1Abb[1663]=7. + Box1Abb[1662]*m_z12;

  Box1Abb[1664]=-2. + Box1Abb[1651]*m_z12 + m_z1k + Box1Abb[1653]*m_z12*m_z1k + Box1Abb[0]*Box1Abb[1656]*m_z1k_2 + Box1Abb[1660]*m_z1k_3 + 2.*Box1Abb[1663]*m_z1k_4 - 3.*Box1Abb[653]*m_z12*m_z1k_5;

  Box1Abb[1665]=8. + Box1Abb[824]*m_z12 - 10.*m_z1k;

  Box1Abb[1666]=6. - Box1Abb[1665]*m_z12;

  Box1Abb[1667]=6. + m_z1k - m_z1k_2;

  Box1Abb[1668]=-1. + Box1Abb[785]*m_z1k;

  Box1Abb[1669]=-1. + Box1Abb[68]*m_z1k;

  Box1Abb[1670]=1. + Box1Abb[1669]*Box1Abb[210]*m_z1k;

  Box1Abb[1671]=4. + 3.*m_z1k;

  Box1Abb[1672]=6. + Box1Abb[1446]*Box1Abb[1671]*m_z1k;

  Box1Abb[1673]=4. + Box1Abb[1672]*m_z1k;

  Box1Abb[1674]=-2.*Box1Abb[15]*Box1Abb[1668]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[1673]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[1670]*m_z12_3 + 2.*pow(Box1Abb[68],3.)*m_z1k + Box1Abb[1667]*m_z12_4*m_z1k;

  Box1Abb[1675]=27. + 44.*m_z1k;

  Box1Abb[1676]=3. - 22.*m_z1k;

  Box1Abb[1677]=-14. + Box1Abb[1676]*m_z1k;

  Box1Abb[1678]=-7. + 45.*m_z1k;

  Box1Abb[1679]=11. + Box1Abb[1678]*m_z1k;

  Box1Abb[1680]=-47. + 5.*m_z1k;

  Box1Abb[1681]=-13. + Box1Abb[1680]*m_z1k;

  Box1Abb[1682]=-17. + 2.*Box1Abb[1681]*m_z1k;

  Box1Abb[1683]=16. + Box1Abb[1682]*m_z12 - Box1Abb[1679]*Box1Abb[68]*m_z12_2 + Box1Abb[1677]*m_z12_3 + Box1Abb[1675]*m_z1k;

  Box1Abb[1684]=23. - 9.*m_z1k;

  Box1Abb[1685]=4. + 6.*m_z1k;

  Box1Abb[1686]=-10. + 39.*m_z1k;

  Box1Abb[1687]=-5. + Box1Abb[1686]*m_z1k;

  Box1Abb[1688]=28. + Box1Abb[1687]*m_z12 + Box1Abb[1685]*m_z12_2 + 2.*Box1Abb[1684]*m_z1k;

  Box1Abb[1689]=-21. + Box1Abb[1688]*m_z12 - 26.*m_z1k;

  Box1Abb[1690]=-Box1Abb[1674]*Box1Abb[68]*m_x + Box1Abb[1664]*m_x_2 + Box1Abb[1650]*m_x_3 + Box1Abb[1683]*m_x_4 + Box1Abb[1689]*m_x_5 + Box1Abb[1666]*m_x_6 + Box1Abb[133]*m_x_7*m_z12 + Box1Abb[1641]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[1691]=2. + Box1Abb[358]*m_z12 - 2.*m_z1k;

  Box1Abb[1692]=3. + 7.*m_z1k;

  Box1Abb[1693]=-1. + 9.*m_z1k_2;

  Box1Abb[1694]=-11. + 9.*m_z1k;

  Box1Abb[1695]=-3. + Box1Abb[1694]*m_z1k;

  Box1Abb[1696]=-36. + 11.*m_z1k;

  Box1Abb[1697]=-4. + Box1Abb[1696]*m_z1k;

  Box1Abb[1698]=2.*pow(Box1Abb[68],3.) - 2.*Box1Abb[1693]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[1697]*Box1Abb[68]*m_z12_2*m_z1k + 2.*Box1Abb[1695]*m_z12_3*m_z1k + Box1Abb[1692]*m_z12_4*m_z1k;

  Box1Abb[1699]=10. + 31.*m_z1k;

  Box1Abb[1700]=12. - Box1Abb[1699]*m_z12 + 38.*m_z1k;

  Box1Abb[1701]=2. + Box1Abb[1700]*m_z12;

  Box1Abb[1702]=6. + m_z1k + 42.*m_z1k_2;

  Box1Abb[1703]=52. - 115.*m_z1k;

  Box1Abb[1704]=-4. + Box1Abb[1703]*m_z1k;

  Box1Abb[1705]=5. + Box1Abb[1704]*m_z1k;

  Box1Abb[1706]=3. + m_z1k + 25.*m_z1k_2;

  Box1Abb[1707]=1. + 2.*Box1Abb[1706]*m_z1k;

  Box1Abb[1708]=2. + 3.*Box1Abb[1707]*m_z12 + Box1Abb[1705]*m_z12_2 - Box1Abb[1702]*m_z12_3 + m_z1k + 4.*m_z1k_2;

  Box1Abb[1709]=5. + Box1Abb[1515]*m_z1k;

  Box1Abb[1710]=-6. + Box1Abb[1709]*m_z12 + 2.*Box1Abb[532]*m_z12_2 - 34.*Box1Abb[265]*m_z1k;

  Box1Abb[1711]=-7. + Box1Abb[1710]*m_z12 - 6.*m_z1k;

  Box1Abb[1712]=1. + Box1Abb[834]*m_z1k;

  Box1Abb[1713]=-15. + 17.*m_z1k;

  Box1Abb[1714]=3. + Box1Abb[1713]*m_z1k;

  Box1Abb[1715]=-38. + 65.*m_z1k;

  Box1Abb[1716]=15. + Box1Abb[1715]*m_z1k;

  Box1Abb[1717]=17. + Box1Abb[1716]*m_z1k;

  Box1Abb[1718]=-148. + 95.*m_z1k;

  Box1Abb[1719]=28. + Box1Abb[1718]*m_z1k;

  Box1Abb[1720]=46. + Box1Abb[1719]*m_z1k;

  Box1Abb[1721]=-10. + Box1Abb[1720]*m_z1k;

  Box1Abb[1722]=Box1Abb[1712]*Box1Abb[210] + Box1Abb[1721]*m_z12_2 + 2.*Box1Abb[1714]*Box1Abb[431]*m_z12_3 - 2.*Box1Abb[1717]*m_z12*m_z1k + 5.*m_z12_4*m_z1k_2;

  Box1Abb[1723]=1. - 2.*m_z1k;

  Box1Abb[1724]=-1. + 13.*m_z1k;

  Box1Abb[1725]=2. + Box1Abb[1724]*m_z1k;

  Box1Abb[1726]=-19. + 33.*m_z1k;

  Box1Abb[1727]=6. + Box1Abb[1726]*m_z1k;

  Box1Abb[1728]=-29. + 2.*Box1Abb[1727]*m_z1k;

  Box1Abb[1729]=3. + Box1Abb[1728]*m_z1k;

  Box1Abb[1730]=142. - 45.*m_z1k;

  Box1Abb[1731]=-86. + Box1Abb[1730]*m_z1k;

  Box1Abb[1732]=45. + Box1Abb[1731]*m_z1k;

  Box1Abb[1733]=-29. + Box1Abb[1732]*m_z1k;

  Box1Abb[1734]=5. + Box1Abb[1733]*m_z1k;

  Box1Abb[1735]=Box1Abb[1729]*Box1Abb[68]*m_z12 + Box1Abb[1734]*m_z12_2 - pow(Box1Abb[1723],2.)*Box1Abb[1725]*m_z12_3 - 3.*pow(Box1Abb[68],2.)*Box1Abb[873]*m_z1k + Box1Abb[1519]*m_z12_4*m_z1k_2;

  Box1Abb[1736]=Box1Abb[1735]*m_x_2 + Box1Abb[1722]*m_x_3 + Box1Abb[1708]*m_x_4 + Box1Abb[1711]*m_x_5 + Box1Abb[1701]*m_x_6 + Box1Abb[699]*m_x_7*m_z12 + Box1Abb[1698]*Box1Abb[68]*m_x*m_z1k - Box1Abb[1691]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[1737]=15. + m_z12;

  Box1Abb[1738]=28. + 5.*m_z12;

  Box1Abb[1739]=-16. + Box1Abb[1738]*m_z12;

  Box1Abb[1740]=-10. + Box1Abb[1739]*m_z12;

  Box1Abb[1741]=5. + Box1Abb[1740]*m_z12;

  Box1Abb[1742]=4. + 17.*m_z12_2;

  Box1Abb[1743]=-15. + Box1Abb[1737]*m_z12 - 8.*m_z1k + 16.*Box1Abb[8]*m_z12*m_z1k + Box1Abb[1741]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[1742]*m_z1k_3 + 5.*Box1Abb[1509]*m_z12*m_z1k_4;

  Box1Abb[1744]=-3. + 8.*Box1Abb[0]*Box1Abb[158]*m_z12;

  Box1Abb[1745]=11. + 2.*m_z12;

  Box1Abb[1746]=-26. + Box1Abb[1745]*m_z12;

  Box1Abb[1747]=7. + Box1Abb[1746]*m_z12;

  Box1Abb[1748]=4. - 11.*m_z12;

  Box1Abb[1749]=-6. + Box1Abb[1748]*m_z12;

  Box1Abb[1750]=12. + Box1Abb[1749]*m_z12;

  Box1Abb[1751]=-11. + Box1Abb[1750]*m_z12;

  Box1Abb[1752]=23. - 13.*m_z12;

  Box1Abb[1753]=-6. + Box1Abb[1752]*m_z12;

  Box1Abb[1754]=1. + Box1Abb[1753]*m_z12;

  Box1Abb[1755]=8. - 15.*m_z12;

  Box1Abb[1756]=3. - 3.*m_z12 + Box1Abb[1744]*m_z1k - Box1Abb[0]*Box1Abb[1747]*m_z1k_2 + Box1Abb[1751]*m_z1k_3 + 4.*Box1Abb[1754]*m_z1k_4 + 3.*Box1Abb[1755]*m_z12*m_z1k_5;

  Box1Abb[1757]=7. + 2.*Box1Abb[1446]*m_z1k;

  Box1Abb[1758]=5. + Box1Abb[881]*m_z1k;

  Box1Abb[1759]=Box1Abb[1758]*Box1Abb[68] + Box1Abb[1757]*m_z12 + Box1Abb[170]*m_z12_2;

  Box1Abb[1760]=16. + 31.*m_z1k;

  Box1Abb[1761]=8. - Box1Abb[1760]*m_z12 + 24.*m_z1k;

  Box1Abb[1762]=4. + Box1Abb[1761]*m_z12;

  Box1Abb[1763]=46. + 81.*m_z1k;

  Box1Abb[1764]=18. + Box1Abb[1763]*m_z1k;

  Box1Abb[1765]=5. + Box1Abb[1764]*m_z12 - 12.*Box1Abb[451]*m_z1k + 10.*m_z12_2*m_z1k;

  Box1Abb[1766]=-17. + Box1Abb[1765]*m_z12 - 16.*m_z1k;

  Box1Abb[1767]=8. + 21.*m_z1k;

  Box1Abb[1768]=23. + 24.*m_z1k;

  Box1Abb[1769]=-1. + 4.*Box1Abb[532]*m_z1k;

  Box1Abb[1770]=-21. + 4.*Box1Abb[1769]*m_z1k;

  Box1Abb[1771]=18. + 115.*m_z1k;

  Box1Abb[1772]=17. + Box1Abb[1771]*m_z1k;

  Box1Abb[1773]=8. + Box1Abb[1772]*m_z1k;

  Box1Abb[1774]=25. + Box1Abb[1770]*m_z12 - Box1Abb[1773]*m_z12_2 + Box1Abb[1768]*m_z1k - 2.*Box1Abb[1767]*m_z12_3*m_z1k;

  Box1Abb[1775]=8. - m_z1k + 4.*m_z1k_3;

  Box1Abb[1776]=6. + 7.*m_z1k;

  Box1Abb[1777]=-5. + Box1Abb[1776]*m_z1k;

  Box1Abb[1778]=-14. + 9.*m_z1k;

  Box1Abb[1779]=-7. + Box1Abb[1778]*m_z1k;

  Box1Abb[1780]=9. + Box1Abb[1779]*m_z1k;

  Box1Abb[1781]=-1. + Box1Abb[1780]*m_z1k;

  Box1Abb[1782]=-31. + 11.*m_z1k;

  Box1Abb[1783]=7. + Box1Abb[1782]*m_z1k;

  Box1Abb[1784]=17. + Box1Abb[1783]*m_z1k;

  Box1Abb[1785]=-10. + Box1Abb[1784]*m_z1k;

  Box1Abb[1786]=-Box1Abb[1775]*pow(Box1Abb[68],2.) + Box1Abb[1785]*Box1Abb[68]*m_z12 + 2.*Box1Abb[1781]*m_z12_2 + Box1Abb[1777]*m_z12_3*m_z1k;

  Box1Abb[1787]=Box1Abb[1756]*m_x_2 + Box1Abb[1743]*m_x_3 + Box1Abb[1774]*m_x_4 + Box1Abb[1766]*m_x_5 + Box1Abb[1762]*m_x_6 + Box1Abb[1280]*m_x_7*m_z12 + Box1Abb[1786]*m_x*m_z12*m_z1k - Box1Abb[1759]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k_2;

  Box1Abb[1788]=-2. + m_z12 + 2.*m_z12*m_z1k;

  Box1Abb[1789]=-Box1Abb[1788]*m_x + m_x_2*m_z12 + Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[1790]=8. + 3.*m_z12;

  Box1Abb[1791]=-3. + m_z12 + 5.*m_z1k;

  Box1Abb[1792]=-13. - 10.*m_z12 - 20.*m_z1k + 8.*m_z12*m_z1k;

  Box1Abb[1793]=4. + 3.*m_z12;

  Box1Abb[1794]=-35. + Box1Abb[1793]*m_z12;

  Box1Abb[1795]=12. + 7.*m_z12;

  Box1Abb[1796]=-26. + Box1Abb[1795]*m_z12;

  Box1Abb[1797]=-4. + Box1Abb[1796]*m_z12;

  Box1Abb[1798]=5. + 11.*m_z12;

  Box1Abb[1799]=-10. + Box1Abb[1492]*m_z12 + 24.*m_z1k + Box1Abb[1794]*m_z12*m_z1k + 2.*Box1Abb[1797]*m_z1k_2 + 8.*Box1Abb[1798]*m_z12*m_z1k_3;

  Box1Abb[1800]=-12. + m_z12;

  Box1Abb[1801]=7. + Box1Abb[1800]*m_z12;

  Box1Abb[1802]=-23. + 10.*Box1Abb[12]*m_z12;

  Box1Abb[1803]=-3. + Box1Abb[79]*m_z12 + 7.*m_z1k + 4.*m_z12*m_z1k + Box1Abb[1801]*m_z1k_2 + Box1Abb[1802]*m_z1k_3 + 12.*m_z1k_4;

  Box1Abb[1804]=6. - 20.*m_z12 + 3.*m_z12_3;

  Box1Abb[1805]=37. + Box1Abb[158]*m_z12;

  Box1Abb[1806]=-27. + 2.*Box1Abb[222]*m_z12;

  Box1Abb[1807]=-2. + Box1Abb[1806]*m_z12;

  Box1Abb[1808]=40. + 47.*m_z12;

  Box1Abb[1809]=-2. + m_z12 - m_z12_2 + Box1Abb[1804]*m_z1k + Box1Abb[1805]*m_z12*m_z1k_2 + 2.*Box1Abb[1807]*m_z1k_3 + Box1Abb[1808]*m_z12*m_z1k_4;

  Box1Abb[1810]=10. + m_z12 + 57.*m_z1k;

  Box1Abb[1811]=1. + 12.*m_z12 + 34.*m_z1k - Box1Abb[1810]*m_z12*m_z1k;

  Box1Abb[1812]=4.*Box1Abb[178] + Box1Abb[1811]*m_z12;

  Box1Abb[1813]=-Box1Abb[1809]*m_x_2 + Box1Abb[1799]*m_x_3 + Box1Abb[1812]*m_x_4 + Box1Abb[1792]*m_x_5*m_z12 + Box1Abb[1790]*m_x_6*m_z12 + Box1Abb[1803]*m_x*m_z12*m_z1k + Box1Abb[1791]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_2;

  Box1Abb[1814]=-6. + 11.*m_z12 - 28.*m_z1k;

  Box1Abb[1815]=-8. + Box1Abb[1814]*m_z12;

  Box1Abb[1816]=2. - 2.*m_z12 + m_z12_2 + 3.*Box1Abb[79]*m_z1k + 4.*m_z1k_2;

  Box1Abb[1817]=71. - 3.*Box1Abb[772]*m_z12;

  Box1Abb[1818]=21. + m_z12;

  Box1Abb[1819]=15. + Box1Abb[1818]*m_z12;

  Box1Abb[1820]=-107. + Box1Abb[1819]*m_z12;

  Box1Abb[1821]=8. + 45.*m_z12;

  Box1Abb[1822]=-65. + Box1Abb[1821]*m_z12;

  Box1Abb[1823]=-12. + Box1Abb[1822]*m_z12;

  Box1Abb[1824]=2. + 11.*Box1Abb[141]*m_z12;

  Box1Abb[1825]=-32. + Box1Abb[1817]*m_z12 + 80.*m_z1k + Box1Abb[1820]*m_z12*m_z1k + 2.*Box1Abb[1823]*m_z1k_2 + 4.*Box1Abb[1824]*m_z1k_3 - 120.*m_z12*m_z1k_4;

  Box1Abb[1826]=12. + m_z1k;

  Box1Abb[1827]=1. + Box1Abb[1826]*m_z1k_2;

  Box1Abb[1828]=-4. + 7.*m_z1k;

  Box1Abb[1829]=-1. + m_z1k + 6.*Box1Abb[1828]*m_z1k_2;

  Box1Abb[1830]=-5. + 8.*m_z1k;

  Box1Abb[1831]=2. + Box1Abb[174]*Box1Abb[1830]*m_z1k;

  Box1Abb[1832]=3. + Box1Abb[1327]*m_z1k;

  Box1Abb[1833]=2.*Box1Abb[317]*pow(Box1Abb[68],4.) - Box1Abb[1831]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[1829]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[1827]*Box1Abb[68]*m_z12_3 + Box1Abb[1832]*m_z12_4*m_z1k;

  Box1Abb[1834]=-44. + 38.*m_z1k;

  Box1Abb[1835]=25. + Box1Abb[1834]*m_z12 - 3.*m_z12_2 + 4.*Box1Abb[834]*m_z1k;

  Box1Abb[1836]=26. + Box1Abb[1835]*m_z12 + 20.*m_z1k;

  Box1Abb[1837]=5. - 4.*m_z1k;

  Box1Abb[1838]=40. + 159.*m_z1k;

  Box1Abb[1839]=55. - Box1Abb[1838]*m_z1k;

  Box1Abb[1840]=9. + m_z1k + 10.*m_z1k_2;

  Box1Abb[1841]=-83. + 6.*Box1Abb[1840]*m_z1k;

  Box1Abb[1842]=4. + Box1Abb[1841]*m_z12 + Box1Abb[1839]*m_z12_2 + 2.*Box1Abb[691]*m_z12_3 + 4.*Box1Abb[1837]*m_z1k;

  Box1Abb[1843]=3. + 13.*m_z1k;

  Box1Abb[1844]=-13. + 5.*m_z1k - 80.*m_z1k_2;

  Box1Abb[1845]=6. + Box1Abb[1844]*m_z1k;

  Box1Abb[1846]=-8. + 13.*m_z1k;

  Box1Abb[1847]=59. + Box1Abb[1846]*m_z1k;

  Box1Abb[1848]=-2. + Box1Abb[1847]*m_z1k;

  Box1Abb[1849]=-2. + Box1Abb[1848]*m_z1k;

  Box1Abb[1850]=-39. + 46.*m_z1k;

  Box1Abb[1851]=25. + Box1Abb[1850]*m_z1k;

  Box1Abb[1852]=-39. + 2.*Box1Abb[1851]*m_z1k;

  Box1Abb[1853]=17. + Box1Abb[1852]*m_z1k;

  Box1Abb[1854]=-4.*Box1Abb[210]*Box1Abb[317]*pow(Box1Abb[68],2.) + Box1Abb[1853]*Box1Abb[68]*m_z12 + Box1Abb[1849]*m_z12_2 + Box1Abb[1845]*m_z12_3 - Box1Abb[1843]*m_z12_4*m_z1k;

  Box1Abb[1855]=Box1Abb[1833]*m_x + Box1Abb[1854]*m_x_2 + Box1Abb[1825]*m_x_3 + Box1Abb[1842]*m_x_4 + Box1Abb[1836]*m_x_5 + Box1Abb[1815]*m_x_6 + 8.*m_x_7*m_z12 + Box1Abb[1816]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[1856]=4. + m_x;

  Box1Abb[1857]=-7. + Box1Abb[1856]*m_x;

  Box1Abb[1858]=7. - 3.*m_x + Box1Abb[1857]*m_z12;

  Box1Abb[1859]=2.*Box1Abb[0] + Box1Abb[1858]*m_x;

  Box1Abb[1860]=6. + m_x - 8.*m_x_2;

  Box1Abb[1861]=-3. + Box1Abb[1860]*m_x;

  Box1Abb[1862]=6. + m_x;

  Box1Abb[1863]=-6. + Box1Abb[1862]*m_x;

  Box1Abb[1864]=-3. + Box1Abb[1863]*m_x;

  Box1Abb[1865]=1. + Box1Abb[1864]*m_x;

  Box1Abb[1866]=-2. + m_x + 13.*m_x_2;

  Box1Abb[1867]=Box1Abb[1861]*m_x + 2.*Box1Abb[1]*Box1Abb[1865]*m_z12 - pow(Box1Abb[1],2.)*Box1Abb[1866]*m_z12_2;

  Box1Abb[1868]=7. + 2.*m_z12;

  Box1Abb[1869]=2. - 7.*m_z12;

  Box1Abb[1870]=-6. + 4.*Box1Abb[1869]*m_z12;

  Box1Abb[1871]=-7. + m_z12 + 5.*m_z12_2;

  Box1Abb[1872]=13. + 8.*m_z12;

  Box1Abb[1873]=-9. + Box1Abb[1872]*m_z12_2;

  Box1Abb[1874]=m_x + Box1Abb[1873]*m_x_2 + Box1Abb[1870]*m_x_3 + Box1Abb[0]*Box1Abb[1868]*m_z12 - 2.*Box1Abb[1871]*m_x*m_z12 - 5.*m_x_4*m_z12;

  Box1Abb[1875]=3. + m_x + 4.*m_x_2 + 5.*m_x_3;

  Box1Abb[1876]=-11. + 6.*m_x;

  Box1Abb[1877]=-15. + Box1Abb[1876]*m_x;

  Box1Abb[1878]=-3.*m_x + 4.*Box1Abb[1875]*m_z12 + Box1Abb[1877]*m_z12_2 + 2.*Box1Abb[770]*m_z12_3;

  Box1Abb[1879]=26. + 25.*m_x;

  Box1Abb[1880]=14. + Box1Abb[1879]*m_x;

  Box1Abb[1881]=23. + 20.*m_x;

  Box1Abb[1882]=m_x - Box1Abb[1880]*m_z12 + Box1Abb[1881]*m_z12_2 - 4.*m_z12_3;

  Box1Abb[1883]=10. + 14.*m_x - 11.*m_z12;

  Box1Abb[1884]=pow(Box1Abb[1],2.)*Box1Abb[1859]*m_x - Box1Abb[1867]*m_z1k + Box1Abb[1874]*m_z1k_2 + Box1Abb[1878]*m_z1k_3 + Box1Abb[1882]*m_z1k_4 + Box1Abb[1883]*m_z12*m_z1k_5 - 3.*m_z12*m_z1k_6;

  Box1Abb[1885]=11. + m_z12 + 17.*m_z1k;

  Box1Abb[1886]=-13. + 2.*Box1Abb[1885]*m_z12 + 15.*m_z1k;

  Box1Abb[1887]=-4. + m_z1k + 2.*m_z1k_2;

  Box1Abb[1888]=1. + 2.*Box1Abb[1887]*m_z1k;

  Box1Abb[1889]=-14. + m_z1k + 6.*m_z1k_2;

  Box1Abb[1890]=3. + Box1Abb[1889]*m_z1k;

  Box1Abb[1891]=pow(Box1Abb[68],3.) - 3.*Box1Abb[317]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[1890]*Box1Abb[68]*m_z12_2 + Box1Abb[1888]*m_z12_3;

  Box1Abb[1892]=9. + 2.*Box1Abb[1776]*m_z1k;

  Box1Abb[1893]=-3. + 20.*Box1Abb[15]*m_z1k;

  Box1Abb[1894]=-2. + Box1Abb[1893]*m_z1k;

  Box1Abb[1895]=21. + 5.*m_z1k;

  Box1Abb[1896]=5. + Box1Abb[1895]*m_z1k;

  Box1Abb[1897]=13. + Box1Abb[1896]*m_z1k;

  Box1Abb[1898]=8. - Box1Abb[1897]*m_z12 + Box1Abb[1894]*m_z12_2 + Box1Abb[1892]*m_z12_3 + Box1Abb[873]*m_z1k;

  Box1Abb[1899]=-5. + 2.*m_z1k;

  Box1Abb[1900]=7. + 10.*m_z1k;

  Box1Abb[1901]=49. + 50.*m_z1k;

  Box1Abb[1902]=17. + Box1Abb[1901]*m_z1k;

  Box1Abb[1903]=-30. + Box1Abb[1902]*m_z12 + Box1Abb[1900]*m_z12_2 + 8.*Box1Abb[1899]*m_z1k;

  Box1Abb[1904]=13. + Box1Abb[1903]*m_z12 + m_z1k;

  Box1Abb[1905]=2. - 13.*m_z1k;

  Box1Abb[1906]=7. + 42.*m_z1k - 20.*m_z1k_2;

  Box1Abb[1907]=5. + Box1Abb[1906]*m_z1k;

  Box1Abb[1908]=2. + Box1Abb[873]*m_z1k;

  Box1Abb[1909]=5. + Box1Abb[1908]*m_z1k;

  Box1Abb[1910]=15. + 4.*Box1Abb[514]*m_z1k;

  Box1Abb[1911]=4. + Box1Abb[1910]*m_z1k;

  Box1Abb[1912]=7. + Box1Abb[1911]*m_z1k;

  Box1Abb[1913]=1. + m_z12 - Box1Abb[1912]*m_z12_2 + Box1Abb[1909]*m_z12_3 + Box1Abb[1905]*m_z1k + Box1Abb[1907]*m_z12*m_z1k;

  Box1Abb[1914]=10. + m_z1k + 8.*m_z1k_2;

  Box1Abb[1915]=-1. + Box1Abb[1914]*m_z1k;

  Box1Abb[1916]=-1. + Box1Abb[1915]*m_z1k;

  Box1Abb[1917]=13. + Box1Abb[1696]*m_z1k;

  Box1Abb[1918]=6. + Box1Abb[1917]*m_z1k;

  Box1Abb[1919]=1. + Box1Abb[1918]*m_z1k;

  Box1Abb[1920]=-23. - 4.*m_z1k + 22.*m_z1k_2;

  Box1Abb[1921]=-2. + Box1Abb[1920]*m_z1k;

  Box1Abb[1922]=3. + Box1Abb[1921]*m_z1k;

  Box1Abb[1923]=2. + Box1Abb[1922]*m_z1k;

  Box1Abb[1924]=Box1Abb[1919]*Box1Abb[68]*m_z12 + Box1Abb[1923]*m_z12_2 + Box1Abb[1916]*m_z12_3 + pow(Box1Abb[68],2.)*Box1Abb[873]*m_z1k;

  Box1Abb[1925]=Box1Abb[1924]*m_x_2 + Box1Abb[1913]*m_x_3 - Box1Abb[1898]*m_x_4 + Box1Abb[1904]*m_x_5 - Box1Abb[1886]*m_x_6*m_z12 + 4.*Box1Abb[14]*m_x_7*m_z12 - Box1Abb[1891]*Box1Abb[68]*m_x*m_z1k - pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[1926]=-1. + m_z12 + 6.*m_z12_2 - 5.*m_z12*m_z1k;

  Box1Abb[1927]=-13. + 7.*m_z12;

  Box1Abb[1928]=8. + Box1Abb[1927]*m_z12;

  Box1Abb[1929]=-4. + Box1Abb[12]*Box1Abb[69]*m_z12;

  Box1Abb[1930]=-1. + Box1Abb[653]*m_z12;

  Box1Abb[1931]=-21. + 2.*m_z12;

  Box1Abb[1932]=3. + Box1Abb[1931]*m_z12;

  Box1Abb[1933]=pow(Box1Abb[0],2.)*Box1Abb[326] - Box1Abb[0]*Box1Abb[1928]*m_z1k + Box1Abb[1929]*m_z1k_2 - Box1Abb[158]*Box1Abb[1930]*m_z1k_3 + Box1Abb[1932]*m_z1k_4 + 9.*m_z12*m_z1k_5;

  Box1Abb[1934]=13. + 2.*m_z12;

  Box1Abb[1935]=6.*Box1Abb[170] + Box1Abb[1934]*m_z12 + 20.*m_z12*m_z1k - 9.*m_z1k_2;

  Box1Abb[1936]=3.*Box1Abb[15] - Box1Abb[1935]*m_z12;

  Box1Abb[1937]=7. + 8.*m_z1k;

  Box1Abb[1938]=1. + 22.*Box1Abb[15]*m_z1k;

  Box1Abb[1939]=27. + Box1Abb[1103]*m_z1k;

  Box1Abb[1940]=17. + Box1Abb[1939]*m_z1k;

  Box1Abb[1941]=8. - Box1Abb[1940]*m_z12 + Box1Abb[1938]*m_z12_2 + Box1Abb[1937]*m_z12_3 - 2.*m_z1k_2;

  Box1Abb[1942]=1. + Box1Abb[210]*Box1Abb[68]*m_z1k;

  Box1Abb[1943]=5. + 7.*m_z1k;

  Box1Abb[1944]=9. + Box1Abb[1943]*m_z1k;

  Box1Abb[1945]=5. + 4.*m_z1k;

  Box1Abb[1946]=5. - 2.*Box1Abb[1945]*m_z1k;

  Box1Abb[1947]=16. + Box1Abb[1946]*m_z1k;

  Box1Abb[1948]=9. + Box1Abb[711]*m_z1k;

  Box1Abb[1949]=-2. + Box1Abb[1948]*m_z1k;

  Box1Abb[1950]=-5. + Box1Abb[1949]*m_z1k;

  Box1Abb[1951]=-2.*Box1Abb[1942] + Box1Abb[1950]*m_z12 + Box1Abb[1947]*m_z12_2 - Box1Abb[1944]*m_z12_3;

  Box1Abb[1952]=1. - 3.*m_z1k + m_z1k_3;

  Box1Abb[1953]=8. - 5.*m_z1k;

  Box1Abb[1954]=-7. + Box1Abb[1953]*m_z1k;

  Box1Abb[1955]=6. + Box1Abb[1954]*m_z1k;

  Box1Abb[1956]=-2. + Box1Abb[1955]*m_z1k;

  Box1Abb[1957]=-pow(Box1Abb[68],3.) + Box1Abb[1956]*m_z12 + Box1Abb[1952]*m_z12_2;

  Box1Abb[1958]=Box1Abb[1957]*Box1Abb[4]*Box1Abb[68]*m_x + Box1Abb[1933]*m_x_2 + Box1Abb[1951]*m_x_3 + Box1Abb[1941]*m_x_4 + Box1Abb[1936]*m_x_5 + Box1Abb[1926]*m_x_6 + m_x_7*m_z12 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[1959]=10. + 9.*m_z12 + 4.*m_z1k;

  Box1Abb[1960]=-16. + 3.*Box1Abb[1959]*m_z12;

  Box1Abb[1961]=2. + 3.*m_z1k_2;

  Box1Abb[1962]=2. - 3.*m_z1k + 9.*m_z1k_2;

  Box1Abb[1963]=Box1Abb[1961]*Box1Abb[68]*m_z12 + Box1Abb[1962]*m_z12_2 - 2.*Box1Abb[317]*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[1964]=38. + 31.*m_z1k;

  Box1Abb[1965]=57. + 34.*m_z1k;

  Box1Abb[1966]=3. + 3.*Box1Abb[1964]*m_z12 + 7.*m_z12_2 + 2.*Box1Abb[1965]*m_z1k;

  Box1Abb[1967]=-6.*Box1Abb[230] + Box1Abb[1966]*m_z12;

  Box1Abb[1968]=21. + 38.*m_z1k;

  Box1Abb[1969]=6. + Box1Abb[1968]*m_z1k;

  Box1Abb[1970]=304. + 155.*m_z1k;

  Box1Abb[1971]=165. + Box1Abb[1970]*m_z1k;

  Box1Abb[1972]=47. + 78.*m_z1k;

  Box1Abb[1973]=-59. + 2.*Box1Abb[1972]*m_z1k;

  Box1Abb[1974]=-145. + Box1Abb[1973]*m_z1k;

  Box1Abb[1975]=-2.*Box1Abb[1969] + Box1Abb[1974]*m_z12 + Box1Abb[1971]*m_z12_2 + 7.*Box1Abb[1151]*m_z12_3;

  Box1Abb[1976]=-17. + 6.*m_z1k;

  Box1Abb[1977]=7. + Box1Abb[1976]*m_z1k;

  Box1Abb[1978]=113. + 58.*m_z1k;

  Box1Abb[1979]=44. + Box1Abb[1978]*m_z1k;

  Box1Abb[1980]=-58. + 45.*m_z1k;

  Box1Abb[1981]=-82. + Box1Abb[1980]*Box1Abb[431]*m_z1k;

  Box1Abb[1982]=-179. + 2.*Box1Abb[1981]*m_z1k;

  Box1Abb[1983]=274. + 229.*m_z1k;

  Box1Abb[1984]=190. + Box1Abb[1983]*m_z1k;

  Box1Abb[1985]=87. + Box1Abb[1984]*m_z1k;

  Box1Abb[1986]=56. + Box1Abb[1982]*m_z12 + Box1Abb[1985]*m_z12_2 + Box1Abb[1979]*m_z12_3 - 4.*Box1Abb[1977]*m_z1k + 3.*m_z12_4*m_z1k;

  Box1Abb[1987]=-7. + 18.*m_z1k;

  Box1Abb[1988]=9. + Box1Abb[1987]*m_z1k_2;

  Box1Abb[1989]=5. + 6.*m_z1k;

  Box1Abb[1990]=-15. + 2.*Box1Abb[1989]*m_z1k;

  Box1Abb[1991]=10. + Box1Abb[1990]*m_z1k;

  Box1Abb[1992]=2. + Box1Abb[1991]*m_z1k;

  Box1Abb[1993]=-78. + 49.*m_z1k;

  Box1Abb[1994]=34. + Box1Abb[1993]*m_z1k;

  Box1Abb[1995]=-7. + Box1Abb[1994]*m_z1k;

  Box1Abb[1996]=-4. + Box1Abb[1995]*m_z1k;

  Box1Abb[1997]=-77. + 79.*m_z1k;

  Box1Abb[1998]=24. + Box1Abb[1997]*m_z1k;

  Box1Abb[1999]=10. + Box1Abb[1998]*m_z1k;

  Box1Abb[2000]=-2. + Box1Abb[1999]*m_z1k;

  Box1Abb[2001]=-Box1Abb[1992]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[1996]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[2000]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[317]*pow(Box1Abb[68],4.)*m_z1k + Box1Abb[1988]*m_z12_4*m_z1k;

  Box1Abb[2002]=11. + 6.*m_z1k;

  Box1Abb[2003]=80. + 47.*m_z1k;

  Box1Abb[2004]=107. + 2.*Box1Abb[2003]*m_z1k;

  Box1Abb[2005]=34. + Box1Abb[2004]*m_z1k;

  Box1Abb[2006]=7. + 6.*m_z1k;

  Box1Abb[2007]=-5. + Box1Abb[2006]*m_z1k;

  Box1Abb[2008]=-9. + Box1Abb[2007]*m_z1k;

  Box1Abb[2009]=9. + Box1Abb[2008]*m_z1k;

  Box1Abb[2010]=-96. + 281.*m_z1k;

  Box1Abb[2011]=-76. + Box1Abb[2010]*m_z1k;

  Box1Abb[2012]=-75. + Box1Abb[2011]*m_z1k;

  Box1Abb[2013]=-4. + Box1Abb[2012]*m_z1k;

  Box1Abb[2014]=-139. + 50.*m_z1k;

  Box1Abb[2015]=33. + Box1Abb[2014]*m_z1k;

  Box1Abb[2016]=37. + Box1Abb[2015]*m_z1k;

  Box1Abb[2017]=13. + 2.*Box1Abb[2016]*m_z1k;

  Box1Abb[2018]=-65. + Box1Abb[2017]*m_z1k;

  Box1Abb[2019]=4.*Box1Abb[2009] + Box1Abb[2018]*m_z12 + Box1Abb[2013]*m_z12_2 + Box1Abb[2005]*m_z12_3 + Box1Abb[2002]*m_z12_4*m_z1k;

  Box1Abb[2020]=25. + 12.*m_z1k;

  Box1Abb[2021]=15. + Box1Abb[2020]*m_z1k;

  Box1Abb[2022]=3. + Box1Abb[1128]*m_z1k;

  Box1Abb[2023]=42. - 131.*m_z1k;

  Box1Abb[2024]=13. + Box1Abb[2023]*m_z1k;

  Box1Abb[2025]=-21. + Box1Abb[2024]*m_z1k;

  Box1Abb[2026]=-13. + Box1Abb[2025]*m_z1k;

  Box1Abb[2027]=-49. + 6.*m_z1k;

  Box1Abb[2028]=-29. + 2.*Box1Abb[2027]*m_z1k;

  Box1Abb[2029]=-40. + Box1Abb[2028]*Box1Abb[68]*m_z1k;

  Box1Abb[2030]=2. + Box1Abb[2029]*m_z1k;

  Box1Abb[2031]=-362. + 191.*m_z1k;

  Box1Abb[2032]=150. + Box1Abb[2031]*m_z1k;

  Box1Abb[2033]=61. + Box1Abb[2032]*m_z1k;

  Box1Abb[2034]=59. - Box1Abb[2033]*m_z1k;

  Box1Abb[2035]=17. + Box1Abb[2034]*m_z1k;

  Box1Abb[2036]=-2.*Box1Abb[2022]*pow(Box1Abb[68],2.) - Box1Abb[2030]*Box1Abb[68]*m_z12 + Box1Abb[2035]*m_z12_2 + Box1Abb[2026]*m_z12_3 - Box1Abb[2021]*m_z12_4*m_z1k;

  Box1Abb[2037]=Box1Abb[2001]*m_x + Box1Abb[2036]*m_x_2 + Box1Abb[2019]*m_x_3 - Box1Abb[1986]*m_x_4 + Box1Abb[1975]*m_x_5 - Box1Abb[1967]*m_x_6 + Box1Abb[1960]*m_x_7 - Box1Abb[1963]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[2038]=-1. + m_z1k - 2.*m_z12*m_z1k;

  Box1Abb[2039]=-1. + m_z12 + 6.*m_z1k - 2.*Box1Abb[15]*m_z12*m_z1k;

  Box1Abb[2040]=Box1Abb[2039]*m_x_2 + Box1Abb[1408]*m_x_3 + m_x_4*m_z12 + Box1Abb[2038]*m_x*m_z1k + pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[2041]=1. + m_z12 - m_z1k;

  Box1Abb[2042]=15. + Box1Abb[1422]*Box1Abb[488]*m_z12;

  Box1Abb[2043]=64. + 9.*m_z12;

  Box1Abb[2044]=18. + Box1Abb[2043]*m_z12;

  Box1Abb[2045]=-8. + Box1Abb[2044]*m_z12;

  Box1Abb[2046]=120. + 23.*m_z12;

  Box1Abb[2047]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box1Abb[2042]*m_z1k + Box1Abb[2045]*m_z1k_2 + Box1Abb[2046]*m_z12*m_z1k_3;

  Box1Abb[2048]=10. + 3.*m_z12;

  Box1Abb[2049]=-10. + Box1Abb[2048]*m_z12;

  Box1Abb[2050]=12. + Box1Abb[2049]*m_z12;

  Box1Abb[2051]=19. + 10.*m_z12;

  Box1Abb[2052]=32. + Box1Abb[2051]*m_z12;

  Box1Abb[2053]=-15. + Box1Abb[2052]*m_z12;

  Box1Abb[2054]=9. - 4.*Box1Abb[707]*m_z12;

  Box1Abb[2055]=2. + Box1Abb[2054]*m_z12;

  Box1Abb[2056]=-80. + 23.*m_z12;

  Box1Abb[2057]=pow(Box1Abb[0],2.) - Box1Abb[2050]*m_z1k - Box1Abb[2053]*m_z1k_2 + 2.*Box1Abb[2055]*m_z1k_3 + Box1Abb[2056]*m_z12*m_z1k_4;

  Box1Abb[2058]=-4. + 9.*m_z1k;

  Box1Abb[2059]=2. + Box1Abb[2058]*m_z12 + 20.*m_z1k;

  Box1Abb[2060]=28. + m_z12;

  Box1Abb[2061]=20. + Box1Abb[2060]*m_z12;

  Box1Abb[2062]=80. + 33.*m_z12;

  Box1Abb[2063]=6. - 6.*m_z12 + Box1Abb[2061]*m_z1k + Box1Abb[2062]*m_z1k_2;

  Box1Abb[2064]=-1. + Box1Abb[2063]*m_z12 - 4.*m_z1k;

  Box1Abb[2065]=-10. + 9.*m_z1k;

  Box1Abb[2066]=-1. + 9.*m_z1k;

  Box1Abb[2067]=2. + Box1Abb[2066]*m_z1k;

  Box1Abb[2068]=2.*pow(Box1Abb[68],2.) + Box1Abb[2067]*m_z12_2 + Box1Abb[2065]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[2069]=-1. + 2.*m_z1k_2;

  Box1Abb[2070]=64. - 33.*m_z1k;

  Box1Abb[2071]=-19. + Box1Abb[2070]*m_z1k;

  Box1Abb[2072]=12. + Box1Abb[2071]*m_z1k;

  Box1Abb[2073]=7. - 4.*m_z1k;

  Box1Abb[2074]=3. + 2.*Box1Abb[2073]*m_z1k;

  Box1Abb[2075]=1. + Box1Abb[2074]*m_z1k;

  Box1Abb[2076]=pow(Box1Abb[68],2.) + Box1Abb[2075]*m_z12_3 + 10.*Box1Abb[2069]*Box1Abb[68]*m_z12*m_z1k + Box1Abb[2072]*m_z12_2*m_z1k;

  Box1Abb[2077]=Box1Abb[2057]*m_x_3 + Box1Abb[2047]*m_x_4 - Box1Abb[2064]*m_x_5 + Box1Abb[2059]*m_x_6*m_z12 + m_x_7*m_z12_2 + Box1Abb[2076]*m_x_2*m_z1k + Box1Abb[2068]*Box1Abb[68]*m_x*m_z12*m_z1k_2 - Box1Abb[2041]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_3;

  Box1Abb[2078]=-2. + m_z12 + 2.*m_z1k;

  Box1Abb[2079]=-1. + m_z1k + 8.*m_z1k_2;

  Box1Abb[2080]=1. + m_z1k + 8.*m_z1k_2;

  Box1Abb[2081]=1. + 4.*Box1Abb[66]*m_z1k;

  Box1Abb[2082]=-2.*pow(Box1Abb[68],3.) + 2.*Box1Abb[2079]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2081]*Box1Abb[68]*m_z12_2 + Box1Abb[2080]*m_z12_3;

  Box1Abb[2083]=2. - 4.*m_z1k;

  Box1Abb[2084]=3. + 16.*m_z1k;

  Box1Abb[2085]=-1. + Box1Abb[2083]*m_z12 + Box1Abb[2084]*m_z1k;

  Box1Abb[2086]=Box1Abb[2085]*m_z12 + m_z1k;

  Box1Abb[2087]=7. + 2.*m_z1k;

  Box1Abb[2088]=3. + 8.*m_z1k;

  Box1Abb[2089]=1. + Box1Abb[68]*m_z1k;

  Box1Abb[2090]=-3. + 20.*Box1Abb[431]*m_z1k;

  Box1Abb[2091]=23. + 4.*Box1Abb[1837]*m_z1k;

  Box1Abb[2092]=-4. + Box1Abb[2091]*m_z1k;

  Box1Abb[2093]=-2. - 2.*Box1Abb[2089]*Box1Abb[2090]*m_z12 + Box1Abb[2092]*m_z12_2 + 2.*Box1Abb[2087]*m_z1k + Box1Abb[2088]*m_z12_3*m_z1k;

  Box1Abb[2094]=-2. + Box1Abb[228]*m_z1k;

  Box1Abb[2095]=12. + m_z1k + 35.*m_z1k_2;

  Box1Abb[2096]=-3. + Box1Abb[2095]*m_z1k;

  Box1Abb[2097]=1. + 2.*Box1Abb[2096]*m_z12 - 3.*Box1Abb[2094]*m_z12_2 + 4.*Box1Abb[178]*m_z1k - m_z12_3*m_z1k;

  Box1Abb[2098]=-1. + 4.*Box1Abb[210]*m_z1k;

  Box1Abb[2099]=-11. + 25.*m_z1k;

  Box1Abb[2100]=3. + Box1Abb[2099]*m_z1k;

  Box1Abb[2101]=-11. + Box1Abb[2100]*m_z1k;

  Box1Abb[2102]=1. + Box1Abb[2101]*m_z1k;

  Box1Abb[2103]=-24. + 43.*m_z1k;

  Box1Abb[2104]=18. + Box1Abb[2103]*m_z1k;

  Box1Abb[2105]=-6. + Box1Abb[2104]*m_z1k;

  Box1Abb[2106]=1. + Box1Abb[2105]*m_z1k;

  Box1Abb[2107]=-Box1Abb[2098]*pow(Box1Abb[68],2.) + 2.*Box1Abb[2102]*Box1Abb[68]*m_z12 + Box1Abb[2106]*m_z12_2 - Box1Abb[1692]*m_z12_3*m_z1k;

  Box1Abb[2108]=Box1Abb[2107]*m_x_2 + Box1Abb[2093]*m_x_3 + Box1Abb[2097]*m_x_4 - 2.*Box1Abb[2086]*m_x_5 + Box1Abb[1459]*m_x_6*m_z12 - Box1Abb[2082]*Box1Abb[68]*m_x*m_z1k + Box1Abb[2078]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[2109]=-2. - Box1Abb[170]*m_z12 + m_z1k + m_z1k_2;

  Box1Abb[2110]=2. - 21.*m_z12;

  Box1Abb[2111]=5. + m_z12_2;

  Box1Abb[2112]=6. + Box1Abb[1509]*m_z12;

  Box1Abb[2113]=18. + Box1Abb[2110]*m_z12 + 8.*m_z1k - 3.*Box1Abb[2111]*m_z12*m_z1k + 2.*Box1Abb[2112]*m_z1k_2 + 40.*m_z12*m_z1k_3;

  Box1Abb[2114]=2. - 5.*m_z12 + 20.*m_z1k;

  Box1Abb[2115]=4. + Box1Abb[2114]*m_z12;

  Box1Abb[2116]=17. - 3.*m_z1k;

  Box1Abb[2117]=5. + Box1Abb[2116]*m_z12 + 4.*m_z1k - 40.*m_z1k_2;

  Box1Abb[2118]=-6.*Box1Abb[873] + Box1Abb[2117]*m_z12;

  Box1Abb[2119]=-7. + 11.*m_z1k;

  Box1Abb[2120]=7. - 2.*m_z1k + 4.*m_z1k_2;

  Box1Abb[2121]=2. + Box1Abb[2120]*m_z1k;

  Box1Abb[2122]=-2. + 3.*Box1Abb[785]*m_z1k;

  Box1Abb[2123]=2. + Box1Abb[2122]*m_z1k;

  Box1Abb[2124]=Box1Abb[2121]*pow(Box1Abb[68],2.) + Box1Abb[2123]*Box1Abb[68]*m_z12 + Box1Abb[15]*Box1Abb[2119]*m_z12_2*m_z1k;

  Box1Abb[2125]=-8. + 7.*m_z1k;

  Box1Abb[2126]=-1. - 8.*m_z1k + 46.*m_z1k_2;

  Box1Abb[2127]=-11. + Box1Abb[2126]*m_z1k;

  Box1Abb[2128]=5. + 4.*Box1Abb[480]*m_z1k;

  Box1Abb[2129]=4. + Box1Abb[2128]*m_z1k;

  Box1Abb[2130]=7. + Box1Abb[2129]*m_z1k;

  Box1Abb[2131]=2.*Box1Abb[170]*Box1Abb[431]*Box1Abb[68] + Box1Abb[2130]*m_z12 + Box1Abb[2127]*m_z12_2 + Box1Abb[2125]*m_z12_3*m_z1k;

  Box1Abb[2132]=-Box1Abb[2131]*m_x_2 + Box1Abb[2113]*m_x_3 + Box1Abb[2118]*m_x_4 + Box1Abb[2115]*m_x_5 + Box1Abb[2124]*m_x*m_z12 - 4.*m_x_6*m_z12 + Box1Abb[2109]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k;

  Box1Abb[2133]=2. + 7.*m_x;

  Box1Abb[2134]=1. - m_x + 10.*m_x_2 + Box1Abb[1]*Box1Abb[2133]*Box1Abb[766]*m_z12 - 3.*pow(Box1Abb[1],2.)*m_z12_2;

  Box1Abb[2135]=5. - 26.*m_x;

  Box1Abb[2136]=-7. + Box1Abb[2135]*m_x;

  Box1Abb[2137]=17. - 8.*m_x;

  Box1Abb[2138]=2. + Box1Abb[2137]*m_x;

  Box1Abb[2139]=-11. + 2.*Box1Abb[2138]*m_x;

  Box1Abb[2140]=4. + Box1Abb[2139]*m_x;

  Box1Abb[2141]=-59. + 40.*m_x;

  Box1Abb[2142]=7. + Box1Abb[2141]*m_x;

  Box1Abb[2143]=3. + m_x + Box1Abb[2142]*m_x_2;

  Box1Abb[2144]=8. + Box1Abb[2143]*m_x;

  Box1Abb[2145]=1. + 2.*m_x + 4.*m_x_2;

  Box1Abb[2146]=2. + Box1Abb[2136]*m_x - 7.*m_z12 + Box1Abb[2140]*m_x*m_z12 + Box1Abb[2144]*m_z12_2 - 3.*pow(Box1Abb[1],2.)*Box1Abb[2145]*m_z12_3;

  Box1Abb[2147]=38. - 15.*m_z12;

  Box1Abb[2148]=-7. + 5.*m_z12;

  Box1Abb[2149]=5. + 2.*Box1Abb[2148]*m_z12;

  Box1Abb[2150]=8. + Box1Abb[676]*m_z12;

  Box1Abb[2151]=19. + 7.*m_z12;

  Box1Abb[2152]=23. - Box1Abb[2151]*m_z12;

  Box1Abb[2153]=1. + Box1Abb[2152]*m_z12;

  Box1Abb[2154]=-37. + 9.*m_z12;

  Box1Abb[2155]=74. + Box1Abb[2154]*m_z12;

  Box1Abb[2156]=-4. + Box1Abb[2155]*m_z12;

  Box1Abb[2157]=2.*Box1Abb[0]*Box1Abb[2149]*m_x + Box1Abb[2153]*m_x_2 + 4.*Box1Abb[2150]*m_x_3 + Box1Abb[2156]*m_x_4 - 2.*pow(Box1Abb[0],2.)*m_z12 + 2.*Box1Abb[2147]*m_x_5*m_z12;

  Box1Abb[2158]=17. + 12.*m_x;

  Box1Abb[2159]=18. + Box1Abb[2158]*m_x;

  Box1Abb[2160]=99. + 70.*m_x;

  Box1Abb[2161]=70. + Box1Abb[2160]*m_x;

  Box1Abb[2162]=61. + 2.*Box1Abb[2161]*m_x;

  Box1Abb[2163]=97. - 30.*m_x;

  Box1Abb[2164]=116. + Box1Abb[2163]*m_x;

  Box1Abb[2165]=85. + Box1Abb[2164]*m_x;

  Box1Abb[2166]=-20. + Box1Abb[2165]*m_x;

  Box1Abb[2167]=-4. + 3.*m_x;

  Box1Abb[2168]=-33. + 4.*Box1Abb[2167]*m_x;

  Box1Abb[2169]=8. + Box1Abb[2168]*m_x;

  Box1Abb[2170]=Box1Abb[2159]*m_x + 12.*m_z12 - Box1Abb[2162]*m_x*m_z12 + Box1Abb[2166]*m_z12_2 + Box1Abb[2169]*m_z12_3;

  Box1Abb[2171]=7. + 6.*m_x;

  Box1Abb[2172]=2. + 5.*m_x;

  Box1Abb[2173]=-7. + 2.*Box1Abb[1269]*Box1Abb[2172]*m_x;

  Box1Abb[2174]=-49. + 54.*m_x;

  Box1Abb[2175]=-50. + Box1Abb[2174]*m_x;

  Box1Abb[2176]=31. + Box1Abb[2175]*m_x;

  Box1Abb[2177]=8. - 15.*m_x;

  Box1Abb[2178]=-7. + Box1Abb[2177]*m_x;

  Box1Abb[2179]=-2.*Box1Abb[2171]*m_x + 4.*Box1Abb[2173]*m_z12 + Box1Abb[2176]*m_z12_2 + Box1Abb[2178]*m_z12_3;

  Box1Abb[2180]=32. - Box1Abb[1934]*m_z12;

  Box1Abb[2181]=2. + 15.*m_z12;

  Box1Abb[2182]=4. + Box1Abb[2181]*m_z12;

  Box1Abb[2183]=Box1Abb[2182]*m_x + Box1Abb[2180]*m_z12 - 8.*Box1Abb[1032]*m_x_2*m_z12;

  Box1Abb[2184]=-18. - 2.*Box1Abb[1649]*m_x + Box1Abb[472]*m_z12;

  Box1Abb[2185]=pow(Box1Abb[1],3.)*pow(Box1Abb[185],2.)*m_x_3 - Box1Abb[1]*Box1Abb[185]*Box1Abb[2134]*m_x_2*m_z1k + Box1Abb[2146]*m_x*m_z1k_2 + Box1Abb[2157]*m_z1k_3 + Box1Abb[2170]*m_z1k_4 + Box1Abb[2179]*m_z1k_5 + Box1Abb[2183]*m_z1k_6 + Box1Abb[2184]*m_z12*m_z1k_7 + Box1Abb[261]*m_z12*m_z1k_8;

  Box1Abb[2186]=m_x + Box1Abb[1]*m_x*m_z12;

  Box1Abb[2187]=-4. + m_x + 12.*m_x_2;

  Box1Abb[2188]=-1. + 3.*m_x + 8.*m_x_2 + Box1Abb[1]*Box1Abb[2187]*m_z12 - 3.*pow(Box1Abb[1],2.)*m_z12_2;

  Box1Abb[2189]=5. - 2.*m_x;

  Box1Abb[2190]=10. + Box1Abb[2189]*m_x;

  Box1Abb[2191]=9. + 2.*Box1Abb[2190]*m_x;

  Box1Abb[2192]=-9. + Box1Abb[2191]*m_x;

  Box1Abb[2193]=-15. + 6.*m_x + 4.*m_x_2;

  Box1Abb[2194]=-4. + Box1Abb[2193]*m_x;

  Box1Abb[2195]=-32. + Box1Abb[2194]*m_x;

  Box1Abb[2196]=32. + Box1Abb[2195]*m_x;

  Box1Abb[2197]=-11. + 31.*m_x;

  Box1Abb[2198]=5. + Box1Abb[2197]*m_x;

  Box1Abb[2199]=-33. + Box1Abb[2198]*m_x;

  Box1Abb[2200]=-2. + Box1Abb[2199]*m_x;

  Box1Abb[2201]=1. + 10.*Box1Abb[183]*m_x;

  Box1Abb[2202]=Box1Abb[2192]*m_x + m_z12 + Box1Abb[2196]*m_x*m_z12 - Box1Abb[1]*Box1Abb[2200]*m_z12_2 + pow(Box1Abb[1],2.)*Box1Abb[2201]*m_z12_3;

  Box1Abb[2203]=-1. + 8.*m_x;

  Box1Abb[2204]=2. + Box1Abb[2203]*m_x;

  Box1Abb[2205]=25. + 2.*Box1Abb[2204]*m_x;

  Box1Abb[2206]=7. + 4.*m_x;

  Box1Abb[2207]=34. + 3.*Box1Abb[2206]*m_x;

  Box1Abb[2208]=22. + Box1Abb[2207]*m_x;

  Box1Abb[2209]=77. + 2.*Box1Abb[2208]*m_x;

  Box1Abb[2210]=6. + Box1Abb[2209]*m_x;

  Box1Abb[2211]=-13. + 17.*m_x;

  Box1Abb[2212]=32. + Box1Abb[2211]*m_x;

  Box1Abb[2213]=5. + Box1Abb[183]*Box1Abb[2212]*m_x;

  Box1Abb[2214]=-1. + 3.*m_x;

  Box1Abb[2215]=8. + Box1Abb[2214]*m_x;

  Box1Abb[2216]=2. + Box1Abb[2215]*m_x;

  Box1Abb[2217]=-Box1Abb[2205]*m_x + Box1Abb[2210]*m_z12 - 2.*Box1Abb[2213]*m_z12_2 + 2.*Box1Abb[2216]*m_z12_3;

  Box1Abb[2218]=31. + 26.*m_x + 24.*m_x_2;

  Box1Abb[2219]=9. + 5.*m_x;

  Box1Abb[2220]=12. + Box1Abb[2219]*m_x;

  Box1Abb[2221]=9. + Box1Abb[2220]*m_x;

  Box1Abb[2222]=3. + 2.*Box1Abb[2221]*m_x;

  Box1Abb[2223]=61. + 31.*m_x;

  Box1Abb[2224]=82. + Box1Abb[2223]*m_x;

  Box1Abb[2225]=28. + Box1Abb[2224]*m_x;

  Box1Abb[2226]=11. + 10.*m_x;

  Box1Abb[2227]=Box1Abb[2218]*m_x - 6.*Box1Abb[2222]*m_z12 + Box1Abb[2225]*m_z12_2 - Box1Abb[2226]*m_z12_3;

  Box1Abb[2228]=9. + 8.*m_x;

  Box1Abb[2229]=33. + 20.*m_x;

  Box1Abb[2230]=28. + Box1Abb[2229]*m_x;

  Box1Abb[2231]=8. + Box1Abb[2230]*m_x;

  Box1Abb[2232]=15. + 8.*m_x;

  Box1Abb[2233]=44. + 5.*Box1Abb[2232]*m_x;

  Box1Abb[2234]=2. + m_x;

  Box1Abb[2235]=-2.*Box1Abb[2228]*m_x + 4.*Box1Abb[2231]*m_z12 - Box1Abb[2233]*m_z12_2 + 7.*Box1Abb[2234]*m_z12_3;

  Box1Abb[2236]=11. + 26.*m_x + 20.*m_x_2;

  Box1Abb[2237]=34. + 33.*m_x;

  Box1Abb[2238]=4.*m_x - 3.*Box1Abb[2236]*m_z12 + Box1Abb[2237]*m_z12_2 - 6.*m_z12_3;

  Box1Abb[2239]=9. + 12.*m_x - 5.*m_z12;

  Box1Abb[2240]=pow(Box1Abb[1],3.)*pow(Box1Abb[2186],2.) - Box1Abb[1]*Box1Abb[185]*Box1Abb[2188]*m_x*m_z1k - Box1Abb[2202]*m_z1k_2 + Box1Abb[2217]*m_z1k_3 + Box1Abb[2227]*m_z1k_4 + Box1Abb[2235]*m_z1k_5 + Box1Abb[2238]*m_z1k_6 + 2.*Box1Abb[2239]*m_z12*m_z1k_7 - 4.*m_z12*m_z1k_8;

  Box1Abb[2241]=14. + m_z12;

  Box1Abb[2242]=-16. + Box1Abb[2241]*m_z12;

  Box1Abb[2243]=20. + Box1Abb[2242]*m_z12;

  Box1Abb[2244]=52. + m_z12;

  Box1Abb[2245]=-26. + Box1Abb[2244]*m_z12;

  Box1Abb[2246]=3. - 22.*m_z12 + 26.*m_z12_2 + Box1Abb[2243]*m_z1k - Box1Abb[2245]*m_z1k_2 + 18.*m_z12*m_z1k_3;

  Box1Abb[2247]=-14. + Box1Abb[2060]*m_z12;

  Box1Abb[2248]=-24. + Box1Abb[2247]*m_z12;

  Box1Abb[2249]=-98. + m_z12;

  Box1Abb[2250]=44. + Box1Abb[2249]*m_z12;

  Box1Abb[2251]=8. + 6.*Box1Abb[858]*m_z12 + m_z1k + Box1Abb[12]*Box1Abb[261]*m_z12*m_z1k - Box1Abb[2248]*m_z1k_2 + Box1Abb[2250]*m_z1k_3 + 10.*m_z12*m_z1k_4;

  Box1Abb[2252]=14. + 5.*m_z1k;

  Box1Abb[2253]=-6. + Box1Abb[2252]*m_z12 + 10.*Box1Abb[68]*m_z1k;

  Box1Abb[2254]=Box1Abb[2253]*m_z12 + 6.*m_z1k;

  Box1Abb[2255]=-5. + m_z1k + 5.*m_z1k_2;

  Box1Abb[2256]=1. + Box1Abb[2255]*m_z1k;

  Box1Abb[2257]=7. + m_z1k;

  Box1Abb[2258]=-4. + Box1Abb[2257]*m_z1k;

  Box1Abb[2259]=2. + Box1Abb[2258]*m_z1k;

  Box1Abb[2260]=8. + 11.*m_z1k;

  Box1Abb[2261]=-12. + Box1Abb[2260]*m_z1k;

  Box1Abb[2262]=4. + Box1Abb[2261]*m_z1k;

  Box1Abb[2263]=2.*Box1Abb[2256]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2262]*Box1Abb[68]*m_z12_2 + Box1Abb[2259]*m_z12_3 - 2.*pow(Box1Abb[68],3.)*m_z1k;

  Box1Abb[2264]=-2. + 6.*m_z1k - 26.*m_z1k_3 + 30.*m_z1k_4 + m_z1k_5 - 9.*m_z1k_6;

  Box1Abb[2265]=-9. + 2.*m_z1k;

  Box1Abb[2266]=-3. + 2.*Box1Abb[2265]*m_z1k;

  Box1Abb[2267]=7. + Box1Abb[2266]*m_z1k;

  Box1Abb[2268]=-3. + 2.*Box1Abb[1041]*m_z1k;

  Box1Abb[2269]=2. + Box1Abb[2268]*m_z1k;

  Box1Abb[2270]=47. + m_z1k;

  Box1Abb[2271]=-2. + Box1Abb[2270]*m_z1k;

  Box1Abb[2272]=-10. + Box1Abb[2271]*m_z1k;

  Box1Abb[2273]=2. + Box1Abb[2272]*m_z1k;

  Box1Abb[2274]=Box1Abb[2269]*pow(Box1Abb[68],2.) + 2.*Box1Abb[2264]*m_z12 - Box1Abb[2273]*Box1Abb[68]*m_z12_2 + Box1Abb[2267]*m_z12_3*m_z1k;

  Box1Abb[2275]=5. - 9.*m_z1k;

  Box1Abb[2276]=-5. + 4.*Box1Abb[2275]*m_z1k;

  Box1Abb[2277]=12. + Box1Abb[2276]*m_z1k;

  Box1Abb[2278]=68. - 9.*m_z1k;

  Box1Abb[2279]=15. + Box1Abb[2278]*m_z1k;

  Box1Abb[2280]=8. + Box1Abb[2279]*m_z1k;

  Box1Abb[2281]=-11. + Box1Abb[2280]*m_z1k;

  Box1Abb[2282]=36. + 5.*m_z1k;

  Box1Abb[2283]=-39. + Box1Abb[2282]*m_z1k;

  Box1Abb[2284]=-11. + Box1Abb[2283]*m_z1k;

  Box1Abb[2285]=-9. + Box1Abb[2284]*m_z1k;

  Box1Abb[2286]=9. + Box1Abb[2285]*m_z1k;

  Box1Abb[2287]=-7. + 2.*Box1Abb[2286]*m_z12 + Box1Abb[2281]*m_z12_2 + Box1Abb[2277]*m_z1k - Box1Abb[1369]*m_z12_3*m_z1k;

  Box1Abb[2288]=Box1Abb[2274]*m_x_2 + Box1Abb[2287]*m_x_3 + Box1Abb[2251]*m_x_4 - Box1Abb[2246]*m_x_5 + Box1Abb[2254]*m_x_6 - Box1Abb[724]*m_x_7*m_z12 + Box1Abb[2263]*Box1Abb[68]*m_x*m_z1k - Box1Abb[4]*pow(Box1Abb[68],3.)*Box1Abb[739]*m_z12*m_z1k_3;

  Box1Abb[2289]=2. + m_z12 - 5.*m_z1k + 2.*m_z12*m_z1k + 3.*m_z1k_2;

  Box1Abb[2290]=10. + m_z12 + 12.*m_z1k;

  Box1Abb[2291]=-2. + Box1Abb[2290]*m_z12;

  Box1Abb[2292]=-Box1Abb[2291]*m_x_2 + 2.*Box1Abb[2289]*m_x*m_z12 + 6.*m_x_3*m_z12 - pow(Box1Abb[68],2.)*m_z12_2;

  Box1Abb[2293]=-3. + m_z12 - 6.*m_z1k;

  Box1Abb[2294]=-2. + Box1Abb[2293]*m_z12;

  Box1Abb[2295]=7. + 4.*m_z1k;

  Box1Abb[2296]=-1. + Box1Abb[2295]*m_z1k;

  Box1Abb[2297]=2.*Box1Abb[265]*pow(Box1Abb[68],2.) + Box1Abb[2296]*m_z12;

  Box1Abb[2298]=36. - Box1Abb[2060]*m_z12 + 54.*m_z1k + 12.*m_z12*m_z1k + 90.*m_z1k_2;

  Box1Abb[2299]=64. + Box1Abb[2298]*m_z12 + 48.*m_z1k;

  Box1Abb[2300]=-2. + Box1Abb[2299]*m_z12;

  Box1Abb[2301]=-2. + 3.*m_z1k + m_z1k_3;

  Box1Abb[2302]=-5. + 20.*m_z1k - 46.*m_z1k_3 + 31.*m_z1k_4;

  Box1Abb[2303]=-7. + 3.*m_z1k;

  Box1Abb[2304]=12. + Box1Abb[2303]*m_z1k;

  Box1Abb[2305]=1. + Box1Abb[2304]*m_z1k;

  Box1Abb[2306]=6.*Box1Abb[2301]*pow(Box1Abb[68],3.) + 4.*Box1Abb[2305]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2302]*m_z12_2 + 12.*m_z12_3*m_z1k_3;

  Box1Abb[2307]=15. + 2.*m_z1k - 6.*m_z1k_3 + 15.*m_z1k_4;

  Box1Abb[2308]=-5. + 31.*m_z1k;

  Box1Abb[2309]=-10. + Box1Abb[2308]*m_z1k;

  Box1Abb[2310]=29. + 20.*m_z1k + 42.*m_z1k_2;

  Box1Abb[2311]=-10. + Box1Abb[2310]*m_z1k;

  Box1Abb[2312]=-4. + 3.*m_z1k;

  Box1Abb[2313]=-13. + 4.*Box1Abb[2312]*m_z1k;

  Box1Abb[2314]=64. + 4.*Box1Abb[2313]*m_z1k;

  Box1Abb[2315]=-2.*Box1Abb[1446]*Box1Abb[68] + Box1Abb[2314]*m_z12 + 6.*Box1Abb[2307]*m_z12_2 + 4.*Box1Abb[2311]*m_z12_3 + Box1Abb[2309]*m_z12_4;

  Box1Abb[2316]=-25. + 51.*m_z1k;

  Box1Abb[2317]=13. + Box1Abb[725]*m_z1k;

  Box1Abb[2318]=6. + Box1Abb[235]*m_z1k;

  Box1Abb[2319]=10. - 8.*Box1Abb[2317]*m_z12 - 12.*Box1Abb[2318]*m_z12_2 - 2.*Box1Abb[15]*Box1Abb[2316]*m_z12_3 + Box1Abb[1945]*m_z12_4 + 4.*m_z1k;

  Box1Abb[2320]=-9. + m_z1k + 3.*m_z1k_2 - 3.*m_z1k_3 + 6.*m_z1k_4;

  Box1Abb[2321]=23. + 34.*m_z1k;

  Box1Abb[2322]=5. + Box1Abb[2321]*m_z1k;

  Box1Abb[2323]=-5. + Box1Abb[2322]*m_z1k;

  Box1Abb[2324]=-40. + 51.*m_z1k;

  Box1Abb[2325]=-54. + Box1Abb[2324]*m_z1k;

  Box1Abb[2326]=60. + Box1Abb[2325]*m_z1k;

  Box1Abb[2327]=-5. + Box1Abb[2326]*m_z1k;

  Box1Abb[2328]=2.*Box1Abb[457]*pow(Box1Abb[68],2.)*Box1Abb[881] + 3.*Box1Abb[2320]*Box1Abb[68]*m_z12 + Box1Abb[2327]*m_z12_2 + Box1Abb[2323]*m_z12_3;

  Box1Abb[2329]=Box1Abb[2315]*m_x_4 + Box1Abb[2319]*m_x_5 + Box1Abb[2300]*m_x_6 - 2.*Box1Abb[2328]*m_x_3*m_z12 + 6.*Box1Abb[2294]*m_x_7*m_z12 + Box1Abb[2306]*m_x_2*m_z12_2 + 6.*m_x_8*m_z12_2 + Box1Abb[2297]*pow(Box1Abb[68],3.)*m_x*m_z12_3 - pow(Box1Abb[68],5.)*m_z12_4*m_z1k;

  Box1Abb[2330]=-36. + 11.*m_z12 - 84.*m_z1k;

  Box1Abb[2331]=68. + Box1Abb[2330]*m_z12 + 116.*m_z1k;

  Box1Abb[2332]=-78. + 5.*Box1Abb[677]*m_z12;

  Box1Abb[2333]=96. + Box1Abb[2332]*m_z12;

  Box1Abb[2334]=18. + Box1Abb[79]*m_z12;

  Box1Abb[2335]=-14. + Box1Abb[2334]*m_z12;

  Box1Abb[2336]=119. - 57.*m_z12;

  Box1Abb[2337]=162. + Box1Abb[2336]*m_z12;

  Box1Abb[2338]=-172. + Box1Abb[2337]*m_z12;

  Box1Abb[2339]=10. + Box1Abb[2338]*m_z12;

  Box1Abb[2340]=130. - 183.*m_z12;

  Box1Abb[2341]=48. + Box1Abb[2340]*m_z12;

  Box1Abb[2342]=232. + Box1Abb[2341]*m_z12;

  Box1Abb[2343]=12. + Box1Abb[2342]*m_z12;

  Box1Abb[2344]=-78. + 53.*m_z12;

  Box1Abb[2345]=58. + Box1Abb[2344]*m_z12;

  Box1Abb[2346]=4. + Box1Abb[2345]*m_z12;

  Box1Abb[2347]=17. - 9.*m_z12;

  Box1Abb[2348]=-34. + Box1Abb[2333]*m_z12 + 52.*m_z1k + 10.*Box1Abb[2335]*m_z12*m_z1k + Box1Abb[2339]*m_z1k_2 + Box1Abb[2343]*m_z1k_3 - 10.*Box1Abb[2346]*m_z1k_4 + 28.*Box1Abb[2347]*m_z12*m_z1k_5;

  Box1Abb[2349]=1. + Box1Abb[881]*m_z1k;

  Box1Abb[2350]=-2.*Box1Abb[317]*pow(Box1Abb[68],3.) - 3.*pow(Box1Abb[68],2.)*Box1Abb[834]*m_z12 + 2.*Box1Abb[2349]*Box1Abb[68]*m_z12_2 + Box1Abb[15]*Box1Abb[431]*m_z12_3;

  Box1Abb[2351]=25. + m_z12;

  Box1Abb[2352]=69. - 2.*Box1Abb[2351]*m_z12;

  Box1Abb[2353]=16. + m_z12;

  Box1Abb[2354]=-13. + 9.*m_z12;

  Box1Abb[2355]=50. + 137.*m_z1k;

  Box1Abb[2356]=-2.*Box1Abb[2355] + Box1Abb[2352]*m_z12 + 6.*Box1Abb[2353]*m_z12*m_z1k + 28.*Box1Abb[2354]*m_z1k_2;

  Box1Abb[2357]=Box1Abb[2356]*m_z12 + 4.*m_z1k;

  Box1Abb[2358]=10. + 11.*m_z1k;

  Box1Abb[2359]=1. + 32.*m_z1k;

  Box1Abb[2360]=85. - 6.*Box1Abb[2359]*m_z1k;

  Box1Abb[2361]=6. - 35.*m_z1k;

  Box1Abb[2362]=5. + Box1Abb[2361]*m_z1k;

  Box1Abb[2363]=-129. + 12.*Box1Abb[2362]*m_z1k;

  Box1Abb[2364]=81. + 161.*m_z1k;

  Box1Abb[2365]=51. + Box1Abb[2364]*m_z1k;

  Box1Abb[2366]=53. + 2.*Box1Abb[2365]*m_z1k;

  Box1Abb[2367]=2.*Box1Abb[2366]*m_z12 + Box1Abb[2363]*m_z12_2 + Box1Abb[2360]*m_z12_3 + Box1Abb[2358]*m_z12_4 - 2.*Box1Abb[230]*m_z1k;

  Box1Abb[2368]=-3. + Box1Abb[884]*m_z1k;

  Box1Abb[2369]=1. + 2.*Box1Abb[2368]*m_z1k;

  Box1Abb[2370]=-9. + m_z1k;

  Box1Abb[2371]=3. + Box1Abb[2370]*m_z1k;

  Box1Abb[2372]=-1. + 4.*Box1Abb[2371]*m_z1k;

  Box1Abb[2373]=1. + Box1Abb[2372]*m_z1k;

  Box1Abb[2374]=25. - 78.*m_z1k + 56.*m_z1k_2;

  Box1Abb[2375]=-16. + Box1Abb[2374]*m_z1k;

  Box1Abb[2376]=5. + Box1Abb[2375]*m_z1k;

  Box1Abb[2377]=-15. + 67.*m_z1k;

  Box1Abb[2378]=17. + Box1Abb[2377]*m_z1k;

  Box1Abb[2379]=-3. + Box1Abb[2378]*m_z1k;

  Box1Abb[2380]=2. + Box1Abb[2379]*m_z1k;

  Box1Abb[2381]=-2.*Box1Abb[317]*pow(Box1Abb[68],5.) + 2.*Box1Abb[2369]*pow(Box1Abb[68],4.)*m_z12 - 3.*Box1Abb[2373]*pow(Box1Abb[68],3.)*m_z12_2 - Box1Abb[2376]*pow(Box1Abb[68],2.)*m_z12_3 - Box1Abb[2380]*Box1Abb[68]*m_z12_4 - 24.*m_z12_5*m_z1k_4;

  Box1Abb[2382]=15. + 4.*Box1Abb[1151]*m_z1k;

  Box1Abb[2383]=-25. + 46.*m_z1k;

  Box1Abb[2384]=-20. + Box1Abb[2383]*m_z1k;

  Box1Abb[2385]=67. + 243.*m_z1k;

  Box1Abb[2386]=11. + Box1Abb[2385]*m_z1k;

  Box1Abb[2387]=-30. + Box1Abb[2386]*m_z1k;

  Box1Abb[2388]=-9. + 7.*m_z1k;

  Box1Abb[2389]=-129. + 20.*Box1Abb[2388]*m_z1k;

  Box1Abb[2390]=-53. + Box1Abb[2389]*m_z1k;

  Box1Abb[2391]=50. + Box1Abb[2390]*m_z1k;

  Box1Abb[2392]=-11. + 70.*m_z1k;

  Box1Abb[2393]=64. + 5.*Box1Abb[2392]*m_z1k;

  Box1Abb[2394]=13. + Box1Abb[2393]*m_z1k;

  Box1Abb[2395]=64. + Box1Abb[2394]*m_z1k;

  Box1Abb[2396]=22. - 2.*Box1Abb[2395]*m_z12 + 3.*Box1Abb[2391]*m_z12_2 + 2.*Box1Abb[2387]*m_z12_3 + Box1Abb[2384]*m_z12_4 + 2.*Box1Abb[2382]*m_z1k;

  Box1Abb[2397]=-7. + 2.*Box1Abb[514]*m_z1k;

  Box1Abb[2398]=-6. + 97.*m_z1k;

  Box1Abb[2399]=-5. + Box1Abb[2398]*m_z1k;

  Box1Abb[2400]=5. + Box1Abb[2399]*m_z1k;

  Box1Abb[2401]=-5. + Box1Abb[2400]*m_z1k;

  Box1Abb[2402]=-65. + 98.*m_z1k;

  Box1Abb[2403]=6. + Box1Abb[2402]*m_z1k;

  Box1Abb[2404]=-37. + Box1Abb[2403]*m_z1k;

  Box1Abb[2405]=14. + Box1Abb[2404]*m_z1k;

  Box1Abb[2406]=-5. + m_z1k;

  Box1Abb[2407]=41. + 28.*Box1Abb[2406]*m_z1k;

  Box1Abb[2408]=-21. + Box1Abb[2407]*m_z1k;

  Box1Abb[2409]=23. + Box1Abb[2408]*m_z1k;

  Box1Abb[2410]=-3. + Box1Abb[2409]*m_z1k;

  Box1Abb[2411]=-195. + 137.*m_z1k;

  Box1Abb[2412]=-6. + Box1Abb[2411]*m_z1k;

  Box1Abb[2413]=-16. + Box1Abb[2412]*m_z1k;

  Box1Abb[2414]=-15. + Box1Abb[2413]*m_z1k;

  Box1Abb[2415]=7. + Box1Abb[2414]*m_z1k;

  Box1Abb[2416]=2.*Box1Abb[2397]*pow(Box1Abb[68],3.) - 2.*Box1Abb[2405]*pow(Box1Abb[68],2.)*m_z12 + 3.*Box1Abb[2410]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[2415]*m_z12_3 + 2.*Box1Abb[2401]*m_z12_4 + 24.*m_z12_5*m_z1k_3;

  Box1Abb[2417]=Box1Abb[2381]*m_x_2 + Box1Abb[2416]*m_x_3 + Box1Abb[2348]*m_x_4 + Box1Abb[2396]*m_x_5 + Box1Abb[2367]*m_x_6 + Box1Abb[2357]*m_x_7 + Box1Abb[2331]*m_x_8*m_z12 + 4.*Box1Abb[129]*m_x_9*m_z12 + Box1Abb[2350]*pow(Box1Abb[68],3.)*m_x*m_z12*m_z1k - Box1Abb[4]*pow(Box1Abb[68],5.)*m_z12_3*m_z1k_2;

  Box1Abb[2418]=18. - 7.*m_z12;

  Box1Abb[2419]=-7. + Box1Abb[2418]*m_z12;

  Box1Abb[2420]=-9. + 2.*m_z12;

  Box1Abb[2421]=4. + Box1Abb[2420]*Box1Abb[699]*m_z12;

  Box1Abb[2422]=-51. + 32.*m_z12;

  Box1Abb[2423]=12. + Box1Abb[2422]*m_z12;

  Box1Abb[2424]=4. + Box1Abb[2423]*m_z12;

  Box1Abb[2425]=-4. + Box1Abb[653]*m_z12;

  Box1Abb[2426]=-31. + 5.*Box1Abb[2425]*m_z12;

  Box1Abb[2427]=27. + Box1Abb[2426]*m_z12;

  Box1Abb[2428]=-5. + Box1Abb[2427]*m_z12;

  Box1Abb[2429]=-59. + 64.*m_z12;

  Box1Abb[2430]=144. + Box1Abb[2429]*m_z12;

  Box1Abb[2431]=-117. + Box1Abb[2430]*m_z12;

  Box1Abb[2432]=8. + Box1Abb[2431]*m_z12;

  Box1Abb[2433]=-108. + 43.*m_z12;

  Box1Abb[2434]=108. + Box1Abb[2433]*m_z12;

  Box1Abb[2435]=-2. + Box1Abb[2434]*m_z12;

  Box1Abb[2436]=-4. + Box1Abb[2419]*m_z12 + 8.*m_z1k - Box1Abb[2421]*m_z12*m_z1k + 2.*Box1Abb[2424]*m_z12*m_z1k_2 + 2.*Box1Abb[2428]*m_z1k_3 + Box1Abb[2432]*m_z1k_4 + Box1Abb[2435]*m_z1k_5 + 14.*Box1Abb[493]*m_z12*m_z1k_6;

  Box1Abb[2437]=-11. + m_z12 - 17.*m_z1k;

  Box1Abb[2438]=27. + 2.*Box1Abb[2437]*m_z12 + 42.*m_z1k;

  Box1Abb[2439]=2. + Box1Abb[2438]*m_z12;

  Box1Abb[2440]=6. + 5.*m_z1k;

  Box1Abb[2441]=-11. + 17.*m_z1k;

  Box1Abb[2442]=44. + 80.*m_z1k + 98.*m_z1k_2;

  Box1Abb[2443]=-41. + Box1Abb[2442]*m_z12 + Box1Abb[2441]*m_z12_2 - 18.*Box1Abb[1776]*m_z1k;

  Box1Abb[2444]=-2.*Box1Abb[2440] + Box1Abb[2443]*m_z12;

  Box1Abb[2445]=3. + Box1Abb[178]*m_z1k_2;

  Box1Abb[2446]=-8. + 5.*Box1Abb[873]*m_z1k;

  Box1Abb[2447]=1. + Box1Abb[2446]*m_z1k;

  Box1Abb[2448]=5. + 9.*m_z1k;

  Box1Abb[2449]=-10. + Box1Abb[2448]*m_z1k;

  Box1Abb[2450]=4. + Box1Abb[2449]*m_z1k;

  Box1Abb[2451]=Box1Abb[2445]*pow(Box1Abb[68],3.) + Box1Abb[2450]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2447]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[178]*m_z12_3*m_z1k_2;

  Box1Abb[2452]=33. + 53.*m_z1k;

  Box1Abb[2453]=25. - 2.*Box1Abb[2452]*m_z1k;

  Box1Abb[2454]=9. + 14.*m_z1k;

  Box1Abb[2455]=31. + 5.*Box1Abb[2454]*m_z1k;

  Box1Abb[2456]=22. + 3.*Box1Abb[2455]*m_z1k;

  Box1Abb[2457]=36. + 77.*m_z1k;

  Box1Abb[2458]=9. + Box1Abb[2457]*m_z1k;

  Box1Abb[2459]=59. + 2.*Box1Abb[2458]*m_z1k;

  Box1Abb[2460]=30. + Box1Abb[2456]*m_z12 - Box1Abb[2459]*m_z12_2 + Box1Abb[2453]*m_z12_3 + 4.*Box1Abb[818]*m_z1k + m_z12_4*m_z1k;

  Box1Abb[2461]=15. + 41.*m_z1k + 105.*m_z1k_3;

  Box1Abb[2462]=1. + Box1Abb[2461]*m_z1k;

  Box1Abb[2463]=5. + 2.*Box1Abb[514]*m_z1k;

  Box1Abb[2464]=18. + Box1Abb[2463]*m_z1k;

  Box1Abb[2465]=192. + 205.*m_z1k;

  Box1Abb[2466]=106. + Box1Abb[2465]*m_z1k;

  Box1Abb[2467]=-30. + Box1Abb[2466]*m_z1k;

  Box1Abb[2468]=-9. + 5.*Box1Abb[317]*m_z1k;

  Box1Abb[2469]=-95. + 14.*Box1Abb[2468]*m_z1k;

  Box1Abb[2470]=61. + Box1Abb[2469]*m_z1k;

  Box1Abb[2471]=-2.*Box1Abb[2464] - 2.*Box1Abb[2462]*m_z12 + Box1Abb[2470]*m_z12_2 + Box1Abb[2467]*m_z12_3 + Box1Abb[1445]*m_z12_4*m_z1k;

  Box1Abb[2472]=-32. - Box1Abb[1446]*Box1Abb[228]*m_z1k;

  Box1Abb[2473]=5. + Box1Abb[2472]*m_z1k;

  Box1Abb[2474]=4. - 9.*m_z1k + 6.*m_z1k_2;

  Box1Abb[2475]=-5. + Box1Abb[2474]*m_z1k;

  Box1Abb[2476]=-2. + Box1Abb[2475]*m_z1k;

  Box1Abb[2477]=13. + m_z1k;

  Box1Abb[2478]=-9. + Box1Abb[2477]*m_z1k;

  Box1Abb[2479]=-9. + 2.*Box1Abb[2478]*m_z1k;

  Box1Abb[2480]=3. + Box1Abb[2479]*m_z1k;

  Box1Abb[2481]=-35. + 6.*Box1Abb[475]*m_z1k;

  Box1Abb[2482]=-45. + Box1Abb[2481]*m_z1k;

  Box1Abb[2483]=21. + Box1Abb[2482]*m_z1k;

  Box1Abb[2484]=-1. + Box1Abb[2483]*m_z1k;

  Box1Abb[2485]=Box1Abb[2476]*pow(Box1Abb[68],3.) + Box1Abb[2480]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[2484]*Box1Abb[68]*m_z12_2 + Box1Abb[2473]*m_z12_3*m_z1k - 2.*Box1Abb[2295]*m_z12_4*m_z1k_3;

  Box1Abb[2486]=-10. + 77.*m_z1k;

  Box1Abb[2487]=3. + Box1Abb[1151]*m_z1k;

  Box1Abb[2488]=89. + 170.*m_z1k;

  Box1Abb[2489]=29. + Box1Abb[2488]*m_z1k;

  Box1Abb[2490]=96. + Box1Abb[2489]*m_z1k;

  Box1Abb[2491]=-20. + Box1Abb[2490]*m_z1k;

  Box1Abb[2492]=17. - 7.*m_z1k;

  Box1Abb[2493]=22. + 5.*Box1Abb[2492]*m_z1k;

  Box1Abb[2494]=53. + Box1Abb[2493]*m_z1k;

  Box1Abb[2495]=53. + Box1Abb[2494]*m_z1k;

  Box1Abb[2496]=-22. + Box1Abb[2495]*m_z1k;

  Box1Abb[2497]=-15. + 14.*m_z1k;

  Box1Abb[2498]=98. + 9.*Box1Abb[2497]*m_z1k;

  Box1Abb[2499]=4. + Box1Abb[2498]*m_z1k;

  Box1Abb[2500]=8. + Box1Abb[2499]*m_z1k;

  Box1Abb[2501]=5. + Box1Abb[2500]*m_z1k;

  Box1Abb[2502]=20. + Box1Abb[2501]*m_z12 + 2.*Box1Abb[2496]*m_z12_2 - Box1Abb[2491]*m_z12_3 + 2.*Box1Abb[170]*Box1Abb[2487]*m_z1k - Box1Abb[15]*Box1Abb[2486]*m_z12_4*m_z1k;

  Box1Abb[2503]=Box1Abb[2436]*m_x_3 + Box1Abb[2502]*m_x_4 + Box1Abb[2471]*m_x_5 + Box1Abb[2460]*m_x_6 + Box1Abb[2444]*m_x_7 + Box1Abb[2439]*m_x_8 + Box1Abb[2485]*m_x_2*m_z12 + Box1Abb[699]*m_x_9*m_z12 - Box1Abb[2451]*Box1Abb[68]*m_x*m_z12_2*m_z1k + Box1Abb[0]*pow(Box1Abb[68],5.)*m_z12_3*m_z1k_2;

  Box1Abb[2504]=3. + 7.*m_z12 - 15.*m_z1k;

  Box1Abb[2505]=-6. + Box1Abb[2504]*m_z12;

  Box1Abb[2506]=-26. + 5.*Box1Abb[14]*m_z12;

  Box1Abb[2507]=23. + Box1Abb[2506]*m_z12;

  Box1Abb[2508]=18. - Box1Abb[2507]*m_z12;

  Box1Abb[2509]=32. + 5.*Box1Abb[2420]*m_z12;

  Box1Abb[2510]=-5. + Box1Abb[2509]*m_z12;

  Box1Abb[2511]=9. + 2.*m_z12;

  Box1Abb[2512]=3. + 5.*Box1Abb[2511]*m_z12;

  Box1Abb[2513]=-14. + 5.*m_z12;

  Box1Abb[2514]=86. + 5.*Box1Abb[2513]*m_z12;

  Box1Abb[2515]=-54. + Box1Abb[2514]*m_z12;

  Box1Abb[2516]=3. + Box1Abb[2515]*m_z12;

  Box1Abb[2517]=-402. + 83.*m_z12;

  Box1Abb[2518]=433. + Box1Abb[2517]*m_z12;

  Box1Abb[2519]=-212. + Box1Abb[2518]*m_z12;

  Box1Abb[2520]=8. + Box1Abb[2519]*m_z12;

  Box1Abb[2521]=-88. + 73.*m_z12;

  Box1Abb[2522]=36. + Box1Abb[2521]*m_z12;

  Box1Abb[2523]=-6. + Box1Abb[2508]*m_z12 + 16.*m_z1k + 2.*Box1Abb[2510]*m_z12*m_z1k + 2.*Box1Abb[0]*Box1Abb[2512]*m_z1k_2 - 4.*Box1Abb[2516]*m_z1k_3 + Box1Abb[2520]*m_z1k_4 + 2.*Box1Abb[2522]*m_z12*m_z1k_5;

  Box1Abb[2524]=1. - 2.*m_z1k + 4.*m_z1k_2;

  Box1Abb[2525]=1. + Box1Abb[1041]*m_z1k;

  Box1Abb[2526]=Box1Abb[534]*pow(Box1Abb[68],3.) + 2.*Box1Abb[2524]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2525]*Box1Abb[68]*m_z12_2 + 2.*m_z12_3*m_z1k_2;

  Box1Abb[2527]=54. + 58.*m_z1k;

  Box1Abb[2528]=-3. + Box1Abb[2527]*m_z12 + 7.*m_z12_2 + 32.*Box1Abb[1367]*m_z1k;

  Box1Abb[2529]=60. - Box1Abb[2528]*m_z12 + 72.*m_z1k;

  Box1Abb[2530]=32. + 39.*m_z1k;

  Box1Abb[2531]=41. + 90.*m_z1k;

  Box1Abb[2532]=82. + 31.*m_z1k;

  Box1Abb[2533]=32. + Box1Abb[2532]*m_z1k;

  Box1Abb[2534]=-61. + 42.*m_z1k;

  Box1Abb[2535]=35. + 4.*Box1Abb[2534]*m_z1k;

  Box1Abb[2536]=60. + Box1Abb[2535]*m_z1k;

  Box1Abb[2537]=2. - Box1Abb[2531]*Box1Abb[431]*m_z12 - Box1Abb[2536]*m_z12_2 + 2.*Box1Abb[2533]*m_z12_3 + Box1Abb[2530]*m_z12_4 + m_z12_5 - 2.*m_z1k;

  Box1Abb[2538]=-119. + 31.*m_z1k;

  Box1Abb[2539]=-33. + Box1Abb[2538]*m_z1k;

  Box1Abb[2540]=76. + 59.*m_z1k;

  Box1Abb[2541]=55. + Box1Abb[2540]*m_z1k;

  Box1Abb[2542]=37. + 60.*m_z1k;

  Box1Abb[2543]=11. + Box1Abb[2542]*m_z1k;

  Box1Abb[2544]=-19. + 7.*m_z1k;

  Box1Abb[2545]=65. + 8.*Box1Abb[2544]*m_z1k;

  Box1Abb[2546]=148. + 3.*Box1Abb[2545]*m_z1k;

  Box1Abb[2547]=53. + Box1Abb[2546]*m_z1k;

  Box1Abb[2548]=18. + Box1Abb[2543]*Box1Abb[834]*m_z12 + Box1Abb[2547]*m_z12_2 - Box1Abb[2541]*m_z12_4 - 5.*Box1Abb[15]*m_z12_5 + 4.*m_z1k + 2.*Box1Abb[2539]*m_z12_3*m_z1k + 8.*m_z1k_2;

  Box1Abb[2549]=6. + 5.*Box1Abb[170]*m_z1k;

  Box1Abb[2550]=-4. + Box1Abb[2549]*m_z1k;

  Box1Abb[2551]=1. + Box1Abb[2550]*m_z1k;

  Box1Abb[2552]=26. + m_z1k;

  Box1Abb[2553]=-32. + Box1Abb[2552]*m_z1k;

  Box1Abb[2554]=18. + Box1Abb[2553]*m_z1k;

  Box1Abb[2555]=-3. + Box1Abb[2554]*m_z1k;

  Box1Abb[2556]=-17. + 4.*Box1Abb[536]*m_z1k;

  Box1Abb[2557]=8. + Box1Abb[2556]*m_z1k;

  Box1Abb[2558]=-3. + Box1Abb[2557]*m_z1k;

  Box1Abb[2559]=6. + 11.*m_z1k;

  Box1Abb[2560]=-20. + Box1Abb[2559]*m_z1k;

  Box1Abb[2561]=11. + Box1Abb[2560]*m_z1k;

  Box1Abb[2562]=-2. + Box1Abb[2561]*m_z1k;

  Box1Abb[2563]=pow(Box1Abb[68],4.) - Box1Abb[2558]*pow(Box1Abb[68],3.)*m_z12 - 2.*Box1Abb[2562]*pow(Box1Abb[68],2.)*m_z12_2 - Box1Abb[2555]*Box1Abb[68]*m_z12_3 + Box1Abb[2551]*m_z12_4;

  Box1Abb[2564]=58. - 7.*m_z1k;

  Box1Abb[2565]=15. + Box1Abb[2564]*m_z1k;

  Box1Abb[2566]=40. + Box1Abb[2565]*m_z1k;

  Box1Abb[2567]=-2. + m_z1k + 3.*m_z1k_2;

  Box1Abb[2568]=1. + Box1Abb[2567]*m_z1k;

  Box1Abb[2569]=14. + 3.*Box1Abb[457]*Box1Abb[785]*m_z1k;

  Box1Abb[2570]=5. + 2.*Box1Abb[2569]*m_z1k;

  Box1Abb[2571]=-31. + 19.*m_z1k;

  Box1Abb[2572]=-3. + Box1Abb[2571]*Box1Abb[438]*m_z1k;

  Box1Abb[2573]=25. + Box1Abb[2572]*m_z1k;

  Box1Abb[2574]=-447. + 440.*m_z1k - 84.*m_z1k_2;

  Box1Abb[2575]=4. + Box1Abb[2574]*m_z1k;

  Box1Abb[2576]=17. + Box1Abb[2575]*m_z1k;

  Box1Abb[2577]=14. + Box1Abb[2576]*m_z1k;

  Box1Abb[2578]=-4.*Box1Abb[2568] - 2.*Box1Abb[2570]*m_z12 + Box1Abb[2577]*m_z12_2 - 2.*Box1Abb[2573]*m_z12_3 + Box1Abb[2566]*m_z12_4 + 10.*Box1Abb[821]*m_z12_5;

  Box1Abb[2579]=-5. + 10.*m_z1k - 16.*m_z1k_3 + 13.*m_z1k_4;

  Box1Abb[2580]=-8. + 3.*m_z1k;

  Box1Abb[2581]=15. + 4.*Box1Abb[2580]*m_z1k;

  Box1Abb[2582]=-3. + Box1Abb[2581]*m_z1k;

  Box1Abb[2583]=-35. + 4.*m_z1k + 8.*m_z1k_2;

  Box1Abb[2584]=86. + 3.*Box1Abb[2583]*m_z1k;

  Box1Abb[2585]=-41. + Box1Abb[2584]*m_z1k;

  Box1Abb[2586]=-75. + 11.*m_z1k;

  Box1Abb[2587]=89. + Box1Abb[2586]*m_z1k;

  Box1Abb[2588]=-41. + Box1Abb[2587]*m_z1k;

  Box1Abb[2589]=8. + Box1Abb[2588]*m_z1k;

  Box1Abb[2590]=2. + Box1Abb[2589]*m_z1k;

  Box1Abb[2591]=-28. + 11.*m_z1k;

  Box1Abb[2592]=94. + 5.*Box1Abb[2591]*m_z1k;

  Box1Abb[2593]=4. + Box1Abb[2592]*m_z1k;

  Box1Abb[2594]=-25. + Box1Abb[2593]*m_z1k;

  Box1Abb[2595]=8. + Box1Abb[2594]*m_z1k;

  Box1Abb[2596]=2.*pow(Box1Abb[68],5.) + Box1Abb[2582]*pow(Box1Abb[68],3.)*m_z12 + 2.*Box1Abb[2590]*Box1Abb[68]*m_z12_3 + Box1Abb[2595]*m_z12_4 + Box1Abb[2579]*m_z12_5 - Box1Abb[2585]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k;

  Box1Abb[2597]=-Box1Abb[2596]*m_x_2 + Box1Abb[2523]*m_x_3 + Box1Abb[2578]*m_x_4 + Box1Abb[2548]*m_x_5 + Box1Abb[2537]*m_x_6 + Box1Abb[2563]*Box1Abb[68]*m_x*m_z12 + Box1Abb[2529]*m_x_7*m_z12 + 2.*Box1Abb[2505]*m_x_8*m_z12 + 4.*m_x_9*m_z12_2 + Box1Abb[2526]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k;

  Box1Abb[2598]=26. + m_z12;

  Box1Abb[2599]=6. + Box1Abb[2598]*m_z12;

  Box1Abb[2600]=-56. + Box1Abb[2599]*m_z12;

  Box1Abb[2601]=29. + Box1Abb[2600]*m_z12;

  Box1Abb[2602]=91. + 34.*m_z12;

  Box1Abb[2603]=-39. + Box1Abb[2602]*m_z12;

  Box1Abb[2604]=-26. + Box1Abb[2603]*m_z12;

  Box1Abb[2605]=240. + 17.*m_z12;

  Box1Abb[2606]=-9. + Box1Abb[2605]*m_z12;

  Box1Abb[2607]=30. + Box1Abb[2606]*m_z12;

  Box1Abb[2608]=10. + Box1Abb[2601]*m_z12 + 20.*m_z1k + Box1Abb[2604]*m_z12*m_z1k + Box1Abb[2607]*m_z1k_2 + 84.*Box1Abb[527]*m_z12*m_z1k_3;

  Box1Abb[2609]=-14. + Box1Abb[707]*m_z12;

  Box1Abb[2610]=3. - 5.*Box1Abb[2609]*m_z12;

  Box1Abb[2611]=-55. + Box1Abb[2610]*m_z12;

  Box1Abb[2612]=12. - 5.*Box1Abb[708]*m_z12;

  Box1Abb[2613]=100. + Box1Abb[2612]*m_z12;

  Box1Abb[2614]=-87. + Box1Abb[2613]*m_z12;

  Box1Abb[2615]=18. + 5.*m_z12;

  Box1Abb[2616]=-27. + 2.*Box1Abb[2615]*m_z12;

  Box1Abb[2617]=7. + Box1Abb[2616]*m_z12;

  Box1Abb[2618]=-614. + 115.*m_z12;

  Box1Abb[2619]=285. + Box1Abb[2618]*m_z12;

  Box1Abb[2620]=-40. + Box1Abb[2619]*m_z12;

  Box1Abb[2621]=20. + Box1Abb[2611]*m_z12 + Box1Abb[2614]*m_z12*m_z1k - 5.*Box1Abb[2617]*m_z12*m_z1k_2 + Box1Abb[2620]*m_z1k_3 + 42.*Box1Abb[1647]*m_z12*m_z1k_4;

  Box1Abb[2622]=3. + m_z12 - m_z12_2;

  Box1Abb[2623]=-11. + 2.*Box1Abb[2622]*m_z12;

  Box1Abb[2624]=29. + 5.*Box1Abb[2623]*m_z12;

  Box1Abb[2625]=-17. + 4.*m_z12;

  Box1Abb[2626]=102. + 5.*Box1Abb[2625]*m_z12;

  Box1Abb[2627]=-53. + Box1Abb[2626]*m_z12;

  Box1Abb[2628]=2. + Box1Abb[69]*Box1Abb[708]*m_z12;

  Box1Abb[2629]=-73. + 26.*m_z12;

  Box1Abb[2630]=116. + Box1Abb[2629]*m_z12;

  Box1Abb[2631]=-91. + Box1Abb[2630]*m_z12;

  Box1Abb[2632]=10. + Box1Abb[2631]*m_z12;

  Box1Abb[2633]=825. - 510.*m_z12 + 86.*m_z12_2;

  Box1Abb[2634]=-595. + Box1Abb[2633]*m_z12;

  Box1Abb[2635]=40. + Box1Abb[2634]*m_z12;

  Box1Abb[2636]=-558. + 185.*m_z12;

  Box1Abb[2637]=495. + Box1Abb[2636]*m_z12;

  Box1Abb[2638]=-12. + Box1Abb[2637]*m_z12;

  Box1Abb[2639]=-4. + Box1Abb[2624]*m_z12 + 12.*m_z1k + Box1Abb[2627]*m_z12*m_z1k + 2.*Box1Abb[2628]*m_z1k_2 - 4.*Box1Abb[2632]*m_z1k_3 + Box1Abb[2635]*m_z1k_4 + Box1Abb[2638]*m_z1k_5 + 84.*Box1Abb[79]*m_z12*m_z1k_6;

  Box1Abb[2640]=-9. + Box1Abb[77]*m_z12;

  Box1Abb[2641]=39. + 5.*Box1Abb[2640]*m_z12;

  Box1Abb[2642]=-7. + 2.*Box1Abb[2641]*m_z12;

  Box1Abb[2643]=23. - 14.*m_z12 + 10.*m_z12_3;

  Box1Abb[2644]=-12. + Box1Abb[2643]*m_z12;

  Box1Abb[2645]=24. + 5.*m_z12;

  Box1Abb[2646]=-27. + Box1Abb[2645]*m_z12;

  Box1Abb[2647]=14. + Box1Abb[2646]*m_z12;

  Box1Abb[2648]=-12. + Box1Abb[2647]*m_z12;

  Box1Abb[2649]=55. - 2.*m_z12;

  Box1Abb[2650]=-515. + 7.*Box1Abb[2649]*m_z12;

  Box1Abb[2651]=340. + Box1Abb[2650]*m_z12;

  Box1Abb[2652]=-40. + Box1Abb[2651]*m_z12;

  Box1Abb[2653]=160. - 47.*m_z12;

  Box1Abb[2654]=-111.+Box1Abb[2653]*m_z12;

  Box1Abb[2655]=6. + Box1Abb[2654]*m_z12;

  Box1Abb[2656]=8. - 5.*m_z12;

  Box1Abb[2657]=2.*Box1Abb[1830] + Box1Abb[2642]*m_z12 + Box1Abb[2644]*m_z12*m_z1k + 2.*Box1Abb[2648]*m_z12*m_z1k_2 + Box1Abb[2652]*m_z1k_3 + 5.*Box1Abb[2655]*m_z1k_4 + 42.*Box1Abb[2656]*m_z12*m_z1k_5;

  Box1Abb[2658]=-2. + 10.*m_z12 - 39.*m_z1k;

  Box1Abb[2659]=15. + Box1Abb[2658]*m_z12 + 48.*m_z1k;

  Box1Abb[2660]=2. + Box1Abb[2659]*m_z12;

  Box1Abb[2661]=5. + m_z12;

  Box1Abb[2662]=-19. + 6.*Box1Abb[2661]*m_z12;

  Box1Abb[2663]=34. + 37.*m_z12;

  Box1Abb[2664]=14. - 11.*m_z12;

  Box1Abb[2665]=-3. + Box1Abb[2662]*m_z12 + 51.*m_z1k + Box1Abb[2663]*m_z12*m_z1k + 12.*Box1Abb[2664]*m_z1k_2;

  Box1Abb[2666]=8. + Box1Abb[2665]*m_z12 + 12.*m_z1k;

  Box1Abb[2667]=1. - 3.*m_z1k + 6.*m_z1k_2;

  Box1Abb[2668]=-15. + m_z1k;

  Box1Abb[2669]=-8. + Box1Abb[2668]*Box1Abb[68]*m_z1k;

  Box1Abb[2670]=2. + Box1Abb[2669]*m_z1k;

  Box1Abb[2671]=22. + 3.*m_z1k;

  Box1Abb[2672]=-39. + Box1Abb[2671]*m_z1k;

  Box1Abb[2673]=16. + Box1Abb[2672]*m_z1k;

  Box1Abb[2674]=-4. + Box1Abb[2673]*m_z1k;

  Box1Abb[2675]=37. + 7.*m_z1k;

  Box1Abb[2676]=-54. + Box1Abb[2675]*m_z1k;

  Box1Abb[2677]=24. + Box1Abb[2676]*m_z1k;

  Box1Abb[2678]=-6. + Box1Abb[2677]*m_z1k;

  Box1Abb[2679]=-Box1Abb[2667]*pow(Box1Abb[68],5.) - Box1Abb[2674]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[2678]*pow(Box1Abb[68],2.)*m_z12_2 + 2.*Box1Abb[2670]*Box1Abb[68]*m_z12_3 + Box1Abb[2551]*m_z12_4;

  Box1Abb[2680]=5. - 10.*m_z1k + 16.*m_z1k_3 - 13.*m_z1k_4;

  Box1Abb[2681]=53. - 75.*m_z1k + 48.*m_z1k_2;

  Box1Abb[2682]=-21. + Box1Abb[2681]*m_z1k;

  Box1Abb[2683]=7. + Box1Abb[2682]*m_z1k;

  Box1Abb[2684]=-46. + 29.*m_z1k;

  Box1Abb[2685]=2. + Box1Abb[2684]*Box1Abb[68]*m_z1k;

  Box1Abb[2686]=-15. + Box1Abb[2685]*m_z1k;

  Box1Abb[2687]=7. + Box1Abb[2686]*m_z1k;

  Box1Abb[2688]=-232. + 53.*m_z1k;

  Box1Abb[2689]=218. + Box1Abb[2688]*m_z1k;

  Box1Abb[2690]=-62. + Box1Abb[2689]*m_z1k;

  Box1Abb[2691]=10. + m_z1k + Box1Abb[2690]*m_z1k_2;

  Box1Abb[2692]=-40. + 3.*m_z1k;

  Box1Abb[2693]=217. + 4.*Box1Abb[2692]*m_z1k;

  Box1Abb[2694]=-110. + Box1Abb[2693]*m_z1k;

  Box1Abb[2695]=33. + Box1Abb[2694]*m_z1k;

  Box1Abb[2696]=-4. + Box1Abb[2695]*m_z1k;

  Box1Abb[2697]=2.*pow(Box1Abb[68],6.) + Box1Abb[2683]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[2696]*pow(Box1Abb[68],2.)*m_z12_2 - Box1Abb[2691]*Box1Abb[68]*m_z12_3 - 2.*Box1Abb[2687]*m_z12_4 + Box1Abb[2680]*m_z12_5;

  Box1Abb[2698]=Box1Abb[2697]*m_x_2 + Box1Abb[2639]*m_x_3 + Box1Abb[2657]*m_x_4 + Box1Abb[2621]*m_x_5 + Box1Abb[2608]*m_x_6 - Box1Abb[2666]*m_x_7 + Box1Abb[2660]*m_x_8 + Box1Abb[2679]*Box1Abb[68]*m_x*m_z12 + Box1Abb[699]*m_x_9*m_z12 + Box1Abb[4]*Box1Abb[541]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_3;

  Box1Abb[2699]=119. + 22.*m_z12;

  Box1Abb[2700]=-81. + Box1Abb[2699]*m_z12;

  Box1Abb[2701]=-88. + Box1Abb[2700]*m_z12;

  Box1Abb[2702]=232. + 31.*m_z12;

  Box1Abb[2703]=81. + Box1Abb[2702]*m_z12;

  Box1Abb[2704]=-208. + Box1Abb[2703]*m_z12;

  Box1Abb[2705]=-5. + Box1Abb[708]*m_z12;

  Box1Abb[2706]=72. + 13.*Box1Abb[2705]*m_z12;

  Box1Abb[2707]=80. + Box1Abb[2701]*m_z12 + 138.*m_z1k + Box1Abb[2704]*m_z12*m_z1k + 2.*Box1Abb[2706]*m_z1k_2 - 336.*Box1Abb[0]*m_z12*m_z1k_3;

  Box1Abb[2708]=7. + 5.*m_z12;

  Box1Abb[2709]=-76. + 5.*Box1Abb[2708]*m_z12;

  Box1Abb[2710]=5. + Box1Abb[2709]*m_z12;

  Box1Abb[2711]=99. + 2.*m_z12;

  Box1Abb[2712]=299. + Box1Abb[2711]*m_z12;

  Box1Abb[2713]=-238. + Box1Abb[2712]*m_z12;

  Box1Abb[2714]=-98. + Box1Abb[2713]*m_z12;

  Box1Abb[2715]=424. + 51.*m_z12;

  Box1Abb[2716]=93. + Box1Abb[2715]*m_z12;

  Box1Abb[2717]=-212. + Box1Abb[2716]*m_z12;

  Box1Abb[2718]=130. + Box1Abb[2717]*m_z12;

  Box1Abb[2719]=418. - 59.*m_z12;

  Box1Abb[2720]=-257. + Box1Abb[2719]*m_z12;

  Box1Abb[2721]=110. + Box1Abb[2720]*m_z12;

  Box1Abb[2722]=50. + 2.*Box1Abb[2710]*m_z12 + 110.*m_z1k + Box1Abb[2714]*m_z12*m_z1k + Box1Abb[2718]*m_z1k_2 + 2.*Box1Abb[2721]*m_z1k_3 - 336.*Box1Abb[0]*m_z12*m_z1k_4;

  Box1Abb[2723]=-7. + 12.*m_z12;

  Box1Abb[2724]=-82. + 5.*Box1Abb[2723]*m_z12;

  Box1Abb[2725]=76. + Box1Abb[2724]*m_z12;

  Box1Abb[2726]=63. + 5.*Box1Abb[488]*m_z12;

  Box1Abb[2727]=-146. + Box1Abb[2726]*m_z12;

  Box1Abb[2728]=44. + Box1Abb[2727]*m_z12;

  Box1Abb[2729]=84. + 5.*m_z12;

  Box1Abb[2730]=335. + 2.*Box1Abb[2729]*m_z12;

  Box1Abb[2731]=-323. + Box1Abb[2730]*m_z12;

  Box1Abb[2732]=30. + Box1Abb[2731]*m_z12;

  Box1Abb[2733]=20. + Box1Abb[2732]*m_z12;

  Box1Abb[2734]=696. + m_z12;

  Box1Abb[2735]=-515. + Box1Abb[2734]*m_z12;

  Box1Abb[2736]=280. + Box1Abb[2735]*m_z12;

  Box1Abb[2737]=-60. + Box1Abb[2736]*m_z12;

  Box1Abb[2738]=65. - 11.*m_z12;

  Box1Abb[2739]=-345. + 8.*Box1Abb[2738]*m_z12;

  Box1Abb[2740]=100. + Box1Abb[2739]*m_z12;

  Box1Abb[2741]=-10. + Box1Abb[2725]*m_z12 + 46.*m_z1k + 2.*Box1Abb[2728]*m_z12*m_z1k + Box1Abb[2733]*m_z1k_2 + Box1Abb[2737]*m_z1k_3 + 2.*Box1Abb[2740]*m_z1k_4 - 168.*Box1Abb[0]*m_z12*m_z1k_5;

  Box1Abb[2742]=-8. + 15.*m_z12 - 60.*m_z1k;

  Box1Abb[2743]=16. + Box1Abb[2742]*m_z12 + 60.*m_z1k;

  Box1Abb[2744]=8. + Box1Abb[2743]*m_z12;

  Box1Abb[2745]=5. + 32.*m_z1k;

  Box1Abb[2746]=70. + 53.*m_z1k;

  Box1Abb[2747]=7. - 48.*m_z1k;

  Box1Abb[2748]=-7. + 4.*Box1Abb[2747]*m_z1k;

  Box1Abb[2749]=-44. + Box1Abb[2748]*m_z12 + Box1Abb[2746]*m_z12_2 + 4.*m_z12_3 + 6.*Box1Abb[2745]*m_z1k;

  Box1Abb[2750]=40. + Box1Abb[2749]*m_z12 + 52.*m_z1k;

  Box1Abb[2751]=3. + 4.*Box1Abb[68]*m_z1k;

  Box1Abb[2752]=2. + Box1Abb[2751]*m_z1k;

  Box1Abb[2753]=-3. + 5.*Box1Abb[15]*m_z1k;

  Box1Abb[2754]=1. + Box1Abb[2753]*m_z1k;

  Box1Abb[2755]=12. + 17.*m_z1k;

  Box1Abb[2756]=-9. + Box1Abb[2755]*m_z1k;

  Box1Abb[2757]=6. + Box1Abb[2756]*m_z1k;

  Box1Abb[2758]=21. + 23.*m_z1k;

  Box1Abb[2759]=-16. + Box1Abb[2758]*m_z1k;

  Box1Abb[2760]=6. + Box1Abb[2759]*m_z1k;

  Box1Abb[2761]=Box1Abb[2752]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[2757]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[2760]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[2754]*m_z12_4 - 2.*Box1Abb[317]*pow(Box1Abb[68],4.)*m_z1k;

  Box1Abb[2762]=1. + Box1Abb[757]*m_z1k;

  Box1Abb[2763]=1. + 3.*Box1Abb[15]*Box1Abb[2762]*m_z1k;

  Box1Abb[2764]=14. + Box1Abb[2440]*m_z1k;

  Box1Abb[2765]=-10. + Box1Abb[2764]*m_z1k;

  Box1Abb[2766]=5. + Box1Abb[2765]*m_z1k;

  Box1Abb[2767]=-38. + 66.*m_z1k + 69.*m_z1k_2;

  Box1Abb[2768]=-30. + Box1Abb[2767]*m_z1k;

  Box1Abb[2769]=27. + Box1Abb[2768]*m_z1k;

  Box1Abb[2770]=-2. + Box1Abb[2769]*m_z1k;

  Box1Abb[2771]=-165. + 84.*m_z1k + 86.*m_z1k_2;

  Box1Abb[2772]=16. + Box1Abb[2771]*m_z1k;

  Box1Abb[2773]=25. + Box1Abb[2772]*m_z1k;

  Box1Abb[2774]=-6. + Box1Abb[2773]*m_z1k;

  Box1Abb[2775]=-109. + 4.*Box1Abb[66]*m_z1k;

  Box1Abb[2776]=31. + Box1Abb[2775]*m_z1k;

  Box1Abb[2777]=4. + Box1Abb[2776]*m_z1k;

  Box1Abb[2778]=-6. + Box1Abb[2777]*m_z1k;

  Box1Abb[2779]=-2.*Box1Abb[2763]*pow(Box1Abb[68],4.)*m_z12 + Box1Abb[2778]*pow(Box1Abb[68],3.)*m_z12_2 + Box1Abb[2774]*pow(Box1Abb[68],2.)*m_z12_3 + Box1Abb[2770]*Box1Abb[68]*m_z12_4 + 2.*Box1Abb[317]*pow(Box1Abb[68],5.)*m_z1k + 2.*Box1Abb[2766]*m_z12_5*m_z1k;

  Box1Abb[2780]=5. + m_z1k;

  Box1Abb[2781]=5. + Box1Abb[2780]*m_z1k;

  Box1Abb[2782]=235. - 27.*m_z1k;

  Box1Abb[2783]=160. + Box1Abb[2782]*m_z1k;

  Box1Abb[2784]=10. + Box1Abb[2783]*m_z1k;

  Box1Abb[2785]=40. + Box1Abb[2784]*m_z1k;

  Box1Abb[2786]=2. + 9.*Box1Abb[536]*m_z1k;

  Box1Abb[2787]=-7. + 2.*Box1Abb[2786]*m_z1k;

  Box1Abb[2788]=7. + Box1Abb[2787]*m_z1k;

  Box1Abb[2789]=310. - 187.*m_z1k;

  Box1Abb[2790]=-106. + Box1Abb[2789]*m_z1k;

  Box1Abb[2791]=12. + Box1Abb[2790]*m_z1k;

  Box1Abb[2792]=9. + Box1Abb[2791]*m_z1k;

  Box1Abb[2793]=28. + 2.*Box1Abb[2792]*m_z1k;

  Box1Abb[2794]=-37. + m_z1k;

  Box1Abb[2795]=133. + 20.*Box1Abb[2794]*m_z1k;

  Box1Abb[2796]=136. + Box1Abb[2795]*m_z1k;

  Box1Abb[2797]=-43. + Box1Abb[2796]*m_z1k;

  Box1Abb[2798]=70. + Box1Abb[2797]*m_z1k;

  Box1Abb[2799]=-125. + 68.*m_z1k;

  Box1Abb[2800]=92. + 3.*Box1Abb[2799]*m_z1k;

  Box1Abb[2801]=-4. + 3.*Box1Abb[2800]*m_z1k;

  Box1Abb[2802]=-88. + Box1Abb[2801]*m_z1k;

  Box1Abb[2803]=17. + Box1Abb[2802]*m_z1k;

  Box1Abb[2804]=2.*Box1Abb[2788]*Box1Abb[68] + Box1Abb[2793]*m_z12 + Box1Abb[2803]*m_z12_2 - Box1Abb[2798]*m_z12_3 + Box1Abb[2785]*m_z12_4 + 4.*Box1Abb[2781]*m_z12_5*m_z1k;

  Box1Abb[2805]=-8. + m_z1k;

  Box1Abb[2806]=-5. + Box1Abb[2805]*m_z1k_2;

  Box1Abb[2807]=1. + m_z1k - 7.*m_z1k_2 + 16.*m_z1k_3;

  Box1Abb[2808]=43. + 24.*m_z1k;

  Box1Abb[2809]=-66. + Box1Abb[2808]*m_z1k;

  Box1Abb[2810]=3. + Box1Abb[2809]*m_z1k;

  Box1Abb[2811]=2. + Box1Abb[2810]*m_z1k;

  Box1Abb[2812]=2. + Box1Abb[2811]*m_z1k;

  Box1Abb[2813]=194. + 41.*m_z1k;

  Box1Abb[2814]=6. + Box1Abb[2813]*m_z1k;

  Box1Abb[2815]=66. + Box1Abb[2814]*m_z1k;

  Box1Abb[2816]=-45. + Box1Abb[2815]*m_z1k;

  Box1Abb[2817]=14. + Box1Abb[2816]*m_z1k;

  Box1Abb[2818]=17. + 6.*m_z1k;

  Box1Abb[2819]=-581. + 8.*Box1Abb[2818]*m_z1k;

  Box1Abb[2820]=300. + Box1Abb[2819]*m_z1k;

  Box1Abb[2821]=-64. + Box1Abb[2820]*m_z1k;

  Box1Abb[2822]=20. + Box1Abb[2821]*m_z1k;

  Box1Abb[2823]=-27. + Box1Abb[2822]*m_z1k;

  Box1Abb[2824]=140. + 59.*m_z1k;

  Box1Abb[2825]=-595. + 2.*Box1Abb[2824]*m_z1k;

  Box1Abb[2826]=212. + Box1Abb[2825]*m_z1k;

  Box1Abb[2827]=-128. + Box1Abb[2826]*m_z1k;

  Box1Abb[2828]=68. + Box1Abb[2827]*m_z1k;

  Box1Abb[2829]=-35. + Box1Abb[2828]*m_z1k;

  Box1Abb[2830]=2.*Box1Abb[2807]*pow(Box1Abb[68],3.) - 2.*Box1Abb[2812]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[2823]*Box1Abb[68]*m_z12_2 + Box1Abb[2829]*m_z12_3 + Box1Abb[2817]*m_z12_4 - 4.*Box1Abb[2806]*m_z12_5*m_z1k;

  Box1Abb[2831]=-Box1Abb[2779]*m_x_2 + Box1Abb[2830]*m_x_3 - Box1Abb[2804]*m_x_4 + Box1Abb[2741]*m_x_5 - Box1Abb[2722]*m_x_6 + Box1Abb[2707]*m_x_7 - Box1Abb[2750]*m_x_8 + Box1Abb[2744]*m_x_9 + 8.*Box1Abb[0]*m_x_10*m_z12 + Box1Abb[2761]*pow(Box1Abb[68],2.)*m_x*m_z12*m_z1k + Box1Abb[4]*pow(Box1Abb[68],5.)*m_z12_3*m_z1k_3;

  Box1Abb[2902]=(4.*m_cL*m_x_2)/(m_s_2*m_z1k);

  Box1Abb[2903]=(-4.*m_cL*m_x)/(m_s_2*m_z1k);

  Box1Abb[2904]=(4.*m_cL)/m_s_2;

  Box1Abb[2905]=(4.*m_cL*m_x)/m_s_2;

  Box1Abb[2906]=(2.*m_cL*m_z12)/m_s;

  Box1Abb[2907]=(-4.*m_cL*m_z12)/m_s_2;

  Box1Abb[2908]=(-4.*m_cL)/m_s_2;

  Box1Abb[2909]=m_cL*m_z12_2;

  Box1Abb[2910]=(4.*m_cL*m_z12_2)/m_s;

  Box1Abb[2911]=(-4.*m_cL*m_z12)/m_s;

  Box1Abb[2912]=(4.*m_cL*m_x*m_z12)/m_s;

  Box1Abb[2913]=(-4.*m_cL*m_x*m_z12)/m_s;

  Box1Abb[2914]=m_cL*m_s*m_z12_2;

  Box1Abb[2915]=4.*m_cL*m_z12_2;

  Box1Abb[2916]=-2.*m_cL;

  Box1Abb[2917]=-4.*m_cL*m_z12;

  Box1Abb[2918]=4.*m_cL*m_x*m_z12;

  Box1Abb[2919]=-4.*m_cL*m_x*m_z12;

  Box1Abb[2920]=4.*m_cL*m_z12;

  Box1Abb[2921]=2.*m_cL;

  Box1Abb[2922]=4.*m_cL;

  Box1Abb[2923]=4.*m_cL*m_x;

  Box1Abb[2924]=(4.*m_cL)/m_s;

  Box1Abb[2925]=(4.*m_cL*m_x)/m_s;

  Box1Abb[2926]=(-4.*m_cL*m_x)/m_s_2;

  Box1Abb[2927]=(-8.*m_cL*m_x)/m_s_2;

  Box1Abb[2928]=(-4.*m_cL*m_z1k)/m_s_2;

  Box1Abb[2929]=(-4.*m_cL)/m_s;

  Box1Abb[2930]=(-4.*m_cL*m_x)/m_s;

  Box1Abb[3001]=(4.*m_cR*m_x_2)/(m_s_2*m_z1k);

  Box1Abb[3002]=(-4.*m_cR*m_x)/(m_s_2*m_z1k);

  Box1Abb[3003]=(4.*m_cR)/m_s_2;

  Box1Abb[3004]=(4.*m_cR*m_x)/m_s_2;

  Box1Abb[3005]=(2.*m_cR*m_z12)/m_s;

  Box1Abb[3006]=(-4.*m_cR*m_z12)/m_s_2;

  Box1Abb[3007]=(-4.*m_cR)/m_s_2;

  Box1Abb[3008]=m_cR*m_z12_2;

  Box1Abb[3009]=(4.*m_cR*m_z12_2)/m_s;

  Box1Abb[3010]=(-4.*m_cR*m_z12)/m_s;

  Box1Abb[3011]=(4.*m_cR*m_x*m_z12)/m_s;

  Box1Abb[3012]=(-4.*m_cR*m_x*m_z12)/m_s;

  Box1Abb[3013]=m_cR*m_s*m_z12_2;

  Box1Abb[3014]=4.*m_cR*m_z12_2;

  Box1Abb[3015]=-2.*m_cR;

  Box1Abb[3016]=-4.*m_cR*m_z12;

  Box1Abb[3017]=4.*m_cR*m_x*m_z12;

  Box1Abb[3018]=-4.*m_cR*m_x*m_z12;

  Box1Abb[3019]=4.*m_cR*m_z12;

  Box1Abb[3020]=2.*m_cR;

  Box1Abb[3021]=4.*m_cR;

  Box1Abb[3022]=4.*m_cR*m_x;

  Box1Abb[3023]=(4.*m_cR)/m_s;

  Box1Abb[3024]=(4.*m_cR*m_x)/m_s;

  Box1Abb[3025]=(-4.*m_cR*m_x)/m_s_2;

  Box1Abb[3026]=(-8.*m_cR*m_x)/m_s_2;

  Box1Abb[3027]=(-4.*m_cR*m_z1k)/m_s_2;

  Box1Abb[3028]=(-4.*m_cR)/m_s;

  Box1Abb[3029]=(-4.*m_cR*m_x)/m_s;

  Box1Abb[3030]=m_z12 + 3.*m_z1k;

  Box1Abb[3031]=-3. + 4.*m_z12 + 3.*m_z1k;

  Box1Abb[3032]=-1. + Box1Abb[119]*m_z12 + 3.*Box1Abb[174]*m_z1k;

  Box1Abb[3033]=1. + Box1Abb[1028]*m_z1k;

  Box1Abb[3034]=-1. + m_z12 + 3.*Box1Abb[3033]*m_z1k + Box1Abb[692]*m_z12*m_z1k;

  Box1Abb[3035]=2. + 3.*Box1Abb[834]*m_z1k;

  Box1Abb[3036]=2. + 9.*m_z1k;

  Box1Abb[3037]=3. + Box1Abb[3036]*m_z1k;

  Box1Abb[3038]=Box1Abb[3035]*Box1Abb[68] + Box1Abb[3037]*m_z12;

  Box1Abb[3039]=-Box1Abb[3034]*m_x_2 + Box1Abb[3032]*m_x_3 - Box1Abb[3030]*m_x_4 + Box1Abb[3038]*m_x*m_z1k - Box1Abb[3031]*pow(Box1Abb[68],2.)*m_z1k_2;

  Box1Abb[3040]=2.*Box1Abb[2]*m_cL + Box1Abb[79]*m_cR;

  Box1Abb[3041]=-Box1Abb[379]*Box1Abb[79]*m_x + Box1Abb[3040]*m_z12*m_z1k + 2.*Box1Abb[56]*m_z12*m_z1k_2;

  Box1Abb[3042]=-5. - 4.*m_x + 2.*m_z12;

  Box1Abb[3043]=4. + Box1Abb[3042]*m_z12;

  Box1Abb[3044]=-1. + Box1Abb[3043]*m_x + m_z12;

  Box1Abb[3045]=4. + Box1Abb[129]*m_z12_2;

  Box1Abb[3046]=Box1Abb[3045]*m_x - 2.*Box1Abb[0]*Box1Abb[79]*m_z12 - 12.*m_x_2*m_z12;

  Box1Abb[3047]=3. - 4.*m_z12;

  Box1Abb[3048]=-2. + Box1Abb[3047]*m_z12;

  Box1Abb[3049]=Box1Abb[3048]*m_x + Box1Abb[0]*Box1Abb[79]*m_z12 + 12.*m_x_2*m_z12;

  Box1Abb[3050]=Box1Abb[988]*m_x + Box1Abb[133]*m_z12;

  Box1Abb[3051]=Box1Abb[3044]*Box1Abb[79]*m_x_2 + Box1Abb[3046]*m_x*m_z1k + Box1Abb[3049]*m_z12*m_z1k_2 + Box1Abb[3050]*m_z12*m_z1k_3;

  Box1Abb[3052]=2. + 2.*Box1Abb[492]*m_z12 - 6.*m_z1k + 3.*Box1Abb[222]*m_z12*m_z1k + 4.*Box1Abb[14]*m_z1k_2;

  Box1Abb[3053]=11. - 12.*m_z12 - 48.*m_z1k;

  Box1Abb[3054]=2.*Box1Abb[1172] + Box1Abb[3053]*m_z12;

  Box1Abb[3055]=8. + Box1Abb[3054]*m_z12;

  Box1Abb[3056]=7. + 2.*Box1Abb[79]*m_z12 - 2.*m_z1k + 23.*m_z12*m_z1k + 36.*m_z1k_2;

  Box1Abb[3057]=-7. + Box1Abb[3056]*m_z12 - 8.*m_z1k;

  Box1Abb[3058]=2. + Box1Abb[3057]*m_z12 - 4.*m_z1k;

  Box1Abb[3059]=Box1Abb[3058]*m_x_2 + Box1Abb[3055]*m_x_3 + 4.*Box1Abb[505]*m_x_4*m_z12 - Box1Abb[3052]*m_x*m_z12*m_z1k + 3.*Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[3060]=Box1Abb[3059]*m_cL + Box1Abb[3051]*m_cR;

  Box1Abb[3061]=10. - 3.*m_z12;

  Box1Abb[3062]=10. + 11.*m_z12;

  Box1Abb[3063]=-2. + 2.*m_z12 - 9.*m_z1k + 3.*Box1Abb[160]*m_z12*m_z1k + Box1Abb[3062]*m_z1k_2;

  Box1Abb[3064]=3. - 7.*m_z1k;

  Box1Abb[3065]=4. + 9.*m_z1k;

  Box1Abb[3066]=4. - Box1Abb[3065]*m_z12 + Box1Abb[3064]*m_z1k;

  Box1Abb[3067]=9. + m_z1k;

  Box1Abb[3068]=3. - Box1Abb[3067]*m_z12 + 2.*m_z12_2 - 20.*m_z1k;

  Box1Abb[3069]=2. + Box1Abb[3068]*m_z12;

  Box1Abb[3070]=Box1Abb[3069]*m_x_3 + Box1Abb[3063]*m_x_2*m_z12 + Box1Abb[3061]*m_x_4*m_z12 + Box1Abb[3066]*m_x*m_z12_2*m_z1k + 2.*Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[3071]=-10. + 13.*m_z12;

  Box1Abb[3072]=-7. + 19.*m_z12;

  Box1Abb[3073]=-10. + 19.*m_z12;

  Box1Abb[3074]=pow(Box1Abb[0],2.)*Box1Abb[14] + Box1Abb[0]*Box1Abb[3072]*m_z1k + Box1Abb[3073]*m_z1k_2;

  Box1Abb[3075]=-14. + 11.*m_z12 + 29.*m_z1k;

  Box1Abb[3076]=1. + Box1Abb[3075]*m_z12 - 20.*m_z1k;

  Box1Abb[3077]=2. + Box1Abb[3076]*m_z12;

  Box1Abb[3078]=-Box1Abb[3077]*m_x_2 + Box1Abb[3074]*m_x*m_z12 + Box1Abb[3071]*m_x_3*m_z12 - 3.*pow(Box1Abb[4],2.)*m_z12_2*m_z1k;

  Box1Abb[3079]=Box1Abb[3070]*m_cR + Box1Abb[3078]*m_cL*m_x;

  Box1Abb[3080]=-10. + 11.*m_z12;

  Box1Abb[3081]=2.*pow(Box1Abb[0],2.) + 5.*Box1Abb[0]*m_z1k + m_z1k_2;

  Box1Abb[3082]=2. + 2.*Box1Abb[79]*m_z12 + m_z1k + Box1Abb[1438]*m_z12*m_z1k + Box1Abb[3071]*m_z1k_2;

  Box1Abb[3083]=5. + m_z12 - 4.*m_z12_2 + 20.*m_z1k - 23.*m_z12*m_z1k;

  Box1Abb[3084]=-2. + Box1Abb[3083]*m_z12;

  Box1Abb[3085]=Box1Abb[3084]*m_x_2 + Box1Abb[3082]*m_x*m_z12 + Box1Abb[3080]*m_x_3*m_z12 - Box1Abb[3081]*m_z12_2*m_z1k;

  Box1Abb[3086]=-3. + 2.*Box1Abb[2661]*m_z12;

  Box1Abb[3087]=10. + 17.*m_z12;

  Box1Abb[3088]=-1. + m_z12 + Box1Abb[3086]*m_z1k + Box1Abb[3087]*m_z1k_2;

  Box1Abb[3089]=2. + 7.*m_z1k;

  Box1Abb[3090]=-2. + Box1Abb[3089]*m_z12 + Box1Abb[3036]*m_z1k;

  Box1Abb[3091]=3. + Box1Abb[3089]*m_z12 + 20.*m_z1k;

  Box1Abb[3092]=-2. + Box1Abb[3091]*m_z12;

  Box1Abb[3093]=-Box1Abb[3092]*m_x_3 + Box1Abb[3088]*m_x_2*m_z12 - Box1Abb[803]*m_x_4*m_z12 - Box1Abb[3090]*m_x*m_z12_2*m_z1k + Box1Abb[367]*m_z12_3*m_z1k_2;

  Box1Abb[3094]=Box1Abb[3093]*m_cR + Box1Abb[3085]*m_cL*m_x;

  Box1Abb[3095]=-27. + 8.*m_z12;

  Box1Abb[3096]=26. + Box1Abb[3095]*m_z12 - 12.*m_z1k;

  Box1Abb[3097]=-8. + Box1Abb[3096]*m_z12;

  Box1Abb[3098]=2. - 2.*Box1Abb[1640]*m_z12 - 6.*m_z1k + Box1Abb[1018]*m_z12*m_z1k + 4.*Box1Abb[1408]*m_z1k_2;

  Box1Abb[3099]=-11. - 2.*Box1Abb[158]*m_z12 + 2.*m_z1k + 5.*m_z12*m_z1k + 12.*m_z1k_2;

  Box1Abb[3100]=7. + Box1Abb[3099]*m_z12 - 8.*m_z1k;

  Box1Abb[3101]=-2. + Box1Abb[3100]*m_z12 + 4.*m_z1k;

  Box1Abb[3102]=Box1Abb[3101]*m_x_2 + Box1Abb[3097]*m_x_3 - 4.*Box1Abb[79]*m_x_4*m_z12 + Box1Abb[3098]*m_x*m_z12*m_z1k + Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[3103]=3. + Box1Abb[69]*m_z12;

  Box1Abb[3104]=-2. + Box1Abb[3103]*m_z12 - 6.*m_z1k + 3.*Box1Abb[677]*m_z12*m_z1k + 4.*Box1Abb[14]*m_z1k_2;

  Box1Abb[3105]=45. - 34.*m_z12 - 48.*m_z1k;

  Box1Abb[3106]=-22. + Box1Abb[3105]*m_z12 + 12.*m_z1k;

  Box1Abb[3107]=8. + Box1Abb[3106]*m_z12;

  Box1Abb[3108]=-29. + 49.*m_z1k;

  Box1Abb[3109]=25. + Box1Abb[3108]*m_z12 + 15.*m_z12_2 + 36.*Box1Abb[68]*m_z1k;

  Box1Abb[3110]=-17. + Box1Abb[3109]*m_z12;

  Box1Abb[3111]=6. + Box1Abb[3110]*m_z12 - 4.*m_z1k;

  Box1Abb[3112]=Box1Abb[3111]*m_x_2 + Box1Abb[3107]*m_x_3 - Box1Abb[3104]*Box1Abb[4]*m_x*m_z12 + 4.*Box1Abb[505]*m_x_4*m_z12 + pow(Box1Abb[4],2.)*Box1Abb[77]*m_z12_2*m_z1k;

  Box1Abb[3113]=Box1Abb[3112]*m_cL + Box1Abb[3102]*m_cR;

  Box1Abb[3114]=-1. + m_z12 + 5.*m_z1k;

  Box1Abb[3115]=-pow(Box1Abb[0],2.)*Box1Abb[69] - 3.*Box1Abb[0]*Box1Abb[856]*m_z1k + 5.*Box1Abb[77]*m_z1k_2;

  Box1Abb[3116]=9. - 18.*m_z12 + 7.*m_z12_2 + 5.*Box1Abb[158]*m_z1k;

  Box1Abb[3117]=2. + Box1Abb[3116]*m_z12;

  Box1Abb[3118]=Box1Abb[3117]*m_x_2 + Box1Abb[3115]*m_x*m_z12 - 5.*Box1Abb[79]*m_x_3*m_z12 - Box1Abb[3114]*Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[3119]=1. + Box1Abb[69]*m_z12 - 8.*m_z1k + 11.*m_z12*m_z1k + 5.*m_z1k_2;

  Box1Abb[3120]=15. + 13.*Box1Abb[79]*m_z12;

  Box1Abb[3121]=-26. + 21.*m_z12;

  Box1Abb[3122]=-2. + Box1Abb[3120]*m_z12 + 13.*m_z1k + 2.*Box1Abb[3121]*m_z12*m_z1k + 5.*Box1Abb[505]*m_z1k_2;

  Box1Abb[3123]=-33. + 25.*m_z12 + 35.*m_z1k;

  Box1Abb[3124]=7. + Box1Abb[3123]*m_z12 - 20.*m_z1k;

  Box1Abb[3125]=2. + Box1Abb[3124]*m_z12;

  Box1Abb[3126]=-Box1Abb[3125]*m_x_3 + Box1Abb[3122]*m_x_2*m_z12 + 5.*Box1Abb[133]*m_x_4*m_z12 - Box1Abb[3119]*Box1Abb[4]*m_x*m_z12_2 + pow(Box1Abb[4],2.)*m_z12_3*m_z1k;

  Box1Abb[3127]=Box1Abb[3126]*m_cL + Box1Abb[3118]*m_cR*m_x;

  Box1Abb[3128]=-10. + m_z12 + 4.*m_z12_2;

  Box1Abb[3129]=-34. + 33.*m_z12;

  Box1Abb[3130]=5. + Box1Abb[3128]*m_z12 + 7.*m_z1k + Box1Abb[3129]*m_z12*m_z1k + Box1Abb[3073]*m_z1k_2;

  Box1Abb[3131]=2. + 13.*m_z1k;

  Box1Abb[3132]=-2. + Box1Abb[3131]*m_z12 + Box1Abb[2303]*m_z1k;

  Box1Abb[3133]=-18. + 17.*m_z12 + 29.*m_z1k;

  Box1Abb[3134]=1. + Box1Abb[3133]*m_z12 - 20.*m_z1k;

  Box1Abb[3135]=2. + Box1Abb[3134]*m_z12;

  Box1Abb[3136]=-Box1Abb[3135]*m_x_3 + Box1Abb[3130]*m_x_2*m_z12 + Box1Abb[3071]*m_x_4*m_z12 - Box1Abb[3132]*Box1Abb[4]*m_x*m_z12_2 + 2.*pow(Box1Abb[4],2.)*m_z12_3*m_z1k;

  Box1Abb[3137]=Box1Abb[3136]*m_cL + Box1Abb[3070]*m_cR;

  Box1Abb[3138]=-Box1Abb[1279]*m_cR + Box1Abb[1277]*m_cL*m_x;

  Box1Abb[3139]=1. - m_z12 - 2.*m_z1k - 3.*Box1Abb[0]*m_z12*m_z1k + m_z1k_2;

  Box1Abb[3140]=Box1Abb[3139]*m_x + Box1Abb[133]*Box1Abb[15]*m_x_2 + Box1Abb[1408]*m_x_3 - Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[3141]=1. + Box1Abb[457]*m_z12 - m_z1k;

  Box1Abb[3142]=Box1Abb[3141]*Box1Abb[4] + 2.*Box1Abb[15]*m_x + Box1Abb[1422]*m_x_2 + Box1Abb[2293]*m_x*m_z12;

  Box1Abb[3143]=Box1Abb[3140]*m_cR + Box1Abb[3142]*m_cL*m_x;

  Box1Abb[3144]=2.*Box1Abb[183]*m_x + m_z12 - 3.*m_x*m_z12 + Box1Abb[1]*m_z12_2;

  Box1Abb[3145]=m_x - 3.*Box1Abb[0]*m_z12 + 3.*m_x*m_z12;

  Box1Abb[3146]=-Box1Abb[0]*pow(Box1Abb[1],2.)*m_x - Box1Abb[3144]*m_z1k + Box1Abb[3145]*m_z1k_2 - 2.*m_z12*m_z1k_3;

  Box1Abb[3147]=2. + 3.*Box1Abb[0]*m_z12;

  Box1Abb[3148]=-1. + m_z12 + Box1Abb[3147]*m_z1k - m_z1k_2;

  Box1Abb[3149]=Box1Abb[3148]*m_x - Box1Abb[133]*Box1Abb[15]*m_x_2 + Box1Abb[8]*m_x_3 + Box1Abb[4]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[3150]=Box1Abb[3149]*m_cL + Box1Abb[3146]*m_cR;

  Box1Abb[3151]=3. + 4.*m_x - 3.*m_z12;

  Box1Abb[3152]=-1. + m_x - 2.*m_x_2 + m_z12 - m_x*m_z12 + Box1Abb[3151]*m_z1k - 2.*m_z1k_2;

  Box1Abb[3153]=Box1Abb[420]*m_cL + Box1Abb[3152]*m_cR;

  Box1Abb[3154]=2. + m_z12 + m_z12_2 + 4.*m_z1k - 5.*m_z12*m_z1k;

  Box1Abb[3155]=1. + m_z12 + m_z12_2;

  Box1Abb[3156]=1. + m_z12 + m_z12_2 + 2.*Box1Abb[3155]*m_z1k + 2.*Box1Abb[208]*m_z1k_2;

  Box1Abb[3157]=-5. + 2.*Box1Abb[474]*m_z1k;

  Box1Abb[3158]=1. + Box1Abb[317]*m_z1k;

  Box1Abb[3159]=m_z12 + 2.*Box1Abb[3158]*m_z1k + Box1Abb[3157]*m_z12*m_z1k + 4.*m_z12_2*m_z1k;

  Box1Abb[3160]=-1. + Box1Abb[873]*m_z1k;

  Box1Abb[3161]=Box1Abb[3160]*m_z12_2 - pow(Box1Abb[68],2.)*m_z1k + Box1Abb[438]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[3162]=Box1Abb[3159]*m_x_2 - Box1Abb[3156]*m_x_3 + Box1Abb[3154]*m_x_4 + Box1Abb[0]*m_x_5 + Box1Abb[3161]*m_x*m_z1k - Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3163]=Box1Abb[3162]*m_cL - Box1Abb[3142]*pow(Box1Abb[7],2.)*m_cR*m_x;

  Box1Abb[3164]=1. + m_x - m_z12;

  Box1Abb[3165]=1. + Box1Abb[3164]*m_z12;

  Box1Abb[3166]=-2. - 4.*m_x + m_z12;

  Box1Abb[3167]=-pow(Box1Abb[79],2.)*m_x_2 + 2.*Box1Abb[3165]*m_x*m_z1k + Box1Abb[3166]*m_z12*m_z1k_2 + 2.*m_z12*m_z1k_3;

  Box1Abb[3168]=-2. + 5.*m_z12 + 2.*m_z1k;

  Box1Abb[3169]=4. - 2.*Box1Abb[175]*m_z12 + 5.*m_z12_2;

  Box1Abb[3170]=12. + 2.*Box1Abb[542]*m_z12 - 11.*m_z12*m_z1k - 6.*m_z1k_2;

  Box1Abb[3171]=2.*Box1Abb[68] + Box1Abb[3170]*m_z12;

  Box1Abb[3172]=8. + 3.*m_z1k;

  Box1Abb[3173]=4. + m_z12 - m_z12_2 + 2.*Box1Abb[3172]*m_z1k - m_z12*m_z1k;

  Box1Abb[3174]=-4. + Box1Abb[3173]*m_z12 - 6.*m_z1k;

  Box1Abb[3175]=Box1Abb[3174]*m_x_2 + Box1Abb[3169]*m_x_3 + Box1Abb[3171]*m_x*m_z1k + Box1Abb[3168]*Box1Abb[4]*m_z12*m_z1k_2;

  Box1Abb[3176]=Box1Abb[3175]*m_cL + Box1Abb[3167]*Box1Abb[7]*m_cR;

  Box1Abb[3177]=2. + Box1Abb[145]*m_z12;

  Box1Abb[3178]=-Box1Abb[3177]*m_x_2 + 2.*m_x_3*m_z12 + Box1Abb[12]*m_x*m_z1k + Box1Abb[68]*m_z12*m_z1k_2;

  Box1Abb[3179]=-15. + 2.*m_z12;

  Box1Abb[3180]=15. + Box1Abb[3179]*m_z12;

  Box1Abb[3181]=-2. + Box1Abb[3180]*m_z12 + 6.*m_z1k - 6.*Box1Abb[77]*m_z12*m_z1k - 4.*m_z1k_2;

  Box1Abb[3182]=4. + Box1Abb[0]*Box1Abb[69]*m_z12 + 30.*m_z1k + 6.*Box1Abb[129]*m_z12*m_z1k + 12.*Box1Abb[1422]*m_z1k_2;

  Box1Abb[3183]=4. - Box1Abb[3182]*m_z12 + 12.*m_z1k;

  Box1Abb[3184]=-24. + 11.*m_z12 + 36.*m_z1k;

  Box1Abb[3185]=24. + Box1Abb[3184]*m_z12 - 20.*m_z1k;

  Box1Abb[3186]=-8. + Box1Abb[3185]*m_z12;

  Box1Abb[3187]=2. + m_z1k + m_z1k_2;

  Box1Abb[3188]=-18. + 12.*Box1Abb[3187]*m_z12 + Box1Abb[761]*m_z12_2 + 3.*m_z12_3 + 4.*m_z1k_2;

  Box1Abb[3189]=2. + Box1Abb[3188]*m_z12 - 4.*m_z1k;

  Box1Abb[3190]=Box1Abb[3183]*m_x_3 + Box1Abb[3186]*m_x_4 + 4.*Box1Abb[1010]*m_x_5*m_z12 + Box1Abb[3189]*m_x_2*m_z1k + Box1Abb[3181]*m_x*m_z12*m_z1k_2 + 3.*Box1Abb[4]*m_z12_3*m_z1k_3;

  Box1Abb[3191]=Box1Abb[3190]*m_cL + 2.*Box1Abb[3178]*Box1Abb[384]*Box1Abb[7]*m_cR;

  Box1Abb[3192]=3. + m_x - m_z12;

  Box1Abb[3193]=-3. + Box1Abb[3192]*m_z12;

  Box1Abb[3194]=-3. + 2.*m_x + m_z12;

  Box1Abb[3195]=m_x - m_z12;

  Box1Abb[3196]=Box1Abb[3193]*m_x_2 - Box1Abb[3194]*m_x*m_z12*m_z1k + Box1Abb[3195]*m_z12*m_z1k_2;

  Box1Abb[3197]=pow(Box1Abb[0],2.)*Box1Abb[505] + 3.*Box1Abb[0]*Box1Abb[1422]*m_z1k + 3.*m_z12*m_z1k_2;

  Box1Abb[3198]=-3. + 2.*m_z12 + m_z1k;

  Box1Abb[3199]=1. + Box1Abb[3198]*m_z12;

  Box1Abb[3200]=Box1Abb[3197]*m_x - 3.*Box1Abb[3199]*m_x_2 - pow(Box1Abb[4],3.)*m_z12 + m_x_3*m_z12;

  Box1Abb[3201]=-Box1Abb[117]*Box1Abb[3196]*m_cR + Box1Abb[3200]*m_cL*m_x;

  Box1Abb[3202]=pow(Box1Abb[0],3.) + 3.*pow(Box1Abb[0],2.)*m_z1k + 2.*Box1Abb[0]*m_z1k_2 + m_z1k_3;

  Box1Abb[3203]=1. + 5.*m_z12;

  Box1Abb[3204]=-pow(Box1Abb[0],3.) - pow(Box1Abb[0],2.)*Box1Abb[3203]*m_z1k - 3.*Box1Abb[0]*Box1Abb[1422]*m_z1k_2 - 4.*m_z12*m_z1k_3;

  Box1Abb[3205]=1. + 3.*Box1Abb[834]*m_z12_2 + 2.*m_z12_3 + 6.*m_z1k + 6.*Box1Abb[881]*m_z12*m_z1k;

  Box1Abb[3206]=-8. + 5.*m_z12 + 4.*m_z1k;

  Box1Abb[3207]=3. + Box1Abb[3206]*m_z12;

  Box1Abb[3208]=Box1Abb[3204]*m_x + Box1Abb[3205]*m_x_2 - Box1Abb[3207]*m_x_3 + m_x_4*m_z12 + Box1Abb[3202]*m_z12*m_z1k;

  Box1Abb[3209]=-Box1Abb[3196]*pow(Box1Abb[7],2.)*m_cR + Box1Abb[3208]*m_cL*m_x;

  Box1Abb[3210]=5. + 2.*Box1Abb[448]*m_x;

  Box1Abb[3211]=5. + 7.*m_x;

  Box1Abb[3212]=-2. - 4.*m_x + Box1Abb[3210]*m_z12 - Box1Abb[3211]*m_z12_2 + 2.*m_z12_3;

  Box1Abb[3213]=2. + m_z12_3;

  Box1Abb[3214]=15. + Box1Abb[1800]*m_z12;

  Box1Abb[3215]=-6. + Box1Abb[3214]*m_z12;

  Box1Abb[3216]=Box1Abb[0]*Box1Abb[3213] + 2.*Box1Abb[3215]*m_x - 4.*Box1Abb[856]*m_x_2*m_z12;

  Box1Abb[3217]=1. - 6.*m_x_2;

  Box1Abb[3218]=3. + 14.*m_x;

  Box1Abb[3219]=-3.*pow(Box1Abb[1400],2.) + 4.*m_x + 2.*Box1Abb[3217]*m_z12 + Box1Abb[3218]*m_z12_3 - 2.*m_z12_4;

  Box1Abb[3220]=3. + Box1Abb[699]*m_z12;

  Box1Abb[3221]=-2.*Box1Abb[3220]*m_x + 4.*Box1Abb[1086]*m_x_2 + Box1Abb[0]*m_z12_2;

  Box1Abb[3222]=Box1Abb[988]*m_x + m_z12_2;

  Box1Abb[3223]=Box1Abb[3212]*Box1Abb[79]*m_x_3 + Box1Abb[3216]*m_x_2*m_z1k + Box1Abb[3219]*m_x*m_z1k_2 + Box1Abb[3221]*m_z12*m_z1k_3 + Box1Abb[3222]*m_z12*m_z1k_4;

  Box1Abb[3224]=-pow(Box1Abb[0],3.) + Box1Abb[386]*m_z1k + m_z12*m_z1k_2;

  Box1Abb[3225]=-6. + 4.*m_z12 + 5.*m_z1k;

  Box1Abb[3226]=2. + Box1Abb[3225]*m_z12;

  Box1Abb[3227]=-8. + 5.*m_z12;

  Box1Abb[3228]=pow(Box1Abb[0],2.) + Box1Abb[3227]*m_z1k + 3.*m_z1k_2;

  Box1Abb[3229]=Box1Abb[3228]*m_z12 + 3.*m_z1k;

  Box1Abb[3230]=Box1Abb[3229]*m_x_2 - Box1Abb[3226]*m_x_3 + 2.*m_x_4*m_z12 + Box1Abb[3224]*m_x*m_z1k - pow(Box1Abb[4],2.)*m_z12*m_z1k_2;

  Box1Abb[3231]=-2.*Box1Abb[3230]*Box1Abb[383]*m_cL + Box1Abb[3223]*m_cR;

  Box1Abb[3232]=-Box1Abb[1640]*m_x + 2.*m_x_2*m_z12 - Box1Abb[2078]*m_z12*m_z1k;

  Box1Abb[3233]=-2. + 7.*m_z12;

  Box1Abb[3234]=2. + m_z12 - 2.*m_z1k;

  Box1Abb[3235]=2. + Box1Abb[0]*m_z12;

  Box1Abb[3236]=-2.*Box1Abb[15] + Box1Abb[3235]*m_z12 + 2.*Box1Abb[1086]*m_z12*m_z1k + 2.*Box1Abb[77]*m_z1k_2;

  Box1Abb[3237]=-2. + 7.*m_z12 + 18.*m_z1k;

  Box1Abb[3238]=2. + Box1Abb[3237]*m_z12;

  Box1Abb[3239]=4. - Box1Abb[3238]*m_z12;

  Box1Abb[3240]=Box1Abb[3239]*m_x_2 + Box1Abb[3236]*m_x*m_z12 + 2.*Box1Abb[3233]*m_x_3*m_z12 - Box1Abb[3234]*Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[3241]=Box1Abb[3240]*m_cL - Box1Abb[3232]*Box1Abb[384]*m_cR;

  Box1Abb[3242]=3. + Box1Abb[156]*m_z1k;

  Box1Abb[3243]=3. + Box1Abb[1003]*m_z12;

  Box1Abb[3244]=3. + 3.*Box1Abb[68]*m_z12 + m_z12_2 + 4.*Box1Abb[881]*m_z1k;

  Box1Abb[3245]=-4. + Box1Abb[3244]*m_z12;

  Box1Abb[3246]=6. - Box1Abb[411]*m_z12 + m_z12_2;

  Box1Abb[3247]=-4. + Box1Abb[3246]*m_z12 + 12.*m_z1k;

  Box1Abb[3248]=2. + Box1Abb[3247]*m_z12;

  Box1Abb[3249]=Box1Abb[3248]*m_x_3 - 4.*Box1Abb[3243]*m_x_4*m_z12 + 2.*m_x_5*m_z12_2 + Box1Abb[3245]*m_x_2*m_z12*m_z1k + Box1Abb[3242]*m_x*m_z12_2*m_z1k_2 + Box1Abb[68]*m_z12_3*m_z1k_3;

  Box1Abb[3250]=12. + 4.*Box1Abb[2805]*m_z12 + 15.*m_z12_2;

  Box1Abb[3251]=-34. + 15.*m_z12;

  Box1Abb[3252]=12. + Box1Abb[3251]*m_z12;

  Box1Abb[3253]=4. + 3.*Box1Abb[158]*m_z12;

  Box1Abb[3254]=Box1Abb[3252]*m_z12 + 6.*Box1Abb[3253]*m_z1k;

  Box1Abb[3255]=4. + Box1Abb[3254]*m_z12;

  Box1Abb[3256]=3. + 2.*Box1Abb[170]*m_z1k;

  Box1Abb[3257]=1. - 6.*Box1Abb[68]*m_z1k;

  Box1Abb[3258]=2. + 3.*Box1Abb[3257]*m_z12 - 3.*Box1Abb[3089]*m_z12_2 + m_z12_3 - 2.*Box1Abb[3256]*m_z1k;

  Box1Abb[3259]=-3. + 6.*Box1Abb[228]*m_z1k;

  Box1Abb[3260]=2. + Box1Abb[3259]*m_z12 + Box1Abb[2275]*m_z12_2 - 2.*m_z12_3 - 30.*m_z1k + 8.*m_z1k_3;

  Box1Abb[3261]=-2. + Box1Abb[3260]*m_z12;

  Box1Abb[3262]=Box1Abb[3255]*m_x_3 + Box1Abb[3261]*m_x_2*m_z12 - 2.*Box1Abb[3250]*m_x_4*m_z12 + 4.*m_x_5*m_z12_2 + Box1Abb[3258]*m_x*m_z12_2*m_z1k + 3.*Box1Abb[4]*m_z12_4*m_z1k_2;

  Box1Abb[3263]=Box1Abb[3262]*m_cL - 2.*Box1Abb[3249]*m_cR;

  Box1Abb[3264]=-1. + m_z12 + 6.*m_z1k;

  Box1Abb[3265]=1. + Box1Abb[3264]*m_z12;

  Box1Abb[3266]=-Box1Abb[3265]*m_x_2 + 3.*m_x_3*m_z12 + Box1Abb[370]*m_x*m_z12*m_z1k - m_z12_2*m_z1k_2;

  Box1Abb[3267]=-6. + 7.*m_z12;

  Box1Abb[3268]=pow(Box1Abb[0],2.)*Box1Abb[220] + 3.*Box1Abb[0]*Box1Abb[326]*m_z1k + 3.*Box1Abb[133]*m_z1k_2;

  Box1Abb[3269]=27. - 14.*m_z12 - 15.*m_z1k;

  Box1Abb[3270]=3.*Box1Abb[665] + Box1Abb[3269]*m_z12;

  Box1Abb[3271]=2. + Box1Abb[3270]*m_z12;

  Box1Abb[3272]=Box1Abb[3271]*m_x_2 + Box1Abb[3268]*m_x*m_z12 + Box1Abb[3267]*m_x_3*m_z12 - pow(Box1Abb[4],3.)*m_z12_2;

  Box1Abb[3273]=Box1Abb[117]*Box1Abb[3266]*Box1Abb[384]*m_cR - Box1Abb[3272]*Box1Abb[7]*m_cL*m_x;

  Box1Abb[3274]=1. + 5.*Box1Abb[8]*m_z12;

  Box1Abb[3275]=pow(Box1Abb[0],2.)*Box1Abb[77]*Box1Abb[8] + 2.*Box1Abb[0]*Box1Abb[3274]*m_z1k + 6.*Box1Abb[60]*m_z12*m_z1k_2;

  Box1Abb[3276]=pow(Box1Abb[0],3.) + 3.*pow(Box1Abb[0],2.)*m_z1k + 4.*Box1Abb[0]*m_z1k_2 + m_z1k_3;

  Box1Abb[3277]=-1. + 7.*m_z12;

  Box1Abb[3278]=-pow(Box1Abb[0],3.) - pow(Box1Abb[0],2.)*Box1Abb[3277]*m_z1k - 3.*Box1Abb[0]*Box1Abb[1086]*m_z1k_2 + 2.*Box1Abb[208]*m_z1k_3;

  Box1Abb[3279]=16. - 9.*m_z12 - 22.*m_z1k;

  Box1Abb[3280]=-9. + Box1Abb[3279]*m_z12 + 18.*m_z1k;

  Box1Abb[3281]=2. + Box1Abb[3280]*m_z12;

  Box1Abb[3282]=Box1Abb[3275]*m_x_2 + Box1Abb[3281]*m_x_3 + Box1Abb[3278]*m_x*m_z12 + Box1Abb[3267]*m_x_4*m_z12 + Box1Abb[3276]*m_z12_2*m_z1k;

  Box1Abb[3283]=-Box1Abb[3266]*Box1Abb[384]*Box1Abb[7]*m_cR + Box1Abb[3282]*m_cL*m_x;

  Box1Abb[3284]=12. + 4.*Box1Abb[1446]*m_z12 + 7.*m_z12_2;

  Box1Abb[3285]=-10. + m_z1k;

  Box1Abb[3286]=9. + 2.*Box1Abb[2406]*m_z1k;

  Box1Abb[3287]=-2. + Box1Abb[3286]*m_z12 + Box1Abb[3285]*m_z12_2 + 3.*m_z12_3 + 6.*m_z1k + 4.*Box1Abb[170]*m_z1k_2;

  Box1Abb[3288]=-16. + 7.*m_z1k;

  Box1Abb[3289]=9. + 11.*m_z1k;

  Box1Abb[3290]=9. + 4.*Box1Abb[1446]*m_z1k;

  Box1Abb[3291]=10. - 15.*m_z12 + Box1Abb[3289]*m_z12_2 - 2.*m_z12_3 + 2.*Box1Abb[3290]*m_z1k + 2.*Box1Abb[3288]*m_z12*m_z1k;

  Box1Abb[3292]=-2. + Box1Abb[3291]*m_z12;

  Box1Abb[3293]=-42. + 11.*m_z12 + 2.*m_z1k;

  Box1Abb[3294]=60. + Box1Abb[3293]*m_z12 - 8.*m_z1k;

  Box1Abb[3295]=8.*Box1Abb[2312] + Box1Abb[3294]*m_z12;

  Box1Abb[3296]=4. + Box1Abb[3295]*m_z12;

  Box1Abb[3297]=Box1Abb[3296]*m_x_3 + Box1Abb[3292]*m_x_2*m_z12 - 2.*Box1Abb[3284]*m_x_4*m_z12 + 4.*m_x_5*m_z12_2 - Box1Abb[3287]*m_x*m_z12_2*m_z1k - Box1Abb[4]*m_z12_4*m_z1k_2;

  Box1Abb[3298]=pow(Box1Abb[0],3.) - 3.*Box1Abb[0]*Box1Abb[8]*m_z1k + 2.*Box1Abb[208]*m_z1k_2 - 2.*m_z1k_3;

  Box1Abb[3299]=4. + Box1Abb[1438]*m_z12;

  Box1Abb[3300]=-8. + m_z12;

  Box1Abb[3301]=4. + Box1Abb[3300]*m_z12;

  Box1Abb[3302]=pow(Box1Abb[0],3.)*m_z12_2 + pow(Box1Abb[0],2.)*Box1Abb[3299]*m_z1k + Box1Abb[0]*Box1Abb[3301]*m_z1k_2 + 6.*Box1Abb[1010]*m_z12*m_z1k_3 - 6.*m_z12*m_z1k_4;

  Box1Abb[3303]=-6. + 4.*m_z12 + m_z1k;

  Box1Abb[3304]=2. + Box1Abb[3303]*m_z12;

  Box1Abb[3305]=89. - 40.*m_z12;

  Box1Abb[3306]=-63. + Box1Abb[3305]*m_z12;

  Box1Abb[3307]=16. + Box1Abb[3306]*m_z12;

  Box1Abb[3308]=3. + Box1Abb[1620]*m_z12;

  Box1Abb[3309]=-3.*pow(Box1Abb[0],2.)*Box1Abb[133]*m_z12 + Box1Abb[3307]*m_z1k - 4.*Box1Abb[3308]*m_z1k_2 + 4.*m_z12*m_z1k_3;

  Box1Abb[3310]=Box1Abb[3309]*m_z12 - 2.*m_z1k;

  Box1Abb[3311]=-19. + m_z1k;

  Box1Abb[3312]=-58. + 51.*m_z1k;

  Box1Abb[3313]=42. + Box1Abb[3312]*m_z12 + 26.*m_z12_2 + 4.*Box1Abb[3311]*m_z1k;

  Box1Abb[3314]=-12. + Box1Abb[3313]*m_z12 + 24.*m_z1k;

  Box1Abb[3315]=2. + Box1Abb[3314]*m_z12;

  Box1Abb[3316]=Box1Abb[3310]*m_x_3 + Box1Abb[3315]*m_x_4 + Box1Abb[3302]*m_x_2*m_z12 - 6.*Box1Abb[3304]*m_x_5*m_z12 + 2.*m_x_6*m_z12_2 - Box1Abb[3298]*Box1Abb[4]*m_x*m_z12_2*m_z1k - pow(Box1Abb[4],3.)*m_z12_3*m_z1k_2;

  Box1Abb[3317]=2.*Box1Abb[3316]*m_cL - Box1Abb[3297]*Box1Abb[7]*m_cR;

  Box1Abb[3318]=6. - 7.*m_z12;

  Box1Abb[3319]=pow(Box1Abb[0],3.) + pow(Box1Abb[0],2.)*Box1Abb[3277]*m_z1k + 3.*Box1Abb[0]*Box1Abb[1086]*m_z1k_2 + 2.*Box1Abb[326]*m_z1k_3;

  Box1Abb[3320]=6. + 5.*Box1Abb[69]*m_z12;

  Box1Abb[3321]=-7. + 6.*m_z12 + m_z12_2 - 2.*m_z12_3 - 2.*Box1Abb[3320]*m_z1k + 6.*Box1Abb[3047]*m_z1k_2;

  Box1Abb[3322]=2.*Box1Abb[15] + Box1Abb[3321]*m_z12;

  Box1Abb[3323]=-16. + 9.*m_z12 + 22.*m_z1k;

  Box1Abb[3324]=9. + Box1Abb[3323]*m_z12 - 18.*m_z1k;

  Box1Abb[3325]=-2. + Box1Abb[3324]*m_z12;

  Box1Abb[3326]=Box1Abb[3322]*m_x_2 + Box1Abb[3325]*m_x_3 + Box1Abb[3319]*m_x*m_z12 + Box1Abb[3318]*m_x_4*m_z12 - Box1Abb[3276]*m_z12_2*m_z1k;

  Box1Abb[3327]=Box1Abb[3266]*Box1Abb[384]*Box1Abb[7]*m_cR + Box1Abb[3326]*m_cL*m_x;

  Box1Abb[3328]=-35. + 9.*m_z12;

  Box1Abb[3329]=9. + Box1Abb[3328]*m_z12;

  Box1Abb[3330]=5. + Box1Abb[3329]*m_z12;

  Box1Abb[3331]=-26. + 5.*m_z12;

  Box1Abb[3332]=60. + Box1Abb[3331]*m_z12;

  Box1Abb[3333]=-14. + Box1Abb[3332]*m_z12;

  Box1Abb[3334]=11. + Box1Abb[3333]*m_z12;

  Box1Abb[3335]=-10. + 17.*m_z12;

  Box1Abb[3336]=9. + Box1Abb[3335]*m_z12;

  Box1Abb[3337]=2. + 19.*m_z12;

  Box1Abb[3338]=-13. + 20.*m_z12 - 10.*m_z12_2 + 6.*m_z12_3 - 2.*Box1Abb[3330]*m_z1k + Box1Abb[3334]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[3336]*m_z1k_3 + 5.*Box1Abb[3337]*m_z12*m_z1k_4;

  Box1Abb[3339]=-7. + Box1Abb[1016]*m_z12;

  Box1Abb[3340]=-16. + 9.*m_z12;

  Box1Abb[3341]=21. + Box1Abb[3340]*m_z12;

  Box1Abb[3342]=-20. + Box1Abb[3341]*m_z12;

  Box1Abb[3343]=56. - 11.*m_z12;

  Box1Abb[3344]=-86. + Box1Abb[3343]*m_z12;

  Box1Abb[3345]=34. + Box1Abb[3344]*m_z12;

  Box1Abb[3346]=-29. + Box1Abb[3345]*m_z12;

  Box1Abb[3347]=51. - 26.*m_z12;

  Box1Abb[3348]=-2. + Box1Abb[3347]*m_z12;

  Box1Abb[3349]=7. + Box1Abb[3348]*m_z12;

  Box1Abb[3350]=4. + Box1Abb[3339]*m_z12 - 9.*m_z1k + 9.*Box1Abb[0]*Box1Abb[158]*m_z12*m_z1k + Box1Abb[0]*Box1Abb[3342]*m_z1k_2 + Box1Abb[3346]*m_z1k_3 + 2.*Box1Abb[3349]*m_z1k_4 - 9.*Box1Abb[1649]*m_z12*m_z1k_5;

  Box1Abb[3351]=2.*Box1Abb[68] + Box1Abb[358]*m_z12;

  Box1Abb[3352]=9. + 5.*m_z12;

  Box1Abb[3353]=-2. + 9.*m_z12;

  Box1Abb[3354]=34. + Box1Abb[1032]*m_z12 + 46.*m_z1k + 2.*Box1Abb[3352]*m_z12*m_z1k + 9.*Box1Abb[3353]*m_z1k_2;

  Box1Abb[3355]=-23. + Box1Abb[3354]*m_z12 - 26.*m_z1k;

  Box1Abb[3356]=8. + Box1Abb[1699]*m_z12 - 10.*m_z1k;

  Box1Abb[3357]=6. - Box1Abb[3356]*m_z12;

  Box1Abb[3358]=32. - 115.*m_z1k;

  Box1Abb[3359]=-20. + Box1Abb[3358]*m_z1k;

  Box1Abb[3360]=5. + Box1Abb[3359]*m_z1k;

  Box1Abb[3361]=-27. + Box1Abb[1680]*m_z1k;

  Box1Abb[3362]=-37. + 2.*Box1Abb[3361]*m_z1k;

  Box1Abb[3363]=26. + Box1Abb[3362]*m_z12 + Box1Abb[3360]*m_z12_2 - Box1Abb[1702]*m_z12_3 + 11.*Box1Abb[233]*m_z1k;

  Box1Abb[3364]=1. + Box1Abb[451]*m_z1k;

  Box1Abb[3365]=-27. + 11.*m_z1k;

  Box1Abb[3366]=20. + Box1Abb[3365]*m_z1k;

  Box1Abb[3367]=-4. + Box1Abb[3366]*m_z1k_2;

  Box1Abb[3368]=-2.*pow(Box1Abb[68],3.) + 2.*Box1Abb[3364]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3367]*m_z12_2 + 2.*Box1Abb[1695]*m_z12_3*m_z1k + Box1Abb[1692]*m_z12_4*m_z1k;

  Box1Abb[3369]=Box1Abb[3350]*m_x_2 + Box1Abb[3338]*m_x_3 + Box1Abb[3363]*m_x_4 + Box1Abb[3355]*m_x_5 + Box1Abb[3357]*m_x_6 + Box1Abb[505]*m_x_7*m_z12 + Box1Abb[3368]*Box1Abb[68]*m_x*m_z1k - Box1Abb[3351]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3370]=9. + 5.*Box1Abb[149]*m_z12;

  Box1Abb[3371]=5. + Box1Abb[3370]*m_z12;

  Box1Abb[3372]=-20. + m_z12;

  Box1Abb[3373]=24. + Box1Abb[3372]*m_z12;

  Box1Abb[3374]=-22. + Box1Abb[3373]*m_z12;

  Box1Abb[3375]=5. + Box1Abb[3374]*m_z12;

  Box1Abb[3376]=5. + 3.*m_z12;

  Box1Abb[3377]=-9. + Box1Abb[3376]*m_z12;

  Box1Abb[3378]=17. - 27.*m_z12 + 11.*m_z12_2 + 2.*Box1Abb[3371]*m_z1k - Box1Abb[3375]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[3377]*m_z1k_3 + 5.*Box1Abb[3353]*m_z12*m_z1k_4;

  Box1Abb[3379]=14. + Box1Abb[803]*m_z12;

  Box1Abb[3380]=-31. + 14.*m_z12 - 2.*Box1Abb[3379]*m_z1k + 3.*Box1Abb[653]*m_z1k_2;

  Box1Abb[3381]=23. + Box1Abb[3380]*m_z12 + 26.*m_z1k;

  Box1Abb[3382]=4. + Box1Abb[536]*m_z12 - 10.*m_z1k;

  Box1Abb[3383]=-6. + Box1Abb[3382]*m_z12;

  Box1Abb[3384]=35. + 44.*m_z1k;

  Box1Abb[3385]=17. + 25.*m_z1k;

  Box1Abb[3386]=20. + Box1Abb[15]*Box1Abb[3385]*m_z1k;

  Box1Abb[3387]=32. - 5.*m_z1k;

  Box1Abb[3388]=24. + Box1Abb[3387]*m_z1k;

  Box1Abb[3389]=47. + 2.*Box1Abb[3388]*m_z1k;

  Box1Abb[3390]=-31. + Box1Abb[3389]*m_z12 - Box1Abb[3386]*m_z12_2 - Box1Abb[3384]*m_z1k + 2.*Box1Abb[2406]*m_z12_3*m_z1k;

  Box1Abb[3391]=-5. + 6.*Box1Abb[68]*m_z1k;

  Box1Abb[3392]=2. + Box1Abb[1446]*m_z1k;

  Box1Abb[3393]=1. + Box1Abb[3392]*m_z1k;

  Box1Abb[3394]=3.*Box1Abb[3393]*m_z12 + Box1Abb[3391]*m_z12_2 + Box1Abb[119]*m_z12_3 - 2.*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[3395]=-6. - 11.*m_z1k + 10.*m_z1k_3;

  Box1Abb[3396]=-6. + 5.*m_z1k;

  Box1Abb[3397]=9. + Box1Abb[3396]*m_z1k;

  Box1Abb[3398]=-18. + 11.*m_z1k;

  Box1Abb[3399]=12. + Box1Abb[3398]*m_z1k;

  Box1Abb[3400]=-11. + Box1Abb[3399]*m_z1k;

  Box1Abb[3401]=4. + 2.*Box1Abb[3400]*m_z1k;

  Box1Abb[3402]=-43. + 17.*m_z1k;

  Box1Abb[3403]=7. + Box1Abb[3402]*m_z1k;

  Box1Abb[3404]=-9. + Box1Abb[3403]*m_z1k;

  Box1Abb[3405]=10. + Box1Abb[3404]*m_z1k;

  Box1Abb[3406]=-Box1Abb[3395]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3405]*Box1Abb[68]*m_z12_2 + Box1Abb[3401]*m_z12_3 + 2.*pow(Box1Abb[68],3.)*m_z1k + Box1Abb[3397]*m_z12_4*m_z1k;

  Box1Abb[3407]=-9. + 14.*m_z1k;

  Box1Abb[3408]=-3. + Box1Abb[3407]*m_z1k_2;

  Box1Abb[3409]=-7. + 4.*Box1Abb[1041]*m_z1k;

  Box1Abb[3410]=22. + Box1Abb[3409]*m_z1k;

  Box1Abb[3411]=-97. + 39.*m_z1k;

  Box1Abb[3412]=87. + Box1Abb[3411]*m_z1k;

  Box1Abb[3413]=-56. + Box1Abb[3412]*m_z1k;

  Box1Abb[3414]=2. + Box1Abb[3413]*m_z1k;

  Box1Abb[3415]=-26. + Box1Abb[3036]*m_z1k;

  Box1Abb[3416]=39. + 2.*Box1Abb[3415]*m_z1k;

  Box1Abb[3417]=-38. + Box1Abb[3416]*m_z1k;

  Box1Abb[3418]=5. + Box1Abb[3417]*m_z1k;

  Box1Abb[3419]=Box1Abb[3408]*Box1Abb[68] - Box1Abb[3418]*m_z12 + Box1Abb[15]*Box1Abb[3414]*m_z12_2 + Box1Abb[3410]*m_z12_3*m_z1k + Box1Abb[175]*m_z12_4*m_z1k_2;

  Box1Abb[3420]=-Box1Abb[3419]*m_x_2 + Box1Abb[3378]*m_x_3 + Box1Abb[3390]*m_x_4 + Box1Abb[3381]*m_x_5 + Box1Abb[3383]*m_x_6 - Box1Abb[79]*m_x_7*m_z12 + Box1Abb[3406]*m_x*m_z1k - Box1Abb[3394]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3421]=Box1Abb[3369]*m_cL + Box1Abb[3420]*m_cR;

  Box1Abb[3422]=-11. - 2.*m_z12 + 3.*Box1Abb[3300]*m_z1k;

  Box1Abb[3423]=-11. + m_z12;

  Box1Abb[3424]=20. + m_z12;

  Box1Abb[3425]=5. + 14.*m_z12 + 27.*m_z1k - 2.*Box1Abb[3423]*m_z12*m_z1k + 3.*Box1Abb[3424]*m_z1k_2;

  Box1Abb[3426]=-1. + 15.*m_z12;

  Box1Abb[3427]=-22. + m_z12;

  Box1Abb[3428]=13. - Box1Abb[3427]*m_z12;

  Box1Abb[3429]=16. + Box1Abb[3428]*m_z12;

  Box1Abb[3430]=-29. + 22.*m_z12 + 6.*m_z12_2;

  Box1Abb[3431]=-9. + 11.*m_z12 + m_z1k + 2.*Box1Abb[3426]*m_z12*m_z1k + Box1Abb[3429]*m_z1k_2 + 2.*Box1Abb[3430]*m_z1k_3 + 15.*Box1Abb[1793]*m_z1k_4;

  Box1Abb[3432]=13. + 3.*Box1Abb[12]*m_z12;

  Box1Abb[3433]=-17. + 2.*Box1Abb[3432]*m_z12;

  Box1Abb[3434]=-23. + m_z12_2;

  Box1Abb[3435]=-57. + 14.*Box1Abb[8]*m_z12;

  Box1Abb[3436]=8. + 13.*m_z12;

  Box1Abb[3437]=-2. + 2.*m_z12 + 6.*m_z1k + 22.*Box1Abb[0]*m_z12*m_z1k + Box1Abb[3433]*m_z1k_2 + Box1Abb[3434]*Box1Abb[79]*m_z1k_3 + Box1Abb[3435]*m_z1k_4 + 3.*Box1Abb[3436]*m_z1k_5;

  Box1Abb[3438]=Box1Abb[119]*m_z12 + 3.*Box1Abb[68]*m_z1k;

  Box1Abb[3439]=-4. - 8.*m_z1k + 17.*m_z1k_3;

  Box1Abb[3440]=-4. + 11.*m_z1k;

  Box1Abb[3441]=2. + Box1Abb[3440]*m_z1k;

  Box1Abb[3442]=-2. + Box1Abb[3441]*m_z1k;

  Box1Abb[3443]=Box1Abb[3439]*pow(Box1Abb[68],2.)*m_z12 + 2.*Box1Abb[3442]*Box1Abb[68]*m_z12_2 + Box1Abb[665]*pow(Box1Abb[68],3.)*m_z1k + Box1Abb[3397]*m_z12_3*m_z1k;

  Box1Abb[3444]=-1. + m_z1k - 40.*m_z1k_2;

  Box1Abb[3445]=56. + 25.*m_z1k;

  Box1Abb[3446]=41. + Box1Abb[3445]*m_z1k;

  Box1Abb[3447]=20. + Box1Abb[3446]*m_z1k;

  Box1Abb[3448]=9. - Box1Abb[3447]*m_z12 + 2.*Box1Abb[3444]*m_z1k + 2.*Box1Abb[2406]*m_z12_2*m_z1k;

  Box1Abb[3449]=-Box1Abb[3437]*m_x_2 + Box1Abb[3431]*m_x_3 + Box1Abb[3448]*m_x_4 + Box1Abb[3425]*m_x_5 + Box1Abb[3422]*m_x_6 - Box1Abb[158]*m_x_7 + Box1Abb[3443]*m_x*m_z1k - Box1Abb[3438]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3450]=-6. + Box1Abb[1793]*m_z12;

  Box1Abb[3451]=9. + 7.*Box1Abb[79]*m_z12;

  Box1Abb[3452]=1. - m_z12 + Box1Abb[3450]*m_z1k + Box1Abb[3451]*m_z1k_2 + Box1Abb[1488]*m_z1k_3;

  Box1Abb[3453]=9. + 31.*m_z1k;

  Box1Abb[3454]=11. - Box1Abb[3453]*m_z12 + 24.*m_z1k;

  Box1Abb[3455]=9. + 20.*m_z1k;

  Box1Abb[3456]=16. + 81.*m_z1k;

  Box1Abb[3457]=-6. + 3.*m_z12 + 2.*Box1Abb[532]*m_z12_2 - 3.*Box1Abb[3455]*m_z1k + Box1Abb[3456]*m_z12*m_z1k;

  Box1Abb[3458]=3. + m_z1k + 21.*m_z1k_2;

  Box1Abb[3459]=-1. + 40.*m_z1k;

  Box1Abb[3460]=3. + Box1Abb[3459]*m_z1k;

  Box1Abb[3461]=31. - 115.*m_z1k;

  Box1Abb[3462]=9. + Box1Abb[3461]*m_z1k;

  Box1Abb[3463]=5. + Box1Abb[3462]*m_z1k;

  Box1Abb[3464]=-1. + Box1Abb[3463]*m_z12 - 2.*Box1Abb[3458]*m_z12_2 + 2.*Box1Abb[3460]*m_z1k;

  Box1Abb[3465]=-15. - 22.*m_z1k + 68.*m_z1k_2;

  Box1Abb[3466]=6. + Box1Abb[3465]*m_z1k;

  Box1Abb[3467]=-29. + 30.*m_z1k;

  Box1Abb[3468]=11. + Box1Abb[3467]*m_z1k;

  Box1Abb[3469]=-1. + m_z1k + Box1Abb[3468]*m_z1k_2;

  Box1Abb[3470]=-9. + 95.*m_z1k;

  Box1Abb[3471]=22. + Box1Abb[3470]*Box1Abb[68]*m_z1k;

  Box1Abb[3472]=-8. + Box1Abb[3471]*m_z1k;

  Box1Abb[3473]=-2.*Box1Abb[3469] + Box1Abb[3472]*m_z12 + Box1Abb[3466]*m_z12_2 + 5.*m_z12_3*m_z1k_2;

  Box1Abb[3474]=8. + 3.*Box1Abb[1163]*m_z1k;

  Box1Abb[3475]=-2. + Box1Abb[3474]*m_z1k;

  Box1Abb[3476]=-11. + 45.*m_z1k;

  Box1Abb[3477]=-14. + Box1Abb[3476]*Box1Abb[68]*m_z1k;

  Box1Abb[3478]=4. + Box1Abb[3477]*m_z1k;

  Box1Abb[3479]=-25. + 26.*m_z1k;

  Box1Abb[3480]=5. + Box1Abb[3479]*m_z1k;

  Box1Abb[3481]=-3. + Box1Abb[3480]*m_z1k;

  Box1Abb[3482]=1. + Box1Abb[3481]*m_z1k;

  Box1Abb[3483]=Box1Abb[3475]*pow(Box1Abb[68],2.) - Box1Abb[3478]*Box1Abb[68]*m_z12 - 2.*Box1Abb[3482]*m_z12_2 + Box1Abb[1519]*m_z12_3*m_z1k_2;

  Box1Abb[3484]=Box1Abb[3483]*m_x_2 + Box1Abb[3473]*m_x_3 + Box1Abb[3464]*m_x_4 + Box1Abb[3457]*m_x_5 + Box1Abb[3454]*m_x_6 + Box1Abb[1280]*m_x_7 + Box1Abb[3452]*Box1Abb[4]*Box1Abb[68]*m_x*m_z1k - pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3485]=Box1Abb[3484]*m_cL + Box1Abb[3449]*m_cR;

  Box1Abb[3486]=-1. + Box1Abb[158]*m_x + m_z12;

  Box1Abb[3487]=23. - 24.*m_x;

  Box1Abb[3488]=-11. + Box1Abb[3487]*m_x;

  Box1Abb[3489]=6. + Box1Abb[3488]*m_x;

  Box1Abb[3490]=-3. + m_x + m_x_2;

  Box1Abb[3491]=2. + m_x + 3.*Box1Abb[3490]*m_x_2;

  Box1Abb[3492]=Box1Abb[3489]*m_x + Box1Abb[3491]*m_z12 - 2.*pow(Box1Abb[1],2.)*Box1Abb[183]*m_z12_2;

  Box1Abb[3493]=11. + m_z12;

  Box1Abb[3494]=14. + Box1Abb[542]*m_z12;

  Box1Abb[3495]=-38. + Box1Abb[1931]*m_z12;

  Box1Abb[3496]=13. + m_z12;

  Box1Abb[3497]=-29. + Box1Abb[3496]*m_z12;

  Box1Abb[3498]=6. + Box1Abb[3497]*m_z12;

  Box1Abb[3499]=Box1Abb[3498]*m_x_2 - Box1Abb[12]*Box1Abb[3493]*Box1Abb[79]*m_x_3 + Box1Abb[3495]*m_x_4 + 3.*Box1Abb[3424]*m_x_5 + Box1Abb[3494]*m_x*m_z12 - Box1Abb[0]*m_z12_2;

  Box1Abb[3500]=9. + 40.*m_x;

  Box1Abb[3501]=11. + Box1Abb[3500]*m_x;

  Box1Abb[3502]=5. + 2.*Box1Abb[3501]*m_x;

  Box1Abb[3503]=-24. + 25.*m_x;

  Box1Abb[3504]=-9. + Box1Abb[3503]*m_x;

  Box1Abb[3505]=6. + Box1Abb[3504]*m_x;

  Box1Abb[3506]=2. + Box1Abb[2171]*m_x;

  Box1Abb[3507]=-4. + Box1Abb[3506]*m_x;

  Box1Abb[3508]=8. + m_x;

  Box1Abb[3509]=-5. + Box1Abb[3508]*m_x;

  Box1Abb[3510]=Box1Abb[3502]*m_x + Box1Abb[3505]*m_x*m_z12 - 2.*Box1Abb[3507]*m_z12_2 + Box1Abb[3509]*m_z12_3;

  Box1Abb[3511]=9. - 28.*m_z12;

  Box1Abb[3512]=37. + Box1Abb[3511]*m_z12;

  Box1Abb[3513]=3. - Box1Abb[1621]*m_z12;

  Box1Abb[3514]=14. + m_z12 - 14.*m_z12_2 + 5.*m_z12_3;

  Box1Abb[3515]=Box1Abb[3514]*m_x + Box1Abb[3512]*m_x_2 + 15.*Box1Abb[1793]*m_x_3 + Box1Abb[3513]*m_z12;

  Box1Abb[3516]=13. + 24.*m_z12 - 22.*m_z12_2;

  Box1Abb[3517]=Box1Abb[3516]*m_x + 3.*Box1Abb[3436]*m_x_2 + 3.*Box1Abb[1436]*m_z12;

  Box1Abb[3518]=4. + 17.*m_z12;

  Box1Abb[3519]=-Box1Abb[3518]*m_x + 3.*Box1Abb[69]*m_z12;

  Box1Abb[3520]=pow(Box1Abb[1],4.)*Box1Abb[3486]*m_x_2 - Box1Abb[1]*Box1Abb[3492]*m_x*m_z1k - Box1Abb[3499]*m_z1k_2 + Box1Abb[3510]*m_z1k_3 - Box1Abb[3515]*m_z1k_4 + Box1Abb[3517]*m_z1k_5 + Box1Abb[3519]*m_z1k_6 + 3.*m_z12*m_z1k_7;

  Box1Abb[3521]=4. - 5.*m_z12;

  Box1Abb[3522]=20. - 27.*m_z12;

  Box1Abb[3523]=19. - 18.*m_z12 + 47.*m_z1k - 10.*Box1Abb[2661]*m_z12*m_z1k + 3.*Box1Abb[3522]*m_z1k_2;

  Box1Abb[3524]=7. - 4.*m_z12;

  Box1Abb[3525]=-13. + Box1Abb[77]*m_z12;

  Box1Abb[3526]=2. + 11.*m_z12;

  Box1Abb[3527]=19. + Box1Abb[3526]*m_z12;

  Box1Abb[3528]=-18. + Box1Abb[3527]*m_z12;

  Box1Abb[3529]=37. - 76.*m_z12 + 52.*m_z12_2;

  Box1Abb[3530]=-8. + 15.*m_z12;

  Box1Abb[3531]=-6. + 2.*Box1Abb[3524]*m_z12 + 11.*m_z1k + 2.*Box1Abb[3525]*m_z12*m_z1k + Box1Abb[3528]*m_z1k_2 + Box1Abb[3529]*m_z1k_3 + 3.*Box1Abb[3530]*m_z1k_4;

  Box1Abb[3532]=-7. + 2.*m_z12;

  Box1Abb[3533]=-17. + 5.*Box1Abb[653]*m_z12;

  Box1Abb[3534]=-16. + Box1Abb[3533]*m_z12;

  Box1Abb[3535]=9. - 22.*m_z12 + 34.*m_z12_2;

  Box1Abb[3536]=-1. + m_z12 + Box1Abb[3532]*m_z1k + Box1Abb[3534]*m_z1k_2 + 2.*Box1Abb[3535]*m_z1k_3 + 5.*Box1Abb[1509]*m_z1k_4;

  Box1Abb[3537]=5. + 8.*m_z1k;

  Box1Abb[3538]=-3.*Box1Abb[3537] + Box1Abb[1760]*m_z12;

  Box1Abb[3539]=-2. + Box1Abb[834]*m_z1k;

  Box1Abb[3540]=5. - Box1Abb[1776]*m_z1k;

  Box1Abb[3541]=-16. + 11.*m_z1k;

  Box1Abb[3542]=4. + Box1Abb[3541]*m_z1k_2;

  Box1Abb[3543]=-4. + m_z1k + 11.*m_z1k_2 - 9.*m_z1k_3;

  Box1Abb[3544]=1. + Box1Abb[3543]*m_z1k;

  Box1Abb[3545]=Box1Abb[3539]*pow(Box1Abb[68],3.) - Box1Abb[3542]*pow(Box1Abb[68],2.)*m_z12 + 2.*Box1Abb[3544]*m_z12_2 + Box1Abb[3540]*m_z12_3*m_z1k;

  Box1Abb[3546]=19. + 40.*m_z1k;

  Box1Abb[3547]=15. + Box1Abb[3546]*m_z1k;

  Box1Abb[3548]=34. + 115.*m_z1k;

  Box1Abb[3549]=11. + Box1Abb[3548]*m_z1k;

  Box1Abb[3550]=8. + Box1Abb[3549]*m_z1k;

  Box1Abb[3551]=-9. + Box1Abb[3550]*m_z12 - 2.*Box1Abb[3547]*m_z1k + 2.*Box1Abb[1767]*m_z12_2*m_z1k;

  Box1Abb[3552]=-Box1Abb[3536]*m_x_3 + Box1Abb[3551]*m_x_4 + Box1Abb[3523]*m_x_5 + Box1Abb[3538]*m_x_6 + Box1Abb[3521]*m_x_7 + Box1Abb[3545]*m_x*m_z1k + Box1Abb[3531]*m_x_2*m_z1k + Box1Abb[170]*pow(Box1Abb[4],2.)*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3553]=Box1Abb[3552]*m_cL + Box1Abb[3520]*m_cR;

  Box1Abb[3554]=17. + 6.*m_z12;

  Box1Abb[3555]=-25. + Box1Abb[3554]*m_z12;

  Box1Abb[3556]=6. + Box1Abb[3555]*m_z12;

  Box1Abb[3557]=3. + 4.*m_z12;

  Box1Abb[3558]=20. + Box1Abb[3557]*m_z12;

  Box1Abb[3559]=-19. + Box1Abb[3558]*m_z12;

  Box1Abb[3560]=172. + 5.*Box1Abb[3300]*m_z12;

  Box1Abb[3561]=-130. + Box1Abb[3560]*m_z12;

  Box1Abb[3562]=29. + Box1Abb[3561]*m_z12;

  Box1Abb[3563]=-15. + 17.*m_z12;

  Box1Abb[3564]=9. + Box1Abb[3563]*m_z12;

  Box1Abb[3565]=-1. + Box1Abb[3556]*m_z12 + 8.*m_z1k + 2.*Box1Abb[3559]*m_z12*m_z1k + Box1Abb[3562]*m_z1k_2 + 4.*Box1Abb[0]*Box1Abb[3564]*m_z1k_3 + 5.*Box1Abb[3337]*m_z12*m_z1k_4;

  Box1Abb[3566]=6. + 31.*m_z1k;

  Box1Abb[3567]=12. + Box1Abb[3566]*m_z12 - 10.*m_z1k;

  Box1Abb[3568]=6. - Box1Abb[3567]*m_z12;

  Box1Abb[3569]=19. + 44.*m_z1k;

  Box1Abb[3570]=-62. + 5.*m_z1k;

  Box1Abb[3571]=9. + Box1Abb[3570]*m_z1k;

  Box1Abb[3572]=29. + 2.*Box1Abb[1767]*m_z1k;

  Box1Abb[3573]=62. - 115.*m_z1k;

  Box1Abb[3574]=-63. + Box1Abb[3573]*m_z1k;

  Box1Abb[3575]=15. + Box1Abb[3574]*m_z1k;

  Box1Abb[3576]=3. + m_z12 + Box1Abb[3575]*m_z12_2 - Box1Abb[3572]*m_z12_3 - 2.*m_z12_4 + Box1Abb[3569]*m_z1k + 2.*Box1Abb[3571]*m_z12*m_z1k;

  Box1Abb[3577]=32. - 9.*m_z1k;

  Box1Abb[3578]=6. + 81.*m_z1k_2;

  Box1Abb[3579]=19. + Box1Abb[3578]*m_z12 + 2.*Box1Abb[2440]*m_z12_2 + 2.*Box1Abb[3577]*m_z1k;

  Box1Abb[3580]=-19. + Box1Abb[3579]*m_z12 - 26.*m_z1k;

  Box1Abb[3581]=2. + 7.*m_z1k_2;

  Box1Abb[3582]=-3. + 2.*Box1Abb[1151]*m_z1k;

  Box1Abb[3583]=-9. + 11.*m_z1k;

  Box1Abb[3584]=16. + Box1Abb[3583]*m_z1k;

  Box1Abb[3585]=-3. + Box1Abb[3584]*m_z1k;

  Box1Abb[3586]=-2.*pow(Box1Abb[68],2.) + Box1Abb[3582]*Box1Abb[68]*m_z12 + Box1Abb[3585]*m_z12_2 + Box1Abb[3581]*m_z12_3;

  Box1Abb[3587]=3. + 14.*m_z1k;

  Box1Abb[3588]=81. - 88.*m_z1k + 52.*m_z1k_2;

  Box1Abb[3589]=-5. + Box1Abb[3588]*m_z1k_2;

  Box1Abb[3590]=20. + 9.*m_z1k;

  Box1Abb[3591]=-19. + Box1Abb[3590]*m_z1k;

  Box1Abb[3592]=-5. + 2.*Box1Abb[3591]*m_z1k;

  Box1Abb[3593]=-12. + 11.*m_z1k;

  Box1Abb[3594]=3. + Box1Abb[3593]*m_z1k;

  Box1Abb[3595]=6. + Box1Abb[3594]*m_z1k;

  Box1Abb[3596]=-19. + 15.*m_z1k;

  Box1Abb[3597]=39. + Box1Abb[3596]*m_z1k;

  Box1Abb[3598]=-10. + Box1Abb[3597]*m_z1k;

  Box1Abb[3599]=-1. + 3.*Box1Abb[3598]*m_z1k;

  Box1Abb[3600]=-Box1Abb[3587]*pow(Box1Abb[68],3.) + Box1Abb[3592]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3599]*Box1Abb[68]*m_z12_2 + Box1Abb[3589]*m_z12_3 + Box1Abb[3595]*m_z12_4;

  Box1Abb[3601]=Box1Abb[3586]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_x - Box1Abb[3600]*m_x_2 + Box1Abb[3565]*m_x_3 + Box1Abb[3576]*m_x_4 + Box1Abb[3580]*m_x_5 + Box1Abb[3568]*m_x_6 + Box1Abb[505]*m_x_7*m_z12 - pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*Box1Abb[77]*m_z12*m_z1k;

  Box1Abb[3602]=-31. + 15.*m_z12;

  Box1Abb[3603]=19. + Box1Abb[3602]*m_z12;

  Box1Abb[3604]=-50. + m_z12;

  Box1Abb[3605]=76. + Box1Abb[3604]*m_z12;

  Box1Abb[3606]=-62. + Box1Abb[3605]*m_z12;

  Box1Abb[3607]=23. + Box1Abb[3606]*m_z12;

  Box1Abb[3608]=-9. + Box1Abb[2048]*m_z12;

  Box1Abb[3609]=2. - 9.*m_z12;

  Box1Abb[3610]=-1. + 2.*Box1Abb[0]*Box1Abb[129]*m_z12 - 2.*m_z1k - 2.*Box1Abb[3603]*m_z12*m_z1k + Box1Abb[3607]*m_z1k_2 - 4.*Box1Abb[0]*Box1Abb[3608]*m_z1k_3 + 5.*Box1Abb[3609]*m_z12*m_z1k_4;

  Box1Abb[3611]=2. + 3.*Box1Abb[358]*m_z12 - 2.*m_z1k;

  Box1Abb[3612]=-8. - 3.*Box1Abb[170]*m_z12 + 10.*m_z1k;

  Box1Abb[3613]=6. + Box1Abb[3612]*m_z12;

  Box1Abb[3614]=-19. + 2.*m_z12;

  Box1Abb[3615]=-19. + m_z12;

  Box1Abb[3616]=30. + Box1Abb[3614]*m_z12 + 46.*m_z1k + 2.*Box1Abb[3615]*m_z12*m_z1k - 3.*Box1Abb[653]*m_z1k_2;

  Box1Abb[3617]=-19. + Box1Abb[3616]*m_z12 - 26.*m_z1k;

  Box1Abb[3618]=-13. + 11.*m_z1k;

  Box1Abb[3619]=-3. + Box1Abb[3618]*m_z1k;

  Box1Abb[3620]=4. + 2.*Box1Abb[3619]*m_z1k;

  Box1Abb[3621]=-32. + 17.*m_z1k;

  Box1Abb[3622]=-8. + Box1Abb[3621]*m_z1k;

  Box1Abb[3623]=12. + Box1Abb[3622]*m_z1k;

  Box1Abb[3624]=2.*pow(Box1Abb[68],3.) - 2.*Box1Abb[1159]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3623]*Box1Abb[68]*m_z12_2 + Box1Abb[3620]*m_z12_3 + Box1Abb[532]*m_z12_4*m_z1k;

  Box1Abb[3625]=21. + 44.*m_z1k;

  Box1Abb[3626]=19. - 2.*m_z1k;

  Box1Abb[3627]=-6. + Box1Abb[3626]*m_z1k;

  Box1Abb[3628]=-15. + Box1Abb[1680]*m_z1k;

  Box1Abb[3629]=-25. + 2.*Box1Abb[3628]*m_z1k;

  Box1Abb[3630]=72. + 25.*m_z1k;

  Box1Abb[3631]=4. + Box1Abb[3630]*m_z1k;

  Box1Abb[3632]=21. + Box1Abb[3631]*m_z1k;

  Box1Abb[3633]=14. + Box1Abb[3629]*m_z12 + Box1Abb[3632]*m_z12_2 + Box1Abb[3627]*m_z12_3 + Box1Abb[3625]*m_z1k;

  Box1Abb[3634]=-5. + 14.*m_z1k;

  Box1Abb[3635]=-37. + 4.*Box1Abb[3089]*m_z1k;

  Box1Abb[3636]=5. + Box1Abb[3635]*m_z1k;

  Box1Abb[3637]=-2. + Box1Abb[3636]*m_z1k;

  Box1Abb[3638]=11. + 9.*m_z1k;

  Box1Abb[3639]=-26. + Box1Abb[3638]*m_z1k;

  Box1Abb[3640]=17. + 2.*Box1Abb[3639]*m_z1k;

  Box1Abb[3641]=-3. + Box1Abb[3640]*m_z1k;

  Box1Abb[3642]=-19. + 39.*m_z1k;

  Box1Abb[3643]=-57. + Box1Abb[3642]*m_z1k;

  Box1Abb[3644]=12. + Box1Abb[3643]*m_z1k;

  Box1Abb[3645]=-5. + Box1Abb[3644]*m_z1k;

  Box1Abb[3646]=-Box1Abb[3641]*Box1Abb[68]*m_z12 + Box1Abb[3645]*Box1Abb[68]*m_z12_2 + Box1Abb[3637]*m_z12_3 + Box1Abb[3634]*pow(Box1Abb[68],2.)*m_z1k + Box1Abb[2477]*m_z12_4*m_z1k_2;

  Box1Abb[3647]=Box1Abb[3646]*m_x_2 + Box1Abb[3610]*m_x_3 + Box1Abb[3633]*m_x_4 + Box1Abb[3617]*m_x_5 + Box1Abb[3613]*m_x_6 + Box1Abb[79]*m_x_7*m_z12 - Box1Abb[3624]*Box1Abb[68]*m_x*m_z1k + Box1Abb[3611]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3648]=Box1Abb[3601]*m_cL - Box1Abb[3647]*m_cR;

  Box1Abb[3649]=-24. + 31.*m_z12;

  Box1Abb[3650]=-7. + Box1Abb[3649]*m_z1k;

  Box1Abb[3651]=23. + 2.*m_z12;

  Box1Abb[3652]=-22. + Box1Abb[3651]*m_z12;

  Box1Abb[3653]=1. + Box1Abb[493]*m_z12;

  Box1Abb[3654]=-16. + 23.*m_z12;

  Box1Abb[3655]=2. + Box1Abb[3652]*m_z12 + 10.*m_z1k + 8.*m_z12_2*m_z1k + 42.*Box1Abb[3653]*m_z1k_2 + 5.*Box1Abb[3654]*m_z1k_3;

  Box1Abb[3656]=11. + 6.*m_z12;

  Box1Abb[3657]=-26. + Box1Abb[3656]*m_z12;

  Box1Abb[3658]=10. + m_z12 + 8.*m_z12_2;

  Box1Abb[3659]=-54. + 5.*m_z12;

  Box1Abb[3660]=101. + Box1Abb[3659]*m_z12;

  Box1Abb[3661]=-40. + Box1Abb[3660]*m_z12;

  Box1Abb[3662]=49. + 34.*Box1Abb[493]*m_z12;

  Box1Abb[3663]=10. + Box1Abb[3657]*m_z12 - 8.*m_z1k + Box1Abb[3658]*m_z12*m_z1k + Box1Abb[3661]*m_z1k_2 + 2.*Box1Abb[3662]*m_z1k_3 + 5.*Box1Abb[1509]*m_z1k_4;

  Box1Abb[3664]=7. + 60.*m_z1k;

  Box1Abb[3665]=7. + 30.*m_z1k - 81.*m_z1k_2;

  Box1Abb[3666]=-5. + Box1Abb[3665]*m_z12 - 10.*Box1Abb[15]*m_z12_2 + Box1Abb[3664]*m_z1k;

  Box1Abb[3667]=-6. + 11.*m_z1k;

  Box1Abb[3668]=3. + Box1Abb[3667]*m_z1k;

  Box1Abb[3669]=-pow(Box1Abb[68],2.)*Box1Abb[834] + Box1Abb[3668]*Box1Abb[68]*m_z12 + Box1Abb[3581]*m_z12_2;

  Box1Abb[3670]=-5. + 24.*m_z1k;

  Box1Abb[3671]=3. + Box1Abb[3670]*m_z1k;

  Box1Abb[3672]=9. - 50.*m_z1k + 52.*m_z1k_2;

  Box1Abb[3673]=7. + Box1Abb[3672]*m_z1k;

  Box1Abb[3674]=-22. + 15.*m_z1k;

  Box1Abb[3675]=4. + Box1Abb[3674]*m_z1k;

  Box1Abb[3676]=-2. + 3.*Box1Abb[3675]*m_z1k;

  Box1Abb[3677]=-Box1Abb[3671]*pow(Box1Abb[68],3.) + Box1Abb[3676]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3673]*Box1Abb[68]*m_z12_2 + Box1Abb[3595]*m_z12_3;

  Box1Abb[3678]=-Box1Abb[3669]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_x + Box1Abb[3677]*m_x_2 - Box1Abb[3663]*m_x_3 + Box1Abb[3655]*m_x_4 + Box1Abb[3666]*m_x_5 + Box1Abb[3650]*m_x_6 + Box1Abb[3521]*m_x_7 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12*m_z1k;

  Box1Abb[3679]=7. - 3.*m_z1k;

  Box1Abb[3680]=7. + Box1Abb[3679]*m_z12 + 24.*m_z1k;

  Box1Abb[3681]=6. + Box1Abb[1931]*m_z12 - 7.*m_z1k + 2.*Box1Abb[3427]*m_z12*m_z1k - 3.*Box1Abb[3424]*m_z1k_2;

  Box1Abb[3682]=26. - 27.*m_z12;

  Box1Abb[3683]=-2. + Box1Abb[3682]*m_z12;

  Box1Abb[3684]=29. + Box1Abb[3604]*m_z12;

  Box1Abb[3685]=-34. + Box1Abb[3684]*m_z12;

  Box1Abb[3686]=-49. + 32.*m_z12 + 6.*m_z12_2;

  Box1Abb[3687]=6.*pow(Box1Abb[0],2.) + Box1Abb[3683]*m_z1k + Box1Abb[3685]*m_z1k_2 - 2.*Box1Abb[3686]*m_z1k_3 - 15.*Box1Abb[1793]*m_z1k_4;

  Box1Abb[3688]=3. + 17.*m_z1k;

  Box1Abb[3689]=-5. + Box1Abb[3688]*m_z1k;

  Box1Abb[3690]=Box1Abb[665]*pow(Box1Abb[68],2.) + Box1Abb[3689]*Box1Abb[68]*m_z12 + Box1Abb[532]*m_z12_2*m_z1k;

  Box1Abb[3691]=-21. + 40.*m_z1k;

  Box1Abb[3692]=-7. + Box1Abb[3691]*m_z1k;

  Box1Abb[3693]=91. + 25.*m_z1k;

  Box1Abb[3694]=37. + Box1Abb[3693]*m_z1k;

  Box1Abb[3695]=21. + Box1Abb[3694]*m_z1k;

  Box1Abb[3696]=-13. + Box1Abb[3695]*m_z12 - 2.*Box1Abb[2371]*m_z12_2 + 2.*Box1Abb[3692]*m_z1k;

  Box1Abb[3697]=-1. + 7.*Box1Abb[15]*m_z1k;

  Box1Abb[3698]=-1. + Box1Abb[317]*Box1Abb[3697]*m_z1k;

  Box1Abb[3699]=-29. + 24.*m_z1k;

  Box1Abb[3700]=4. + Box1Abb[3699]*m_z1k;

  Box1Abb[3701]=-2. + Box1Abb[3700]*m_z1k;

  Box1Abb[3702]=20. + 39.*m_z1k;

  Box1Abb[3703]=-35. + Box1Abb[3702]*m_z1k;

  Box1Abb[3704]=6. + Box1Abb[3703]*m_z1k;

  Box1Abb[3705]=-4. + Box1Abb[3704]*m_z1k;

  Box1Abb[3706]=Box1Abb[3701]*pow(Box1Abb[68],2.) + Box1Abb[3705]*Box1Abb[68]*m_z12 + 2.*Box1Abb[3698]*m_z12_2 + Box1Abb[2477]*m_z12_3*m_z1k_2;

  Box1Abb[3707]=Box1Abb[3706]*m_x_2 + Box1Abb[3687]*m_x_3 + Box1Abb[3696]*m_x_4 + Box1Abb[3681]*m_x_5 + Box1Abb[3680]*m_x_6 + Box1Abb[158]*m_x_7 - Box1Abb[3690]*Box1Abb[4]*Box1Abb[68]*m_x*m_z1k + 3.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3708]=Box1Abb[3678]*m_cL + Box1Abb[3707]*m_cR;

  Box1Abb[3709]=-11. + 9.*m_z12;

  Box1Abb[3710]=-11. + 2.*m_z12;

  Box1Abb[3711]=16. + Box1Abb[3710]*m_z12;

  Box1Abb[3712]=-2. + Box1Abb[3711]*m_z12;

  Box1Abb[3713]=13. + Box1Abb[3331]*m_z12;

  Box1Abb[3714]=-10. + Box1Abb[3713]*m_z12;

  Box1Abb[3715]=29. - 72.*m_z12 + 34.*m_z12_2;

  Box1Abb[3716]=6. + 2.*Box1Abb[3709]*m_z12 + Box1Abb[3712]*m_z1k + Box1Abb[3714]*m_z1k_2 + 2.*Box1Abb[3715]*m_z1k_3 + 5.*Box1Abb[1509]*m_z1k_4;

  Box1Abb[3717]=7. + 31.*m_z1k;

  Box1Abb[3718]=-11. + Box1Abb[3717]*m_z12 - 24.*m_z1k;

  Box1Abb[3719]=4. + 81.*m_z1k;

  Box1Abb[3720]=3. - Box1Abb[3719]*m_z1k;

  Box1Abb[3721]=4. + Box1Abb[3720]*m_z12 - 2.*Box1Abb[451]*m_z12_2 + 3.*Box1Abb[3455]*m_z1k;

  Box1Abb[3722]=-6. + m_z1k + 7.*m_z1k_2;

  Box1Abb[3723]=-15. + 11.*m_z1k;

  Box1Abb[3724]=-1. + Box1Abb[3723]*m_z1k;

  Box1Abb[3725]=2. + Box1Abb[3724]*m_z1k;

  Box1Abb[3726]=-Box1Abb[3539]*pow(Box1Abb[68],2.) + Box1Abb[3725]*Box1Abb[68]*m_z12 + Box1Abb[3722]*m_z12_2*m_z1k;

  Box1Abb[3727]=1. + m_z1k - 40.*m_z1k_2;

  Box1Abb[3728]=1. + 7.*m_z1k;

  Box1Abb[3729]=7. + 3.*Box1Abb[3728]*m_z1k;

  Box1Abb[3730]=-61. + 115.*m_z1k;

  Box1Abb[3731]=-23. + Box1Abb[3730]*m_z1k;

  Box1Abb[3732]=-15. + Box1Abb[3731]*m_z1k;

  Box1Abb[3733]=5. + Box1Abb[3732]*m_z12 + 2.*Box1Abb[3729]*m_z12_2 + 2.*Box1Abb[3727]*m_z1k;

  Box1Abb[3734]=-2. + Box1Abb[1163]*m_z1k_2;

  Box1Abb[3735]=6. + Box1Abb[2119]*m_z1k;

  Box1Abb[3736]=-33. + 26.*m_z1k;

  Box1Abb[3737]=5. + Box1Abb[3736]*m_z1k;

  Box1Abb[3738]=-7. + Box1Abb[3737]*m_z1k;

  Box1Abb[3739]=5. + Box1Abb[3738]*m_z1k;

  Box1Abb[3740]=-86. + 45.*m_z1k;

  Box1Abb[3741]=17. + Box1Abb[3740]*m_z1k;

  Box1Abb[3742]=-10. + Box1Abb[3741]*m_z1k;

  Box1Abb[3743]=16. + Box1Abb[3742]*m_z1k;

  Box1Abb[3744]=-3.*Box1Abb[3734]*pow(Box1Abb[68],2.) + Box1Abb[3743]*Box1Abb[68]*m_z12 + 2.*Box1Abb[3739]*m_z12_2 + Box1Abb[3735]*m_z12_3*m_z1k;

  Box1Abb[3745]=-Box1Abb[3726]*Box1Abb[4]*Box1Abb[68]*m_x + Box1Abb[3744]*m_x_2 - Box1Abb[3716]*m_x_3 + Box1Abb[3733]*m_x_4 + Box1Abb[3721]*m_x_5 + Box1Abb[3718]*m_x_6 + Box1Abb[3521]*m_x_7 + Box1Abb[170]*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[3746]=11. + 2.*m_z12 - 3.*Box1Abb[3300]*m_z1k;

  Box1Abb[3747]=41. + 10.*m_z12;

  Box1Abb[3748]=2. + Box1Abb[3747]*m_z12;

  Box1Abb[3749]=-28. + m_z12;

  Box1Abb[3750]=1. + Box1Abb[3749]*m_z12;

  Box1Abb[3751]=16. + 5.*m_z12;

  Box1Abb[3752]=-9. + 20.*m_z12 + Box1Abb[3748]*m_z1k - 2.*Box1Abb[3750]*m_z1k_2 + 5.*Box1Abb[3751]*m_z1k_3;

  Box1Abb[3753]=15. + Box1Abb[2002]*m_z1k;

  Box1Abb[3754]=8. + Box1Abb[3467]*m_z1k;

  Box1Abb[3755]=1. + 2.*Box1Abb[3754]*m_z1k;

  Box1Abb[3756]=44. + 45.*m_z1k;

  Box1Abb[3757]=13. + Box1Abb[3756]*m_z1k;

  Box1Abb[3758]=-2. + Box1Abb[3757]*m_z1k;

  Box1Abb[3759]=11. + Box1Abb[3758]*m_z1k;

  Box1Abb[3760]=9. - Box1Abb[3759]*m_z12 - Box1Abb[3755]*m_z1k - 2.*Box1Abb[3753]*m_z12_2*m_z1k + m_z12_3*m_z1k_2;

  Box1Abb[3761]=Box1Abb[3437]*m_x_2 + Box1Abb[3760]*m_x_3 + Box1Abb[3752]*m_x_4 - Box1Abb[3425]*m_x_5 + Box1Abb[3746]*m_x_6 + Box1Abb[158]*m_x_7 - Box1Abb[3443]*m_x*m_z1k + Box1Abb[3438]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3762]=Box1Abb[3745]*m_cL + Box1Abb[3761]*m_cR;

  Box1Abb[3763]=Box1Abb[418]*m_cL*m_x - Box1Abb[416]*Box1Abb[7]*m_cR*m_z1k;

  Box1Abb[3764]=2.*pow(Box1Abb[68],2.) + 9.*Box1Abb[68]*m_z12 + 5.*m_z12_2;

  Box1Abb[3765]=12. + Box1Abb[1112]*m_z12;

  Box1Abb[3766]=16. - 9.*m_z12_2 + 6.*m_z1k - pow(Box1Abb[251],2.)*m_z12*m_z1k + 2.*Box1Abb[3765]*m_z1k_2;

  Box1Abb[3767]=15. + 2.*m_z12;

  Box1Abb[3768]=-16. + Box1Abb[3767]*m_z12;

  Box1Abb[3769]=18. - 5.*m_z12;

  Box1Abb[3770]=-11. + Box1Abb[3769]*m_z12;

  Box1Abb[3771]=20. + Box1Abb[3770]*m_z12;

  Box1Abb[3772]=6. + Box1Abb[141]*m_z12;

  Box1Abb[3773]=-2. + 3.*Box1Abb[0]*m_z12 - 6.*m_z1k + Box1Abb[3768]*m_z12*m_z1k + Box1Abb[3771]*m_z1k_2 - 2.*Box1Abb[3772]*m_z1k_3 + 10.*m_z12*m_z1k_4;

  Box1Abb[3774]=-2. - 3.*m_z12 + 8.*m_z1k;

  Box1Abb[3775]=6. + Box1Abb[3774]*m_z12;

  Box1Abb[3776]=6. - 5.*m_z1k;

  Box1Abb[3777]=7. - 3.*Box1Abb[881]*m_z12 + 2.*Box1Abb[3776]*m_z1k;

  Box1Abb[3778]=-20.*Box1Abb[15] + Box1Abb[3777]*m_z12;

  Box1Abb[3779]=-1. + Box1Abb[2805]*m_z1k;

  Box1Abb[3780]=-13. - 2.*m_z1k + 8.*m_z1k_2;

  Box1Abb[3781]=2.*pow(Box1Abb[68],3.) - Box1Abb[3780]*pow(Box1Abb[68],2.)*m_z12 - 3.*Box1Abb[1039]*Box1Abb[15]*Box1Abb[68]*m_z12_2 + Box1Abb[3779]*m_z12_3;

  Box1Abb[3782]=Box1Abb[3773]*m_x_2 + Box1Abb[3766]*m_x_3 + Box1Abb[3778]*m_x_4 + Box1Abb[3775]*m_x_5 - 2.*m_x_6*m_z12 + Box1Abb[3781]*m_x*m_z1k + Box1Abb[3764]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3783]=-2. + 4.*m_z1k;

  Box1Abb[3784]=m_z12 + Box1Abb[3783]*m_z12_2 - 3.*m_z1k + Box1Abb[2073]*m_z12*m_z1k;

  Box1Abb[3785]=1. + 4.*Box1Abb[437]*m_z1k;

  Box1Abb[3786]=-1. + Box1Abb[1945]*m_z1k;

  Box1Abb[3787]=2.*pow(Box1Abb[68],3.) - 2.*Box1Abb[3786]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3785]*Box1Abb[68]*m_z12_2 + Box1Abb[2080]*m_z12_3;

  Box1Abb[3788]=24. + 19.*m_z1k;

  Box1Abb[3789]=-6. + Box1Abb[3788]*m_z1k;

  Box1Abb[3790]=29. - 5.*m_z1k;

  Box1Abb[3791]=8. + Box1Abb[3790]*m_z1k;

  Box1Abb[3792]=3. + Box1Abb[3791]*m_z1k;

  Box1Abb[3793]=-1. + 2.*Box1Abb[3792]*m_z12 + Box1Abb[3789]*m_z12_2 - 4.*Box1Abb[2440]*m_z1k + m_z12_3*m_z1k;

  Box1Abb[3794]=23. + 36.*m_z1k;

  Box1Abb[3795]=-4. + Box1Abb[3794]*m_z1k;

  Box1Abb[3796]=5. + 6.*Box1Abb[431]*m_z1k;

  Box1Abb[3797]=1. + Box1Abb[3796]*m_z1k;

  Box1Abb[3798]=-3. + Box1Abb[3459]*m_z1k;

  Box1Abb[3799]=6. + 2.*Box1Abb[3798]*m_z1k;

  Box1Abb[3800]=-2.*Box1Abb[3797] + Box1Abb[3799]*m_z12 + Box1Abb[3795]*m_z12_2 + Box1Abb[2088]*m_z12_3*m_z1k;

  Box1Abb[3801]=1. + 4.*Box1Abb[265]*m_z1k;

  Box1Abb[3802]=-1. + Box1Abb[1895]*m_z1k;

  Box1Abb[3803]=-1. + m_z1k + Box1Abb[3802]*m_z1k_2;

  Box1Abb[3804]=-40. + 19.*m_z1k;

  Box1Abb[3805]=-6. + Box1Abb[3804]*m_z1k;

  Box1Abb[3806]=-6. + Box1Abb[3805]*m_z1k;

  Box1Abb[3807]=1. + Box1Abb[3806]*m_z1k;

  Box1Abb[3808]=Box1Abb[3801]*pow(Box1Abb[68],2.) - 2.*Box1Abb[3803]*Box1Abb[68]*m_z12 + Box1Abb[3807]*m_z12_2 - Box1Abb[1692]*m_z12_3*m_z1k;

  Box1Abb[3809]=Box1Abb[3808]*m_x_2 + Box1Abb[3800]*m_x_3 - Box1Abb[3793]*m_x_4 + 2.*Box1Abb[3784]*m_x_5 + Box1Abb[63]*m_x_6*m_z12 - Box1Abb[3787]*Box1Abb[68]*m_x*m_z1k + Box1Abb[3234]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3810]=Box1Abb[3809]*m_cL + Box1Abb[3782]*m_cR*m_z1k;

  Box1Abb[3811]=-3. + 5.*m_z12 + 3.*m_z1k;

  Box1Abb[3812]=-4. + 9.*m_z12;

  Box1Abb[3813]=9. + m_z12;

  Box1Abb[3814]=21. + Box1Abb[3813]*m_z12;

  Box1Abb[3815]=6. + Box1Abb[1434]*m_z12;

  Box1Abb[3816]=-3. + Box1Abb[3812]*m_z12 + m_z1k + Box1Abb[3814]*m_z12*m_z1k - 4.*Box1Abb[3815]*m_z1k_2 + 40.*m_z1k_3;

  Box1Abb[3817]=-4. + Box1Abb[3813]*m_z12;

  Box1Abb[3818]=19. - 5.*m_z12;

  Box1Abb[3819]=-11. + Box1Abb[3818]*m_z12;

  Box1Abb[3820]=27. + Box1Abb[3819]*m_z12;

  Box1Abb[3821]=21. + Box1Abb[1475]*m_z12;

  Box1Abb[3822]=1. + Box1Abb[129]*m_z12 - 6.*m_z1k + 2.*Box1Abb[3817]*m_z12*m_z1k + Box1Abb[3820]*m_z1k_2 - 2.*Box1Abb[3821]*m_z1k_3 + 5.*Box1Abb[1793]*m_z1k_4;

  Box1Abb[3823]=9. - 2.*m_z1k;

  Box1Abb[3824]=1. + m_z12 + Box1Abb[3823]*m_z12_2 + 10.*Box1Abb[174]*m_z1k + 3.*Box1Abb[711]*m_z12*m_z1k;

  Box1Abb[3825]=2. - 3.*m_z12 + 12.*m_z1k;

  Box1Abb[3826]=-9. + Box1Abb[3825]*m_z12 - 20.*m_z1k;

  Box1Abb[3827]=1. - Box1Abb[2805]*m_z1k;

  Box1Abb[3828]=-13. + 6.*Box1Abb[431]*m_z1k;

  Box1Abb[3829]=-7. + Box1Abb[732]*m_z1k;

  Box1Abb[3830]=Box1Abb[665]*pow(Box1Abb[68],3.) + Box1Abb[3828]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3829]*Box1Abb[68]*m_z12_2 + Box1Abb[3827]*m_z12_3;

  Box1Abb[3831]=Box1Abb[3822]*m_x_2 - Box1Abb[3816]*m_x_3 + Box1Abb[3824]*m_x_4 + Box1Abb[3826]*m_x_5 + Box1Abb[527]*m_x_6 - Box1Abb[3830]*m_x*m_z1k + Box1Abb[3811]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[3832]=m_z12_2 - 4.*m_z1k + 7.*m_z12*m_z1k;

  Box1Abb[3833]=9. + 16.*m_z1k;

  Box1Abb[3834]=7. + 29.*m_z1k;

  Box1Abb[3835]=2. - Box1Abb[3834]*m_z1k;

  Box1Abb[3836]=Box1Abb[3835]*m_z12 + Box1Abb[1828]*m_z12_2 + Box1Abb[3833]*m_z1k;

  Box1Abb[3837]=20. + 3.*m_z1k;

  Box1Abb[3838]=-6. + Box1Abb[3837]*m_z1k;

  Box1Abb[3839]=8. - 23.*m_z1k;

  Box1Abb[3840]=-9. + Box1Abb[3839]*m_z1k;

  Box1Abb[3841]=3. + Box1Abb[3840]*m_z1k;

  Box1Abb[3842]=-1. + 2.*Box1Abb[3841]*m_z12 + Box1Abb[3838]*m_z12_2 + m_z12_3*m_z1k + m_z1k_2 + 24.*m_z1k_3;

  Box1Abb[3843]=-4. + 17.*m_z1k - 21.*m_z1k_3;

  Box1Abb[3844]=7. + 16.*m_z1k;

  Box1Abb[3845]=-2. + Box1Abb[3844]*m_z1k;

  Box1Abb[3846]=-16. + m_z1k + 29.*m_z1k_2 - 17.*m_z1k_3;

  Box1Abb[3847]=3. + Box1Abb[3846]*m_z1k;

  Box1Abb[3848]=Box1Abb[3845]*pow(Box1Abb[68],2.) + 2.*Box1Abb[3847]*m_z12 + Box1Abb[3843]*m_z12_2 + Box1Abb[1692]*m_z12_3*m_z1k;

  Box1Abb[3849]=Box1Abb[3452]*Box1Abb[4]*Box1Abb[68]*m_x + Box1Abb[3848]*m_x_2 - Box1Abb[3842]*m_x_3 + Box1Abb[3836]*m_x_4 + Box1Abb[3832]*m_x_5 - pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[3850]=Box1Abb[3849]*Box1Abb[7]*m_cL + Box1Abb[3831]*m_cR*m_z1k;

  Box1Abb[3851]=18. - 11.*m_z12;

  Box1Abb[3852]=-3. + Box1Abb[3277]*m_z12;

  Box1Abb[3853]=-9. + Box1Abb[3852]*m_z12;

  Box1Abb[3854]=-34. + 29.*m_z12;

  Box1Abb[3855]=13. + Box1Abb[3854]*m_z12;

  Box1Abb[3856]=-7. + Box1Abb[3851]*m_z12 + 10.*m_z1k - 8.*Box1Abb[3235]*m_z12*m_z1k + Box1Abb[3853]*m_z1k_2 + 2.*Box1Abb[3855]*m_z1k_3 + 5.*Box1Abb[3812]*m_z1k_4;

  Box1Abb[3857]=13. - 4.*Box1Abb[3065]*m_z12 + 5.*m_z12_2 + 20.*m_z1k;

  Box1Abb[3858]=-17. + 4.*m_z1k;

  Box1Abb[3859]=13. + 20.*m_z1k;

  Box1Abb[3860]=26. + 75.*m_z1k;

  Box1Abb[3861]=21. + Box1Abb[3860]*m_z1k;

  Box1Abb[3862]=-11. + Box1Abb[3861]*m_z12 + Box1Abb[3858]*m_z12_2 - 2.*Box1Abb[3859]*m_z1k;

  Box1Abb[3863]=5. + 44.*m_z1k;

  Box1Abb[3864]=21. - Box1Abb[3863]*m_z1k;

  Box1Abb[3865]=3. - 10.*m_z1k;

  Box1Abb[3866]=11. + 8.*Box1Abb[3865]*m_z1k;

  Box1Abb[3867]=-26. + Box1Abb[3866]*m_z1k;

  Box1Abb[3868]=7. + Box1Abb[3867]*m_z12 + Box1Abb[3864]*m_z12_2 + 7.*m_z1k + 3.*m_z12_3*m_z1k + 40.*m_z1k_3;

  Box1Abb[3869]=-11. + m_z1k + 33.*m_z1k_2 - 25.*m_z1k_3;

  Box1Abb[3870]=2. + Box1Abb[3869]*m_z1k;

  Box1Abb[3871]=-1. + 4.*Box1Abb[2312]*m_z1k;

  Box1Abb[3872]=4. + Box1Abb[3871]*m_z1k;

  Box1Abb[3873]=Box1Abb[3539]*pow(Box1Abb[68],3.) - Box1Abb[3872]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3870]*m_z12_2 - Box1Abb[15]*Box1Abb[2119]*m_z12_3*m_z1k;

  Box1Abb[3874]=Box1Abb[3873]*m_x + Box1Abb[3856]*m_x_2 + Box1Abb[3868]*m_x_3 + Box1Abb[3862]*m_x_4 + Box1Abb[3857]*m_x_5 + Box1Abb[220]*m_x_6 + Box1Abb[170]*pow(Box1Abb[4],2.)*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[3875]=-1. + 3.*m_z12 + m_z1k;

  Box1Abb[3876]=-3. + 2.*m_z12 + 3.*m_z1k;

  Box1Abb[3877]=4. + 11.*m_z12;

  Box1Abb[3878]=-13. + Box1Abb[3877]*m_z12;

  Box1Abb[3879]=-14. + 9.*m_z12 + Box1Abb[3878]*m_z1k + 6.*Box1Abb[158]*m_z1k_2;

  Box1Abb[3880]=13. - 9.*Box1Abb[15]*m_z12 + 16.*m_z1k;

  Box1Abb[3881]=-20. + Box1Abb[2448]*m_z1k;

  Box1Abb[3882]=11. - Box1Abb[2260]*m_z1k;

  Box1Abb[3883]=-Box1Abb[665]*pow(Box1Abb[68],2.) - Box1Abb[3881]*Box1Abb[68]*m_z12 + Box1Abb[3882]*m_z12_2 + 6.*m_z12_3*m_z1k;

  Box1Abb[3884]=-13. + 16.*m_z1k;

  Box1Abb[3885]=25. + 2.*Box1Abb[475]*m_z1k;

  Box1Abb[3886]=-3. + Box1Abb[3885]*m_z1k;

  Box1Abb[3887]=5. + Box1Abb[3886]*m_z12 - 11.*Box1Abb[178]*m_z12_2*m_z1k + Box1Abb[3884]*m_z1k_2;

  Box1Abb[3888]=Box1Abb[3887]*m_x_2 + Box1Abb[3879]*m_x_3 + Box1Abb[3880]*m_x_4 + Box1Abb[129]*m_x_5 + Box1Abb[3883]*m_x*m_z1k + Box1Abb[3875]*Box1Abb[3876]*Box1Abb[68]*m_z12*m_z1k_2;

  Box1Abb[3889]=Box1Abb[3874]*m_cL - Box1Abb[3888]*Box1Abb[7]*m_cR;

  Box1Abb[3890]=-1. + 2.*m_z12 - 4.*m_z1k;

  Box1Abb[3891]=Box1Abb[3890]*Box1Abb[68]*m_z12 - 3.*m_z1k;

  Box1Abb[3892]=-6. + m_z12_2;

  Box1Abb[3893]=-11. + 2.*Box1Abb[251]*m_z12;

  Box1Abb[3894]=2. + Box1Abb[3893]*m_z12;

  Box1Abb[3895]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box1Abb[1422]*Box1Abb[3892]*m_z1k + 2.*Box1Abb[3894]*m_z1k_2 - 8.*Box1Abb[0]*Box1Abb[493]*m_z1k_3;

  Box1Abb[3896]=-4. + Box1Abb[707]*m_z12;

  Box1Abb[3897]=26. + 7.*m_z12;

  Box1Abb[3898]=20. - Box1Abb[3897]*m_z12;

  Box1Abb[3899]=1. + 6.*Box1Abb[0]*m_z12 + 16.*m_z1k - Box1Abb[3896]*m_z12*m_z1k + Box1Abb[3898]*m_z1k_2 + 10.*m_z12*m_z1k_3;

  Box1Abb[3900]=11. - 4.*m_z1k;

  Box1Abb[3901]=-1. + Box1Abb[3900]*m_z1k;

  Box1Abb[3902]=-3. + m_z1k + 4.*m_z1k_2;

  Box1Abb[3903]=-5. + 4.*Box1Abb[210]*m_z1k;

  Box1Abb[3904]=-2.*pow(Box1Abb[68],3.) + 2.*Box1Abb[3902]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[3903]*Box1Abb[68]*m_z12_2 + Box1Abb[3901]*m_z12_3;

  Box1Abb[3905]=1. + 4.*Box1Abb[457]*m_z1k;

  Box1Abb[3906]=-3. + Box1Abb[1945]*m_z1k;

  Box1Abb[3907]=6. - 12.*m_z1k + 5.*m_z1k_3;

  Box1Abb[3908]=1. + Box1Abb[3907]*m_z1k_2;

  Box1Abb[3909]=-24. + 11.*m_z1k;

  Box1Abb[3910]=2. + Box1Abb[3909]*m_z1k;

  Box1Abb[3911]=10. + Box1Abb[3910]*m_z1k;

  Box1Abb[3912]=1. + Box1Abb[3911]*m_z1k;

  Box1Abb[3913]=Box1Abb[3905]*pow(Box1Abb[68],2.) - 2.*Box1Abb[3908]*m_z12 + Box1Abb[3912]*m_z12_2 + Box1Abb[3906]*m_z12_3*m_z1k;

  Box1Abb[3914]=Box1Abb[3913]*m_x_2 + Box1Abb[3895]*m_x_3 + Box1Abb[3899]*m_x_4 + 2.*Box1Abb[3891]*m_x_5 + Box1Abb[63]*m_x_6*m_z12 + Box1Abb[3904]*Box1Abb[68]*m_x*m_z1k - Box1Abb[4]*pow(Box1Abb[68],3.)*Box1Abb[739]*m_z12*m_z1k_2;

  Box1Abb[3915]=-2. + 5.*m_z12 + 14.*m_z12*m_z1k + 2.*m_z1k_2;

  Box1Abb[3916]=-8. + Box1Abb[1737]*m_z12;

  Box1Abb[3917]=-18. + 37.*m_z12;

  Box1Abb[3918]=6. + Box1Abb[3917]*m_z12;

  Box1Abb[3919]=1. + 2.*Box1Abb[1280]*m_z12 + 3.*Box1Abb[3916]*m_z12*m_z1k + Box1Abb[3918]*m_z1k_2 + 10.*m_z12*m_z1k_3;

  Box1Abb[3920]=-1. + 6.*m_z12_2;

  Box1Abb[3921]=-14. + Box1Abb[527]*m_z12;

  Box1Abb[3922]=3. + Box1Abb[3921]*m_z12;

  Box1Abb[3923]=5. + 6.*m_z12 + 3.*m_z12_3;

  Box1Abb[3924]=121. + 8.*m_z12;

  Box1Abb[3925]=-200. + Box1Abb[3924]*m_z12;

  Box1Abb[3926]=26. + Box1Abb[3925]*m_z12;

  Box1Abb[3927]=22. + 5.*m_z12;

  Box1Abb[3928]=-6. + Box1Abb[3927]*m_z12;

  Box1Abb[3929]=-pow(Box1Abb[0],2.) - 2.*Box1Abb[0]*Box1Abb[3920]*m_z1k + Box1Abb[3922]*m_z1k_2 + 2.*Box1Abb[3923]*m_z1k_3 + Box1Abb[3926]*m_z1k_4 + 6.*Box1Abb[3928]*m_z1k_5 + 10.*m_z12*m_z1k_6;

  Box1Abb[3930]=1. - 3.*m_z1k + 2.*m_z1k_3;

  Box1Abb[3931]=-3. + 2.*Box1Abb[2087]*m_z1k;

  Box1Abb[3932]=-2.*pow(Box1Abb[68],2.) + 5.*Box1Abb[3930]*m_z12 + Box1Abb[3931]*m_z12_2;

  Box1Abb[3933]=10. + 13.*m_z1k;

  Box1Abb[3934]=5. + Box1Abb[3933]*m_z1k;

  Box1Abb[3935]=23. + 5.*Box1Abb[475]*m_z1k;

  Box1Abb[3936]=5. + Box1Abb[3935]*m_z1k;

  Box1Abb[3937]=-41. + 18.*Box1Abb[2406]*m_z1k;

  Box1Abb[3938]=-52. + Box1Abb[3937]*m_z1k;

  Box1Abb[3939]=-12. + Box1Abb[3938]*m_z1k;

  Box1Abb[3940]=3. + Box1Abb[3939]*m_z12 + 2.*Box1Abb[3936]*m_z12_2 + 2.*Box1Abb[3934]*m_z1k + 12.*Box1Abb[15]*m_z12_3*m_z1k;

  Box1Abb[3941]=-1. + 2.*Box1Abb[3089]*m_z1k;

  Box1Abb[3942]=26. + 9.*m_z1k;

  Box1Abb[3943]=8. + Box1Abb[3942]*m_z1k;

  Box1Abb[3944]=-3. + Box1Abb[3943]*m_z1k_2;

  Box1Abb[3945]=-25. + 9.*Box1Abb[210]*m_z1k;

  Box1Abb[3946]=-3. + Box1Abb[3945]*m_z1k;

  Box1Abb[3947]=-1. + 2.*Box1Abb[3946]*m_z1k;

  Box1Abb[3948]=86. + 25.*m_z1k;

  Box1Abb[3949]=-30. + Box1Abb[3948]*m_z1k;

  Box1Abb[3950]=-2. + Box1Abb[3949]*m_z1k;

  Box1Abb[3951]=-5. + Box1Abb[3950]*m_z1k;

  Box1Abb[3952]=Box1Abb[3941]*pow(Box1Abb[68],3.) - Box1Abb[3947]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[3951]*Box1Abb[68]*m_z12_2 - Box1Abb[3944]*m_z12_3;

  Box1Abb[3953]=9. + 4.*m_z1k;

  Box1Abb[3954]=9. + Box1Abb[3953]*m_z1k;

  Box1Abb[3955]=-5. + 11.*Box1Abb[431]*m_z1k;

  Box1Abb[3956]=11. + 2.*Box1Abb[3955]*m_z1k;

  Box1Abb[3957]=52. + 9.*m_z1k;

  Box1Abb[3958]=2. + Box1Abb[3957]*m_z1k;

  Box1Abb[3959]=-8. + Box1Abb[3958]*m_z1k;

  Box1Abb[3960]=5. - Box1Abb[3959]*m_z1k;

  Box1Abb[3961]=-84. + 5.*m_z1k;

  Box1Abb[3962]=29. + Box1Abb[3961]*m_z1k;

  Box1Abb[3963]=2. + Box1Abb[3962]*m_z1k;

  Box1Abb[3964]=-31. + 2.*Box1Abb[3963]*m_z1k;

  Box1Abb[3965]=-8. + Box1Abb[3964]*m_z1k;

  Box1Abb[3966]=3. + Box1Abb[3965]*m_z12 + Box1Abb[3960]*m_z12_2 + Box1Abb[3956]*m_z1k + 2.*Box1Abb[3954]*m_z12_3*m_z1k;

  Box1Abb[3967]=Box1Abb[3929]*m_x_3 + Box1Abb[3966]*m_x_4 - Box1Abb[3940]*m_x_5 + Box1Abb[3919]*m_x_6 - Box1Abb[3915]*m_x_7*m_z12 + m_x_8*m_z12_2 + Box1Abb[3952]*m_x_2*m_z1k + Box1Abb[3932]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_x*m_z1k_2 - 2.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12*m_z1k_3;

  Box1Abb[3968]=Box1Abb[3967]*m_cL + Box1Abb[3914]*Box1Abb[7]*m_cR*m_z1k;

  Box1Abb[3969]=m_z12_2 - 4.*m_z1k + 3.*m_z12*m_z1k;

  Box1Abb[3970]=-8. + 4.*m_z12 + m_z12_3;

  Box1Abb[3971]=37. + 2.*m_z12;

  Box1Abb[3972]=-6. + Box1Abb[3971]*m_z12;

  Box1Abb[3973]=8. - 3.*m_z12;

  Box1Abb[3974]=-1. + 6.*m_z12 - 6.*m_z12_2 + Box1Abb[3970]*m_z1k + Box1Abb[3972]*m_z1k_2 + 5.*Box1Abb[3973]*m_z1k_3;

  Box1Abb[3975]=-13. + 3.*m_z12;

  Box1Abb[3976]=2. + Box1Abb[3975]*m_z12;

  Box1Abb[3977]=2. + Box1Abb[3976]*m_z12;

  Box1Abb[3978]=-1. + Box1Abb[251]*m_z12;

  Box1Abb[3979]=-33. + 5.*m_z12;

  Box1Abb[3980]=24. + Box1Abb[3979]*m_z12;

  Box1Abb[3981]=-2. + 6.*m_z12 - 4.*m_z12_2 + Box1Abb[3977]*m_z1k + 4.*Box1Abb[3978]*m_z12*m_z1k_2 - 2.*Box1Abb[3980]*m_z1k_3 + 40.*m_z1k_4;

  Box1Abb[3982]=2. + 3.*Box1Abb[856]*m_z12;

  Box1Abb[3983]=10. + Box1Abb[2513]*m_z12;

  Box1Abb[3984]=11. + Box1Abb[3983]*m_z12;

  Box1Abb[3985]=-9. + m_z12;

  Box1Abb[3986]=25. + 2.*Box1Abb[3985]*m_z12;

  Box1Abb[3987]=-26. + Box1Abb[3986]*m_z12;

  Box1Abb[3988]=58. + Box1Abb[1482]*m_z12;

  Box1Abb[3989]=pow(Box1Abb[0],2.) - Box1Abb[0]*Box1Abb[3982]*m_z1k + Box1Abb[3984]*m_z1k_2 + 2.*Box1Abb[3987]*m_z1k_3 + Box1Abb[3988]*m_z1k_4 - 5.*Box1Abb[1793]*m_z1k_5;

  Box1Abb[3990]=5. - 12.*m_z1k;

  Box1Abb[3991]=2. + Box1Abb[3990]*m_z1k;

  Box1Abb[3992]=Box1Abb[3991]*m_z12 + 2.*Box1Abb[170]*m_z12_2 + 5.*Box1Abb[174]*m_z1k;

  Box1Abb[3993]=-11. + 4.*m_z1k;

  Box1Abb[3994]=1. + Box1Abb[3993]*m_z1k;

  Box1Abb[3995]=5. + m_z1k - 4.*m_z1k_2;

  Box1Abb[3996]=-2. + Box1Abb[3995]*m_z1k;

  Box1Abb[3997]=-Box1Abb[665]*pow(Box1Abb[68],2.) + 3.*Box1Abb[3996]*m_z12 + Box1Abb[3994]*m_z12_2;

  Box1Abb[3998]=Box1Abb[3989]*m_x_2 + Box1Abb[3981]*m_x_3 - Box1Abb[3974]*m_x_4 + Box1Abb[3992]*m_x_5 + Box1Abb[3969]*m_x_6 - Box1Abb[3997]*Box1Abb[4]*Box1Abb[68]*m_x*m_z1k - 3.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[3999]=13. + 6.*m_z12;

  Box1Abb[4000]=-18. + Box1Abb[3999]*m_z12;

  Box1Abb[4001]=33. + 10.*m_z12;

  Box1Abb[4002]=-16. + Box1Abb[4001]*m_z12;

  Box1Abb[4003]=9. + Box1Abb[4002]*m_z12;

  Box1Abb[4004]=1. - 9.*m_z12 + 6.*m_z12_2;

  Box1Abb[4005]=3. + 2.*Box1Abb[699]*m_z12 + 8.*m_z1k + 2.*Box1Abb[4000]*m_z12*m_z1k + Box1Abb[4003]*m_z1k_2 + 6.*Box1Abb[4004]*m_z1k_3 + 5.*Box1Abb[3530]*m_z1k_4;

  Box1Abb[4006]=5. + 12.*m_z1k;

  Box1Abb[4007]=2. - 7.*m_z1k_2;

  Box1Abb[4008]=Box1Abb[4007]*m_z12 - Box1Abb[4006]*m_z12_2 + 4.*m_z1k_2;

  Box1Abb[4009]=3. - 2.*m_z1k + 8.*m_z1k_2;

  Box1Abb[4010]=1. + Box1Abb[536]*m_z1k;

  Box1Abb[4011]=-pow(Box1Abb[68],2.)*Box1Abb[834] + 4.*Box1Abb[4010]*Box1Abb[68]*m_z12 + Box1Abb[4009]*m_z12_2;

  Box1Abb[4012]=10. + 7.*Box1Abb[1945]*m_z1k;

  Box1Abb[4013]=-5. + Box1Abb[2066]*m_z1k;

  Box1Abb[4014]=-2. + Box1Abb[4013]*m_z1k;

  Box1Abb[4015]=1. + 4.*Box1Abb[4014]*m_z12 + Box1Abb[4012]*m_z12_2 + 3.*m_z12_3*m_z1k - 5.*Box1Abb[174]*m_z1k_2;

  Box1Abb[4016]=12. - 9.*m_z1k + 5.*m_z1k_3;

  Box1Abb[4017]=-1. + 2.*m_z1k + 20.*m_z1k_2;

  Box1Abb[4018]=1. + Box1Abb[4017]*m_z1k;

  Box1Abb[4019]=-17. + 23.*m_z1k;

  Box1Abb[4020]=-7. + Box1Abb[4019]*m_z1k;

  Box1Abb[4021]=21. + 2.*Box1Abb[4020]*m_z1k;

  Box1Abb[4022]=-1. + Box1Abb[4021]*m_z1k;

  Box1Abb[4023]=-6. + Box1Abb[1980]*m_z1k;

  Box1Abb[4024]=10. + Box1Abb[4023]*m_z1k;

  Box1Abb[4025]=-2. + Box1Abb[4024]*m_z1k;

  Box1Abb[4026]=-Box1Abb[4018]*pow(Box1Abb[68],3.) + Box1Abb[4025]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4022]*Box1Abb[68]*m_z12_2 + Box1Abb[4016]*m_z12_3*m_z1k;

  Box1Abb[4027]=18. + 5.*Box1Abb[178]*m_z1k;

  Box1Abb[4028]=-17. + m_z1k + 48.*m_z1k_2 - 40.*m_z1k_3;

  Box1Abb[4029]=5. + Box1Abb[4028]*m_z1k;

  Box1Abb[4030]=-17. + 47.*m_z1k;

  Box1Abb[4031]=7. + Box1Abb[4030]*m_z1k;

  Box1Abb[4032]=-12. + Box1Abb[4031]*m_z1k;

  Box1Abb[4033]=5. + Box1Abb[4032]*m_z1k;

  Box1Abb[4034]=47. + 16.*Box1Abb[1103]*m_z1k;

  Box1Abb[4035]=8. + Box1Abb[4034]*m_z1k;

  Box1Abb[4036]=-7. + Box1Abb[4035]*m_z1k;

  Box1Abb[4037]=-8. + Box1Abb[4036]*m_z1k;

  Box1Abb[4038]=3. + Box1Abb[4037]*m_z12 + Box1Abb[4033]*m_z12_2 + Box1Abb[4029]*m_z1k + Box1Abb[4027]*m_z12_3*m_z1k;

  Box1Abb[4039]=-Box1Abb[4026]*m_x_2 + Box1Abb[4038]*m_x_3 - Box1Abb[4005]*m_x_4 + Box1Abb[4015]*m_x_5 + Box1Abb[4008]*m_x_6 + m_x_7*m_z12_2 + Box1Abb[4]*Box1Abb[4011]*pow(Box1Abb[68],2.)*m_x*m_z1k - pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12*m_z1k_2;

  Box1Abb[4040]=Box1Abb[4039]*m_cL + Box1Abb[3998]*m_cR*m_z1k;

  Box1Abb[4041]=3.*m_z12_2 - 4.*m_z1k + 7.*m_z12*m_z1k;

  Box1Abb[4042]=5. - 4.*m_z12;

  Box1Abb[4043]=3. + 5.*m_z12;

  Box1Abb[4044]=12. + Box1Abb[4043]*m_z12;

  Box1Abb[4045]=-47. + 19.*m_z12;

  Box1Abb[4046]=12. + Box1Abb[4045]*m_z12;

  Box1Abb[4047]=-8. + 6.*Box1Abb[4042]*m_z12 + 6.*m_z1k - Box1Abb[4044]*m_z12*m_z1k + 2.*Box1Abb[0]*Box1Abb[12]*Box1Abb[251]*m_z1k_2 - 2.*Box1Abb[4046]*m_z1k_3 + 40.*Box1Abb[1408]*m_z1k_4;

  Box1Abb[4048]=5. + 36.*m_z1k;

  Box1Abb[4049]=6. - Box1Abb[4048]*m_z1k;

  Box1Abb[4050]=Box1Abb[4049]*m_z12 - 2.*Box1Abb[2087]*m_z12_2 + Box1Abb[3455]*m_z1k;

  Box1Abb[4051]=13. + 5.*Box1Abb[15]*m_z1k;

  Box1Abb[4052]=-1. + 5.*Box1Abb[174]*m_z1k;

  Box1Abb[4053]=-19. + 75.*m_z1k;

  Box1Abb[4054]=-22. + Box1Abb[4053]*m_z1k_2;

  Box1Abb[4055]=3. + Box1Abb[4054]*m_z12 + 2.*Box1Abb[4051]*m_z12_2 - 2.*Box1Abb[4052]*m_z1k + m_z12_3*m_z1k;

  Box1Abb[4056]=-1. + 10.*m_z1k;

  Box1Abb[4057]=-7. + Box1Abb[4056]*m_z1k;

  Box1Abb[4058]=-2. + 3.*Box1Abb[665]*m_z1k;

  Box1Abb[4059]=2. + Box1Abb[4058]*m_z1k;

  Box1Abb[4060]=-Box1Abb[3539]*pow(Box1Abb[68],2.) + Box1Abb[4059]*Box1Abb[68]*m_z12 + Box1Abb[4057]*m_z12_2*m_z1k;

  Box1Abb[4061]=9. + Box1Abb[1028]*m_z1k;

  Box1Abb[4062]=-1. + Box1Abb[4056]*m_z1k;

  Box1Abb[4063]=-7. + 2.*Box1Abb[4062]*m_z1k;

  Box1Abb[4064]=-73. + 45.*m_z1k;

  Box1Abb[4065]=-5. + Box1Abb[4064]*m_z1k;

  Box1Abb[4066]=-3. + Box1Abb[4065]*m_z1k;

  Box1Abb[4067]=18. + Box1Abb[4066]*m_z1k;

  Box1Abb[4068]=-50. + 49.*m_z1k;

  Box1Abb[4069]=-6. + Box1Abb[4068]*m_z1k;

  Box1Abb[4070]=-12. + Box1Abb[4069]*m_z1k;

  Box1Abb[4071]=11. + Box1Abb[4070]*m_z1k;

  Box1Abb[4072]=-Box1Abb[4063]*pow(Box1Abb[68],2.) + Box1Abb[4067]*Box1Abb[68]*m_z12 + Box1Abb[4071]*m_z12_2 + Box1Abb[4061]*m_z12_3*m_z1k;

  Box1Abb[4073]=-Box1Abb[4]*Box1Abb[4060]*Box1Abb[68]*m_x + Box1Abb[4072]*m_x_2 + Box1Abb[4047]*m_x_3 + Box1Abb[4055]*m_x_4 + Box1Abb[4050]*m_x_5 + Box1Abb[4041]*m_x_6 + Box1Abb[170]*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[4074]=Box1Abb[4073]*m_cL + Box1Abb[3831]*m_cR*m_z1k;

  Box1Abb[4075]=4. - 7.*m_z12 - 6.*m_z12*m_z1k - 4.*Box1Abb[0]*m_z1k_2;

  Box1Abb[4076]=-2. + m_z1k + 8.*m_z1k_2;

  Box1Abb[4077]=-2.*pow(Box1Abb[68],2.) + Box1Abb[4076]*Box1Abb[68]*m_z12 + 2.*Box1Abb[178]*m_z12_2*m_z1k;

  Box1Abb[4078]=Box1Abb[4077]*m_x + Box1Abb[4075]*m_x_2 + Box1Abb[743]*m_x_3 - 2.*Box1Abb[77]*m_x_4 - 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[4079]=4. + 6.*m_z12;

  Box1Abb[4080]=6. - Box1Abb[1868]*m_z12;

  Box1Abb[4081]=2. + Box1Abb[1989]*m_z1k;

  Box1Abb[4082]=-2.*Box1Abb[4081]*m_z12 + 3.*m_z12_2 - 4.*Box1Abb[881]*m_z1k;

  Box1Abb[4083]=-3. + Box1Abb[2088]*m_z1k;

  Box1Abb[4084]=-2.*pow(Box1Abb[68],2.) + Box1Abb[4083]*Box1Abb[68]*m_z12 + Box1Abb[3786]*m_z12_2;

  Box1Abb[4085]=Box1Abb[4084]*m_x + Box1Abb[4082]*m_x_2 + Box1Abb[4080]*m_x_3 + Box1Abb[4079]*m_x_4 - 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[4086]=Box1Abb[4085]*m_cL + Box1Abb[4078]*m_cR;

  Box1Abb[4087]=-2. - 2.*m_z12 + m_z12_2 + 3.*Box1Abb[77]*m_z1k - 4.*m_z1k_2;

  Box1Abb[4088]=-49. + Box1Abb[1818]*m_z12;

  Box1Abb[4089]=-11. + Box1Abb[4088]*m_z12;

  Box1Abb[4090]=16. + Box1Abb[4089]*m_z12;

  Box1Abb[4091]=-32. - 37.*m_z12 + 45.*m_z12_3;

  Box1Abb[4092]=42. + 29.*m_z12;

  Box1Abb[4093]=-2. + Box1Abb[4092]*m_z12;

  Box1Abb[4094]=24. - 3.*Box1Abb[12]*Box1Abb[3557]*m_z12 + Box1Abb[4090]*m_z1k + 2.*Box1Abb[4091]*m_z1k_2 + 4.*Box1Abb[4093]*m_z1k_3 - 80.*m_z12*m_z1k_4;

  Box1Abb[4095]=38. + 11.*m_z12 + 12.*m_z1k;

  Box1Abb[4096]=-16. + Box1Abb[4095]*m_z12;

  Box1Abb[4097]=44. - 30.*m_z1k;

  Box1Abb[4098]=27. + 14.*m_z1k;

  Box1Abb[4099]=75. + Box1Abb[4097]*m_z12 + 3.*m_z12_2 + 4.*Box1Abb[4098]*m_z1k;

  Box1Abb[4100]=-62. + Box1Abb[4099]*m_z12 - 44.*m_z1k;

  Box1Abb[4101]=-23. + 4.*Box1Abb[2087]*m_z1k;

  Box1Abb[4102]=2. + Box1Abb[4101]*m_z1k;

  Box1Abb[4103]=-28. + 9.*m_z1k;

  Box1Abb[4104]=25. + 2.*Box1Abb[4103]*m_z1k;

  Box1Abb[4105]=-1. + Box1Abb[4104]*m_z1k;

  Box1Abb[4106]=-2.*Box1Abb[317]*pow(Box1Abb[68],4.) + Box1Abb[4102]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[4105]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[1827]*Box1Abb[68]*m_z12_3 + Box1Abb[1832]*m_z12_4*m_z1k;

  Box1Abb[4107]=-15. + m_z1k - 8.*m_z1k_2;

  Box1Abb[4108]=4. + 135.*m_z1k;

  Box1Abb[4109]=55. - Box1Abb[4108]*m_z1k;

  Box1Abb[4110]=13. + 50.*m_z1k;

  Box1Abb[4111]=49. + Box1Abb[4110]*m_z1k;

  Box1Abb[4112]=37. + 2.*Box1Abb[4111]*m_z1k;

  Box1Abb[4113]=4.*Box1Abb[4107] + Box1Abb[4112]*m_z12 + Box1Abb[4109]*m_z12_2 + 2.*Box1Abb[691]*m_z12_3;

  Box1Abb[4114]=-3. + Box1Abb[233]*m_z1k;

  Box1Abb[4115]=16. + m_z1k;

  Box1Abb[4116]=5. + 3.*Box1Abb[4115]*m_z1k;

  Box1Abb[4117]=-54. + Box1Abb[4116]*m_z1k;

  Box1Abb[4118]=2. + Box1Abb[4117]*m_z1k;

  Box1Abb[4119]=-13. + 2.*m_z1k;

  Box1Abb[4120]=27. + 5.*Box1Abb[4119]*m_z1k;

  Box1Abb[4121]=45. + 2.*Box1Abb[4120]*m_z1k;

  Box1Abb[4122]=-11. + Box1Abb[4121]*m_z1k;

  Box1Abb[4123]=4.*Box1Abb[4114]*pow(Box1Abb[68],2.) + Box1Abb[4122]*Box1Abb[68]*m_z12 - Box1Abb[4118]*m_z12_2 + Box1Abb[1845]*m_z12_3 - Box1Abb[1843]*m_z12_4*m_z1k;

  Box1Abb[4124]=Box1Abb[4106]*m_x + Box1Abb[4123]*m_x_2 + Box1Abb[4094]*m_x_3 + Box1Abb[4113]*m_x_4 - Box1Abb[4100]*m_x_5 + Box1Abb[4096]*m_x_6 + Box1Abb[4]*Box1Abb[4087]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[4125]=18. + m_z12 + 12.*m_z1k;

  Box1Abb[4126]=-16. + Box1Abb[4125]*m_z12;

  Box1Abb[4127]=15. + 14.*m_z1k;

  Box1Abb[4128]=-49. + 2.*Box1Abb[1830]*m_z12 - 4.*Box1Abb[4127]*m_z1k;

  Box1Abb[4129]=58. + Box1Abb[4128]*m_z12 + 44.*m_z1k;

  Box1Abb[4130]=4. - 3.*m_z1k;

  Box1Abb[4131]=-2. + Box1Abb[317]*m_z1k;

  Box1Abb[4132]=-4. + Box1Abb[818]*m_z1k;

  Box1Abb[4133]=-5. + 19.*m_z1k;

  Box1Abb[4134]=-8. + Box1Abb[4133]*m_z1k;

  Box1Abb[4135]=2.*Box1Abb[4131]*pow(Box1Abb[68],3.) + Box1Abb[4134]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4132]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[4130]*m_z12_3*m_z1k;

  Box1Abb[4136]=17. + m_z1k + 8.*m_z1k_2;

  Box1Abb[4137]=14. + 11.*m_z1k;

  Box1Abb[4138]=26. + Box1Abb[4137]*m_z1k;

  Box1Abb[4139]=18. + 5.*Box1Abb[234]*m_z1k;

  Box1Abb[4140]=33. + 2.*Box1Abb[4139]*m_z1k;

  Box1Abb[4141]=-4.*Box1Abb[4136] + Box1Abb[4140]*m_z12 + Box1Abb[4138]*m_z12_2 - 5.*m_z12_3*m_z1k;

  Box1Abb[4142]=3. + Box1Abb[2780]*m_z1k;

  Box1Abb[4143]=19. + 63.*m_z1k;

  Box1Abb[4144]=35. + Box1Abb[4143]*m_z1k;

  Box1Abb[4145]=14. + Box1Abb[4144]*m_z1k;

  Box1Abb[4146]=9. - 10.*m_z1k;

  Box1Abb[4147]=35. + 4.*Box1Abb[4146]*m_z1k;

  Box1Abb[4148]=51. + 2.*Box1Abb[4147]*m_z1k;

  Box1Abb[4149]=13. + Box1Abb[4148]*m_z1k;

  Box1Abb[4150]=24. + Box1Abb[4149]*m_z12 - 2.*Box1Abb[4145]*m_z12_2 - 8.*Box1Abb[4142]*m_z1k + Box1Abb[1378]*m_z12_3*m_z1k;

  Box1Abb[4151]=-85. + 52.*m_z1k;

  Box1Abb[4152]=10. + Box1Abb[4151]*m_z1k;

  Box1Abb[4153]=19. + Box1Abb[4152]*m_z1k;

  Box1Abb[4154]=-29. + 4.*m_z1k + 8.*m_z1k_2;

  Box1Abb[4155]=3. + Box1Abb[4154]*m_z1k;

  Box1Abb[4156]=4. + Box1Abb[4155]*m_z1k;

  Box1Abb[4157]=-43. + 10.*m_z1k;

  Box1Abb[4158]=18. + Box1Abb[4157]*m_z1k;

  Box1Abb[4159]=11. + Box1Abb[4158]*m_z1k;

  Box1Abb[4160]=-1. + Box1Abb[4159]*m_z1k;

  Box1Abb[4161]=2.*Box1Abb[431]*pow(Box1Abb[68],4.) - Box1Abb[4156]*pow(Box1Abb[68],2.)*m_z12 + 2.*Box1Abb[4160]*Box1Abb[68]*m_z12_2 + Box1Abb[4153]*m_z12_3*m_z1k + 2.*Box1Abb[119]*m_z12_4*m_z1k_2;

  Box1Abb[4162]=1. + m_z1k + 4.*m_z1k_2;

  Box1Abb[4163]=21. + 25.*Box1Abb[431]*m_z1k;

  Box1Abb[4164]=-42. + 127.*m_z1k - 66.*m_z1k_3 + 20.*m_z1k_4;

  Box1Abb[4165]=-19. + Box1Abb[4164]*m_z1k;

  Box1Abb[4166]=-124. + 135.*m_z1k;

  Box1Abb[4167]=-78. + Box1Abb[4166]*m_z1k;

  Box1Abb[4168]=24. + Box1Abb[4167]*m_z1k;

  Box1Abb[4169]=13. + Box1Abb[4168]*m_z1k;

  Box1Abb[4170]=4.*Box1Abb[4162]*pow(Box1Abb[68],2.) + Box1Abb[4165]*m_z12 + Box1Abb[4169]*m_z12_2 + Box1Abb[4163]*m_z12_3*m_z1k - 2.*m_z12_4*m_z1k_2;

  Box1Abb[4171]=-Box1Abb[4161]*m_x_2 + Box1Abb[4170]*m_x_3 + Box1Abb[4150]*m_x_4 + Box1Abb[4141]*m_x_5 + Box1Abb[4129]*m_x_6 + Box1Abb[4126]*m_x_7 - Box1Abb[4135]*Box1Abb[68]*m_x*m_z12*m_z1k + 2.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_2;

  Box1Abb[4172]=-Box1Abb[4171]*m_cR + Box1Abb[4124]*m_cL*m_x;

  Box1Abb[4173]=-2.*Box1Abb[317]*pow(Box1Abb[68],2.) + Box1Abb[119]*Box1Abb[68]*m_z12 + Box1Abb[2066]*m_z12_2;

  Box1Abb[4174]=22. + 11.*m_z12 + 12.*m_z1k;

  Box1Abb[4175]=-16. + Box1Abb[4174]*m_z12;

  Box1Abb[4176]=20. + m_z1k;

  Box1Abb[4177]=11. + 2.*Box1Abb[4176]*m_z12 + 3.*m_z12_2 + 60.*m_z1k + 56.*m_z1k_2;

  Box1Abb[4178]=-46. + Box1Abb[4177]*m_z12 - 44.*m_z1k;

  Box1Abb[4179]=-12. + m_z1k + 26.*m_z1k_2;

  Box1Abb[4180]=-17. + 10.*m_z1k;

  Box1Abb[4181]=1. + 8.*m_z1k_2;

  Box1Abb[4182]=43. - Box1Abb[4180]*Box1Abb[4181]*m_z1k;

  Box1Abb[4183]=73. + 72.*m_z1k + 76.*m_z1k_2;

  Box1Abb[4184]=-3. + Box1Abb[4183]*m_z1k;

  Box1Abb[4185]=-6. + Box1Abb[175]*m_z1k;

  Box1Abb[4186]=3. + Box1Abb[4185]*m_z1k;

  Box1Abb[4187]=-8.*Box1Abb[4186] + Box1Abb[4182]*m_z12 - Box1Abb[4184]*m_z12_2 + Box1Abb[4179]*m_z12_3 + m_z12_4*m_z1k;

  Box1Abb[4188]=18. + 5.*m_z1k;

  Box1Abb[4189]=-3. + Box1Abb[4188]*m_z1k;

  Box1Abb[4190]=-2. + 9.*m_z1k;

  Box1Abb[4191]=-1. + Box1Abb[1041]*Box1Abb[4190]*m_z1k;

  Box1Abb[4192]=-7. + 4.*Box1Abb[873]*m_z1k;

  Box1Abb[4193]=6. + Box1Abb[4192]*m_z1k;

  Box1Abb[4194]=9. + 10.*Box1Abb[3396]*m_z1k;

  Box1Abb[4195]=-5. + Box1Abb[4194]*m_z1k;

  Box1Abb[4196]=2.*Box1Abb[317]*pow(Box1Abb[68],4.) - Box1Abb[4193]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[4195]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[4191]*Box1Abb[68]*m_z12_3 + Box1Abb[4189]*m_z12_4*m_z1k;

  Box1Abb[4197]=10. - 4.*m_z1k;

  Box1Abb[4198]=64. - 7.*m_z1k;

  Box1Abb[4199]=39. + Box1Abb[4198]*m_z1k;

  Box1Abb[4200]=3. + Box1Abb[1830]*m_z1k;

  Box1Abb[4201]=-3. + 50.*m_z1k;

  Box1Abb[4202]=1. + Box1Abb[4201]*m_z1k;

  Box1Abb[4203]=-55. + 2.*Box1Abb[4202]*m_z1k;

  Box1Abb[4204]=-4.*Box1Abb[4200] + Box1Abb[4203]*m_z12 + Box1Abb[4199]*m_z12_2 + Box1Abb[4197]*m_z12_3;

  Box1Abb[4205]=89. + 32.*m_z1k;

  Box1Abb[4206]=-1. + Box1Abb[4205]*m_z1k;

  Box1Abb[4207]=6. + Box1Abb[4206]*m_z1k;

  Box1Abb[4208]=-41. + 10.*m_z1k;

  Box1Abb[4209]=3. + Box1Abb[4208]*m_z1k;

  Box1Abb[4210]=25. + 2.*Box1Abb[4209]*m_z1k;

  Box1Abb[4211]=-7. + Box1Abb[4210]*m_z1k;

  Box1Abb[4212]=-24. + 25.*m_z1k;

  Box1Abb[4213]=-53. + 5.*Box1Abb[4212]*m_z1k;

  Box1Abb[4214]=30. + Box1Abb[4213]*m_z1k;

  Box1Abb[4215]=-18. + Box1Abb[4214]*m_z1k;

  Box1Abb[4216]=4.*Box1Abb[1712]*pow(Box1Abb[68],2.) + Box1Abb[4211]*Box1Abb[68]*m_z12 + Box1Abb[4215]*m_z12_2 + Box1Abb[4207]*m_z12_3 - Box1Abb[514]*m_z12_4*m_z1k;

  Box1Abb[4217]=-Box1Abb[4196]*m_x + Box1Abb[4216]*m_x_2 + Box1Abb[4187]*m_x_3 + Box1Abb[4204]*m_x_4 - Box1Abb[4178]*m_x_5 + Box1Abb[4175]*m_x_6 + Box1Abb[4]*Box1Abb[4173]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[4218]=17. + 32.*m_z12 - 8.*m_z1k;

  Box1Abb[4219]=-6. + Box1Abb[4218]*m_z12;

  Box1Abb[4220]=-2. + m_z1k + 9.*m_z1k_2;

  Box1Abb[4221]=-1. + 3.*Box1Abb[178]*m_z1k;

  Box1Abb[4222]=-1. + Box1Abb[437]*m_z1k;

  Box1Abb[4223]=Box1Abb[4222]*pow(Box1Abb[68],2.) + Box1Abb[4220]*Box1Abb[68]*m_z12 + Box1Abb[4221]*m_z12_2;

  Box1Abb[4224]=67. - 8.*m_z1k;

  Box1Abb[4225]=95. + 119.*m_z1k;

  Box1Abb[4226]=-39. + Box1Abb[4225]*m_z12 + 16.*m_z12_2 + Box1Abb[4224]*m_z1k;

  Box1Abb[4227]=23. - Box1Abb[4226]*m_z12 + 20.*m_z1k;

  Box1Abb[4228]=197. + 157.*m_z1k;

  Box1Abb[4229]=65. + Box1Abb[4228]*m_z1k;

  Box1Abb[4230]=69. + 8.*m_z1k;

  Box1Abb[4231]=-131. + Box1Abb[4230]*m_z1k;

  Box1Abb[4232]=-110. + Box1Abb[4231]*m_z1k;

  Box1Abb[4233]=24. + Box1Abb[4232]*m_z12 + Box1Abb[4229]*m_z12_2 + 10.*Box1Abb[1943]*m_z12_3 + 2.*m_z12_4 - 3.*Box1Abb[1989]*m_z1k;

  Box1Abb[4234]=-2. + m_z1k - 31.*m_z1k_2 + 38.*m_z1k_3 + 26.*m_z1k_4;

  Box1Abb[4235]=-34. + m_z1k + 8.*m_z1k_2;

  Box1Abb[4236]=9. + Box1Abb[4235]*m_z1k;

  Box1Abb[4237]=7. + Box1Abb[2440]*m_z1k;

  Box1Abb[4238]=1. + m_z1k - Box1Abb[4237]*m_z1k_2;

  Box1Abb[4239]=62. + 33.*m_z1k;

  Box1Abb[4240]=-7. + Box1Abb[4239]*m_z1k;

  Box1Abb[4241]=1. + Box1Abb[4240]*m_z1k;

  Box1Abb[4242]=-Box1Abb[4241]*pow(Box1Abb[68],3.)*m_z12_2 - Box1Abb[4234]*Box1Abb[68]*m_z12_3 + Box1Abb[4238]*m_z12_4 + Box1Abb[437]*pow(Box1Abb[68],4.)*m_z1k - Box1Abb[4236]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[4243]=13. + 4.*m_z1k;

  Box1Abb[4244]=-4. + Box1Abb[4243]*m_z1k;

  Box1Abb[4245]=93. + 88.*m_z1k;

  Box1Abb[4246]=53. + Box1Abb[4245]*m_z1k;

  Box1Abb[4247]=76. + 71.*m_z1k;

  Box1Abb[4248]=-26. + Box1Abb[4247]*m_z1k;

  Box1Abb[4249]=-26. + Box1Abb[4248]*m_z1k;

  Box1Abb[4250]=-41. + 20.*m_z1k;

  Box1Abb[4251]=-84. + Box1Abb[4250]*m_z1k;

  Box1Abb[4252]=-17. + Box1Abb[4251]*m_z1k;

  Box1Abb[4253]=-46. + Box1Abb[4252]*m_z1k;

  Box1Abb[4254]=22. + Box1Abb[4253]*m_z12 + Box1Abb[4249]*m_z12_2 + Box1Abb[4246]*m_z12_3 + Box1Abb[1900]*m_z12_4 + 2.*Box1Abb[4244]*m_z1k;

  Box1Abb[4255]=3. + 4.*Box1Abb[15]*m_z1k;

  Box1Abb[4256]=11. + 28.*m_z1k + 22.*m_z1k_2;

  Box1Abb[4257]=18. + Box1Abb[4256]*m_z1k;

  Box1Abb[4258]=32. + 9.*m_z1k;

  Box1Abb[4259]=88. + 3.*Box1Abb[4258]*m_z1k;

  Box1Abb[4260]=-35. + m_z1k - Box1Abb[4259]*m_z1k_2;

  Box1Abb[4261]=-3. + 11.*m_z1k;

  Box1Abb[4262]=-16. + Box1Abb[4261]*m_z1k;

  Box1Abb[4263]=-3. + Box1Abb[4262]*m_z1k;

  Box1Abb[4264]=-109. + 8.*m_z1k;

  Box1Abb[4265]=150. + Box1Abb[4264]*m_z1k;

  Box1Abb[4266]=68. + Box1Abb[4265]*m_z1k;

  Box1Abb[4267]=2. + Box1Abb[4266]*m_z1k;

  Box1Abb[4268]=3. + Box1Abb[4267]*m_z1k;

  Box1Abb[4269]=6. + Box1Abb[4268]*m_z12 + Box1Abb[4260]*m_z12_2 + Box1Abb[4257]*m_z12_3 + 3.*Box1Abb[4255]*m_z12_4 + 2.*Box1Abb[4263]*m_z1k;

  Box1Abb[4270]=5. + Box1Abb[15]*Box1Abb[178]*m_z1k;

  Box1Abb[4271]=-11. + 12.*m_z1k;

  Box1Abb[4272]=-10. + Box1Abb[4271]*m_z1k;

  Box1Abb[4273]=1. + Box1Abb[4272]*m_z1k;

  Box1Abb[4274]=45. + 46.*m_z1k + 28.*m_z1k_2;

  Box1Abb[4275]=-8. + Box1Abb[4274]*m_z1k;

  Box1Abb[4276]=3. + Box1Abb[4275]*m_z1k;

  Box1Abb[4277]=47. + 8.*m_z1k;

  Box1Abb[4278]=-142. + Box1Abb[4277]*m_z1k;

  Box1Abb[4279]=34. + Box1Abb[4278]*m_z1k;

  Box1Abb[4280]=18. + Box1Abb[4279]*m_z1k;

  Box1Abb[4281]=3. + Box1Abb[4280]*m_z1k;

  Box1Abb[4282]=73. + 51.*m_z1k;

  Box1Abb[4283]=-153. + Box1Abb[4282]*m_z1k;

  Box1Abb[4284]=-20. + Box1Abb[4283]*m_z1k;

  Box1Abb[4285]=11. + Box1Abb[4284]*m_z1k;

  Box1Abb[4286]=6. + Box1Abb[4285]*m_z1k;

  Box1Abb[4287]=-Box1Abb[4273]*pow(Box1Abb[68],2.) + Box1Abb[4281]*Box1Abb[68]*m_z12 + Box1Abb[4286]*m_z12_2 + Box1Abb[4276]*m_z12_3 - Box1Abb[4270]*m_z12_4;

  Box1Abb[4288]=Box1Abb[4242]*m_x_2 + Box1Abb[4287]*m_x_3 + Box1Abb[4269]*m_x_4 - Box1Abb[4254]*m_x_5 + Box1Abb[4233]*m_x_6 + Box1Abb[4227]*m_x_7 + Box1Abb[4219]*m_x_8 + 2.*m_x_9*m_z12 + Box1Abb[4]*Box1Abb[4223]*pow(Box1Abb[68],2.)*m_x*m_z12*m_z1k - pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12_2*m_z1k_2;

  Box1Abb[4289]=2.*Box1Abb[4288]*m_cL - Box1Abb[4217]*Box1Abb[7]*m_cR*m_x;

  Box1Abb[4290]=1. + m_z12 - 6.*m_z1k - 9.*m_z12*m_z1k + 5.*m_z1k_2;

  Box1Abb[4291]=49. + 3.*Box1Abb[708]*m_z12;

  Box1Abb[4292]=42. + 5.*m_z12;

  Box1Abb[4293]=2. + 21.*m_z12;

  Box1Abb[4294]=-51. + Box1Abb[4291]*m_z12 - 61.*m_z1k + 3.*Box1Abb[4292]*m_z12*m_z1k + 5.*Box1Abb[4293]*m_z1k_2;

  Box1Abb[4295]=8.*Box1Abb[15] + Box1Abb[4294]*m_z12;

  Box1Abb[4296]=34. + 10.*m_z12 + 25.*m_z1k;

  Box1Abb[4297]=-3. + Box1Abb[4296]*m_z12 + 20.*m_z1k;

  Box1Abb[4298]=2. + Box1Abb[4297]*m_z12;

  Box1Abb[4299]=3. + 5.*Box1Abb[437]*m_z1k;

  Box1Abb[4300]=39. + 5.*m_z1k;

  Box1Abb[4301]=-12. + Box1Abb[4300]*m_z1k;

  Box1Abb[4302]=4. + Box1Abb[4301]*m_z1k;

  Box1Abb[4303]=-115. + 64.*m_z1k;

  Box1Abb[4304]=49. + Box1Abb[4303]*m_z1k;

  Box1Abb[4305]=1. + m_z1k + Box1Abb[4304]*m_z1k_2;

  Box1Abb[4306]=-Box1Abb[4299]*pow(Box1Abb[68],3.) - Box1Abb[4302]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4305]*m_z12_2 + Box1Abb[4189]*m_z12_3*m_z1k;

  Box1Abb[4307]=125. + 47.*m_z1k;

  Box1Abb[4308]=33. + Box1Abb[4307]*m_z1k;

  Box1Abb[4309]=11. + 4.*m_z1k;

  Box1Abb[4310]=2. + Box1Abb[4309]*m_z1k;

  Box1Abb[4311]=23. + 5.*Box1Abb[4310]*m_z1k;

  Box1Abb[4312]=62. + 185.*m_z1k;

  Box1Abb[4313]=31. + Box1Abb[4312]*m_z1k;

  Box1Abb[4314]=-2. + Box1Abb[4313]*m_z1k;

  Box1Abb[4315]=-4.*Box1Abb[823] + 2.*Box1Abb[4311]*m_z12 - Box1Abb[4314]*m_z12_2 - Box1Abb[4308]*m_z12_3 + Box1Abb[3285]*m_z12_4;

  Box1Abb[4316]=49. + 22.*m_z1k;

  Box1Abb[4317]=7. + Box1Abb[4316]*m_z1k;

  Box1Abb[4318]=3. + Box1Abb[4317]*m_z1k;

  Box1Abb[4319]=-25. + 44.*m_z1k;

  Box1Abb[4320]=-42. + Box1Abb[4319]*m_z1k;

  Box1Abb[4321]=7. + Box1Abb[4320]*m_z1k;

  Box1Abb[4322]=41. + 84.*m_z1k;

  Box1Abb[4323]=-76. - Box1Abb[4322]*Box1Abb[437]*m_z1k;

  Box1Abb[4324]=19. + Box1Abb[4323]*m_z1k;

  Box1Abb[4325]=-61. + 17.*m_z1k;

  Box1Abb[4326]=-13. + Box1Abb[4325]*m_z1k;

  Box1Abb[4327]=41. + Box1Abb[4326]*m_z1k;

  Box1Abb[4328]=-6. + Box1Abb[4327]*m_z1k;

  Box1Abb[4329]=-2.*pow(Box1Abb[68],4.) + Box1Abb[4321]*pow(Box1Abb[68],2.)*m_z12 - 3.*Box1Abb[4328]*Box1Abb[68]*m_z12_2 + Box1Abb[4324]*m_z12_3 - 2.*Box1Abb[4318]*m_z12_4 + Box1Abb[514]*m_z12_5*m_z1k;

  Box1Abb[4330]=5. - 8.*m_z1k;

  Box1Abb[4331]=6. + Box1Abb[4330]*m_z1k;

  Box1Abb[4332]=-4. + m_z1k + 51.*m_z1k_2 + 3.*m_z1k_3 - 35.*m_z1k_4;

  Box1Abb[4333]=22. + 25.*m_z1k;

  Box1Abb[4334]=163. + 6.*Box1Abb[4333]*m_z1k;

  Box1Abb[4335]=-7. + Box1Abb[4334]*m_z1k;

  Box1Abb[4336]=-204. + 155.*m_z1k;

  Box1Abb[4337]=-88. + Box1Abb[4336]*m_z1k;

  Box1Abb[4338]=-126. + Box1Abb[4337]*m_z1k;

  Box1Abb[4339]=-7. + Box1Abb[4338]*m_z1k;

  Box1Abb[4340]=8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4332]*m_z12 + Box1Abb[4339]*m_z12_2 + Box1Abb[4335]*m_z12_3 + 2.*Box1Abb[4331]*m_z12_4 - m_z12_5*m_z1k;

  Box1Abb[4341]=Box1Abb[4329]*m_x_2 + Box1Abb[4340]*m_x_3 + Box1Abb[4315]*m_x_4 + Box1Abb[4295]*m_x_5 - Box1Abb[4298]*m_x_6 + Box1Abb[4]*Box1Abb[4306]*m_x*m_z12 + Box1Abb[653]*m_x_7*m_z12 + pow(Box1Abb[4],2.)*Box1Abb[4290]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k;

  Box1Abb[4342]=-6. + 17.*m_z12;

  Box1Abb[4343]=-17. + 45.*m_z12 - 107.*m_z1k;

  Box1Abb[4344]=7. + Box1Abb[4343]*m_z12 + 32.*m_z1k;

  Box1Abb[4345]=2. + Box1Abb[4344]*m_z12;

  Box1Abb[4346]=53. + 14.*m_z12;

  Box1Abb[4347]=133. - 2.*Box1Abb[4346]*m_z12;

  Box1Abb[4348]=7. - 97.*m_z12;

  Box1Abb[4349]=-14. + 67.*m_z12;

  Box1Abb[4350]=9.*Box1Abb[178] + Box1Abb[4347]*m_z12 + 2.*Box1Abb[4348]*m_z12*m_z1k + 4.*Box1Abb[4349]*m_z1k_2;

  Box1Abb[4351]=-4.*Box1Abb[119] + Box1Abb[4350]*m_z12;

  Box1Abb[4352]=3. + 5.*Box1Abb[536]*m_z1k;

  Box1Abb[4353]=1. + Box1Abb[797]*m_z1k;

  Box1Abb[4354]=Box1Abb[4353]*pow(Box1Abb[68],2.) + Box1Abb[4352]*Box1Abb[68]*m_z12 + Box1Abb[3364]*m_z12_2;

  Box1Abb[4355]=41. + 75.*m_z1k;

  Box1Abb[4356]=-11. + 318.*m_z1k + 321.*m_z1k_2;

  Box1Abb[4357]=2. + Box1Abb[1151]*m_z1k;

  Box1Abb[4358]=36. + 97.*m_z1k;

  Box1Abb[4359]=66. - Box1Abb[4358]*m_z1k;

  Box1Abb[4360]=425. + 28.*Box1Abb[3440]*m_z1k;

  Box1Abb[4361]=195. + Box1Abb[4360]*m_z1k;

  Box1Abb[4362]=6.*Box1Abb[4357] + Box1Abb[4359]*m_z12 - Box1Abb[4361]*m_z12_2 + Box1Abb[4356]*m_z12_3 + 2.*Box1Abb[4355]*m_z12_4 + 4.*m_z12_5;

  Box1Abb[4363]=7. + 12.*m_z1k;

  Box1Abb[4364]=2. + Box1Abb[451]*m_z1k;

  Box1Abb[4365]=74. + 46.*Box1Abb[1989]*m_z1k;

  Box1Abb[4366]=388. + 285.*m_z1k;

  Box1Abb[4367]=-154. + Box1Abb[4366]*m_z1k;

  Box1Abb[4368]=-158. + Box1Abb[4367]*m_z1k;

  Box1Abb[4369]=34. - 35.*Box1Abb[233]*m_z1k;

  Box1Abb[4370]=46. + Box1Abb[4369]*m_z1k;

  Box1Abb[4371]=38. + Box1Abb[4370]*m_z1k;

  Box1Abb[4372]=-71. + 35.*m_z1k;

  Box1Abb[4373]=519. + 2.*Box1Abb[4372]*m_z1k;

  Box1Abb[4374]=244. + Box1Abb[4373]*m_z1k;

  Box1Abb[4375]=5. + Box1Abb[4374]*m_z1k;

  Box1Abb[4376]=8. + Box1Abb[4371]*m_z12 - Box1Abb[4375]*m_z12_2 + Box1Abb[4368]*m_z12_3 + Box1Abb[4365]*m_z12_4 + 2.*Box1Abb[4363]*m_z12_5 + 8.*Box1Abb[4364]*m_z1k;

  Box1Abb[4377]=-13. + 10.*m_z1k;

  Box1Abb[4378]=2. + Box1Abb[4377]*m_z1k;

  Box1Abb[4379]=-31. + 111.*m_z1k;

  Box1Abb[4380]=38. + Box1Abb[4379]*m_z1k;

  Box1Abb[4381]=-14. + Box1Abb[4380]*m_z1k;

  Box1Abb[4382]=10. + Box1Abb[536]*m_z1k;

  Box1Abb[4383]=-5. + Box1Abb[4382]*m_z1k;

  Box1Abb[4384]=2. + Box1Abb[4383]*m_z1k;

  Box1Abb[4385]=-25. + 32.*m_z1k;

  Box1Abb[4386]=21. + Box1Abb[4385]*m_z1k;

  Box1Abb[4387]=-19. + Box1Abb[4386]*m_z1k;

  Box1Abb[4388]=7. + Box1Abb[4387]*m_z1k;

  Box1Abb[4389]=-97. + 47.*m_z1k;

  Box1Abb[4390]=55. + Box1Abb[4389]*m_z1k;

  Box1Abb[4391]=-25. + Box1Abb[4390]*m_z1k;

  Box1Abb[4392]=4. + Box1Abb[4391]*m_z1k;

  Box1Abb[4393]=-Box1Abb[4392]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[4381]*pow(Box1Abb[68],3.)*m_z12_2 - 2.*Box1Abb[4388]*Box1Abb[68]*m_z12_3 - 2.*Box1Abb[4384]*m_z12_4 + Box1Abb[4378]*pow(Box1Abb[68],4.)*m_z1k;

  Box1Abb[4394]=38. + 45.*m_z1k;

  Box1Abb[4395]=18. + Box1Abb[4394]*m_z1k;

  Box1Abb[4396]=29. + 56.*m_z1k;

  Box1Abb[4397]=4. + Box1Abb[4396]*m_z1k;

  Box1Abb[4398]=4. + Box1Abb[1118]*m_z1k;

  Box1Abb[4399]=4. + Box1Abb[4398]*m_z1k;

  Box1Abb[4400]=210. + 251.*m_z1k;

  Box1Abb[4401]=-206. + Box1Abb[4400]*m_z1k;

  Box1Abb[4402]=-170. + Box1Abb[4401]*m_z1k;

  Box1Abb[4403]=-97. + Box1Abb[4402]*m_z1k;

  Box1Abb[4404]=125. - 224.*m_z1k;

  Box1Abb[4405]=24. + Box1Abb[4404]*m_z1k;

  Box1Abb[4406]=-100. + Box1Abb[4405]*m_z1k;

  Box1Abb[4407]=56. + Box1Abb[4406]*m_z1k;

  Box1Abb[4408]=-1. + Box1Abb[4407]*m_z1k;

  Box1Abb[4409]=-123. - 270.*m_z1k + 238.*m_z1k_2;

  Box1Abb[4410]=27. + Box1Abb[4409]*m_z1k;

  Box1Abb[4411]=-67. + Box1Abb[4410]*m_z1k;

  Box1Abb[4412]=61. + Box1Abb[4411]*m_z1k;

  Box1Abb[4413]=2. + Box1Abb[4408]*m_z12 + Box1Abb[4412]*m_z12_2 + Box1Abb[4403]*m_z12_3 + 2.*Box1Abb[431]*Box1Abb[4397]*m_z12_4 + Box1Abb[4395]*m_z12_5 + 2.*Box1Abb[4399]*m_z1k;

  Box1Abb[4414]=1. + 3.*m_z1k_2;

  Box1Abb[4415]=34. + 33.*m_z1k;

  Box1Abb[4416]=22. + Box1Abb[4415]*m_z1k;

  Box1Abb[4417]=10. + Box1Abb[4416]*m_z1k;

  Box1Abb[4418]=43. + 55.*m_z1k;

  Box1Abb[4419]=-17. + Box1Abb[4418]*m_z1k;

  Box1Abb[4420]=10. + Box1Abb[4419]*m_z1k;

  Box1Abb[4421]=-9. + Box1Abb[4420]*m_z1k;

  Box1Abb[4422]=-188. + 305.*m_z1k;

  Box1Abb[4423]=123. + Box1Abb[4422]*m_z1k;

  Box1Abb[4424]=150. + Box1Abb[4423]*m_z1k;

  Box1Abb[4425]=-146. + Box1Abb[4424]*m_z1k;

  Box1Abb[4426]=-6. + Box1Abb[4425]*m_z1k;

  Box1Abb[4427]=317. - 168.*m_z1k;

  Box1Abb[4428]=-222. + Box1Abb[4427]*m_z1k;

  Box1Abb[4429]=52. + Box1Abb[4428]*m_z1k;

  Box1Abb[4430]=52. + Box1Abb[4429]*m_z1k;

  Box1Abb[4431]=-13. + Box1Abb[4430]*m_z1k;

  Box1Abb[4432]=-2. + Box1Abb[4431]*m_z1k;

  Box1Abb[4433]=581. - 722.*m_z1k + 308.*m_z1k_2;

  Box1Abb[4434]=-252. + Box1Abb[4433]*m_z1k;

  Box1Abb[4435]=-200. + Box1Abb[4434]*m_z1k;

  Box1Abb[4436]=90. + Box1Abb[4435]*m_z1k;

  Box1Abb[4437]=15. + Box1Abb[4436]*m_z1k;

  Box1Abb[4438]=Box1Abb[4432]*m_z12 + Box1Abb[4437]*m_z12_2 + Box1Abb[4426]*m_z12_3 + 2.*Box1Abb[4421]*m_z12_4 + Box1Abb[4417]*m_z12_5 + 4.*Box1Abb[4414]*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[4439]=-25. + 43.*m_z1k;

  Box1Abb[4440]=59. + Box1Abb[4439]*m_z1k;

  Box1Abb[4441]=-3. - 8.*m_z1k + Box1Abb[4440]*m_z1k_3;

  Box1Abb[4442]=-2. + 5.*Box1Abb[15]*m_z1k;

  Box1Abb[4443]=6. + Box1Abb[4442]*m_z1k;

  Box1Abb[4444]=1. + Box1Abb[4443]*m_z1k;

  Box1Abb[4445]=-93. + 64.*m_z1k;

  Box1Abb[4446]=58. + Box1Abb[4445]*m_z1k;

  Box1Abb[4447]=-49. + Box1Abb[4446]*m_z1k;

  Box1Abb[4448]=12. + Box1Abb[4447]*m_z1k;

  Box1Abb[4449]=-113. + 43.*m_z1k;

  Box1Abb[4450]=421. + 4.*Box1Abb[4449]*m_z1k;

  Box1Abb[4451]=-300. + Box1Abb[4450]*m_z1k;

  Box1Abb[4452]=172. + Box1Abb[4451]*m_z1k;

  Box1Abb[4453]=-46. + Box1Abb[4452]*m_z1k;

  Box1Abb[4454]=1. + Box1Abb[4453]*m_z1k;

  Box1Abb[4455]=-494. + 259.*m_z1k;

  Box1Abb[4456]=469. + Box1Abb[4455]*m_z1k;

  Box1Abb[4457]=-370. + Box1Abb[4456]*m_z1k;

  Box1Abb[4458]=133. + Box1Abb[4457]*m_z1k;

  Box1Abb[4459]=-34. + Box1Abb[4458]*m_z1k;

  Box1Abb[4460]=5. + Box1Abb[4459]*m_z1k;

  Box1Abb[4461]=Box1Abb[4454]*Box1Abb[68]*m_z12_2 + Box1Abb[4460]*m_z12_3 + 2.*Box1Abb[4441]*m_z12_4 + 2.*Box1Abb[4444]*m_z12_5 - Box1Abb[4448]*pow(Box1Abb[68],2.)*m_z12*m_z1k + 2.*pow(Box1Abb[68],4.)*m_z1k_2;

  Box1Abb[4462]=Box1Abb[4461]*m_x_3 - Box1Abb[4438]*m_x_4 + Box1Abb[4413]*m_x_5 - Box1Abb[4376]*m_x_6 + Box1Abb[4362]*m_x_7 + Box1Abb[4351]*m_x_8 + Box1Abb[4345]*m_x_9 + Box1Abb[4342]*m_x_10*m_z12 + Box1Abb[4393]*m_x_2*m_z12*m_z1k + Box1Abb[4]*Box1Abb[4354]*pow(Box1Abb[68],2.)*m_x*m_z12_2*m_z1k_2 - pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12_3*m_z1k_3;

  Box1Abb[4463]=Box1Abb[4462]*m_cL + Box1Abb[4341]*pow(Box1Abb[7],2.)*m_cR*m_x;

  Box1Abb[4464]=35. + 26.*m_z12;

  Box1Abb[4465]=-50. + Box1Abb[4464]*m_z12;

  Box1Abb[4466]=72. - 5.*Box1Abb[158]*m_z12;

  Box1Abb[4467]=-20. + Box1Abb[4466]*m_z12;

  Box1Abb[4468]=-74. + 57.*m_z12 + 6.*m_z12_2;

  Box1Abb[4469]=12. + Box1Abb[4468]*m_z12;

  Box1Abb[4470]=-8. + 23.*m_z12;

  Box1Abb[4471]=12. + Box1Abb[4465]*m_z12 + 8.*m_z1k + Box1Abb[4467]*m_z12*m_z1k + Box1Abb[4469]*m_z1k_2 + 5.*Box1Abb[4470]*m_z12*m_z1k_3;

  Box1Abb[4472]=15. + m_z12 + 11.*m_z1k;

  Box1Abb[4473]=7. + Box1Abb[4472]*m_z12 + 20.*m_z1k;

  Box1Abb[4474]=2. + Box1Abb[4473]*m_z12;

  Box1Abb[4475]=14. - 13.*m_z1k;

  Box1Abb[4476]=3. + Box1Abb[4475]*m_z1k;

  Box1Abb[4477]=-4. + Box1Abb[4476]*m_z1k;

  Box1Abb[4478]=-pow(Box1Abb[68],3.)*Box1Abb[692] + Box1Abb[4477]*m_z12 + 2.*Box1Abb[2312]*m_z12_2*m_z1k;

  Box1Abb[4479]=27. - 10.*m_z1k;

  Box1Abb[4480]=76. + 63.*m_z1k;

  Box1Abb[4481]=47. + Box1Abb[4480]*m_z1k;

  Box1Abb[4482]=26. - Box1Abb[4481]*m_z12 + 5.*Box1Abb[536]*m_z12_2 + Box1Abb[4479]*m_z1k;

  Box1Abb[4483]=-8.*Box1Abb[15] + Box1Abb[4482]*m_z12;

  Box1Abb[4484]=41. + 13.*Box1Abb[174]*m_z1k;

  Box1Abb[4485]=14. + Box1Abb[4484]*m_z1k;

  Box1Abb[4486]=2. + 5.*Box1Abb[1776]*m_z1k;

  Box1Abb[4487]=-11. + Box1Abb[4486]*m_z1k;

  Box1Abb[4488]=104. - 85.*m_z1k;

  Box1Abb[4489]=84. + Box1Abb[4488]*m_z1k;

  Box1Abb[4490]=35. + Box1Abb[4489]*m_z1k;

  Box1Abb[4491]=10. + Box1Abb[4490]*m_z1k;

  Box1Abb[4492]=-8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4487]*Box1Abb[68]*m_z12 + Box1Abb[4491]*m_z12_2 - 2.*Box1Abb[4485]*m_z12_3 + Box1Abb[4056]*m_z12_4*m_z1k;

  Box1Abb[4493]=-9. + 10.*m_z1k;

  Box1Abb[4494]=-2. + Box1Abb[4493]*m_z1k;

  Box1Abb[4495]=45. + 19.*m_z1k;

  Box1Abb[4496]=-7. + Box1Abb[4495]*m_z1k;

  Box1Abb[4497]=4. + Box1Abb[4496]*Box1Abb[68]*m_z1k;

  Box1Abb[4498]=-85. + 46.*m_z1k;

  Box1Abb[4499]=8. + Box1Abb[4498]*m_z1k;

  Box1Abb[4500]=19. + Box1Abb[4499]*m_z1k;

  Box1Abb[4501]=-85. + 9.*m_z1k;

  Box1Abb[4502]=41. + Box1Abb[4501]*m_z1k;

  Box1Abb[4503]=19. + Box1Abb[4502]*m_z1k;

  Box1Abb[4504]=-2. + Box1Abb[4503]*m_z1k;

  Box1Abb[4505]=-Box1Abb[4494]*pow(Box1Abb[68],4.) - Box1Abb[4497]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4504]*Box1Abb[68]*m_z12_2 + Box1Abb[4500]*m_z12_3*m_z1k + 2.*Box1Abb[119]*m_z12_4*m_z1k_2;

  Box1Abb[4506]=21. + Box1Abb[3625]*m_z1k;

  Box1Abb[4507]=-3. + 11.*Box1Abb[3539]*m_z1k;

  Box1Abb[4508]=-13. + m_z1k;

  Box1Abb[4509]=147. + Box1Abb[4190]*Box1Abb[4508]*m_z1k;

  Box1Abb[4510]=-49. + Box1Abb[4509]*m_z1k;

  Box1Abb[4511]=-18. + Box1Abb[4510]*m_z1k;

  Box1Abb[4512]=-132. + 109.*m_z1k;

  Box1Abb[4513]=-94. + Box1Abb[4512]*m_z1k;

  Box1Abb[4514]=34. + Box1Abb[4513]*m_z1k;

  Box1Abb[4515]=13. + Box1Abb[4514]*m_z1k;

  Box1Abb[4516]=2.*pow(Box1Abb[68],4.) - Box1Abb[4507]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4511]*m_z12_2 + Box1Abb[4515]*m_z12_3 + Box1Abb[4506]*m_z12_4*m_z1k - 2.*m_z12_5*m_z1k_2;

  Box1Abb[4517]=Box1Abb[4516]*m_x_3 + Box1Abb[4492]*m_x_4 + Box1Abb[4471]*m_x_5 + Box1Abb[4483]*m_x_6 + Box1Abb[4474]*m_x_7 - Box1Abb[4505]*m_x_2*m_z12 + Box1Abb[542]*m_x_8*m_z12 + Box1Abb[4]*Box1Abb[4478]*Box1Abb[68]*m_x*m_z12_2*m_z1k + 2.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12_3*m_z1k_2;

  Box1Abb[4518]=10. + 24.*m_z12 - 15.*m_z1k;

  Box1Abb[4519]=5. + Box1Abb[4518]*m_z12 + 20.*m_z1k;

  Box1Abb[4520]=2. + Box1Abb[4519]*m_z12;

  Box1Abb[4521]=102. + 7.*m_z12;

  Box1Abb[4522]=-5. + Box1Abb[4521]*m_z12;

  Box1Abb[4523]=86. + 53.*m_z12;

  Box1Abb[4524]=10. + 33.*m_z12;

  Box1Abb[4525]=-35. + Box1Abb[4522]*m_z12 - 37.*m_z1k + Box1Abb[4523]*m_z12*m_z1k + Box1Abb[4524]*m_z1k_2;

  Box1Abb[4526]=8.*Box1Abb[15] + Box1Abb[4525]*m_z12;

  Box1Abb[4527]=2. - 24.*m_z1k_2 + 22.*m_z1k_3;

  Box1Abb[4528]=-2. + Box1Abb[2303]*m_z1k;

  Box1Abb[4529]=Box1Abb[4528]*pow(Box1Abb[68],2.) + Box1Abb[4527]*m_z12 + Box1Abb[3289]*m_z12_2*m_z1k;

  Box1Abb[4530]=15. + 7.*Box1Abb[178]*m_z1k;

  Box1Abb[4531]=-5. + Box1Abb[269]*m_z1k;

  Box1Abb[4532]=-13. + 7.*m_z1k;

  Box1Abb[4533]=-19. + 5.*Box1Abb[4532]*m_z1k;

  Box1Abb[4534]=1. + Box1Abb[4533]*m_z1k;

  Box1Abb[4535]=-12. + 25.*m_z1k;

  Box1Abb[4536]=-62. + 5.*Box1Abb[4535]*m_z1k;

  Box1Abb[4537]=-17. + Box1Abb[4536]*m_z1k;

  Box1Abb[4538]=17. + 89.*m_z1k;

  Box1Abb[4539]=-23. + Box1Abb[4538]*m_z1k;

  Box1Abb[4540]=-13. + Box1Abb[4539]*m_z1k;

  Box1Abb[4541]=-Box1Abb[4531]*pow(Box1Abb[68],4.) + Box1Abb[4534]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[4537]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[4540]*Box1Abb[68]*m_z12_3 + Box1Abb[4530]*m_z12_4*m_z1k;

  Box1Abb[4542]=28. + 27.*m_z1k;

  Box1Abb[4543]=29. + 11.*m_z1k;

  Box1Abb[4544]=21. + Box1Abb[4543]*m_z1k;

  Box1Abb[4545]=7. + 20.*m_z1k;

  Box1Abb[4546]=11. + Box1Abb[178]*Box1Abb[4545]*m_z1k;

  Box1Abb[4547]=94. + 145.*m_z1k;

  Box1Abb[4548]=-11. + Box1Abb[4547]*m_z1k;

  Box1Abb[4549]=-118. + Box1Abb[4548]*m_z1k;

  Box1Abb[4550]=4.*Box1Abb[823] - 2.*Box1Abb[4546]*m_z12 + Box1Abb[4549]*m_z12_2 + 7.*Box1Abb[4544]*m_z12_3 + Box1Abb[4542]*m_z12_4;

  Box1Abb[4551]=22. + 5.*m_z1k;

  Box1Abb[4552]=11. + Box1Abb[4551]*m_z1k;

  Box1Abb[4553]=1. + 5.*Box1Abb[662]*m_z1k;

  Box1Abb[4554]=-16. + Box1Abb[4553]*Box1Abb[68]*m_z1k;

  Box1Abb[4555]=45. + 43.*m_z1k;

  Box1Abb[4556]=127. + 4.*Box1Abb[4555]*m_z1k;

  Box1Abb[4557]=75. + Box1Abb[4556]*m_z1k;

  Box1Abb[4558]=124. + 39.*Box1Abb[711]*m_z1k;

  Box1Abb[4559]=96. + Box1Abb[4558]*m_z1k;

  Box1Abb[4560]=145. + Box1Abb[4559]*m_z1k;

  Box1Abb[4561]=-8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4554]*m_z12 + Box1Abb[4560]*m_z12_2 - Box1Abb[4557]*m_z12_3 - 4.*Box1Abb[4552]*m_z12_4 - 3.*m_z12_5*m_z1k;

  Box1Abb[4562]=11. + m_z1k;

  Box1Abb[4563]=23. + 15.*m_z1k;

  Box1Abb[4564]=20. + Box1Abb[4563]*m_z1k;

  Box1Abb[4565]=34. + 4.*Box1Abb[4564]*m_z1k;

  Box1Abb[4566]=-1. + 44.*m_z1k;

  Box1Abb[4567]=-34. + Box1Abb[4566]*m_z1k;

  Box1Abb[4568]=-25. + Box1Abb[4567]*m_z1k;

  Box1Abb[4569]=-67. + 113.*m_z1k;

  Box1Abb[4570]=-117. + 2.*Box1Abb[4569]*m_z1k;

  Box1Abb[4571]=-52. + Box1Abb[4570]*m_z1k;

  Box1Abb[4572]=-7. + Box1Abb[4571]*m_z1k;

  Box1Abb[4573]=-199. + 123.*m_z1k;

  Box1Abb[4574]=-99. + Box1Abb[4573]*m_z1k;

  Box1Abb[4575]=33. + Box1Abb[4574]*m_z1k;

  Box1Abb[4576]=52. + Box1Abb[4575]*m_z1k;

  Box1Abb[4577]=2.*pow(Box1Abb[68],4.) - Box1Abb[4568]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4576]*Box1Abb[68]*m_z12_2 + Box1Abb[4572]*m_z12_3 + Box1Abb[4565]*m_z12_4 + Box1Abb[4562]*m_z12_5*m_z1k;

  Box1Abb[4578]=Box1Abb[4577]*m_x_3 + Box1Abb[4561]*m_x_4 + Box1Abb[4550]*m_x_5 - Box1Abb[4526]*m_x_6 + Box1Abb[4520]*m_x_7 - Box1Abb[4541]*m_x_2*m_z12 + Box1Abb[3267]*m_x_8*m_z12 + Box1Abb[4]*Box1Abb[4529]*pow(Box1Abb[68],2.)*m_x*m_z12_2 - 2.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12_3*m_z1k;

  Box1Abb[4579]=Box1Abb[4578]*m_cL - Box1Abb[4517]*m_cR;

  Box1Abb[4580]=-3. + m_z12 + 3.*m_z1k;

  Box1Abb[4581]=17. + 12.*m_z12;

  Box1Abb[4582]=35. - Box1Abb[4581]*m_z12;

  Box1Abb[4583]=2. + Box1Abb[856]*m_z12_2;

  Box1Abb[4584]=-11. + 20.*m_z12;

  Box1Abb[4585]=-8. + Box1Abb[4584]*m_z12;

  Box1Abb[4586]=-39. + 2.*Box1Abb[4585]*m_z12;

  Box1Abb[4587]=4. + Box1Abb[4586]*m_z12;

  Box1Abb[4588]=78. + 29.*m_z12;

  Box1Abb[4589]=5. + Box1Abb[4588]*m_z12;

  Box1Abb[4590]=-4. + Box1Abb[4589]*m_z12;

  Box1Abb[4591]=14. - 39.*m_z12;

  Box1Abb[4592]=8.*Box1Abb[68] + Box1Abb[4582]*m_z12_2 + Box1Abb[1737]*Box1Abb[4583]*m_z12*m_z1k + 2.*Box1Abb[4587]*m_z1k_2 + 2.*Box1Abb[4590]*m_z1k_3 + 5.*Box1Abb[4591]*m_z12*m_z1k_4;

  Box1Abb[4593]=6. + 10.*m_z12 - 15.*m_z1k;

  Box1Abb[4594]=5. + Box1Abb[4593]*m_z12 + 20.*m_z1k;

  Box1Abb[4595]=2. + Box1Abb[4594]*m_z12;

  Box1Abb[4596]=1. + 3.*m_z12;

  Box1Abb[4597]=14. - 5.*m_z12;

  Box1Abb[4598]=-31. + Box1Abb[3496]*Box1Abb[4596]*m_z12 - 37.*m_z1k + 5.*Box1Abb[4597]*m_z12*m_z1k + Box1Abb[4524]*m_z1k_2;

  Box1Abb[4599]=8.*Box1Abb[15] + Box1Abb[4598]*m_z12;

  Box1Abb[4600]=10. - 9.*m_z1k;

  Box1Abb[4601]=1. - 5.*m_z1k;

  Box1Abb[4602]=7. + 3.*Box1Abb[4601]*m_z1k;

  Box1Abb[4603]=47. + 20.*m_z1k;

  Box1Abb[4604]=10. + Box1Abb[4603]*m_z1k;

  Box1Abb[4605]=15. + Box1Abb[4604]*m_z1k;

  Box1Abb[4606]=74. + 145.*m_z1k;

  Box1Abb[4607]=83. + Box1Abb[4606]*m_z1k;

  Box1Abb[4608]=-30. + Box1Abb[4607]*m_z1k;

  Box1Abb[4609]=4.*Box1Abb[823] - 2.*Box1Abb[4605]*m_z12 + Box1Abb[4608]*m_z12_2 + 7.*Box1Abb[4602]*m_z12_3 + Box1Abb[4600]*m_z12_4;

  Box1Abb[4610]=-1. + Box1Abb[269]*m_z1k;

  Box1Abb[4611]=-39. + 35.*m_z1k;

  Box1Abb[4612]=8. + Box1Abb[4611]*m_z1k;

  Box1Abb[4613]=-31. + 16.*m_z1k;

  Box1Abb[4614]=9. + Box1Abb[4613]*m_z1k;

  Box1Abb[4615]=5. + Box1Abb[4614]*m_z1k;

  Box1Abb[4616]=1. + Box1Abb[4615]*m_z1k;

  Box1Abb[4617]=Box1Abb[4610]*pow(Box1Abb[68],3.) - Box1Abb[4616]*m_z12_2 - Box1Abb[4612]*pow(Box1Abb[68],2.)*m_z12*m_z1k + Box1Abb[1832]*m_z12_3*m_z1k;

  Box1Abb[4618]=3. + m_z1k + 7.*m_z1k_2 - 34.*m_z1k_3;

  Box1Abb[4619]=-25. + 14.*m_z1k;

  Box1Abb[4620]=16. + Box1Abb[4619]*m_z1k;

  Box1Abb[4621]=-1. + Box1Abb[4620]*m_z1k;

  Box1Abb[4622]=-26. + Box1Abb[4566]*m_z1k;

  Box1Abb[4623]=-1. + Box1Abb[4622]*m_z1k;

  Box1Abb[4624]=-179. + 123.*m_z1k;

  Box1Abb[4625]=5. + Box1Abb[4624]*m_z1k;

  Box1Abb[4626]=59. + Box1Abb[4625]*m_z1k;

  Box1Abb[4627]=6. + Box1Abb[4626]*m_z1k;

  Box1Abb[4628]=2.*pow(Box1Abb[68],4.) - Box1Abb[4623]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4627]*Box1Abb[68]*m_z12_2 + Box1Abb[233]*Box1Abb[4621]*m_z12_3 + 2.*Box1Abb[4618]*m_z12_4 - Box1Abb[1843]*m_z12_5*m_z1k;

  Box1Abb[4629]=Box1Abb[4628]*m_x_2 + Box1Abb[4592]*m_x_3 + Box1Abb[4609]*m_x_4 - Box1Abb[4599]*m_x_5 + Box1Abb[4595]*m_x_6 + Box1Abb[4]*Box1Abb[4617]*m_x*m_z12 + Box1Abb[3267]*m_x_7*m_z12 + pow(Box1Abb[4],2.)*Box1Abb[4580]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k;

  Box1Abb[4630]=-Box1Abb[4517]*m_cR + Box1Abb[4629]*m_cL*m_x;

  Box1Abb[4631]=-6. + Box1Abb[0]*Box1Abb[1621]*m_z12;

  Box1Abb[4632]=17. + 2.*m_z12;

  Box1Abb[4633]=-22. + Box1Abb[4632]*m_z12;

  Box1Abb[4634]=-4. + Box1Abb[4633]*m_z12;

  Box1Abb[4635]=-8. + Box1Abb[4634]*m_z12;

  Box1Abb[4636]=11. + 18.*m_z12;

  Box1Abb[4637]=-19. + Box1Abb[4636]*m_z12;

  Box1Abb[4638]=6. + Box1Abb[4637]*m_z12;

  Box1Abb[4639]=8. - 9.*m_z12;

  Box1Abb[4640]=2.*Box1Abb[4631] + Box1Abb[4635]*m_z1k - 2.*Box1Abb[4638]*m_z1k_2 + 5.*Box1Abb[4639]*m_z12*m_z1k_3;

  Box1Abb[4641]=6. + m_z12 + 3.*m_z1k;

  Box1Abb[4642]=-17. + Box1Abb[4641]*m_z12 - 20.*m_z1k;

  Box1Abb[4643]=-2. + Box1Abb[4642]*m_z12;

  Box1Abb[4644]=-1. + Box1Abb[474]*m_z1k;

  Box1Abb[4645]=-2. + Box1Abb[3036]*m_z1k;

  Box1Abb[4646]=7. + 34.*m_z1k;

  Box1Abb[4647]=4. + Box1Abb[4646]*Box1Abb[68]*m_z1k;

  Box1Abb[4648]=-1. + 4.*Box1Abb[68]*m_z1k;

  Box1Abb[4649]=1. + 2.*Box1Abb[4648]*m_z1k;

  Box1Abb[4650]=Box1Abb[4645]*pow(Box1Abb[68],4.) + Box1Abb[4647]*pow(Box1Abb[68],2.)*m_z12 + 2.*Box1Abb[4649]*Box1Abb[68]*m_z12_2 + Box1Abb[4644]*m_z12_3*m_z1k;

  Box1Abb[4651]=4. + 7.*Box1Abb[119]*m_z1k;

  Box1Abb[4652]=9. + Box1Abb[4651]*m_z12 - Box1Abb[234]*m_z12_2 + Box1Abb[1900]*m_z1k;

  Box1Abb[4653]=8.*Box1Abb[15] + Box1Abb[4652]*m_z12;

  Box1Abb[4654]=2. + m_z1k + 39.*m_z1k_2 + 133.*m_z1k_3;

  Box1Abb[4655]=11. + 7.*Box1Abb[4601]*m_z1k;

  Box1Abb[4656]=21. + Box1Abb[4655]*m_z1k;

  Box1Abb[4657]=-8. + Box1Abb[4656]*m_z1k;

  Box1Abb[4658]=-44. + 15.*m_z1k;

  Box1Abb[4659]=-106. + Box1Abb[4658]*m_z1k;

  Box1Abb[4660]=-8. + Box1Abb[4659]*m_z1k;

  Box1Abb[4661]=9. + Box1Abb[4660]*m_z1k;

  Box1Abb[4662]=8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4657]*m_z12 + Box1Abb[4661]*m_z12_2 + Box1Abb[4654]*m_z12_3 + 4.*Box1Abb[68]*m_z12_4*m_z1k;

  Box1Abb[4663]=-2. + 3.*Box1Abb[1041]*m_z1k;

  Box1Abb[4664]=2. + Box1Abb[4663]*m_z1k;

  Box1Abb[4665]=-31. + 24.*m_z1k;

  Box1Abb[4666]=14. + Box1Abb[4665]*m_z1k;

  Box1Abb[4667]=4. + Box1Abb[4666]*m_z1k;

  Box1Abb[4668]=-1. + Box1Abb[4667]*m_z1k;

  Box1Abb[4669]=46. - 33.*m_z1k;

  Box1Abb[4670]=46. + Box1Abb[4669]*m_z1k;

  Box1Abb[4671]=-108. + Box1Abb[4670]*m_z1k;

  Box1Abb[4672]=51. + Box1Abb[4671]*m_z1k;

  Box1Abb[4673]=-2. + Box1Abb[4672]*m_z1k_2;

  Box1Abb[4674]=-Box1Abb[317]*Box1Abb[532]*pow(Box1Abb[68],4.) + Box1Abb[4673]*m_z12 - Box1Abb[15]*Box1Abb[4668]*m_z12_2 + 2.*Box1Abb[4664]*m_z12_3*m_z1k + Box1Abb[400]*m_z12_4*m_z1k_2;

  Box1Abb[4675]=2. + 29.*m_z1k;

  Box1Abb[4676]=22. - 95.*m_z1k;

  Box1Abb[4677]=49. + Box1Abb[4676]*m_z1k;

  Box1Abb[4678]=-9. + Box1Abb[4677]*m_z1k;

  Box1Abb[4679]=-3. + Box1Abb[4678]*m_z1k;

  Box1Abb[4680]=-41. + 44.*m_z1k;

  Box1Abb[4681]=3. + Box1Abb[4680]*m_z1k;

  Box1Abb[4682]=21. + Box1Abb[4681]*m_z1k;

  Box1Abb[4683]=-3. + Box1Abb[4682]*m_z1k;

  Box1Abb[4684]=34. + Box1Abb[4415]*m_z1k;

  Box1Abb[4685]=-44. + Box1Abb[4684]*m_z1k;

  Box1Abb[4686]=15. + Box1Abb[4685]*m_z1k;

  Box1Abb[4687]=2. + Box1Abb[4686]*m_z1k;

  Box1Abb[4688]=-2.*pow(Box1Abb[68],4.) + Box1Abb[4683]*Box1Abb[68]*m_z12 + Box1Abb[4687]*m_z12_2 + Box1Abb[4679]*m_z12_3 - 2.*Box1Abb[4675]*m_z12_4*m_z1k_2 + m_z12_5*m_z1k_2;

  Box1Abb[4689]=Box1Abb[4688]*m_x_3 + Box1Abb[4662]*m_x_4 + Box1Abb[4640]*m_x_5 + Box1Abb[4653]*m_x_6 + Box1Abb[4643]*m_x_7 + Box1Abb[4674]*m_x_2*m_z12 - 3.*Box1Abb[79]*m_x_8*m_z12 + Box1Abb[4650]*m_x*m_z12_2*m_z1k - pow(Box1Abb[4],2.)*Box1Abb[457]*pow(Box1Abb[68],2.)*m_z12_3*m_z1k_2;

  Box1Abb[4690]=7. + Box1Abb[167]*m_z12 + Box1Abb[396]*m_z1k;

  Box1Abb[4691]=-2. + Box1Abb[4690]*m_z1k;

  Box1Abb[4692]=17. + m_z12 + 29.*m_z1k;

  Box1Abb[4693]=-5.*Box1Abb[233] + Box1Abb[4692]*m_z12;

  Box1Abb[4694]=-2. + Box1Abb[4693]*m_z12;

  Box1Abb[4695]=-2. + 57.*m_z1k;

  Box1Abb[4696]=-1. + Box1Abb[2058]*m_z1k;

  Box1Abb[4697]=2. + Box1Abb[4696]*m_z12 + Box1Abb[4695]*m_z12_2 + Box1Abb[3865]*m_z1k;

  Box1Abb[4698]=-8.*Box1Abb[15] + Box1Abb[4697]*m_z12;

  Box1Abb[4699]=13. + Box1Abb[174]*Box1Abb[2440]*m_z1k;

  Box1Abb[4700]=51. + 62.*m_z1k;

  Box1Abb[4701]=-7. + Box1Abb[4700]*m_z1k;

  Box1Abb[4702]=49. + 75.*m_z1k;

  Box1Abb[4703]=102. + Box1Abb[4702]*m_z1k;

  Box1Abb[4704]=9. + Box1Abb[4703]*m_z1k;

  Box1Abb[4705]=4.*Box1Abb[823] - 2.*Box1Abb[4699]*m_z12 + Box1Abb[4704]*m_z12_2 - 2.*Box1Abb[4701]*m_z12_3 - 13.*m_z12_4*m_z1k;

  Box1Abb[4706]=5. + 2.*m_z1k;

  Box1Abb[4707]=19. - 30.*m_z1k;

  Box1Abb[4708]=12. + Box1Abb[4707]*m_z1k;

  Box1Abb[4709]=-5. + Box1Abb[4708]*m_z1k;

  Box1Abb[4710]=-50. + 21.*m_z1k;

  Box1Abb[4711]=14. + Box1Abb[4710]*m_z1k;

  Box1Abb[4712]=3. + Box1Abb[4711]*m_z1k;

  Box1Abb[4713]=-4. + Box1Abb[4712]*m_z1k;

  Box1Abb[4714]=-55. + 59.*m_z1k;

  Box1Abb[4715]=-19. + Box1Abb[4714]*m_z1k;

  Box1Abb[4716]=7. + Box1Abb[4715]*m_z1k;

  Box1Abb[4717]=-2. + Box1Abb[4716]*m_z1k;

  Box1Abb[4718]=Box1Abb[317]*Box1Abb[451]*pow(Box1Abb[68],4.) - Box1Abb[4713]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[4717]*Box1Abb[68]*m_z12_2 + Box1Abb[4709]*m_z12_3*m_z1k + 2.*Box1Abb[4706]*m_z12_4*m_z1k_2;

  Box1Abb[4719]=25. + 62.*m_z1k;

  Box1Abb[4720]=32. + 17.*m_z1k;

  Box1Abb[4721]=16. + Box1Abb[4720]*m_z1k;

  Box1Abb[4722]=-20. + 2.*Box1Abb[4721]*m_z1k;

  Box1Abb[4723]=3. + 35.*m_z1k;

  Box1Abb[4724]=-16. + Box1Abb[4723]*m_z1k;

  Box1Abb[4725]=3. + Box1Abb[4724]*m_z1k;

  Box1Abb[4726]=7. + Box1Abb[4725]*m_z1k;

  Box1Abb[4727]=-56. + 125.*m_z1k;

  Box1Abb[4728]=8. + Box1Abb[4727]*m_z1k;

  Box1Abb[4729]=87. + Box1Abb[4728]*m_z1k;

  Box1Abb[4730]=12. - Box1Abb[4729]*m_z1k;

  Box1Abb[4731]=-8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4726]*m_z12 + Box1Abb[4730]*m_z12_2 + Box1Abb[4722]*m_z12_3 + Box1Abb[4719]*m_z12_4*m_z1k;

  Box1Abb[4732]=57. + 28.*m_z1k;

  Box1Abb[4733]=9. + Box1Abb[4732]*m_z1k;

  Box1Abb[4734]=13. + 44.*m_z1k;

  Box1Abb[4735]=2. + Box1Abb[4734]*m_z1k;

  Box1Abb[4736]=-3. + Box1Abb[4735]*m_z1k;

  Box1Abb[4737]=8. - 48.*m_z1k + 87.*m_z1k_2;

  Box1Abb[4738]=8. + Box1Abb[4737]*m_z1k;

  Box1Abb[4739]=11. + Box1Abb[4738]*m_z1k;

  Box1Abb[4740]=-167. + 81.*m_z1k;

  Box1Abb[4741]=14. + Box1Abb[4740]*m_z1k;

  Box1Abb[4742]=45. + Box1Abb[4741]*m_z1k;

  Box1Abb[4743]=23. + Box1Abb[4742]*m_z1k;

  Box1Abb[4744]=-16. + Box1Abb[4743]*m_z1k;

  Box1Abb[4745]=2.*pow(Box1Abb[68],4.) - Box1Abb[4736]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4744]*m_z12_2 + Box1Abb[4739]*m_z12_3 - Box1Abb[4733]*m_z12_4*m_z1k - 8.*m_z12_5*m_z1k_2;

  Box1Abb[4746]=Box1Abb[4745]*m_x_2 + Box1Abb[4731]*m_x_3 + Box1Abb[4705]*m_x_4 + Box1Abb[4698]*m_x_5 - Box1Abb[4694]*m_x_6 + Box1Abb[4718]*m_x*m_z12 + 3.*Box1Abb[133]*m_x_7*m_z12 + pow(Box1Abb[4],2.)*Box1Abb[4691]*Box1Abb[68]*m_z12_2*m_z1k;

  Box1Abb[4747]=Box1Abb[4689]*m_cR + Box1Abb[4746]*m_cL*m_x;

  Box1Abb[4748]=-4. + 2.*Box1Abb[175]*m_z12 + 5.*m_z12_2;

  Box1Abb[4749]=-9. + 2.*Box1Abb[2241]*m_z12;

  Box1Abb[4750]=-2.*Box1Abb[2006] + 3.*pow(Box1Abb[77],2.)*m_z12 + 2.*Box1Abb[14]*Box1Abb[2661]*m_z12*m_z1k + 2.*Box1Abb[4749]*m_z1k_2 + 20.*m_z12*m_z1k_3;

  Box1Abb[4751]=-12. + m_z1k;

  Box1Abb[4752]=1. + Box1Abb[4751]*m_z1k;

  Box1Abb[4753]=-1. + 5.*Box1Abb[15]*m_z1k;

  Box1Abb[4754]=2.*pow(Box1Abb[68],3.) - 2.*Box1Abb[4753]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[4752]*Box1Abb[68]*m_z12_2 + Box1Abb[4162]*m_z12_3;

  Box1Abb[4755]=23. + 5.*m_z1k;

  Box1Abb[4756]=20. + Box1Abb[4137]*m_z12 + m_z12_2 + 2.*Box1Abb[4755]*m_z1k;

  Box1Abb[4757]=-2.*Box1Abb[662] + Box1Abb[4756]*m_z12;

  Box1Abb[4758]=-3. + 2.*Box1Abb[228]*m_z1k_2;

  Box1Abb[4759]=13. - 2.*m_z1k;

  Box1Abb[4760]=6. + Box1Abb[4759]*m_z1k;

  Box1Abb[4761]=2. + 2.*Box1Abb[4760]*m_z1k;

  Box1Abb[4762]=-2.*Box1Abb[1151]*pow(Box1Abb[68],2.) + 2.*Box1Abb[4758]*Box1Abb[68]*m_z12 + Box1Abb[4761]*m_z12_2 + 3.*Box1Abb[15]*m_z12_3;

  Box1Abb[4763]=-Box1Abb[4754]*Box1Abb[68]*m_x - Box1Abb[4762]*m_x_2 + Box1Abb[4750]*m_x_3 - Box1Abb[4757]*m_x_4 + Box1Abb[4748]*m_x_5 + Box1Abb[3234]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[4764]=2.*pow(Box1Abb[68],2.) + 7.*Box1Abb[68]*m_z12 + 4.*m_z12_2;

  Box1Abb[4765]=4. + m_z12 + 2.*m_z1k;

  Box1Abb[4766]=-4. + Box1Abb[4765]*m_z12;

  Box1Abb[4767]=5. + 3.*m_z1k + 9.*m_z1k_2;

  Box1Abb[4768]=1. + m_z1k + 4.*m_z1k_2 + 5.*m_z1k_3;

  Box1Abb[4769]=3. + 2.*Box1Abb[4706]*m_z1k;

  Box1Abb[4770]=-2.*Box1Abb[4767] + 4.*Box1Abb[4768]*m_z12 + Box1Abb[4769]*m_z12_2;

  Box1Abb[4771]=8. + Box1Abb[210]*m_z12 + 2.*Box1Abb[228]*m_z1k;

  Box1Abb[4772]=-2.*Box1Abb[1776] + Box1Abb[4771]*m_z12;

  Box1Abb[4773]=9. + 19.*m_z1k;

  Box1Abb[4774]=-6. + Box1Abb[4773]*m_z1k;

  Box1Abb[4775]=-2.*pow(Box1Abb[68],3.) + 10.*Box1Abb[15]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[4774]*Box1Abb[68]*m_z12_2 + 4.*Box1Abb[178]*m_z12_3*m_z1k;

  Box1Abb[4776]=-4. + Box1Abb[1039]*m_z1k;

  Box1Abb[4777]=-4. + Box1Abb[4776]*m_z1k;

  Box1Abb[4778]=9. + 8.*m_z1k;

  Box1Abb[4779]=15. + 2.*Box1Abb[4778]*m_z1k;

  Box1Abb[4780]=1. + Box1Abb[4779]*m_z1k;

  Box1Abb[4781]=-2. + Box1Abb[4780]*m_z12_2 + 4.*Box1Abb[4777]*m_z12*m_z1k + 2.*Box1Abb[3776]*m_z1k_2;

  Box1Abb[4782]=Box1Abb[4781]*m_x_2 - Box1Abb[4770]*m_x_3 + Box1Abb[4772]*m_x_4 - Box1Abb[4766]*m_x_5 - Box1Abb[4775]*m_x*m_z1k + Box1Abb[4764]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[4783]=Box1Abb[4763]*m_cL + Box1Abb[4782]*m_cR;

  Box1Abb[4784]=pow(Box1Abb[68],3.) + 6.*pow(Box1Abb[68],2.)*m_z12 + 6.*Box1Abb[68]*m_z12_2 + 2.*m_z12_3;

  Box1Abb[4785]=39. - 19.*m_z12;

  Box1Abb[4786]=15. + Box1Abb[4785]*m_z12;

  Box1Abb[4787]=50. + Box1Abb[4786]*m_z12;

  Box1Abb[4788]=30. - 37.*m_z12;

  Box1Abb[4789]=18. + Box1Abb[4788]*m_z12;

  Box1Abb[4790]=50. + Box1Abb[4789]*m_z12;

  Box1Abb[4791]=40. + 2.*Box1Abb[12]*Box1Abb[2148]*m_z12 + Box1Abb[4787]*m_z1k + Box1Abb[4790]*m_z1k_2 - 14.*Box1Abb[1621]*m_z12*m_z1k_3;

  Box1Abb[4792]=12. - 5.*m_z12;

  Box1Abb[4793]=-45. + 4.*Box1Abb[4792]*m_z12;

  Box1Abb[4794]=2. + Box1Abb[4793]*m_z12;

  Box1Abb[4795]=25. + Box1Abb[3353]*m_z12;

  Box1Abb[4796]=7. + Box1Abb[4795]*m_z12;

  Box1Abb[4797]=8. + m_z12 + m_z12_2;

  Box1Abb[4798]=-27. + Box1Abb[4797]*m_z12;

  Box1Abb[4799]=16. + Box1Abb[4798]*m_z12;

  Box1Abb[4800]=-10. + Box1Abb[4799]*m_z12;

  Box1Abb[4801]=43. - 12.*m_z12;

  Box1Abb[4802]=66. + Box1Abb[4801]*m_z12;

  Box1Abb[4803]=-60. + Box1Abb[4802]*m_z12;

  Box1Abb[4804]=40. + Box1Abb[4803]*m_z12;

  Box1Abb[4805]=-9. + 14.*m_z12;

  Box1Abb[4806]=9. + 4.*Box1Abb[4805]*m_z12;

  Box1Abb[4807]=14. + Box1Abb[4806]*m_z12;

  Box1Abb[4808]=-2. + 2.*m_z12_2 - m_z12_3 + Box1Abb[4794]*m_z1k - 2.*Box1Abb[0]*Box1Abb[4796]*m_z1k_2 + 4.*Box1Abb[4800]*m_z1k_3 + Box1Abb[4804]*m_z1k_4 - Box1Abb[4807]*m_z1k_5 - 14.*Box1Abb[0]*m_z12*m_z1k_6;

  Box1Abb[4809]=-6. + m_z12 - 11.*m_z1k;

  Box1Abb[4810]=14. + Box1Abb[4809]*m_z12 + 26.*m_z1k;

  Box1Abb[4811]=4. + Box1Abb[4810]*m_z12;

  Box1Abb[4812]=-35. + 11.*m_z12;

  Box1Abb[4813]=-10. + Box1Abb[1051]*m_z12 - 43.*m_z1k + 8.*Box1Abb[12]*m_z12*m_z1k + 2.*Box1Abb[4812]*m_z1k_2;

  Box1Abb[4814]=-20. + Box1Abb[4813]*m_z12 - 22.*m_z1k;

  Box1Abb[4815]=-3. + Box1Abb[873]*m_z1k;

  Box1Abb[4816]=Box1Abb[437]*pow(Box1Abb[68],4.) + 2.*Box1Abb[210]*Box1Abb[437]*pow(Box1Abb[68],3.)*m_z12 + 6.*Box1Abb[4815]*pow(Box1Abb[68],2.)*m_z12_2 + 2.*Box1Abb[1163]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[318]*m_z12_4*m_z1k;

  Box1Abb[4817]=7. + 9.*m_z1k;

  Box1Abb[4818]=1. + Box1Abb[4817]*m_z1k;

  Box1Abb[4819]=-5. + 2.*Box1Abb[4818]*m_z1k;

  Box1Abb[4820]=33. + 7.*m_z1k;

  Box1Abb[4821]=27. + Box1Abb[4820]*m_z1k;

  Box1Abb[4822]=20. + Box1Abb[4821]*m_z1k;

  Box1Abb[4823]=-3. + Box1Abb[4822]*m_z1k;

  Box1Abb[4824]=-13. + 14.*m_z1k;

  Box1Abb[4825]=9. + Box1Abb[4824]*m_z1k;

  Box1Abb[4826]=14. + 5.*Box1Abb[4825]*m_z1k;

  Box1Abb[4827]=22. - Box1Abb[4826]*m_z1k;

  Box1Abb[4828]=-38. + Box1Abb[4827]*m_z12 - 2.*Box1Abb[4823]*m_z12_2 + 2.*Box1Abb[4819]*m_z12_3 - 20.*Box1Abb[461]*m_z1k + 12.*m_z12_4*m_z1k_2;

  Box1Abb[4829]=-5. - 8.*m_z1k + 6.*m_z1k_2;

  Box1Abb[4830]=6. + Box1Abb[4829]*m_z1k;

  Box1Abb[4831]=-9. + 2.*Box1Abb[267]*m_z1k;

  Box1Abb[4832]=15. + Box1Abb[4831]*m_z1k;

  Box1Abb[4833]=-3. + 2.*Box1Abb[4562]*m_z1k;

  Box1Abb[4834]=-60. + Box1Abb[4833]*m_z1k;

  Box1Abb[4835]=15. + Box1Abb[4834]*m_z1k;

  Box1Abb[4836]=-50. + 17.*m_z1k;

  Box1Abb[4837]=-36. + Box1Abb[4836]*m_z1k;

  Box1Abb[4838]=62. + Box1Abb[4837]*m_z1k;

  Box1Abb[4839]=-5. + Box1Abb[4838]*m_z1k;

  Box1Abb[4840]=2.*pow(Box1Abb[68],5.) - Box1Abb[4832]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[4835]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[4839]*Box1Abb[68]*m_z12_3 + 4.*Box1Abb[4830]*m_z12_4*m_z1k + 2.*m_z12_5*m_z1k_3;

  Box1Abb[4841]=8. + 25.*m_z1k;

  Box1Abb[4842]=36. + Box1Abb[4841]*m_z1k;

  Box1Abb[4843]=22. + Box1Abb[4842]*m_z1k;

  Box1Abb[4844]=5. + Box1Abb[4843]*m_z1k;

  Box1Abb[4845]=55. + 7.*Box1Abb[2406]*m_z1k;

  Box1Abb[4846]=11. + Box1Abb[4845]*m_z1k;

  Box1Abb[4847]=23. + Box1Abb[4846]*m_z1k;

  Box1Abb[4848]=-4. + Box1Abb[4847]*m_z1k;

  Box1Abb[4849]=-69. + 4.*Box1Abb[1943]*m_z1k;

  Box1Abb[4850]=8. + Box1Abb[4849]*m_z1k;

  Box1Abb[4851]=-29. + Box1Abb[4850]*m_z1k;

  Box1Abb[4852]=-6. + Box1Abb[4851]*m_z1k;

  Box1Abb[4853]=-4.*Box1Abb[1446] + 2.*Box1Abb[4848]*m_z12 + Box1Abb[4852]*m_z12_2 + Box1Abb[4844]*m_z12_3 - 6.*Box1Abb[119]*m_z12_4*m_z1k_2 + 40.*Box1Abb[68]*m_z1k_3;

  Box1Abb[4854]=Box1Abb[4808]*m_x_3 + Box1Abb[4853]*m_x_4 + Box1Abb[4828]*m_x_5 + Box1Abb[4791]*m_x_6 + Box1Abb[4814]*m_x_7 + Box1Abb[4811]*m_x_8 + 2.*Box1Abb[79]*m_x_9*m_z12 + Box1Abb[4840]*m_x_2*m_z1k + Box1Abb[4816]*Box1Abb[68]*m_x*m_z12*m_z1k_2 - Box1Abb[4784]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_3;

  Box1Abb[4855]=78. - 5.*Box1Abb[677]*m_z12;

  Box1Abb[4856]=12. + Box1Abb[4855]*m_z12;

  Box1Abb[4857]=-6. + Box1Abb[803]*m_z12;

  Box1Abb[4858]=-42. - 5.*Box1Abb[4857]*m_z12;

  Box1Abb[4859]=-3. + 77.*m_z12;

  Box1Abb[4860]=30. + Box1Abb[4859]*m_z12;

  Box1Abb[4861]=36. + Box1Abb[4860]*m_z12;

  Box1Abb[4862]=30. + Box1Abb[4861]*m_z12;

  Box1Abb[4863]=-2. + 15.*m_z12;

  Box1Abb[4864]=192. + 13.*Box1Abb[4863]*m_z12;

  Box1Abb[4865]=-424. + Box1Abb[4864]*m_z12;

  Box1Abb[4866]=76. + Box1Abb[4865]*m_z12;

  Box1Abb[4867]=32. + 5.*Box1Abb[1602]*m_z12;

  Box1Abb[4868]=-8. + Box1Abb[4867]*m_z12;

  Box1Abb[4869]=-1. + 9.*m_z12;

  Box1Abb[4870]=-46. + Box1Abb[4856]*m_z12 + 20.*m_z1k + 2.*Box1Abb[4858]*m_z12*m_z1k + Box1Abb[4862]*m_z1k_2 + Box1Abb[4866]*m_z1k_3 + 10.*Box1Abb[4868]*m_z1k_4 + 28.*Box1Abb[4869]*m_z12*m_z1k_5;

  Box1Abb[4871]=16. + Box1Abb[2330]*m_z12 + 52.*m_z1k;

  Box1Abb[4872]=8. + Box1Abb[4871]*m_z12;

  Box1Abb[4873]=-69. + 2.*Box1Abb[2351]*m_z12;

  Box1Abb[4874]=5. - 9.*m_z12;

  Box1Abb[4875]=-44. + Box1Abb[4873]*m_z12 + 14.*m_z1k - 10.*Box1Abb[708]*m_z12*m_z1k + 28.*Box1Abb[4874]*m_z1k_2;

  Box1Abb[4876]=44.*Box1Abb[15] + Box1Abb[4875]*m_z12;

  Box1Abb[4877]=-1. + m_z1k + 3.*m_z1k_2;

  Box1Abb[4878]=17. + 10.*m_z1k;

  Box1Abb[4879]=1. - Box1Abb[4878]*m_z1k;

  Box1Abb[4880]=2.*Box1Abb[317]*pow(Box1Abb[68],3.) + 3.*Box1Abb[174]*pow(Box1Abb[68],2.)*m_z12 - 6.*Box1Abb[4877]*Box1Abb[68]*m_z12_2 + Box1Abb[4879]*m_z12_3;

  Box1Abb[4881]=-85. + 30.*m_z1k + 212.*m_z1k_2;

  Box1Abb[4882]=57. + 50.*m_z1k;

  Box1Abb[4883]=50. + Box1Abb[4882]*m_z1k;

  Box1Abb[4884]=36. - 49.*m_z1k;

  Box1Abb[4885]=37. + Box1Abb[4884]*m_z1k;

  Box1Abb[4886]=57. + 2.*Box1Abb[4885]*m_z1k;

  Box1Abb[4887]=6. + 35.*m_z1k;

  Box1Abb[4888]=5. + Box1Abb[4887]*m_z1k;

  Box1Abb[4889]=43. + 4.*Box1Abb[4888]*m_z1k;

  Box1Abb[4890]=-2.*Box1Abb[4883] + 2.*Box1Abb[4886]*m_z12 + 3.*Box1Abb[4889]*m_z12_2 + Box1Abb[4881]*m_z12_3 - Box1Abb[2358]*m_z12_4;

  Box1Abb[4891]=4. + 5.*Box1Abb[1723]*m_z1k;

  Box1Abb[4892]=25. + 4.*Box1Abb[511]*m_z1k;

  Box1Abb[4893]=107. + 261.*m_z1k;

  Box1Abb[4894]=41. + Box1Abb[4893]*m_z1k;

  Box1Abb[4895]=-30. + Box1Abb[4894]*m_z1k;

  Box1Abb[4896]=-9. + 20.*Box1Abb[515]*m_z1k;

  Box1Abb[4897]=27. + Box1Abb[4896]*m_z1k;

  Box1Abb[4898]=50. + Box1Abb[4897]*m_z1k;

  Box1Abb[4899]=-37. + 14.*m_z1k;

  Box1Abb[4900]=8. + 5.*Box1Abb[4899]*m_z1k;

  Box1Abb[4901]=-29. + Box1Abb[4900]*m_z1k;

  Box1Abb[4902]=-36. + Box1Abb[4901]*m_z1k;

  Box1Abb[4903]=98. + 2.*Box1Abb[4902]*m_z12 - 3.*Box1Abb[4898]*m_z12_2 - 2.*Box1Abb[4895]*m_z12_3 + 5.*Box1Abb[4891]*m_z12_4 + 2.*Box1Abb[4892]*m_z1k;

  Box1Abb[4904]=1. + 2.*Box1Abb[240]*Box1Abb[68]*m_z1k;

  Box1Abb[4905]=69. - 46.*m_z1k + 20.*m_z1k_2;

  Box1Abb[4906]=8. + Box1Abb[4905]*m_z1k;

  Box1Abb[4907]=5. + Box1Abb[4906]*m_z1k;

  Box1Abb[4908]=9. + Box1Abb[210]*m_z1k;

  Box1Abb[4909]=-9. + 4.*Box1Abb[4908]*m_z1k;

  Box1Abb[4910]=1. + Box1Abb[4909]*m_z1k;

  Box1Abb[4911]=-90. + 59.*m_z1k;

  Box1Abb[4912]=8. + Box1Abb[4911]*m_z1k;

  Box1Abb[4913]=20. + Box1Abb[4912]*m_z1k;

  Box1Abb[4914]=5. + Box1Abb[4913]*m_z1k;

  Box1Abb[4915]=-2. + Box1Abb[4914]*m_z1k;

  Box1Abb[4916]=-2.*Box1Abb[317]*pow(Box1Abb[68],5.) + 2.*Box1Abb[4904]*pow(Box1Abb[68],4.)*m_z12 + 3.*Box1Abb[4910]*pow(Box1Abb[68],3.)*m_z12_2 + Box1Abb[4907]*pow(Box1Abb[68],2.)*m_z12_3 + Box1Abb[4915]*m_z12_4 + 24.*m_z12_5*m_z1k_4;

  Box1Abb[4917]=55. + 14.*m_z1k;

  Box1Abb[4918]=-82. + Box1Abb[4917]*m_z1k;

  Box1Abb[4919]=-37. + Box1Abb[4918]*m_z1k;

  Box1Abb[4920]=2. + Box1Abb[4919]*m_z1k;

  Box1Abb[4921]=2. + 101.*m_z1k;

  Box1Abb[4922]=15. + Box1Abb[4921]*m_z1k;

  Box1Abb[4923]=5. + Box1Abb[4922]*m_z1k;

  Box1Abb[4924]=-5. + Box1Abb[4923]*m_z1k;

  Box1Abb[4925]=-5. + 7.*m_z1k;

  Box1Abb[4926]=81. + 4.*Box1Abb[4925]*m_z1k;

  Box1Abb[4927]=-13. + Box1Abb[4926]*m_z1k;

  Box1Abb[4928]=-17. + Box1Abb[4927]*m_z1k;

  Box1Abb[4929]=-3. + Box1Abb[4928]*m_z1k;

  Box1Abb[4930]=-147. + 127.*m_z1k;

  Box1Abb[4931]=22. + Box1Abb[4930]*m_z1k;

  Box1Abb[4932]=16. + Box1Abb[4931]*m_z1k;

  Box1Abb[4933]=15. + Box1Abb[4932]*m_z1k;

  Box1Abb[4934]=7. + Box1Abb[4933]*m_z1k;

  Box1Abb[4935]=-14.*Box1Abb[2069]*pow(Box1Abb[68],3.) + 2.*Box1Abb[4920]*pow(Box1Abb[68],2.)*m_z12 + 3.*Box1Abb[4929]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[4934]*m_z12_3 + 2.*Box1Abb[4924]*m_z12_4 + 24.*m_z12_5*m_z1k_3;

  Box1Abb[4936]=Box1Abb[4916]*m_x_2 - Box1Abb[4935]*m_x_3 + Box1Abb[4870]*m_x_4 + Box1Abb[4903]*m_x_5 + Box1Abb[4890]*m_x_6 + Box1Abb[4876]*m_x_7 - Box1Abb[4872]*m_x_8 + 4.*Box1Abb[1010]*m_x_9*m_z12 - Box1Abb[4880]*pow(Box1Abb[68],3.)*m_x*m_z12*m_z1k - 3.*Box1Abb[4]*pow(Box1Abb[68],5.)*m_z12_3*m_z1k_2;

  Box1Abb[4937]=Box1Abb[4936]*m_cL + 2.*Box1Abb[4854]*m_cR;

  Box1Abb[4938]=-2. + 3.*m_z12 + 4.*m_z1k;

  Box1Abb[4939]=3. + 2.*Box1Abb[4938]*m_z12;

  Box1Abb[4940]=1. + 8.*m_z1k;

  Box1Abb[4941]=1. - 2.*m_z1k + 8.*m_z1k_3 - 7.*m_z1k_4;

  Box1Abb[4942]=1. + Box1Abb[757]*m_z1k_2;

  Box1Abb[4943]=-Box1Abb[536]*pow(Box1Abb[68],5.) - 4.*Box1Abb[317]*pow(Box1Abb[68],6.)*m_z12 - 3.*Box1Abb[4940]*pow(Box1Abb[68],5.)*m_z12_2 - 10.*Box1Abb[4942]*pow(Box1Abb[68],2.)*m_z12_3 + 5.*Box1Abb[4941]*m_z12_4 - 12.*m_z12_5*m_z1k_3;

  Box1Abb[4944]=6. + 13.*m_z1k;

  Box1Abb[4945]=-26. + 3.*Box1Abb[4944]*m_z12 - 5.*m_z12_2 + 4.*Box1Abb[2125]*m_z1k;

  Box1Abb[4946]=13. + Box1Abb[4945]*m_z12 + 18.*m_z1k;

  Box1Abb[4947]=22. - 8.*m_z1k;

  Box1Abb[4948]=4. + Box1Abb[4817]*m_z1k;

  Box1Abb[4949]=13. + 36.*m_z1k;

  Box1Abb[4950]=11. + Box1Abb[4949]*m_z1k;

  Box1Abb[4951]=27. - 14.*m_z1k;

  Box1Abb[4952]=18. + Box1Abb[4951]*m_z1k;

  Box1Abb[4953]=9. + Box1Abb[4952]*m_z1k;

  Box1Abb[4954]=-5.*Box1Abb[4948] + 4.*Box1Abb[4953]*m_z12 - 3.*Box1Abb[4950]*m_z12_2 + Box1Abb[4947]*m_z12_3 + m_z12_4;

  Box1Abb[4955]=1. + m_z1k + m_z1k_2 + 6.*m_z1k_3;

  Box1Abb[4956]=-20. + 7.*m_z1k;

  Box1Abb[4957]=-3. + Box1Abb[4956]*m_z1k_3;

  Box1Abb[4958]=22. + 85.*m_z1k;

  Box1Abb[4959]=-35. + Box1Abb[4958]*m_z1k;

  Box1Abb[4960]=-3. + Box1Abb[4261]*m_z1k;

  Box1Abb[4961]=19. + 5.*Box1Abb[4960]*m_z1k;

  Box1Abb[4962]=10.*Box1Abb[4955] + 10.*Box1Abb[4957]*m_z12 + 3.*Box1Abb[4961]*m_z12_2 + Box1Abb[4959]*m_z12_3 - 5.*Box1Abb[15]*m_z12_4;

  Box1Abb[4963]=-10. + 7.*m_z1k;

  Box1Abb[4964]=-9. + 15.*m_z1k + Box1Abb[1899]*Box1Abb[4963]*m_z1k_3;

  Box1Abb[4965]=-5. + Box1Abb[1330]*m_z1k;

  Box1Abb[4966]=-1. + Box1Abb[1669]*m_z1k;

  Box1Abb[4967]=9. + 40.*m_z1k;

  Box1Abb[4968]=9. + Box1Abb[4967]*m_z1k;

  Box1Abb[4969]=-5. + Box1Abb[4968]*m_z1k;

  Box1Abb[4970]=5. + 4.*Box1Abb[4964]*m_z12 + 30.*Box1Abb[1039]*Box1Abb[4966]*m_z12_2 + 4.*Box1Abb[4969]*m_z12_3 + 2.*Box1Abb[4965]*m_z12_4 + 5.*Box1Abb[2065]*m_z1k_3;

  Box1Abb[4971]=-13. + 2.*Box1Abb[2125]*m_z1k;

  Box1Abb[4972]=-5. + 2.*Box1Abb[3833]*m_z1k_2;

  Box1Abb[4973]=-1. + 18.*m_z1k;

  Box1Abb[4974]=-7. + Box1Abb[4973]*m_z1k;

  Box1Abb[4975]=-18. + 5.*Box1Abb[267]*m_z1k;

  Box1Abb[4976]=4. + Box1Abb[4975]*m_z1k;

  Box1Abb[4977]=1. + Box1Abb[4976]*m_z1k;

  Box1Abb[4978]=-53. + 27.*m_z1k;

  Box1Abb[4979]=-3. + Box1Abb[4978]*m_z1k;

  Box1Abb[4980]=17. + Box1Abb[4979]*m_z1k;

  Box1Abb[4981]=-8. + Box1Abb[4980]*m_z1k;

  Box1Abb[4982]=Box1Abb[4974]*pow(Box1Abb[68],3.) + 2.*Box1Abb[4971]*pow(Box1Abb[68],4.)*m_z12 + 3.*Box1Abb[4981]*Box1Abb[68]*m_z12_2 + 5.*Box1Abb[4977]*m_z12_3 + 2.*Box1Abb[4972]*m_z12_4;

  Box1Abb[4983]=Box1Abb[4943]*m_x + Box1Abb[4982]*m_x_2 - Box1Abb[4970]*m_x_3 + Box1Abb[4962]*m_x_4 + Box1Abb[4954]*m_x_5 + Box1Abb[4946]*m_x_6 - Box1Abb[4939]*m_x_7 + pow(Box1Abb[4],3.)*pow(Box1Abb[68],5.)*m_z12 + m_x_8*m_z12;

  Box1Abb[4984]=3. + m_z12 + 2.*m_z12_2 + 9.*m_z12*m_z1k;

  Box1Abb[4985]=-1. + m_z1k + 9.*m_z1k_2;

  Box1Abb[4986]=3.*Box1Abb[692] + 4.*Box1Abb[4985]*m_z12 + 2.*Box1Abb[1692]*m_z12_2 - m_z12_3;

  Box1Abb[4987]=-18. + 13.*Box1Abb[15]*m_z1k;

  Box1Abb[4988]=8. + m_z1k;

  Box1Abb[4989]=-3. + Box1Abb[4988]*m_z1k;

  Box1Abb[4990]=-15. + Box1Abb[1776]*m_z1k;

  Box1Abb[4991]=-Box1Abb[210]*pow(Box1Abb[68],5.) - Box1Abb[4990]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[4987]*pow(Box1Abb[68],2.)*m_z12_2 - 2.*Box1Abb[4989]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[437]*m_z12_4*m_z1k;

  Box1Abb[4992]=m_z1k + 28.*m_z1k_3;

  Box1Abb[4993]=20. + 43.*m_z1k;

  Box1Abb[4994]=6. + Box1Abb[4993]*m_z1k;

  Box1Abb[4995]=15. + 3.*Box1Abb[4992]*m_z12 + Box1Abb[4994]*m_z12_2 + Box1Abb[4925]*m_z12_3 + 21.*Box1Abb[119]*m_z1k;

  Box1Abb[4996]=1. + Box1Abb[3089]*m_z1k;

  Box1Abb[4997]=7. + 19.*m_z1k;

  Box1Abb[4998]=-5. + Box1Abb[4997]*m_z1k;

  Box1Abb[4999]=18. + 77.*m_z1k;

  Box1Abb[5000]=-18. + Box1Abb[4999]*m_z1k;

  Box1Abb[5001]=4. + Box1Abb[5000]*m_z1k;

  Box1Abb[5002]=36. + 7.*Box1Abb[4190]*m_z1k;

  Box1Abb[5003]=45. + 2.*Box1Abb[5002]*m_z1k;

  Box1Abb[5004]=19. + Box1Abb[5003]*m_z1k;

  Box1Abb[5005]=Box1Abb[5004]*m_z12 + Box1Abb[5001]*m_z12_2 + 2.*Box1Abb[4998]*m_z12_3 + 15.*Box1Abb[4996]*m_z1k;

  Box1Abb[5006]=16. - 9.*m_z1k;

  Box1Abb[5007]=11. + Box1Abb[5006]*m_z1k;

  Box1Abb[5008]=-12. + Box1Abb[5007]*m_z1k;

  Box1Abb[5009]=19. + m_z1k;

  Box1Abb[5010]=16. + Box1Abb[5009]*m_z1k;

  Box1Abb[5011]=-33. + Box1Abb[5010]*m_z1k;

  Box1Abb[5012]=3. + Box1Abb[5011]*m_z1k;

  Box1Abb[5013]=2. + Box1Abb[4817]*m_z1k;

  Box1Abb[5014]=-21. + Box1Abb[5013]*m_z1k;

  Box1Abb[5015]=15. + Box1Abb[5014]*m_z1k;

  Box1Abb[5016]=16. + 23.*m_z1k;

  Box1Abb[5017]=9. + Box1Abb[5016]*m_z1k;

  Box1Abb[5018]=-54. + Box1Abb[5017]*m_z1k;

  Box1Abb[5019]=18. + Box1Abb[5018]*m_z1k;

  Box1Abb[5020]=3.*pow(Box1Abb[68],6.) + Box1Abb[5015]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[5019]*pow(Box1Abb[68],2.)*m_z12_2 + 2.*Box1Abb[5012]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[5008]*m_z12_4*m_z1k - 2.*m_z12_5*m_z1k_3;

  Box1Abb[5021]=4. - 7.*m_z1k;

  Box1Abb[5022]=36. + 55.*m_z1k;

  Box1Abb[5023]=6. - Box1Abb[5022]*m_z1k;

  Box1Abb[5024]=10. + Box1Abb[5023]*m_z1k;

  Box1Abb[5025]=6. - 91.*m_z1k;

  Box1Abb[5026]=9. + Box1Abb[5025]*m_z1k;

  Box1Abb[5027]=16. + Box1Abb[5026]*m_z1k;

  Box1Abb[5028]=-6. + Box1Abb[5027]*m_z1k;

  Box1Abb[5029]=-5. + 9.*m_z1k;

  Box1Abb[5030]=25. + 2.*Box1Abb[5029]*m_z1k;

  Box1Abb[5031]=45. + 7.*Box1Abb[5030]*m_z1k;

  Box1Abb[5032]=55. + Box1Abb[5031]*m_z1k;

  Box1Abb[5033]=23. + Box1Abb[5032]*m_z1k;

  Box1Abb[5034]=13. - Box1Abb[5033]*m_z12 + Box1Abb[5028]*m_z12_2 + Box1Abb[5024]*m_z12_3 - 12.*m_z12_4*m_z1k_2 + 15.*Box1Abb[5021]*m_z1k_3;

  Box1Abb[5035]=-40. + 21.*m_z1k;

  Box1Abb[5036]=20. + Box1Abb[5035]*m_z1k;

  Box1Abb[5037]=15. + 7.*Box1Abb[68]*m_z1k;

  Box1Abb[5038]=-55. + 6.*Box1Abb[5037]*m_z1k;

  Box1Abb[5039]=4. + Box1Abb[5038]*m_z1k_3;

  Box1Abb[5040]=23. + 20.*m_z1k;

  Box1Abb[5041]=-18. + Box1Abb[5040]*m_z1k;

  Box1Abb[5042]=-32. + Box1Abb[5041]*m_z1k;

  Box1Abb[5043]=-5. + Box1Abb[5042]*m_z1k;

  Box1Abb[5044]=-20. + 77.*m_z1k;

  Box1Abb[5045]=99. + Box1Abb[5044]*m_z1k;

  Box1Abb[5046]=16. + Box1Abb[5045]*m_z1k;

  Box1Abb[5047]=50. + Box1Abb[5046]*m_z1k;

  Box1Abb[5048]=6. + Box1Abb[5047]*m_z1k;

  Box1Abb[5049]=-8. + 2.*Box1Abb[5039]*m_z12 + Box1Abb[5048]*m_z12_2 + Box1Abb[5043]*m_z12_3 + m_z1k + 4.*Box1Abb[514]*m_z12_4*m_z1k_2 + 3.*Box1Abb[5036]*m_z1k_3;

  Box1Abb[5050]=9. + Box1Abb[2312]*m_z1k;

  Box1Abb[5051]=-12. + 7.*m_z1k;

  Box1Abb[5052]=4. + Box1Abb[5051]*m_z1k;

  Box1Abb[5053]=-1. + 3.*Box1Abb[5052]*m_z1k_2;

  Box1Abb[5054]=61. + 4.*Box1Abb[5029]*m_z1k;

  Box1Abb[5055]=-89. + Box1Abb[5054]*m_z1k;

  Box1Abb[5056]=-3. + Box1Abb[5055]*m_z1k;

  Box1Abb[5057]=33. + Box1Abb[5056]*m_z1k;

  Box1Abb[5058]=-43. + 11.*m_z1k;

  Box1Abb[5059]=-40. + Box1Abb[5058]*m_z1k;

  Box1Abb[5060]=-38. + Box1Abb[5059]*m_z1k;

  Box1Abb[5061]=25. + Box1Abb[5060]*m_z1k;

  Box1Abb[5062]=1. + Box1Abb[5061]*m_z1k;

  Box1Abb[5063]=24. - 49.*m_z1k;

  Box1Abb[5064]=-96. + Box1Abb[5063]*m_z1k;

  Box1Abb[5065]=108. + Box1Abb[5064]*m_z1k;

  Box1Abb[5066]=51. + Box1Abb[5065]*m_z1k;

  Box1Abb[5067]=-60. + Box1Abb[5066]*m_z1k;

  Box1Abb[5068]=-2. + Box1Abb[5067]*m_z1k;

  Box1Abb[5069]=-Box1Abb[5053]*pow(Box1Abb[68],2.) + Box1Abb[5068]*m_z12_2 + Box1Abb[5062]*m_z12_3 - Box1Abb[5057]*Box1Abb[68]*m_z12*m_z1k + 2.*Box1Abb[5050]*m_z12_4*m_z1k_2 - 4.*m_z12_5*m_z1k_3;

  Box1Abb[5070]=Box1Abb[5069]*m_x_3 + Box1Abb[5049]*m_x_4 + Box1Abb[5034]*m_x_5 + Box1Abb[5005]*m_x_6 - Box1Abb[4995]*m_x_7 + Box1Abb[4986]*m_x_8 - Box1Abb[4984]*m_x_9 + m_x_10*m_z12 + Box1Abb[5020]*m_x_2*m_z1k + Box1Abb[4991]*Box1Abb[68]*m_x*m_z12*m_z1k_2 + Box1Abb[4]*Box1Abb[541]*pow(Box1Abb[68],3.)*m_z12_2*m_z1k_3;

  Box1Abb[5071]=-Box1Abb[5070]*m_cR + Box1Abb[4983]*Box1Abb[7]*m_cL*m_x;

  Box1Abb[5072]=29. + 2.*Box1Abb[3340]*m_z12;

  Box1Abb[5073]=20. + 9.*m_z12 + 3.*m_z12_3;

  Box1Abb[5074]=-16. + m_z12 + 6.*m_z12_2;

  Box1Abb[5075]=6. + Box1Abb[5074]*m_z12;

  Box1Abb[5076]=9. + 2.*Box1Abb[0]*m_z12;

  Box1Abb[5077]=15. + Box1Abb[2513]*m_z12 + 2.*Box1Abb[5072]*m_z12*m_z1k + 2.*Box1Abb[5073]*m_z12*m_z1k_2 - 5.*Box1Abb[5075]*m_z1k_3 + 5.*Box1Abb[5076]*m_z1k_4 + 56.*m_z12*m_z1k_5;

  Box1Abb[5078]=-59. + 24.*m_z12;

  Box1Abb[5079]=54. + Box1Abb[5078]*m_z12;

  Box1Abb[5080]=41. + 9.*Box1Abb[79]*m_z12;

  Box1Abb[5081]=-27. + Box1Abb[5080]*m_z12;

  Box1Abb[5082]=12. + 5.*m_z12;

  Box1Abb[5083]=-46. + Box1Abb[5082]*m_z12;

  Box1Abb[5084]=20. + Box1Abb[5083]*m_z12;

  Box1Abb[5085]=-15. + Box1Abb[5084]*m_z12;

  Box1Abb[5086]=9. + m_z12 - 12.*m_z12_2;

  Box1Abb[5087]=-9. + Box1Abb[5086]*m_z12;

  Box1Abb[5088]=-18. + 11.*m_z12;

  Box1Abb[5089]=18. + Box1Abb[5088]*m_z12;

  Box1Abb[5090]=3. + Box1Abb[493]*m_z12 - 6.*m_z1k + Box1Abb[5079]*m_z12*m_z1k + 2.*Box1Abb[5081]*m_z12*m_z1k_2 - 2.*Box1Abb[5085]*m_z1k_3 + 5.*Box1Abb[5087]*m_z1k_4 + Box1Abb[5089]*m_z1k_5 + 28.*m_z12*m_z1k_6;

  Box1Abb[5091]=2. + m_z12 + 8.*m_z1k;

  Box1Abb[5092]=3. + Box1Abb[5091]*m_z12;

  Box1Abb[5093]=3. + 5.*Box1Abb[68]*m_z1k;

  Box1Abb[5094]=-18. + 5.*Box1Abb[15]*m_z1k;

  Box1Abb[5095]=Box1Abb[210]*pow(Box1Abb[68],5.) + Box1Abb[5094]*pow(Box1Abb[68],3.)*m_z12 + 6.*Box1Abb[210]*pow(Box1Abb[68],3.)*m_z12_2 - 2.*Box1Abb[5093]*Box1Abb[68]*m_z12_3 + 4.*Box1Abb[1723]*m_z12_4*m_z1k;

  Box1Abb[5096]=5. + 14.*m_z1k;

  Box1Abb[5097]=-5. + 5.*Box1Abb[15]*m_z12 + 2.*Box1Abb[5096]*m_z1k;

  Box1Abb[5098]=3.*Box1Abb[1989] + Box1Abb[5097]*m_z12;

  Box1Abb[5099]=9. + 28.*m_z1k;

  Box1Abb[5100]=-10. + Box1Abb[5099]*m_z1k_2;

  Box1Abb[5101]=30. + 2.*Box1Abb[5100]*m_z12 + 10.*Box1Abb[821]*m_z12_2 + 45.*Box1Abb[15]*m_z1k + 6.*m_z12_3*m_z1k;

  Box1Abb[5102]=9. + Box1Abb[2295]*m_z1k;

  Box1Abb[5103]=9. + Box1Abb[5102]*Box1Abb[68]*m_z1k;

  Box1Abb[5104]=-8. - Box1Abb[170]*Box1Abb[2440]*m_z1k;

  Box1Abb[5105]=1. + Box1Abb[5104]*m_z1k;

  Box1Abb[5106]=-41. + 9.*m_z1k + 5.*m_z1k_3;

  Box1Abb[5107]=36. + Box1Abb[5106]*m_z1k;

  Box1Abb[5108]=-9. + Box1Abb[5107]*m_z1k;

  Box1Abb[5109]=3.*pow(Box1Abb[68],5.) + 2.*Box1Abb[5103]*pow(Box1Abb[68],2.)*m_z12 + 2.*Box1Abb[5108]*m_z12_2 + 6.*Box1Abb[5105]*m_z12_3 - 2.*Box1Abb[317]*Box1Abb[4778]*m_z12_4*m_z1k + 2.*m_z12_5*m_z1k_2;

  Box1Abb[5110]=1. + m_z1k + m_z1k_2 + 2.*m_z1k_3;

  Box1Abb[5111]=10. - 18.*m_z1k + 11.*m_z1k_3;

  Box1Abb[5112]=9. + 2.*Box1Abb[3728]*m_z1k;

  Box1Abb[5113]=4. + Box1Abb[5112]*m_z1k;

  Box1Abb[5114]=-5. + Box1Abb[5113]*m_z1k;

  Box1Abb[5115]=30.*Box1Abb[5110] + 5.*Box1Abb[5114]*m_z12 + Box1Abb[5111]*m_z12_2 + 6.*Box1Abb[411]*m_z12_3*m_z1k;

  Box1Abb[5116]=Box1Abb[5090]*m_x_3 - Box1Abb[5077]*m_x_4 + Box1Abb[5115]*m_x_5 - Box1Abb[5101]*m_x_6 + Box1Abb[5098]*m_x_7 - Box1Abb[5092]*m_x_8 + m_x_9*m_z12 - Box1Abb[5109]*m_x_2*m_z1k + Box1Abb[5095]*m_x*m_z12*m_z1k_2 - Box1Abb[4784]*pow(Box1Abb[68],2.)*m_z12_2*m_z1k_3;

  Box1Abb[5117]=-67. - 11.*Box1Abb[158]*m_z12;

  Box1Abb[5118]=86. + 15.*m_z12;

  Box1Abb[5119]=-127. + Box1Abb[5118]*m_z12;

  Box1Abb[5120]=-16. + 17.*m_z12;

  Box1Abb[5121]=9. + Box1Abb[5120]*m_z12;

  Box1Abb[5122]=34. + Box1Abb[5117]*m_z12 + 66.*m_z1k + Box1Abb[5119]*m_z12*m_z1k + 7.*Box1Abb[5121]*m_z1k_2 + 84.*m_z12*m_z1k_3;

  Box1Abb[5123]=61. - 30.*m_z12;

  Box1Abb[5124]=-51. + Box1Abb[5123]*m_z12;

  Box1Abb[5125]=81. - 5.*m_z12;

  Box1Abb[5126]=-131. + Box1Abb[5125]*m_z12;

  Box1Abb[5127]=18. + Box1Abb[5126]*m_z12;

  Box1Abb[5128]=36. + 5.*m_z12;

  Box1Abb[5129]=-171. + 5.*Box1Abb[5128]*m_z12;

  Box1Abb[5130]=-20. + Box1Abb[5129]*m_z12;

  Box1Abb[5131]=30. + Box1Abb[5130]*m_z12;

  Box1Abb[5132]=-20. + 43.*m_z12;

  Box1Abb[5133]=9. + Box1Abb[5132]*m_z12;

  Box1Abb[5134]=4. + Box1Abb[5133]*m_z12;

  Box1Abb[5135]=19. + Box1Abb[5124]*m_z12 + 30.*m_z1k + Box1Abb[5127]*m_z12*m_z1k + Box1Abb[5131]*m_z1k_2 + 5.*Box1Abb[5134]*m_z1k_3 + 35.*Box1Abb[0]*Box1Abb[1417]*m_z1k_4 + 126.*m_z12*m_z1k_5;

  Box1Abb[5136]=-3. + 5.*m_z12 + 9.*m_z1k;

  Box1Abb[5137]=3. + Box1Abb[5136]*m_z12;

  Box1Abb[5138]=-7. + 9.*m_z1k;

  Box1Abb[5139]=22. + 37.*m_z1k;

  Box1Abb[5140]=-30. + Box1Abb[5139]*m_z12 - 2.*m_z12_2 + 4.*Box1Abb[5138]*m_z1k;

  Box1Abb[5141]=16. + Box1Abb[5140]*m_z12 + 21.*m_z1k;

  Box1Abb[5142]=3. + 2.*Box1Abb[881]*m_z1k;

  Box1Abb[5143]=-1. + Box1Abb[871]*m_z1k;

  Box1Abb[5144]=-3. + Box1Abb[2780]*m_z1k;

  Box1Abb[5145]=1. + Box1Abb[5144]*m_z1k;

  Box1Abb[5146]=13. + 5.*m_z1k;

  Box1Abb[5147]=-3. + Box1Abb[5146]*m_z1k;

  Box1Abb[5148]=1. + Box1Abb[5147]*m_z1k;

  Box1Abb[5149]=-Box1Abb[5143]*pow(Box1Abb[68],4.) - Box1Abb[5142]*pow(Box1Abb[68],4.)*m_z12 + 3.*Box1Abb[5145]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[5148]*Box1Abb[68]*m_z12_3 + 2.*Box1Abb[178]*m_z12_4*m_z1k_2;

  Box1Abb[5150]=-5. + 17.*m_z1k;

  Box1Abb[5151]=14. + 3.*Box1Abb[1776]*m_z1k;

  Box1Abb[5152]=7. + Box1Abb[5151]*m_z1k;

  Box1Abb[5153]=12. + 31.*m_z1k;

  Box1Abb[5154]=9. + 7.*Box1Abb[5153]*m_z1k;

  Box1Abb[5155]=59. + Box1Abb[5154]*m_z1k;

  Box1Abb[5156]=-9. + 7.*Box1Abb[170]*m_z1k;

  Box1Abb[5157]=-121. + 18.*Box1Abb[5156]*m_z1k;

  Box1Abb[5158]=-72. + Box1Abb[5157]*m_z1k;

  Box1Abb[5159]=5.*Box1Abb[5152] + Box1Abb[5158]*m_z12 + Box1Abb[5155]*m_z12_2 + Box1Abb[1989]*Box1Abb[5150]*m_z12_3 - m_z12_4*m_z1k;

  Box1Abb[5160]=-10. + 13.*m_z1k;

  Box1Abb[5161]=-24. + Box1Abb[5160]*m_z1k;

  Box1Abb[5162]=22. + Box1Abb[5161]*m_z1k;

  Box1Abb[5163]=-5. + Box1Abb[5162]*m_z1k;

  Box1Abb[5164]=-25. + 9.*m_z1k;

  Box1Abb[5165]=36. + Box1Abb[5164]*m_z1k;

  Box1Abb[5166]=-5. + Box1Abb[5165]*m_z1k;

  Box1Abb[5167]=-8. + Box1Abb[5166]*m_z1k;

  Box1Abb[5168]=3. + Box1Abb[5167]*m_z1k;

  Box1Abb[5169]=-59. + 10.*m_z1k;

  Box1Abb[5170]=28. + Box1Abb[5169]*m_z1k;

  Box1Abb[5171]=46. + Box1Abb[5170]*m_z1k;

  Box1Abb[5172]=-14. + Box1Abb[5171]*m_z1k;

  Box1Abb[5173]=1. + Box1Abb[5172]*m_z1k;

  Box1Abb[5174]=-66. + 19.*m_z1k;

  Box1Abb[5175]=74. + Box1Abb[5174]*m_z1k;

  Box1Abb[5176]=21. + Box1Abb[5175]*m_z1k;

  Box1Abb[5177]=-15. + Box1Abb[5176]*m_z1k;

  Box1Abb[5178]=3. + Box1Abb[5177]*m_z1k;

  Box1Abb[5179]=Box1Abb[4877]*pow(Box1Abb[68],5.) + Box1Abb[5168]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[5178]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[5173]*Box1Abb[68]*m_z12_3 + Box1Abb[5163]*m_z12_4*m_z1k + 2.*Box1Abb[2295]*m_z12_5*m_z1k_3;

  Box1Abb[5180]=10. - 60.*m_z1k + 63.*m_z1k_2;

  Box1Abb[5181]=31. + 38.*m_z1k;

  Box1Abb[5182]=-5. + Box1Abb[5181]*m_z1k;

  Box1Abb[5183]=47. + 105.*m_z1k;

  Box1Abb[5184]=9. + Box1Abb[5183]*m_z1k;

  Box1Abb[5185]=33. + Box1Abb[5184]*m_z1k;

  Box1Abb[5186]=-10. + Box1Abb[5185]*m_z1k;

  Box1Abb[5187]=-5. + 3.*m_z1k;

  Box1Abb[5188]=5. + 14.*Box1Abb[170]*m_z1k;

  Box1Abb[5189]=26. + Box1Abb[5187]*Box1Abb[5188]*m_z1k;

  Box1Abb[5190]=47. + Box1Abb[5189]*m_z1k;

  Box1Abb[5191]=-32. + 2.*Box1Abb[5190]*m_z1k;

  Box1Abb[5192]=-58. + 35.*m_z1k;

  Box1Abb[5193]=-74. + 5.*Box1Abb[5192]*m_z1k;

  Box1Abb[5194]=-130. + Box1Abb[5193]*m_z1k;

  Box1Abb[5195]=-145. + Box1Abb[5194]*m_z1k;

  Box1Abb[5196]=44. + Box1Abb[5195]*m_z1k;

  Box1Abb[5197]=8. + Box1Abb[5191]*m_z12 + Box1Abb[5196]*m_z12_2 + 2.*Box1Abb[5186]*m_z12_3 - 9.*m_z1k + 2.*Box1Abb[5182]*m_z12_4*m_z1k + Box1Abb[5180]*m_z1k_3;

  Box1Abb[5198]=27. + 30.*m_z1k + 34.*m_z1k_2;

  Box1Abb[5199]=-5. + Box1Abb[5198]*m_z1k;

  Box1Abb[5200]=4. + 3.*Box1Abb[1692]*m_z1k;

  Box1Abb[5201]=-4. + Box1Abb[5200]*m_z1k;

  Box1Abb[5202]=-83. + 93.*m_z1k;

  Box1Abb[5203]=-46. + Box1Abb[5202]*m_z1k;

  Box1Abb[5204]=-102. + Box1Abb[5203]*m_z1k;

  Box1Abb[5205]=37. + Box1Abb[5204]*m_z1k;

  Box1Abb[5206]=-7. + Box1Abb[5205]*m_z1k;

  Box1Abb[5207]=-11. + 3.*m_z1k;

  Box1Abb[5208]=61. + 4.*Box1Abb[5207]*m_z1k;

  Box1Abb[5209]=-52. + 3.*Box1Abb[5208]*m_z1k;

  Box1Abb[5210]=-10. + Box1Abb[5209]*m_z1k;

  Box1Abb[5211]=-46. + Box1Abb[5210]*m_z1k;

  Box1Abb[5212]=15. + Box1Abb[5211]*m_z1k;

  Box1Abb[5213]=-258. + 77.*m_z1k;

  Box1Abb[5214]=234. + Box1Abb[5213]*m_z1k;

  Box1Abb[5215]=-38. + Box1Abb[5214]*m_z1k;

  Box1Abb[5216]=63. + Box1Abb[5215]*m_z1k;

  Box1Abb[5217]=-72. + Box1Abb[5216]*m_z1k;

  Box1Abb[5218]=18. + Box1Abb[5217]*m_z1k;

  Box1Abb[5219]=Box1Abb[5201]*pow(Box1Abb[68],3.) + Box1Abb[5212]*Box1Abb[68]*m_z12 + Box1Abb[5218]*m_z12_2 + Box1Abb[5206]*m_z12_3 + 2.*Box1Abb[5199]*m_z12_4*m_z1k + 10.*m_z12_5*m_z1k_3;

  Box1Abb[5220]=Box1Abb[5179]*m_x - Box1Abb[5219]*m_x_2 + Box1Abb[5197]*m_x_3 - Box1Abb[5135]*m_x_4 + Box1Abb[5159]*m_x_5 - Box1Abb[5122]*m_x_6 + Box1Abb[5141]*m_x_7 - Box1Abb[5137]*m_x_8 + m_x_9*m_z12 + Box1Abb[5149]*Box1Abb[68]*m_z12*m_z1k;

  Box1Abb[5221]=-Box1Abb[5116]*Box1Abb[7]*m_cR + Box1Abb[5220]*m_cL*m_x;

  Box1Abb[5222]=11. + 5.*Box1Abb[1032]*m_z12;

  Box1Abb[5223]=-58. + Box1Abb[5222]*m_z12;

  Box1Abb[5224]=58. + 7.*m_z12;

  Box1Abb[5225]=156. + Box1Abb[5224]*m_z12;

  Box1Abb[5226]=-100. + Box1Abb[5225]*m_z12;

  Box1Abb[5227]=-25. + 9.*Box1Abb[3301]*m_z12;

  Box1Abb[5228]=60. + Box1Abb[5223]*m_z12 + 82.*m_z1k + Box1Abb[5226]*m_z12*m_z1k - 4.*Box1Abb[5227]*m_z1k_2 + 28.*Box1Abb[522]*m_z12*m_z1k_3;

  Box1Abb[5229]=4. + 7.*m_z12 - 20.*m_z1k;

  Box1Abb[5230]=16. + Box1Abb[5229]*m_z12 + 52.*m_z1k;

  Box1Abb[5231]=8. + Box1Abb[5230]*m_z12;

  Box1Abb[5232]=23. + 2.*Box1Abb[3496]*m_z12;

  Box1Abb[5233]=40. + 3.*m_z12;

  Box1Abb[5234]=14.*Box1Abb[170] + Box1Abb[5232]*m_z12 + 2.*Box1Abb[5233]*m_z12*m_z1k - 28.*Box1Abb[856]*m_z1k_2;

  Box1Abb[5235]=36. + Box1Abb[5234]*m_z12 + 44.*m_z1k;

  Box1Abb[5236]=9. + 4.*Box1Abb[437]*m_z1k;

  Box1Abb[5237]=-9. + 17.*m_z1k;

  Box1Abb[5238]=5. + Box1Abb[5237]*m_z1k;

  Box1Abb[5239]=-11. + 42.*m_z1k;

  Box1Abb[5240]=3. + Box1Abb[5239]*m_z1k;

  Box1Abb[5241]=-2.*Box1Abb[317]*pow(Box1Abb[68],4.) + Box1Abb[5236]*pow(Box1Abb[68],3.)*m_z12 + 2.*Box1Abb[5238]*pow(Box1Abb[68],2.)*m_z12_2 + Box1Abb[5240]*Box1Abb[68]*m_z12_3 + 16.*m_z12_4*m_z1k_2;

  Box1Abb[5242]=-20. + Box1Abb[1095]*m_z1k;

  Box1Abb[5243]=5. + 4.*m_z1k + 60.*m_z1k_2;

  Box1Abb[5244]=13. + Box1Abb[5243]*m_z1k;

  Box1Abb[5245]=-97. + 13.*m_z1k;

  Box1Abb[5246]=-63. + Box1Abb[5245]*m_z1k;

  Box1Abb[5247]=10. + Box1Abb[5246]*m_z1k;

  Box1Abb[5248]=37. - 14.*m_z1k;

  Box1Abb[5249]=-24. + 5.*Box1Abb[5248]*m_z1k;

  Box1Abb[5250]=5. + Box1Abb[5249]*m_z1k;

  Box1Abb[5251]=4. + Box1Abb[5250]*m_z1k;

  Box1Abb[5252]=19. + 7.*m_z1k;

  Box1Abb[5253]=63. + 20.*Box1Abb[5252]*m_z1k;

  Box1Abb[5254]=-29. + Box1Abb[5253]*m_z1k;

  Box1Abb[5255]=10. + Box1Abb[5254]*m_z1k;

  Box1Abb[5256]=-2.*Box1Abb[5244] + 2.*Box1Abb[5251]*m_z12 - Box1Abb[5255]*m_z12_2 + 2.*Box1Abb[5247]*m_z12_3 + Box1Abb[5242]*m_z12_4;

  Box1Abb[5257]=-5. + Box1Abb[1151]*m_z1k;

  Box1Abb[5258]=5. + 2.*Box1Abb[5257]*m_z1k;

  Box1Abb[5259]=-11. + 13.*m_z1k;

  Box1Abb[5260]=-27. + Box1Abb[5259]*m_z1k;

  Box1Abb[5261]=63. + 4.*Box1Abb[5260]*m_z1k;

  Box1Abb[5262]=-15. + Box1Abb[5261]*m_z1k;

  Box1Abb[5263]=17. + 35.*m_z1k;

  Box1Abb[5264]=-85. + 3.*Box1Abb[5263]*m_z1k;

  Box1Abb[5265]=23. + Box1Abb[5264]*m_z1k;

  Box1Abb[5266]=-2. + Box1Abb[5265]*m_z1k;

  Box1Abb[5267]=-5. + 78.*m_z1k;

  Box1Abb[5268]=-169. + 2.*Box1Abb[5267]*m_z1k;

  Box1Abb[5269]=72. + Box1Abb[5268]*m_z1k;

  Box1Abb[5270]=-9. + Box1Abb[5269]*m_z1k;

  Box1Abb[5271]=2.*Box1Abb[317]*pow(Box1Abb[68],5.) - 2.*Box1Abb[5258]*pow(Box1Abb[68],4.)*m_z12 + Box1Abb[5262]*pow(Box1Abb[68],3.)*m_z12_2 + Box1Abb[5270]*pow(Box1Abb[68],2.)*m_z12_3 + Box1Abb[5266]*Box1Abb[68]*m_z12_4 + 8.*Box1Abb[411]*m_z12_5*m_z1k_3;

  Box1Abb[5272]=1. + 2.*Box1Abb[1041]*m_z1k;

  Box1Abb[5273]=-58. + Box1Abb[4917]*m_z1k;

  Box1Abb[5274]=-13. + Box1Abb[5273]*m_z1k;

  Box1Abb[5275]=10. + Box1Abb[5274]*m_z1k;

  Box1Abb[5276]=70. + 27.*m_z1k;

  Box1Abb[5277]=77. + Box1Abb[5276]*m_z1k;

  Box1Abb[5278]=-25. + Box1Abb[5277]*m_z1k;

  Box1Abb[5279]=5. + Box1Abb[5278]*m_z1k;

  Box1Abb[5280]=102. + m_z1k - 121.*m_z1k_2;

  Box1Abb[5281]=128. + Box1Abb[5280]*m_z1k;

  Box1Abb[5282]=-89. + Box1Abb[5281]*m_z1k;

  Box1Abb[5283]=19. + Box1Abb[5282]*m_z1k;

  Box1Abb[5284]=-19. + 35.*m_z1k;

  Box1Abb[5285]=-343. + 4.*Box1Abb[5284]*m_z1k;

  Box1Abb[5286]=27. + Box1Abb[5285]*m_z1k;

  Box1Abb[5287]=135. + Box1Abb[5286]*m_z1k;

  Box1Abb[5288]=-51. + Box1Abb[5287]*m_z1k;

  Box1Abb[5289]=-2.*Box1Abb[5272]*pow(Box1Abb[68],3.) + 2.*Box1Abb[5275]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[5288]*Box1Abb[68]*m_z12_2 + 2.*Box1Abb[5283]*m_z12_3 - 2.*Box1Abb[5279]*m_z12_4 + 8.*m_z12_5*m_z1k_3;

  Box1Abb[5290]=47. - 31.*m_z1k;

  Box1Abb[5291]=-30. + Box1Abb[5290]*m_z1k;

  Box1Abb[5292]=20. + Box1Abb[5291]*m_z1k;

  Box1Abb[5293]=-7. + 20.*m_z1k;

  Box1Abb[5294]=-9. + 2.*Box1Abb[5293]*m_z1k;

  Box1Abb[5295]=5. + Box1Abb[5294]*m_z1k;

  Box1Abb[5296]=161. + 65.*m_z1k;

  Box1Abb[5297]=191. + 2.*Box1Abb[5296]*m_z1k;

  Box1Abb[5298]=172. + Box1Abb[5297]*m_z1k;

  Box1Abb[5299]=-55. + Box1Abb[5298]*m_z1k;

  Box1Abb[5300]=-80. + 7.*m_z1k;

  Box1Abb[5301]=98. + Box1Abb[5300]*m_z1k;

  Box1Abb[5302]=-13. + Box1Abb[5301]*m_z1k;

  Box1Abb[5303]=13. + Box1Abb[5302]*m_z1k;

  Box1Abb[5304]=1. + Box1Abb[5303]*m_z1k;

  Box1Abb[5305]=25. + 49.*m_z1k;

  Box1Abb[5306]=-92. + Box1Abb[5305]*m_z1k;

  Box1Abb[5307]=-59. + 2.*Box1Abb[5306]*m_z1k;

  Box1Abb[5308]=-86. + Box1Abb[5307]*m_z1k;

  Box1Abb[5309]=25. + Box1Abb[5308]*m_z1k;

  Box1Abb[5310]=2.*Box1Abb[5295]*Box1Abb[68] + 4.*Box1Abb[5304]*m_z12 + 2.*Box1Abb[5309]*m_z12_2 + Box1Abb[5299]*m_z12_3 + Box1Abb[5292]*m_z12_4;

  Box1Abb[5311]=Box1Abb[5271]*m_x_2 + Box1Abb[5289]*m_x_3 + Box1Abb[5310]*m_x_4 + Box1Abb[5256]*m_x_5 + Box1Abb[5228]*m_x_6 - Box1Abb[5235]*m_x_7 + Box1Abb[5231]*m_x_8 + 4.*Box1Abb[79]*m_x_9*m_z12 - Box1Abb[5241]*pow(Box1Abb[68],2.)*m_x*m_z12*m_z1k - Box1Abb[4]*pow(Box1Abb[68],5.)*m_z12_3*m_z1k_2;

  Box1Abb[5312]=32. + m_z12;

  Box1Abb[5313]=53. + Box1Abb[5312]*m_z12;

  Box1Abb[5314]=-82. + Box1Abb[5313]*m_z12;

  Box1Abb[5315]=-36. + Box1Abb[5314]*m_z12;

  Box1Abb[5316]=42. + 11.*m_z12;

  Box1Abb[5317]=-55. + 4.*Box1Abb[5316]*m_z12;

  Box1Abb[5318]=-183. + Box1Abb[5317]*m_z12;

  Box1Abb[5319]=274. + 55.*m_z12;

  Box1Abb[5320]=-191. + Box1Abb[5319]*m_z12;

  Box1Abb[5321]=72. + Box1Abb[5320]*m_z12;

  Box1Abb[5322]=14. - 25.*m_z12;

  Box1Abb[5323]=40. + Box1Abb[5315]*m_z12 + 68.*m_z1k + Box1Abb[5318]*m_z12*m_z1k + Box1Abb[5321]*m_z1k_2 + 12.*Box1Abb[5322]*m_z12*m_z1k_3;

  Box1Abb[5324]=-3. + Box1Abb[2661]*m_z12;

  Box1Abb[5325]=-6. + 5.*Box1Abb[5324]*m_z12;

  Box1Abb[5326]=-22. + 5.*Box1Abb[1790]*m_z12;

  Box1Abb[5327]=-53. + Box1Abb[5326]*m_z12;

  Box1Abb[5328]=7. + Box1Abb[5327]*m_z12;

  Box1Abb[5329]=16. + Box1Abb[5328]*m_z12;

  Box1Abb[5330]=37. + 7.*m_z12;

  Box1Abb[5331]=-125. + 2.*Box1Abb[5330]*m_z12;

  Box1Abb[5332]=31. + Box1Abb[5331]*m_z12;

  Box1Abb[5333]=-707. + 474.*m_z12 + 4.*m_z12_2;

  Box1Abb[5334]=315. + Box1Abb[5333]*m_z12;

  Box1Abb[5335]=-40. + Box1Abb[5334]*m_z12;

  Box1Abb[5336]=842. - 277.*m_z12;

  Box1Abb[5337]=-555. + Box1Abb[5336]*m_z12;

  Box1Abb[5338]=100. + Box1Abb[5337]*m_z12;

  Box1Abb[5339]=2.*Box1Abb[0]*Box1Abb[5325]*m_z12 + Box1Abb[5329]*m_z1k + Box1Abb[12]*Box1Abb[5332]*m_z12*m_z1k_2 + Box1Abb[5335]*m_z1k_3 + Box1Abb[5338]*m_z1k_4 + 42.*Box1Abb[810]*m_z12*m_z1k_5;

  Box1Abb[5340]=-5. + Box1Abb[3493]*m_z12;

  Box1Abb[5341]=69. - 5.*Box1Abb[5340]*m_z12;

  Box1Abb[5342]=-24. + Box1Abb[5341]*m_z12;

  Box1Abb[5343]=65. + 6.*Box1Abb[2353]*m_z12;

  Box1Abb[5344]=-212. + Box1Abb[5343]*m_z12;

  Box1Abb[5345]=9. + Box1Abb[5344]*m_z12;

  Box1Abb[5346]=279. + 76.*m_z12;

  Box1Abb[5347]=227. - Box1Abb[5346]*m_z12;

  Box1Abb[5348]=140. + Box1Abb[5347]*m_z12;

  Box1Abb[5349]=-60. + Box1Abb[5348]*m_z12;

  Box1Abb[5350]=-684. + 101.*m_z12;

  Box1Abb[5351]=467. + Box1Abb[5350]*m_z12;

  Box1Abb[5352]=-110. + Box1Abb[5351]*m_z12;

  Box1Abb[5353]=6. + 25.*m_z1k;

  Box1Abb[5354]=-2.*Box1Abb[5353] + Box1Abb[5342]*m_z12 - Box1Abb[5345]*m_z12*m_z1k + Box1Abb[5349]*m_z1k_2 + Box1Abb[5352]*m_z1k_3 + 168.*Box1Abb[8]*m_z12*m_z1k_4;

  Box1Abb[5355]=-2. + 14.*m_z12 - 47.*m_z1k;

  Box1Abb[5356]=2. + Box1Abb[5355]*m_z12 + 30.*m_z1k;

  Box1Abb[5357]=4. + Box1Abb[5356]*m_z12;

  Box1Abb[5358]=-9. + 32.*m_z1k;

  Box1Abb[5359]=52. + 58.*m_z1k;

  Box1Abb[5360]=40. - 159.*m_z1k;

  Box1Abb[5361]=-19. + Box1Abb[5360]*m_z1k;

  Box1Abb[5362]=-46. + Box1Abb[5361]*m_z12 + Box1Abb[5359]*m_z12_2 + 7.*m_z12_3 + 3.*Box1Abb[5358]*m_z1k;

  Box1Abb[5363]=20. + Box1Abb[5362]*m_z12 + 26.*m_z1k;

  Box1Abb[5364]=-1. + Box1Abb[2087]*m_z1k;

  Box1Abb[5365]=-4. + Box1Abb[4190]*m_z1k;

  Box1Abb[5366]=-5. + m_z1k + 3.*m_z1k_2;

  Box1Abb[5367]=1. + Box1Abb[5366]*m_z1k;

  Box1Abb[5368]=-Box1Abb[431]*pow(Box1Abb[68],4.) + Box1Abb[5365]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[4774]*pow(Box1Abb[68],2.)*m_z12_2 + 4.*Box1Abb[5367]*m_z12_3 + Box1Abb[5364]*m_z12_4;

  Box1Abb[5369]=-7. + 3.*Box1Abb[2295]*m_z1k;

  Box1Abb[5370]=1. + m_z1k + Box1Abb[1103]*Box1Abb[431]*m_z1k_2;

  Box1Abb[5371]=-69. + 4.*m_z1k + 30.*m_z1k_2;

  Box1Abb[5372]=8. + Box1Abb[5371]*m_z1k;

  Box1Abb[5373]=-1. + Box1Abb[5372]*m_z1k;

  Box1Abb[5374]=-16. + 3.*m_z1k;

  Box1Abb[5375]=15. + 4.*Box1Abb[5374]*m_z1k;

  Box1Abb[5376]=2. + Box1Abb[5375]*m_z1k;

  Box1Abb[5377]=3. + Box1Abb[5376]*m_z1k;

  Box1Abb[5378]=68. + 25.*m_z1k;

  Box1Abb[5379]=-72. + Box1Abb[5378]*m_z1k;

  Box1Abb[5380]=2. + Box1Abb[5379]*m_z1k;

  Box1Abb[5381]=-3. + Box1Abb[5380]*m_z1k;

  Box1Abb[5382]=-Box1Abb[5373]*pow(Box1Abb[68],3.)*m_z12_2 - Box1Abb[5381]*pow(Box1Abb[68],2.)*m_z12_3 + Box1Abb[5377]*Box1Abb[68]*m_z12_4 + Box1Abb[5370]*m_z12_5 - 2.*pow(Box1Abb[68],5.)*m_z1k + Box1Abb[5369]*pow(Box1Abb[68],4.)*m_z12*m_z1k;

  Box1Abb[5383]=-2. + m_z1k_2 - 25.*m_z1k_3 + 48.*m_z1k_4 - 23.*m_z1k_5;

  Box1Abb[5384]=12. - 17.*m_z1k;

  Box1Abb[5385]=-10. + Box1Abb[5384]*m_z1k;

  Box1Abb[5386]=5. + Box1Abb[5385]*m_z1k_2;

  Box1Abb[5387]=109. + 24.*m_z1k;

  Box1Abb[5388]=-64. + Box1Abb[5387]*m_z1k;

  Box1Abb[5389]=-3. + Box1Abb[5388]*m_z1k;

  Box1Abb[5390]=-4. + Box1Abb[5389]*m_z1k;

  Box1Abb[5391]=-301. + 6.*Box1Abb[2818]*m_z1k;

  Box1Abb[5392]=90. + Box1Abb[5391]*m_z1k;

  Box1Abb[5393]=-13. + Box1Abb[5392]*m_z1k;

  Box1Abb[5394]=10. + Box1Abb[5393]*m_z1k;

  Box1Abb[5395]=-305. + 67.*m_z1k;

  Box1Abb[5396]=284. + Box1Abb[5395]*m_z1k;

  Box1Abb[5397]=-40. + Box1Abb[5396]*m_z1k;

  Box1Abb[5398]=17. + Box1Abb[5397]*m_z1k;

  Box1Abb[5399]=-3. + Box1Abb[5398]*m_z1k;

  Box1Abb[5400]=-Box1Abb[5390]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[5394]*pow(Box1Abb[68],2.)*m_z12_2 - Box1Abb[5399]*Box1Abb[68]*m_z12_3 + 4.*Box1Abb[5383]*m_z12_4 + Box1Abb[5386]*m_z12_5 + 4.*Box1Abb[174]*pow(Box1Abb[68],4.)*m_z1k;

  Box1Abb[5401]=4. + 27.*m_z1k;

  Box1Abb[5402]=1. + Box1Abb[5401]*m_z1k;

  Box1Abb[5403]=10. + Box1Abb[475]*m_z1k;

  Box1Abb[5404]=10. + Box1Abb[5403]*m_z1k;

  Box1Abb[5405]=-57. + 53.*m_z1k;

  Box1Abb[5406]=20. + m_z1k + 2.*Box1Abb[5405]*m_z1k_2;

  Box1Abb[5407]=-5. + Box1Abb[5406]*m_z1k;

  Box1Abb[5408]=965. - 488.*m_z1k + 42.*m_z1k_2;

  Box1Abb[5409]=-512. + Box1Abb[5408]*m_z1k;

  Box1Abb[5410]=-39. + 72.*m_z1k + Box1Abb[5409]*m_z1k_3;

  Box1Abb[5411]=-297. + 313.*m_z1k;

  Box1Abb[5412]=35. + Box1Abb[5411]*m_z1k;

  Box1Abb[5413]=11. + Box1Abb[5412]*m_z1k;

  Box1Abb[5414]=-8. + Box1Abb[5413]*m_z1k;

  Box1Abb[5415]=-625. + 241.*m_z1k;

  Box1Abb[5416]=375. + Box1Abb[5415]*m_z1k;

  Box1Abb[5417]=35. + Box1Abb[5416]*m_z1k;

  Box1Abb[5418]=-60. + Box1Abb[5417]*m_z1k;

  Box1Abb[5419]=46. + Box1Abb[5418]*m_z1k;

  Box1Abb[5420]=Box1Abb[5414]*Box1Abb[68]*m_z12 + Box1Abb[5410]*m_z12_2 + Box1Abb[5419]*m_z12_3 + Box1Abb[5407]*m_z12_4 - Box1Abb[5404]*m_z12_5 - 2.*Box1Abb[5402]*pow(Box1Abb[68],2.)*m_z1k;

  Box1Abb[5421]=Box1Abb[5382]*Box1Abb[68]*m_x_2 + Box1Abb[5400]*m_x_3 + Box1Abb[5420]*m_x_4 + Box1Abb[5339]*m_x_5 + Box1Abb[5354]*m_x_6 + Box1Abb[5323]*m_x_7 - Box1Abb[5363]*m_x_8 + Box1Abb[5357]*m_x_9 + 2.*Box1Abb[133]*m_x_10*m_z12 + Box1Abb[5368]*pow(Box1Abb[68],3.)*m_x*m_z12*m_z1k - pow(Box1Abb[4],3.)*pow(Box1Abb[68],5.)*m_z12_2*m_z1k_2;

  Box1Abb[5422]=-2.*Box1Abb[5421]*m_cL + Box1Abb[5311]*Box1Abb[7]*m_cR;

  Box1Abb[5423]=-12. + Box1Abb[707]*m_z12;

  Box1Abb[5424]=-9. + 5.*Box1Abb[5423]*m_z12;

  Box1Abb[5425]=32. + Box1Abb[5424]*m_z12;

  Box1Abb[5426]=-6. + 5.*Box1Abb[708]*m_z12;

  Box1Abb[5427]=-84. + Box1Abb[5426]*m_z12;

  Box1Abb[5428]=32. + Box1Abb[5427]*m_z12;

  Box1Abb[5429]=72. + 19.*m_z12;

  Box1Abb[5430]=-63. + Box1Abb[5429]*m_z12;

  Box1Abb[5431]=-5. + Box1Abb[5430]*m_z12;

  Box1Abb[5432]=62. - 17.*m_z12;

  Box1Abb[5433]=-46. + Box1Abb[5432]*m_z12;

  Box1Abb[5434]=10. + Box1Abb[5433]*m_z12;

  Box1Abb[5435]=-7. + Box1Abb[5425]*m_z12 + Box1Abb[5428]*m_z12*m_z1k + 2.*Box1Abb[5431]*m_z12*m_z1k_2 + 10.*Box1Abb[5434]*m_z1k_3 - 35.*Box1Abb[0]*Box1Abb[860]*m_z1k_4 - 126.*m_z12*m_z1k_5;

  Box1Abb[5436]=16. + Box1Abb[2598]*m_z12;

  Box1Abb[5437]=-60. + Box1Abb[5436]*m_z12;

  Box1Abb[5438]=10. + Box1Abb[5437]*m_z12;

  Box1Abb[5439]=77. + 34.*m_z12;

  Box1Abb[5440]=21. - Box1Abb[5439]*m_z12;

  Box1Abb[5441]=71. + Box1Abb[5440]*m_z12;

  Box1Abb[5442]=-74. + 7.*m_z12;

  Box1Abb[5443]=27. + Box1Abb[5442]*m_z12;

  Box1Abb[5444]=-52. + 47.*m_z12;

  Box1Abb[5445]=15. + Box1Abb[5444]*m_z12;

  Box1Abb[5446]=-5.*Box1Abb[178] - Box1Abb[5438]*m_z12 + Box1Abb[5441]*m_z12*m_z1k + 3.*Box1Abb[5443]*m_z12*m_z1k_2 + 7.*Box1Abb[5445]*m_z1k_3 + 126.*m_z12*m_z1k_4;

  Box1Abb[5447]=-5. + 7.*m_z12 + 9.*m_z1k;

  Box1Abb[5448]=3. + Box1Abb[5447]*m_z12;

  Box1Abb[5449]=8. + 53.*m_z1k;

  Box1Abb[5450]=-19. + Box1Abb[5449]*m_z12 - 11.*m_z12_2 + 4.*Box1Abb[1694]*m_z1k;

  Box1Abb[5451]=10. + Box1Abb[5450]*m_z12 + 21.*m_z1k;

  Box1Abb[5452]=2. - 25.*m_z1k;

  Box1Abb[5453]=-25. + 7.*Box1Abb[5452]*m_z1k;

  Box1Abb[5454]=4. - 7.*Box1Abb[170]*m_z1k;

  Box1Abb[5455]=-3. + 12.*Box1Abb[5454]*m_z1k;

  Box1Abb[5456]=-7. + Box1Abb[5455]*m_z12 + Box1Abb[5453]*m_z12_2 + 5.*Box1Abb[2006]*m_z12_3 + 6.*m_z12_4 - 3.*Box1Abb[819]*m_z1k;

  Box1Abb[5457]=1. + 26.*m_z1k;

  Box1Abb[5458]=-5. + Box1Abb[884]*m_z1k;

  Box1Abb[5459]=5. + Box1Abb[5458]*m_z1k;

  Box1Abb[5460]=-33. + 29.*m_z1k;

  Box1Abb[5461]=3. + Box1Abb[5460]*m_z1k;

  Box1Abb[5462]=9. + Box1Abb[5461]*m_z1k;

  Box1Abb[5463]=Box1Abb[536]*pow(Box1Abb[68],4.) + Box1Abb[5029]*pow(Box1Abb[68],5.)*m_z12 + Box1Abb[5457]*pow(Box1Abb[68],4.)*m_z12_2 + Box1Abb[5462]*Box1Abb[68]*m_z12_3 + Box1Abb[5459]*m_z12_4;

  Box1Abb[5464]=10. - 4.*m_z1k_3;

  Box1Abb[5465]=-5. + 3.*Box1Abb[1041]*m_z1k;

  Box1Abb[5466]=-9. + 8.*Box1Abb[3396]*m_z1k;

  Box1Abb[5467]=10. + Box1Abb[5466]*m_z1k;

  Box1Abb[5468]=5. + Box1Abb[5467]*m_z1k;

  Box1Abb[5469]=-29. + 9.*m_z1k;

  Box1Abb[5470]=30. + Box1Abb[5469]*m_z1k;

  Box1Abb[5471]=-1. + 4.*Box1Abb[5470]*m_z1k;

  Box1Abb[5472]=-29. + Box1Abb[5471]*m_z1k;

  Box1Abb[5473]=-316. + 133.*m_z1k;

  Box1Abb[5474]=156. + Box1Abb[5473]*m_z1k;

  Box1Abb[5475]=56. + Box1Abb[5474]*m_z1k;

  Box1Abb[5476]=-53. + Box1Abb[5475]*m_z1k;

  Box1Abb[5477]=-293. + 174.*m_z1k;

  Box1Abb[5478]=39. + Box1Abb[5477]*m_z1k;

  Box1Abb[5479]=81. + Box1Abb[5478]*m_z1k;

  Box1Abb[5480]=-29. + Box1Abb[5479]*m_z1k;

  Box1Abb[5481]=-Box1Abb[5465]*pow(Box1Abb[68],4.) - Box1Abb[5472]*pow(Box1Abb[68],3.)*m_z12 - Box1Abb[5476]*pow(Box1Abb[68],2.)*m_z12_2 - Box1Abb[5480]*Box1Abb[68]*m_z12_3 - 2.*Box1Abb[5468]*m_z12_4 + Box1Abb[5464]*m_z12_5;

  Box1Abb[5482]=-20. - 36.*m_z1k_2 + 34.*m_z1k_3;

  Box1Abb[5483]=-8. + 21.*m_z1k;

  Box1Abb[5484]=-11. + 3.*Box1Abb[5483]*m_z1k;

  Box1Abb[5485]=2. + Box1Abb[5484]*m_z1k;

  Box1Abb[5486]=-362. + 255.*m_z1k;

  Box1Abb[5487]=36. + Box1Abb[5486]*m_z1k;

  Box1Abb[5488]=-18. + Box1Abb[5487]*m_z1k;

  Box1Abb[5489]=85. + Box1Abb[5488]*m_z1k;

  Box1Abb[5490]=614. + 41.*Box1Abb[4956]*m_z1k;

  Box1Abb[5491]=-12. + Box1Abb[5490]*m_z1k;

  Box1Abb[5492]=27. + Box1Abb[5491]*m_z1k;

  Box1Abb[5493]=-72. + Box1Abb[5492]*m_z1k;

  Box1Abb[5494]=35. + 4.*Box1Abb[2406]*m_z1k;

  Box1Abb[5495]=-150. + 7.*Box1Abb[5494]*m_z1k;

  Box1Abb[5496]=8. + Box1Abb[5495]*m_z1k;

  Box1Abb[5497]=4. + Box1Abb[5496]*m_z1k;

  Box1Abb[5498]=5. + Box1Abb[5497]*m_z1k;

  Box1Abb[5499]=Box1Abb[5485]*pow(Box1Abb[68],2.) + 3.*Box1Abb[5498]*m_z12 + Box1Abb[5493]*m_z12_2 + Box1Abb[5489]*m_z12_3 + Box1Abb[5482]*m_z12_4 - 10.*Box1Abb[821]*m_z12_5;

  Box1Abb[5500]=Box1Abb[4]*Box1Abb[5463]*Box1Abb[68]*m_x + Box1Abb[5481]*m_x_2 + Box1Abb[5499]*m_x_3 + Box1Abb[5435]*m_x_4 + Box1Abb[5446]*m_x_5 + Box1Abb[5456]*m_x_6 + Box1Abb[5451]*m_x_7 - Box1Abb[5448]*m_x_8 - pow(Box1Abb[4],4.)*pow(Box1Abb[68],5.)*m_z12 + m_x_9*m_z12;

  Box1Abb[5501]=3. + 3.*m_z12_2 + 10.*m_z12*m_z1k;

  Box1Abb[5502]=4. - 45.*m_z1k;

  Box1Abb[5503]=2. - 5.*Box1Abb[532]*m_z12 + 3.*m_z12_2 + Box1Abb[5502]*m_z1k;

  Box1Abb[5504]=-3.*Box1Abb[2088] + Box1Abb[5503]*m_z12;

  Box1Abb[5505]=2. - 5.*m_z1k + 3.*m_z1k_3;

  Box1Abb[5506]=-9. + Box1Abb[3537]*m_z1k;

  Box1Abb[5507]=Box1Abb[210]*pow(Box1Abb[68],4.) + Box1Abb[5506]*pow(Box1Abb[68],2.)*m_z12 + 3.*Box1Abb[5505]*m_z12_2 + 8.*m_z12_3*m_z1k;

  Box1Abb[5508]=16. + 93.*m_z1k;

  Box1Abb[5509]=8. + Box1Abb[5508]*m_z1k;

  Box1Abb[5510]=-4. + 15.*m_z1k;

  Box1Abb[5511]=2. + Box1Abb[5510]*m_z1k;

  Box1Abb[5512]=2. + Box1Abb[5511]*m_z1k;

  Box1Abb[5513]=3. + 8.*Box1Abb[5512]*m_z12 + Box1Abb[5509]*m_z12_2 - 2.*Box1Abb[2780]*m_z12_3 - m_z12_4 + 33.*m_z1k + 84.*m_z1k_2;

  Box1Abb[5514]=9. - 22.*m_z1k;

  Box1Abb[5515]=5. + 2.*Box1Abb[5514]*m_z1k;

  Box1Abb[5516]=3. + 56.*m_z1k;

  Box1Abb[5517]=4. - Box1Abb[5516]*m_z1k;

  Box1Abb[5518]=5. + Box1Abb[5517]*m_z1k;

  Box1Abb[5519]=-5. + 204.*m_z1k;

  Box1Abb[5520]=13. + Box1Abb[5519]*m_z1k;

  Box1Abb[5521]=22. + Box1Abb[5520]*m_z1k;

  Box1Abb[5522]=69. + 7.*Box1Abb[1118]*m_z1k;

  Box1Abb[5523]=42. + Box1Abb[5522]*m_z1k;

  Box1Abb[5524]=11. + Box1Abb[5523]*m_z1k;

  Box1Abb[5525]=3.*Box1Abb[5518] - 2.*Box1Abb[5524]*m_z12 - Box1Abb[5521]*m_z12_2 + Box1Abb[5515]*m_z12_3 + 5.*Box1Abb[15]*m_z12_4;

  Box1Abb[5526]=6. - 12.*m_z1k + 4.*m_z1k_3 + 5.*m_z1k_4;

  Box1Abb[5527]=21. + 2.*m_z1k;

  Box1Abb[5528]=7. + Box1Abb[5527]*m_z1k;

  Box1Abb[5529]=-15. + Box1Abb[5528]*m_z1k;

  Box1Abb[5530]=23. + 17.*m_z1k;

  Box1Abb[5531]=8. + Box1Abb[5530]*m_z1k;

  Box1Abb[5532]=-33. + Box1Abb[5531]*m_z1k;

  Box1Abb[5533]=3. + Box1Abb[5532]*m_z1k;

  Box1Abb[5534]=22. + 39.*m_z1k;

  Box1Abb[5535]=2. + Box1Abb[5534]*m_z1k;

  Box1Abb[5536]=-54. + Box1Abb[5535]*m_z1k;

  Box1Abb[5537]=15. + Box1Abb[5536]*m_z1k;

  Box1Abb[5538]=3.*pow(Box1Abb[68],7.) + 2.*Box1Abb[5526]*pow(Box1Abb[68],4.)*m_z12 + Box1Abb[5537]*pow(Box1Abb[68],3.)*m_z12_2 + 2.*Box1Abb[5533]*pow(Box1Abb[68],2.)*m_z12_3 + 2.*Box1Abb[5529]*Box1Abb[68]*m_z12_4*m_z1k + 4.*Box1Abb[431]*m_z12_5*m_z1k_2;

  Box1Abb[5539]=-5. + Box1Abb[665]*m_z1k;

  Box1Abb[5540]=-1. + Box1Abb[3407]*m_z1k;

  Box1Abb[5541]=-1. + Box1Abb[5540]*m_z1k;

  Box1Abb[5542]=-9. + 71.*m_z1k;

  Box1Abb[5543]=-9. + Box1Abb[5542]*m_z1k;

  Box1Abb[5544]=10. + Box1Abb[5543]*m_z1k;

  Box1Abb[5545]=-31. + 98.*m_z1k;

  Box1Abb[5546]=21. + Box1Abb[5545]*m_z1k;

  Box1Abb[5547]=31. + Box1Abb[5546]*m_z1k;

  Box1Abb[5548]=5. + Box1Abb[5547]*m_z1k;

  Box1Abb[5549]=-8. + 9.*m_z1k;

  Box1Abb[5550]=95. + 7.*Box1Abb[5549]*m_z1k;

  Box1Abb[5551]=12. + Box1Abb[5550]*m_z1k;

  Box1Abb[5552]=13. + 2.*Box1Abb[5551]*m_z1k;

  Box1Abb[5553]=-8. + 2.*Box1Abb[5552]*m_z1k;

  Box1Abb[5554]=-9. + Box1Abb[5553]*m_z12 + 3.*Box1Abb[5548]*m_z12_2 + 2.*Box1Abb[5544]*m_z12_3 + 2.*Box1Abb[5539]*m_z12_4 + 15.*Box1Abb[5541]*m_z1k;

  Box1Abb[5555]=10. + 6.*m_z1k_2 - 51.*m_z1k_3;

  Box1Abb[5556]=-95. + 56.*m_z1k;

  Box1Abb[5557]=40. + Box1Abb[5556]*m_z1k;

  Box1Abb[5558]=25. - 97.*m_z1k;

  Box1Abb[5559]=15. + 2.*Box1Abb[5558]*m_z1k;

  Box1Abb[5560]=-4. + Box1Abb[5559]*m_z1k;

  Box1Abb[5561]=-35. + Box1Abb[5560]*m_z1k;

  Box1Abb[5562]=-187. + 294.*m_z1k;

  Box1Abb[5563]=241. + Box1Abb[5562]*m_z1k;

  Box1Abb[5564]=96. + Box1Abb[5563]*m_z1k;

  Box1Abb[5565]=91. + Box1Abb[5564]*m_z1k;

  Box1Abb[5566]=-19. + Box1Abb[5565]*m_z1k;

  Box1Abb[5567]=53. + 7.*Box1Abb[2312]*m_z1k;

  Box1Abb[5568]=-30. + Box1Abb[5567]*m_z1k;

  Box1Abb[5569]=-11. - 5.*Box1Abb[5568]*m_z1k;

  Box1Abb[5570]=22. + Box1Abb[5569]*m_z1k;

  Box1Abb[5571]=10. + 2.*Box1Abb[5570]*m_z1k;

  Box1Abb[5572]=-5. + Box1Abb[5571]*m_z12 - Box1Abb[5566]*m_z12_2 + Box1Abb[5561]*m_z12_3 + Box1Abb[5555]*m_z12_4 + 14.*m_z1k - 3.*Box1Abb[5557]*m_z1k_3;

  Box1Abb[5573]=-19. + 64.*m_z1k;

  Box1Abb[5574]=-12. + Box1Abb[5573]*m_z1k;

  Box1Abb[5575]=10. + Box1Abb[5574]*m_z1k;

  Box1Abb[5576]=-5. + Box1Abb[5575]*m_z1k;

  Box1Abb[5577]=-28. + 15.*m_z1k;

  Box1Abb[5578]=51. + Box1Abb[5577]*m_z1k;

  Box1Abb[5579]=-70. + Box1Abb[5578]*m_z1k;

  Box1Abb[5580]=35. + Box1Abb[5579]*m_z1k;

  Box1Abb[5581]=2. - 9.*m_z1k + 2.*Box1Abb[5580]*m_z1k_3;

  Box1Abb[5582]=-87. + 28.*m_z1k;

  Box1Abb[5583]=95. + Box1Abb[5582]*m_z1k;

  Box1Abb[5584]=-40. + Box1Abb[5583]*m_z1k;

  Box1Abb[5585]=8. + 3.*Box1Abb[5584]*m_z1k;

  Box1Abb[5586]=3. + Box1Abb[5585]*m_z1k;

  Box1Abb[5587]=-27. + 71.*m_z1k;

  Box1Abb[5588]=40. + Box1Abb[5587]*m_z1k;

  Box1Abb[5589]=40. + Box1Abb[5588]*m_z1k;

  Box1Abb[5590]=3. + Box1Abb[5589]*m_z1k;

  Box1Abb[5591]=11. + Box1Abb[5590]*m_z1k;

  Box1Abb[5592]=-41. + 42.*m_z1k;

  Box1Abb[5593]=67. + Box1Abb[5592]*m_z1k;

  Box1Abb[5594]=-217. + 5.*Box1Abb[5593]*m_z1k;

  Box1Abb[5595]=-11. + Box1Abb[5594]*m_z1k;

  Box1Abb[5596]=22. + Box1Abb[5595]*m_z1k;

  Box1Abb[5597]=-26. + Box1Abb[5596]*m_z1k;

  Box1Abb[5598]=1. + 4.*Box1Abb[5581]*m_z12 + Box1Abb[5597]*m_z12_2 + 2.*Box1Abb[5591]*m_z12_3 + Box1Abb[5576]*m_z12_4 + Box1Abb[5586]*m_z1k + 6.*m_z12_5*m_z1k_3;

  Box1Abb[5599]=-7. + 8.*m_z1k;

  Box1Abb[5600]=-1. + Box1Abb[5599]*Box1Abb[68]*m_z1k;

  Box1Abb[5601]=1. + 3.*Box1Abb[5600]*m_z1k;

  Box1Abb[5602]=6. + m_z1k + 5.*m_z1k_2;

  Box1Abb[5603]=26. + 5.*Box1Abb[5602]*m_z1k;

  Box1Abb[5604]=5. + Box1Abb[5603]*m_z1k;

  Box1Abb[5605]=-1. + Box1Abb[5604]*m_z1k;

  Box1Abb[5606]=-69. + 6.*m_z1k - 68.*m_z1k_2;

  Box1Abb[5607]=112. + Box1Abb[5606]*m_z1k;

  Box1Abb[5608]=6. + Box1Abb[5607]*m_z1k;

  Box1Abb[5609]=-6. + Box1Abb[5608]*m_z1k;

  Box1Abb[5610]=-5. + Box1Abb[5609]*m_z1k;

  Box1Abb[5611]=-17. + 36.*m_z1k;

  Box1Abb[5612]=44. + Box1Abb[5611]*m_z1k;

  Box1Abb[5613]=-218. + 3.*Box1Abb[5612]*m_z1k;

  Box1Abb[5614]=36. + Box1Abb[5613]*m_z1k;

  Box1Abb[5615]=21. + Box1Abb[5614]*m_z1k;

  Box1Abb[5616]=8. + Box1Abb[5615]*m_z1k;

  Box1Abb[5617]=-22. + 45.*m_z1k;

  Box1Abb[5618]=77. + Box1Abb[5617]*m_z1k;

  Box1Abb[5619]=-148. + Box1Abb[5618]*m_z1k;

  Box1Abb[5620]=45. + Box1Abb[5619]*m_z1k;

  Box1Abb[5621]=6. + Box1Abb[5620]*m_z1k;

  Box1Abb[5622]=5. + Box1Abb[5621]*m_z1k;

  Box1Abb[5623]=-Box1Abb[5601]*pow(Box1Abb[68],3.) - Box1Abb[5622]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[5616]*Box1Abb[68]*m_z12_2 + Box1Abb[5610]*m_z12_3 - Box1Abb[5605]*m_z12_4 - 8.*Box1Abb[68]*m_z12_5*m_z1k_3;

  Box1Abb[5624]=Box1Abb[5623]*m_x_3 + Box1Abb[5598]*m_x_4 + Box1Abb[5572]*m_x_5 + Box1Abb[5554]*m_x_6 + Box1Abb[5525]*m_x_7 + Box1Abb[5513]*m_x_8 + Box1Abb[5504]*m_x_9 + Box1Abb[5501]*m_x_10 - m_x_11*m_z12 + Box1Abb[5538]*m_x_2*m_z1k - Box1Abb[4]*Box1Abb[5507]*pow(Box1Abb[68],2.)*m_x*m_z12*m_z1k_2 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*Box1Abb[72]*m_z12_2*m_z1k_3;

  Box1Abb[5625]=Box1Abb[5624]*m_cR + Box1Abb[5500]*Box1Abb[7]*m_cL*m_x;

  Box1Abb[5626]=2. + 15.*Box1Abb[69]*m_z12;

  Box1Abb[5627]=16. + Box1Abb[5626]*m_z12;

  Box1Abb[5628]=8. + Box1Abb[488]*m_z12;

  Box1Abb[5629]=-119. + 5.*Box1Abb[5628]*m_z12;

  Box1Abb[5630]=17. + Box1Abb[5629]*m_z12;

  Box1Abb[5631]=-94. + Box1Abb[2615]*Box1Abb[707]*m_z12;

  Box1Abb[5632]=48. - Box1Abb[5631]*m_z12;

  Box1Abb[5633]=-261. + 46.*m_z12;

  Box1Abb[5634]=140. + Box1Abb[5633]*m_z12;

  Box1Abb[5635]=-95. + Box1Abb[5634]*m_z12;

  Box1Abb[5636]=20. + Box1Abb[5635]*m_z12;

  Box1Abb[5637]=-197. + 78.*m_z12;

  Box1Abb[5638]=162. + Box1Abb[5637]*m_z12;

  Box1Abb[5639]=-29. + Box1Abb[5638]*m_z12;

  Box1Abb[5640]=-56. + 39.*m_z12;

  Box1Abb[5641]=12. + Box1Abb[5640]*m_z12;

  Box1Abb[5642]=3. - Box1Abb[5627]*m_z12 - 12.*m_z1k - Box1Abb[5630]*m_z12*m_z1k + Box1Abb[5632]*m_z12*m_z1k_2 + Box1Abb[5636]*m_z1k_3 + 5.*Box1Abb[5639]*m_z1k_4 + 14.*Box1Abb[5641]*m_z1k_5 + 210.*m_z12*m_z1k_6;

  Box1Abb[5643]=-13. + 20.*m_z12;

  Box1Abb[5644]=-103. + 42.*m_z12 + 10.*m_z12_3;

  Box1Abb[5645]=50. + Box1Abb[5644]*m_z12;

  Box1Abb[5646]=-8. + Box1Abb[221]*m_z12;

  Box1Abb[5647]=-14. + Box1Abb[5646]*m_z12;

  Box1Abb[5648]=-3. + Box1Abb[5647]*m_z12;

  Box1Abb[5649]=-2. + Box1Abb[5648]*m_z12;

  Box1Abb[5650]=-32. + m_z12;

  Box1Abb[5651]=6. + Box1Abb[5650]*m_z12;

  Box1Abb[5652]=-66. + Box1Abb[5651]*m_z12;

  Box1Abb[5653]=71. + Box1Abb[5652]*m_z12;

  Box1Abb[5654]=-10. + Box1Abb[5653]*m_z12;

  Box1Abb[5655]=509. - 114.*m_z12;

  Box1Abb[5656]=-850. + Box1Abb[5655]*m_z12;

  Box1Abb[5657]=655. + Box1Abb[5656]*m_z12;

  Box1Abb[5658]=-110. + Box1Abb[5657]*m_z12;

  Box1Abb[5659]=1031. - 338.*m_z12;

  Box1Abb[5660]=-996. + Box1Abb[5659]*m_z12;

  Box1Abb[5661]=177. + Box1Abb[5660]*m_z12;

  Box1Abb[5662]=pow(Box1Abb[0],2.)*Box1Abb[5643]*m_z12 + m_z1k + Box1Abb[5645]*m_z12*m_z1k + 2.*Box1Abb[5649]*m_z1k_2 - 2.*Box1Abb[5654]*m_z1k_3 + Box1Abb[5658]*m_z1k_4 + Box1Abb[5661]*m_z1k_5 - 28.*Box1Abb[222]*Box1Abb[69]*m_z1k_6 - 120.*m_z12*m_z1k_7;

  Box1Abb[5663]=-2. + 3.*m_z12 + 5.*m_z1k;

  Box1Abb[5664]=3. + 2.*Box1Abb[5663]*m_z12;

  Box1Abb[5665]=14. + 51.*m_z1k;

  Box1Abb[5666]=-24. + Box1Abb[5665]*m_z12 - 7.*m_z12_2 + 5.*Box1Abb[5549]*m_z1k;

  Box1Abb[5667]=13. + Box1Abb[5666]*m_z12 + 24.*m_z1k;

  Box1Abb[5668]=31. + 18.*m_z1k;

  Box1Abb[5669]=61. + 84.*m_z1k;

  Box1Abb[5670]=43. + 192.*m_z1k;

  Box1Abb[5671]=22. + Box1Abb[5670]*m_z1k;

  Box1Abb[5672]=27. + 44.*m_z1k - 30.*m_z1k_2;

  Box1Abb[5673]=21. + 4.*Box1Abb[5672]*m_z1k;

  Box1Abb[5674]=-18. + Box1Abb[5673]*m_z12 - Box1Abb[5671]*m_z12_2 + Box1Abb[5668]*m_z12_3 + 2.*m_z12_4 - Box1Abb[5669]*m_z1k;

  Box1Abb[5675]=-9. + Box1Abb[4562]*m_z1k;

  Box1Abb[5676]=3. + Box1Abb[5675]*m_z1k;

  Box1Abb[5677]=-Box1Abb[5143]*pow(Box1Abb[68],3.) - Box1Abb[5142]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[5676]*Box1Abb[68]*m_z12_2 + Box1Abb[5145]*m_z12_3;

  Box1Abb[5678]=11. + 16.*m_z1k;

  Box1Abb[5679]=-81. + 40.*m_z1k;

  Box1Abb[5680]=-46. + Box1Abb[5679]*m_z1k;

  Box1Abb[5681]=31. + 56.*m_z1k;

  Box1Abb[5682]=38. + 3.*Box1Abb[5681]*m_z1k;

  Box1Abb[5683]=-61. + 420.*m_z1k;

  Box1Abb[5684]=-16. + Box1Abb[5683]*m_z1k;

  Box1Abb[5685]=50. + Box1Abb[5684]*m_z1k;

  Box1Abb[5686]=-32. + 15.*m_z1k;

  Box1Abb[5687]=-57. + 7.*Box1Abb[5686]*m_z1k;

  Box1Abb[5688]=5. + 2.*Box1Abb[5687]*m_z1k;

  Box1Abb[5689]=-1. + Box1Abb[5688]*m_z1k;

  Box1Abb[5690]=1. + Box1Abb[5689]*m_z12 + Box1Abb[5685]*m_z12_2 + Box1Abb[5680]*m_z12_3 - Box1Abb[5678]*m_z12_4 + Box1Abb[5682]*m_z1k;

  Box1Abb[5691]=25. + 3.*Box1Abb[2818]*m_z1k;

  Box1Abb[5692]=107. - 226.*m_z1k;

  Box1Abb[5693]=98. + Box1Abb[5692]*m_z1k;

  Box1Abb[5694]=10. + Box1Abb[5693]*m_z1k;

  Box1Abb[5695]=2. + m_z1k + 42.*m_z1k_2;

  Box1Abb[5696]=-1. + Box1Abb[5695]*m_z1k;

  Box1Abb[5697]=485. - 588.*m_z1k;

  Box1Abb[5698]=138. + Box1Abb[5697]*m_z1k;

  Box1Abb[5699]=-68. + Box1Abb[5698]*m_z1k;

  Box1Abb[5700]=-56. + Box1Abb[5699]*m_z1k;

  Box1Abb[5701]=-26. + 9.*m_z1k;

  Box1Abb[5702]=60. + 7.*Box1Abb[5701]*m_z1k;

  Box1Abb[5703]=111. + 4.*Box1Abb[5702]*m_z1k;

  Box1Abb[5704]=53. + Box1Abb[5703]*m_z1k;

  Box1Abb[5705]=-16. + Box1Abb[5704]*m_z1k;

  Box1Abb[5706]=6. - Box1Abb[5705]*m_z12 + Box1Abb[5700]*m_z12_2 + Box1Abb[5694]*m_z12_3 + Box1Abb[5691]*m_z12_4 - 5.*Box1Abb[5696]*m_z1k + m_z12_5*m_z1k;

  Box1Abb[5707]=6. - 7.*m_z1k;

  Box1Abb[5708]=10. + Box1Abb[5707]*m_z1k;

  Box1Abb[5709]=-10. + Box1Abb[5708]*m_z1k;

  Box1Abb[5710]=5. + Box1Abb[5709]*m_z1k;

  Box1Abb[5711]=-14. + 5.*m_z1k;

  Box1Abb[5712]=22. + Box1Abb[5711]*m_z1k;

  Box1Abb[5713]=-11. + 2.*Box1Abb[5712]*m_z1k;

  Box1Abb[5714]=-9. + Box1Abb[5713]*m_z1k;

  Box1Abb[5715]=4. + Box1Abb[5714]*m_z1k;

  Box1Abb[5716]=-25. + 6.*m_z1k;

  Box1Abb[5717]=22. + Box1Abb[5716]*m_z1k;

  Box1Abb[5718]=19. + 5.*Box1Abb[5717]*m_z1k;

  Box1Abb[5719]=-26. + Box1Abb[5718]*m_z1k;

  Box1Abb[5720]=4. + Box1Abb[5719]*m_z1k;

  Box1Abb[5721]=-59. + 18.*m_z1k;

  Box1Abb[5722]=20. + Box1Abb[5721]*m_z1k;

  Box1Abb[5723]=28. + Box1Abb[5722]*m_z1k;

  Box1Abb[5724]=-18. + Box1Abb[5723]*m_z1k;

  Box1Abb[5725]=1. + Box1Abb[5724]*m_z1k;

  Box1Abb[5726]=-101. + 30.*m_z1k;

  Box1Abb[5727]=121. + Box1Abb[5726]*m_z1k;

  Box1Abb[5728]=-8. + Box1Abb[5727]*m_z1k;

  Box1Abb[5729]=-20. + Box1Abb[5728]*m_z1k;

  Box1Abb[5730]=6. + Box1Abb[5729]*m_z1k;

  Box1Abb[5731]=-Box1Abb[4877]*pow(Box1Abb[68],6.) - Box1Abb[5715]*pow(Box1Abb[68],4.)*m_z12 - Box1Abb[5730]*pow(Box1Abb[68],3.)*m_z12_2 - Box1Abb[5720]*pow(Box1Abb[68],2.)*m_z12_3 - Box1Abb[5725]*Box1Abb[68]*m_z12_4 + Box1Abb[5710]*m_z12_5*m_z1k;

  Box1Abb[5732]=1. + Box1Abb[119]*m_z1k;

  Box1Abb[5733]=-5. + Box1Abb[1828]*m_z1k_2;

  Box1Abb[5734]=-113. + 84.*m_z1k;

  Box1Abb[5735]=-4. + Box1Abb[5734]*m_z1k;

  Box1Abb[5736]=-36. + Box1Abb[5735]*m_z1k;

  Box1Abb[5737]=30. + Box1Abb[5736]*m_z1k;

  Box1Abb[5738]=-7. + Box1Abb[5737]*m_z1k;

  Box1Abb[5739]=-315. + 152.*m_z1k;

  Box1Abb[5740]=79. + Box1Abb[5739]*m_z1k;

  Box1Abb[5741]=-5. + Box1Abb[5740]*m_z1k;

  Box1Abb[5742]=65. + Box1Abb[5741]*m_z1k;

  Box1Abb[5743]=24. + Box1Abb[5742]*Box1Abb[68]*m_z1k;

  Box1Abb[5744]=-166. + 45.*m_z1k;

  Box1Abb[5745]=265. + Box1Abb[5744]*m_z1k;

  Box1Abb[5746]=-105. + Box1Abb[5745]*m_z1k;

  Box1Abb[5747]=-12. + Box1Abb[5746]*m_z1k;

  Box1Abb[5748]=-31. + Box1Abb[5747]*m_z1k;

  Box1Abb[5749]=16. + Box1Abb[5748]*m_z1k;

  Box1Abb[5750]=-475. + 132.*m_z1k;

  Box1Abb[5751]=557. + Box1Abb[5750]*m_z1k;

  Box1Abb[5752]=-161. + Box1Abb[5751]*m_z1k;

  Box1Abb[5753]=17. + Box1Abb[5752]*m_z1k;

  Box1Abb[5754]=-64. + Box1Abb[5753]*m_z1k;

  Box1Abb[5755]=30. + Box1Abb[5754]*m_z1k;

  Box1Abb[5756]=Box1Abb[1163]*Box1Abb[5732]*pow(Box1Abb[68],4.) + Box1Abb[5749]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[5755]*Box1Abb[68]*m_z12_2 + Box1Abb[5743]*m_z12_3 + Box1Abb[5738]*m_z12_4 + 2.*Box1Abb[5733]*m_z12_5*m_z1k;

  Box1Abb[5757]=Box1Abb[5731]*m_x + Box1Abb[5756]*m_x_2 + Box1Abb[5662]*m_x_3 + Box1Abb[5642]*m_x_4 + Box1Abb[5706]*m_x_5 + Box1Abb[5690]*m_x_6 + Box1Abb[5674]*m_x_7 + Box1Abb[5667]*m_x_8 - Box1Abb[5664]*m_x_9 + m_x_10*m_z12 - Box1Abb[4]*Box1Abb[5677]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[5758]=-Box1Abb[5070]*Box1Abb[7]*m_cR + Box1Abb[5757]*m_cL*m_x;

  Box1Abb[5820]=(4.*m_cL*m_m*m_x)/(m_s_2*m_z1k);

  Box1Abb[5821]=(-4.*m_cL*m_m*m_x)/(m_s_2*m_z1k);

  Box1Abb[5822]=(4.*m_m)/(m_s_2*m_z12);

  Box1Abb[5823]=(-4.*m_m)/(m_s_2*m_z12);

  Box1Abb[5824]=-m_m;

  Box1Abb[5825]=(-4.*m_m)/m_s;

  Box1Abb[5826]=(4.*m_m)/m_s;

  Box1Abb[5827]=4.*m_m;

  Box1Abb[5828]=-4.*m_m;

  Box1Abb[5829]=(4.*m_m*m_z1k)/m_s_2;

  Box1Abb[5830]=(-4.*m_m)/(m_s_2*m_z1k);

  Box1Abb[5831]=Box1Abb[766]*m_cL + Box1Abb[2]*m_cR;

  Box1Abb[5832]=Box1Abb[79]*m_cL + 2.*Box1Abb[2]*m_cR;

  Box1Abb[5833]=Box1Abb[5831]*Box1Abb[79]*m_x + Box1Abb[5832]*m_z12*m_z1k + 2.*Box1Abb[56]*m_z12*m_z1k_2;

  Box1Abb[5834]=2.*Box1Abb[178] + 2.*Box1Abb[493]*m_z12 + Box1Abb[60]*m_z12*m_z1k + 4.*Box1Abb[8]*m_z1k_2;

  Box1Abb[5835]=2. + Box1Abb[4580]*m_z12 - 2.*m_z1k;

  Box1Abb[5836]=-14. + 9.*m_z12 - 2.*m_z12_2 + 12.*m_z1k;

  Box1Abb[5837]=8. + Box1Abb[5836]*m_z12;

  Box1Abb[5838]=-4. + 3.*m_z12 + 12.*m_z1k;

  Box1Abb[5839]=-3. + m_z12 + Box1Abb[5838]*m_z12*m_z1k;

  Box1Abb[5840]=2. + Box1Abb[5839]*m_z12 + 4.*m_z1k;

  Box1Abb[5841]=-Box1Abb[5840]*m_x_2 + Box1Abb[5837]*m_x_3 + 4.*Box1Abb[79]*m_x_4*m_z12 + Box1Abb[5834]*m_x*m_z12*m_z1k - Box1Abb[5835]*m_z12_2*m_z1k_2;

  Box1Abb[5842]=-11. + 12.*m_z12 + 48.*m_z1k;

  Box1Abb[5843]=10. + Box1Abb[5842]*m_z12 - 12.*m_z1k;

  Box1Abb[5844]=-8. + Box1Abb[5843]*m_z12;

  Box1Abb[5845]=7. - Box1Abb[3056]*m_z12 + 8.*m_z1k;

  Box1Abb[5846]=-2. + Box1Abb[5845]*m_z12 + 4.*m_z1k;

  Box1Abb[5847]=Box1Abb[5846]*m_x_2 + Box1Abb[5844]*m_x_3 + 4.*Box1Abb[810]*m_x_4*m_z12 + Box1Abb[3052]*m_x*m_z12*m_z1k - 3.*Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[5848]=Box1Abb[5841]*m_cL + Box1Abb[5847]*m_cR;

  Box1Abb[5849]=-10. + 3.*m_z12;

  Box1Abb[5850]=Box1Abb[68]*Box1Abb[692] + Box1Abb[3065]*m_z12;

  Box1Abb[5851]=9. - 2.*m_z12 + m_z1k;

  Box1Abb[5852]=-3. + Box1Abb[5851]*m_z12 + 20.*m_z1k;

  Box1Abb[5853]=-2. + Box1Abb[5852]*m_z12;

  Box1Abb[5854]=Box1Abb[5853]*m_x_3 - Box1Abb[3063]*m_x_2*m_z12 + Box1Abb[5849]*m_x_4*m_z12 + Box1Abb[5850]*m_x*m_z12_2*m_z1k - 2.*Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[5855]=10. - 13.*m_z12;

  Box1Abb[5856]=Box1Abb[3077]*m_x_2 - Box1Abb[3074]*m_x*m_z12 + Box1Abb[5855]*m_x_3*m_z12 + 3.*pow(Box1Abb[4],2.)*m_z12_2*m_z1k;

  Box1Abb[5857]=Box1Abb[5854]*m_cL + Box1Abb[5856]*m_cR*m_x;

  Box1Abb[5858]=Box1Abb[3092]*m_x_3 - Box1Abb[3088]*m_x_2*m_z12 + Box1Abb[803]*m_x_4*m_z12 + Box1Abb[3090]*m_x*m_z12_2*m_z1k - Box1Abb[367]*m_z12_3*m_z1k_2;

  Box1Abb[5859]=10. - 11.*m_z12;

  Box1Abb[5860]=-1. + 4.*m_z12 + 23.*m_z1k;

  Box1Abb[5861]=-5.*Box1Abb[174] + Box1Abb[5860]*m_z12;

  Box1Abb[5862]=2. + Box1Abb[5861]*m_z12;

  Box1Abb[5863]=Box1Abb[5862]*m_x_2 - Box1Abb[3082]*m_x*m_z12 + Box1Abb[5859]*m_x_3*m_z12 + Box1Abb[3081]*m_z12_2*m_z1k;

  Box1Abb[5864]=Box1Abb[5858]*m_cL + Box1Abb[5863]*m_cR*m_x;

  Box1Abb[5865]=-13. + 10.*m_z12;

  Box1Abb[5866]=-2. + 2.*Box1Abb[1640]*m_z12 + 6.*m_z1k + Box1Abb[5865]*m_z12*m_z1k + 4.*Box1Abb[8]*m_z1k_2;

  Box1Abb[5867]=27. - 8.*m_z12;

  Box1Abb[5868]=-26. + Box1Abb[5867]*m_z12 + 12.*m_z1k;

  Box1Abb[5869]=8. + Box1Abb[5868]*m_z12;

  Box1Abb[5870]=11. + 2.*Box1Abb[158]*m_z12 - 2.*m_z1k - 5.*m_z12*m_z1k - 12.*m_z1k_2;

  Box1Abb[5871]=-7. + Box1Abb[5870]*m_z12 + 8.*m_z1k;

  Box1Abb[5872]=2. + Box1Abb[5871]*m_z12 - 4.*m_z1k;

  Box1Abb[5873]=Box1Abb[5872]*m_x_2 + Box1Abb[5869]*m_x_3 + 4.*Box1Abb[79]*m_x_4*m_z12 + Box1Abb[5866]*m_x*m_z12*m_z1k - Box1Abb[4]*m_z12_3*m_z1k_2;

  Box1Abb[5874]=-45. + 34.*m_z12 + 48.*m_z1k;

  Box1Abb[5875]=22. + Box1Abb[5874]*m_z12 - 12.*m_z1k;

  Box1Abb[5876]=-8. + Box1Abb[5875]*m_z12;

  Box1Abb[5877]=17. - Box1Abb[3109]*m_z12;

  Box1Abb[5878]=-6. + Box1Abb[5877]*m_z12 + 4.*m_z1k;

  Box1Abb[5879]=Box1Abb[5878]*m_x_2 + Box1Abb[5876]*m_x_3 + Box1Abb[3104]*Box1Abb[4]*m_x*m_z12 + 4.*Box1Abb[810]*m_x_4*m_z12 - pow(Box1Abb[4],2.)*Box1Abb[77]*m_z12_2*m_z1k;

  Box1Abb[5880]=Box1Abb[5873]*m_cL + Box1Abb[5879]*m_cR;

  Box1Abb[5881]=pow(Box1Abb[0],2.)*Box1Abb[69] + 3.*Box1Abb[0]*Box1Abb[856]*m_z1k - 5.*Box1Abb[77]*m_z1k_2;

  Box1Abb[5882]=-Box1Abb[3117]*m_x_2 + Box1Abb[5881]*m_x*m_z12 + 5.*Box1Abb[79]*m_x_3*m_z12 + Box1Abb[3114]*Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[5883]=-15. - 13.*Box1Abb[79]*m_z12;

  Box1Abb[5884]=26. - 21.*m_z12;

  Box1Abb[5885]=2. + Box1Abb[5883]*m_z12 - 13.*m_z1k + 2.*Box1Abb[5884]*m_z12*m_z1k + 5.*Box1Abb[810]*m_z1k_2;

  Box1Abb[5886]=Box1Abb[3125]*m_x_3 + Box1Abb[5885]*m_x_2*m_z12 + 5.*Box1Abb[1010]*m_x_4*m_z12 + Box1Abb[3119]*Box1Abb[4]*m_x*m_z12_2 - pow(Box1Abb[4],2.)*m_z12_3*m_z1k;

  Box1Abb[5887]=Box1Abb[5886]*m_cR + Box1Abb[5882]*m_cL*m_x;

  Box1Abb[5888]=Box1Abb[3135]*m_x_3 - Box1Abb[3130]*m_x_2*m_z12 + Box1Abb[5855]*m_x_4*m_z12 + Box1Abb[3132]*Box1Abb[4]*m_x*m_z12_2 - 2.*pow(Box1Abb[4],2.)*m_z12_3*m_z1k;

  Box1Abb[5889]=Box1Abb[5854]*m_cL + Box1Abb[5888]*m_cR;

  Box1Abb[5890]=Box1Abb[1279]*m_cL + Box1Abb[3152]*m_cR*m_x;

  Box1Abb[5891]=1. - 3.*m_z12;

  Box1Abb[5892]=1. - m_z12 - m_z1k + 3.*m_z12*m_z1k;

  Box1Abb[5893]=2.*Box1Abb[15] + Box1Abb[2293]*m_z12;

  Box1Abb[5894]=-Box1Abb[4]*Box1Abb[5892] - Box1Abb[5893]*m_x + Box1Abb[5891]*m_x_2;

  Box1Abb[5895]=Box1Abb[3149]*m_cL + Box1Abb[5894]*m_cR*m_x;

  Box1Abb[5896]=Box1Abb[0]*pow(Box1Abb[1],2.)*m_x + Box1Abb[3144]*m_z1k - Box1Abb[3145]*m_z1k_2 + 2.*m_z12*m_z1k_3;

  Box1Abb[5897]=Box1Abb[5896]*m_cL + Box1Abb[3140]*m_cR;

  Box1Abb[5898]=Box1Abb[1277]*m_cL - Box1Abb[420]*m_cR;

  Box1Abb[5899]=-1. + Box1Abb[3776]*m_z1k;

  Box1Abb[5900]=1. - Box1Abb[873]*m_z1k;

  Box1Abb[5901]=Box1Abb[5900]*m_z12_2 + pow(Box1Abb[68],2.)*m_z1k + Box1Abb[5899]*m_z12*m_z1k;

  Box1Abb[5902]=-Box1Abb[3159]*m_x_2 + Box1Abb[3156]*m_x_3 - Box1Abb[3154]*m_x_4 - Box1Abb[0]*m_x_5 + Box1Abb[5901]*m_x*m_z1k + Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[5903]=Box1Abb[5902]*m_cR + Box1Abb[3142]*pow(Box1Abb[7],2.)*m_cL*m_x;

  Box1Abb[5904]=Box1Abb[3167]*Box1Abb[7]*m_cL + Box1Abb[3175]*m_cR;

  Box1Abb[5905]=15. - 2.*m_z12;

  Box1Abb[5906]=-15. + Box1Abb[5905]*m_z12;

  Box1Abb[5907]=2. + Box1Abb[5906]*m_z12 - 6.*m_z1k + 6.*Box1Abb[77]*m_z12*m_z1k + 4.*m_z1k_2;

  Box1Abb[5908]=-4.*Box1Abb[265] + Box1Abb[3182]*m_z12;

  Box1Abb[5909]=24. - 11.*m_z12 - 36.*m_z1k;

  Box1Abb[5910]=4.*Box1Abb[3396] + Box1Abb[5909]*m_z12;

  Box1Abb[5911]=8. + Box1Abb[5910]*m_z12;

  Box1Abb[5912]=Box1Abb[5908]*m_x_3 + Box1Abb[5911]*m_x_4 + 4.*Box1Abb[133]*m_x_5*m_z12 - Box1Abb[3189]*m_x_2*m_z1k + Box1Abb[5907]*m_x*m_z12*m_z1k_2 - 3.*Box1Abb[4]*m_z12_3*m_z1k_3;

  Box1Abb[5913]=-2.*Box1Abb[3178]*Box1Abb[384]*Box1Abb[7]*m_cL + Box1Abb[5912]*m_cR;

  Box1Abb[5914]=9. + Box1Abb[1615]*m_z12 - 12.*m_z1k + 9.*m_z12*m_z1k + 3.*m_z1k_2;

  Box1Abb[5915]=2. - Box1Abb[5914]*m_z12 - 3.*m_z1k;

  Box1Abb[5916]=Box1Abb[5915]*m_x + 3.*Box1Abb[3199]*m_x_2 + pow(Box1Abb[4],3.)*m_z12 - m_x_3*m_z12;

  Box1Abb[5917]=Box1Abb[117]*Box1Abb[3196]*m_cL + Box1Abb[5916]*m_cR*m_x;

  Box1Abb[5918]=Box1Abb[3196]*pow(Box1Abb[7],2.)*m_cL - Box1Abb[3208]*m_cR*m_x;

  Box1Abb[5919]=Box1Abb[3223]*m_cL - 2.*Box1Abb[3230]*Box1Abb[383]*m_cR;

  Box1Abb[5920]=1. + m_z12 - 5.*m_z12_2;

  Box1Abb[5921]=2. - 2.*m_z12 + m_z12_2 - m_z12_3 + 2.*Box1Abb[5920]*m_z1k - 2.*Box1Abb[77]*m_z1k_2;

  Box1Abb[5922]=-4. + Box1Abb[3238]*m_z12;

  Box1Abb[5923]=Box1Abb[5922]*m_x_2 + Box1Abb[5921]*m_x*m_z12 + 2.*Box1Abb[1869]*m_x_3*m_z12 + Box1Abb[3234]*Box1Abb[4]*m_z12_2*m_z1k;

  Box1Abb[5924]=Box1Abb[3232]*Box1Abb[384]*m_cL + Box1Abb[5923]*m_cR;

  Box1Abb[5925]=2.*Box1Abb[3249]*m_cL - Box1Abb[3262]*m_cR;

  Box1Abb[5926]=Box1Abb[117]*Box1Abb[3266]*Box1Abb[384]*m_cL - Box1Abb[3272]*Box1Abb[7]*m_cR*m_x;

  Box1Abb[5927]=Box1Abb[3266]*Box1Abb[384]*Box1Abb[7]*m_cL + Box1Abb[3326]*m_cR*m_x;

  Box1Abb[5928]=-pow(Box1Abb[0],3.)*m_z12_2 - pow(Box1Abb[0],2.)*Box1Abb[3299]*m_z1k - Box1Abb[0]*Box1Abb[3301]*m_z1k_2 + 6.*Box1Abb[133]*m_z12*m_z1k_3 + 6.*m_z12*m_z1k_4;

  Box1Abb[5929]=-89. + 40.*m_z12;

  Box1Abb[5930]=63. + Box1Abb[5929]*m_z12;

  Box1Abb[5931]=-16. + Box1Abb[5930]*m_z12;

  Box1Abb[5932]=3.*pow(Box1Abb[0],2.)*Box1Abb[133]*m_z12 + Box1Abb[5931]*m_z1k + 4.*Box1Abb[3308]*m_z1k_2 - 4.*m_z12*m_z1k_3;

  Box1Abb[5933]=Box1Abb[5932]*m_z12 + 2.*m_z1k;

  Box1Abb[5934]=Box1Abb[5933]*m_x_3 - Box1Abb[3315]*m_x_4 + Box1Abb[5928]*m_x_2*m_z12 + 6.*Box1Abb[3304]*m_x_5*m_z12 - 2.*m_x_6*m_z12_2 + Box1Abb[3298]*Box1Abb[4]*m_x*m_z12_2*m_z1k + pow(Box1Abb[4],3.)*m_z12_3*m_z1k_2;

  Box1Abb[5935]=Box1Abb[3297]*Box1Abb[7]*m_cL + 2.*Box1Abb[5934]*m_cR;

  Box1Abb[5936]=-17. + 27.*m_z12 - 11.*m_z12_2 - 2.*Box1Abb[3371]*m_z1k + Box1Abb[3375]*m_z1k_2 - 4.*Box1Abb[0]*Box1Abb[3377]*m_z1k_3 + 5.*Box1Abb[3609]*m_z12*m_z1k_4;

  Box1Abb[5937]=-4. + Box1Abb[1367]*m_z12 + 10.*m_z1k;

  Box1Abb[5938]=6. + Box1Abb[5937]*m_z12;

  Box1Abb[5939]=31. - 14.*m_z12 + 2.*Box1Abb[3379]*m_z1k - 3.*Box1Abb[653]*m_z1k_2;

  Box1Abb[5940]=-23. + Box1Abb[5939]*m_z12 - 26.*m_z1k;

  Box1Abb[5941]=-32. + 5.*m_z1k;

  Box1Abb[5942]=-24. + Box1Abb[5941]*m_z1k;

  Box1Abb[5943]=-47. + 2.*Box1Abb[5942]*m_z1k;

  Box1Abb[5944]=31. + Box1Abb[5943]*m_z12 + Box1Abb[3386]*m_z12_2 + Box1Abb[3384]*m_z1k - 2.*Box1Abb[2406]*m_z12_3*m_z1k;

  Box1Abb[5945]=Box1Abb[3419]*m_x_2 + Box1Abb[5936]*m_x_3 + Box1Abb[5944]*m_x_4 + Box1Abb[5940]*m_x_5 + Box1Abb[5938]*m_x_6 + Box1Abb[79]*m_x_7*m_z12 - Box1Abb[3406]*m_x*m_z1k + Box1Abb[3394]*pow(Box1Abb[68],2.)*m_z12*m_z1k_2;

  Box1Abb[5946]=10. + Box1Abb[149]*m_z12;

  Box1Abb[5947]=26. - 5.*m_z12;

  Box1Abb[5948]=-60. + Box1Abb[5947]*m_z12;

  Box1Abb[5949]=14. + Box1Abb[5948]*m_z12;

  Box1Abb[5950]=-11. + Box1Abb[5949]*m_z12;

  Box1Abb[5951]=13. - 2.*Box1Abb[5946]*m_z12 + 10.*m_z1k + 2.*Box1Abb[3329]*m_z12*m_z1k + Box1Abb[5950]*m_z1k_2 - 4.*Box1Abb[0]*Box1Abb[3336]*m_z1k_3 - 5.*Box1Abb[3337]*m_z12*m_z1k_4;

  Box1Abb[5952]=-6. + Box1Abb[3356]*m_z12;

  Box1Abb[5953]=47. - 5.*m_z1k;

  Box1Abb[5954]=27. + Box1Abb[5953]*m_z1k;

  Box1Abb[5955]=37. + 2.*Box1Abb[5954]*m_z1k;

  Box1Abb[5956]=-32. + 115.*m_z1k;

  Box1Abb[5957]=20. + Box1Abb[5956]*m_z1k;

  Box1Abb[5958]=-5. + Box1Abb[5957]*m_z1k;

  Box1Abb[5959]=-26. + Box1Abb[5955]*m_z12 + Box1Abb[5958]*m_z12_2 + Box1Abb[1702]*m_z12_3 - 11.*Box1Abb[233]*m_z1k;

  Box1Abb[5960]=-1. + 14.*m_z1k;

  Box1Abb[5961]=4. + Box1Abb[5960]*m_z1k;

  Box1Abb[5962]=-6. + Box1Abb[3638]*m_z1k;

  Box1Abb[5963]=29. + 2.*Box1Abb[5962]*m_z1k;

  Box1Abb[5964]=-7. + Box1Abb[5963]*m_z1k;

  Box1Abb[5965]=-34. + 15.*m_z1k;

  Box1Abb[5966]=86. + 3.*Box1Abb[5965]*m_z1k;

  Box1Abb[5967]=-37. + Box1Abb[5966]*m_z1k;

  Box1Abb[5968]=45. + Box1Abb[5967]*m_z1k;

  Box1Abb[5969]=-5. + Box1Abb[5968]*m_z1k;

  Box1Abb[5970]=-Box1Abb[5961]*pow(Box1Abb[68],2.) + Box1Abb[5964]*Box1Abb[68]*m_z12 + Box1Abb[5969]*m_z12_2 + pow(Box1Abb[1723],2.)*Box1Abb[1725]*m_z12_3 + Box1Abb[3583]*m_z12_4*m_z1k_2;

  Box1Abb[5971]=Box1Abb[5970]*m_x_2 + Box1Abb[5951]*m_x_3 + Box1Abb[5959]*m_x_4 - Box1Abb[3355]*m_x_5 + Box1Abb[5952]*m_x_6 + Box1Abb[810]*m_x_7*m_z12 - Box1Abb[3368]*Box1Abb[68]*m_x*m_z1k + Box1Abb[3351]*Box1Abb[4]*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[5972]=Box1Abb[5945]*m_cL + Box1Abb[5971]*m_cR;

  Box1Abb[5973]=-11. + Box1Abb[3453]*m_z12 - 24.*m_z1k;

  Box1Abb[5974]=15. + 22.*m_z1k - 68.*m_z1k_2;

  Box1Abb[5975]=-6. + Box1Abb[5974]*m_z1k;

  Box1Abb[5976]=8. - Box1Abb[3471]*m_z1k;

  Box1Abb[5977]=2.*Box1Abb[3469] + Box1Abb[5976]*m_z12 + Box1Abb[5975]*m_z12_2 - 5.*m_z12_3*m_z1k_2;

  Box1Abb[5978]=-3. + m_z1k - 40.*m_z1k_2;

  Box1Abb[5979]=-31. + 115.*m_z1k;

  Box1Abb[5980]=-9. + Box1Abb[5979]*m_z1k;

  Box1Abb[5981]=-5. + Box1Abb[5980]*m_z1k;

  Box1Abb[5982]=1. + Box1Abb[5981]*m_z12 + 2.*Box1Abb[3458]*m_z12_2 + 2.*Box1Abb[5978]*m_z1k;

  Box1Abb[5983]=-Box1Abb[3475]*pow(Box1Abb[68],2.) + Box1Abb[3478]*Box1Abb[68]*m_z12 + 2.*Box1Abb[3482]*m_z12_2 + Box1Abb[3583]*m_z12_3*m_z1k_2;

  Box1Abb[5984]=Box1Abb[5983]*m_x_2 + Box1Abb[5977]*m_x_3 + Box1Abb[5982]*m_x_4 - Box1Abb[3457]*m_x_5 + Box1Abb[5973]*m_x_6 + Box1Abb[3521]*m_x_7 - Box1Abb[3452]*Box1Abb[4]*Box1Abb[68]*m_x*m_z1k + pow(Box1Abb[4],2.)*pow(Box1Abb[68],3.)*m_z12*m_z1k_2;

  Box1Abb[5985]=Box1Abb[3761]*m_cL + Box1Abb[5984]*m_cR;

  Box1Abb[5986]=Box1Abb[3520]*m_cL + Box1Abb[3552]*m_cR;

  Box1Abb[5987]=-1. + Box1Abb[1737]*Box1Abb[8]*m_z12;

  Box1Abb[5988]=63. + 16.*m_z12;

  Box1Abb[5989]=-18. + Box1Abb[5988]*m_z12;

  Box1Abb[5990]=-31. + 21.*m_z12;

  Box1Abb[5991]=62. + Box1Abb[5990]*m_z12;

  Box1Abb[5992]=-22. + Box1Abb[5991]*m_z12;

  Box1Abb[5993]=-2. + 23.*m_z12;

  Box1Abb[5994]=-3. + Box1Abb[5987]*m_z12 - 19.*m_z1k + Box1Abb[5989]*m_z12*m_z1k + 2.*Box1Abb[5992]*m_z1k_2 + 5.*Box1Abb[5993]*m_z12*m_z1k_3;

  Box1Abb[5995]=-6. + Box1Abb[3567]*m_z12;

  Box1Abb[5996]=-Box1Abb[3586]*Box1Abb[4]*pow(Box1Abb[68],2.)*m_x + Box1Abb[3600]*m_x_2 - Box1Abb[3565]*m_x_3 + Box1Abb[5994]*m_x_4 - Box1Abb[3580]*m_x_5 + Box1Abb[5995]*m_x_6 + Box1Abb[810]*m_x_7*m_z12 + pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*Box1Abb[77]*m_z12*m_z1k;

  Box1Abb[5997]=Box1Abb[3647]*m_cL + Box1Abb[5996]*m_cR;

  Box1Abb[5998]=Box1Abb[3707]*m_cL + Box1Abb[3678]*m_cR;

  Box1Abb[5999]=Box1Abb[3761]*m_cL + Box1Abb[3745]*m_cR;

  Box1Abb[6000]=Box1Abb[418]*m_cR*m_x - Box1Abb[416]*Box1Abb[7]*m_cL*m_z1k;

  Box1Abb[6001]=Box1Abb[3809]*m_cR + Box1Abb[3782]*m_cL*m_z1k;

  Box1Abb[6002]=Box1Abb[3849]*Box1Abb[7]*m_cR + Box1Abb[3831]*m_cL*m_z1k;

  Box1Abb[6003]=-Box1Abb[3888]*Box1Abb[7]*m_cL + Box1Abb[3874]*m_cR;

  Box1Abb[6004]=Box1Abb[3967]*m_cR + Box1Abb[3914]*Box1Abb[7]*m_cL*m_z1k;

  Box1Abb[6005]=Box1Abb[4039]*m_cR + Box1Abb[3998]*m_cL*m_z1k;

  Box1Abb[6006]=Box1Abb[4073]*m_cR + Box1Abb[3831]*m_cL*m_z1k;

  Box1Abb[6007]=-4. + 7.*m_z12 + 6.*m_z12*m_z1k + 4.*Box1Abb[0]*m_z1k_2;

  Box1Abb[6008]=2.*pow(Box1Abb[68],2.) - Box1Abb[4076]*Box1Abb[68]*m_z12 - 2.*Box1Abb[178]*m_z12_2*m_z1k;

  Box1Abb[6009]=Box1Abb[6008]*m_x + Box1Abb[6007]*m_x_2 - Box1Abb[743]*m_x_3 + 2.*Box1Abb[77]*m_x_4 + 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[6010]=-6. + Box1Abb[1868]*m_z12;

  Box1Abb[6011]=2.*Box1Abb[4081]*m_z12 - 3.*m_z12_2 + 4.*Box1Abb[881]*m_z1k;

  Box1Abb[6012]=1. - Box1Abb[1945]*m_z1k;

  Box1Abb[6013]=-3. + Box1Abb[4331]*m_z1k;

  Box1Abb[6014]=2.*pow(Box1Abb[68],2.) + Box1Abb[6013]*m_z12 + Box1Abb[6012]*m_z12_2;

  Box1Abb[6015]=Box1Abb[6014]*m_x + Box1Abb[6011]*m_x_2 + Box1Abb[6010]*m_x_3 - 2.*Box1Abb[141]*m_x_4 + 2.*Box1Abb[4]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[6016]=Box1Abb[6009]*m_cL + Box1Abb[6015]*m_cR;

  Box1Abb[6017]=2. - 11.*m_z1k;

  Box1Abb[6018]=-3. + Box1Abb[6017]*m_z1k;

  Box1Abb[6019]=2.*Box1Abb[317]*pow(Box1Abb[68],4.) - Box1Abb[4102]*pow(Box1Abb[68],3.)*m_z12 + Box1Abb[4105]*pow(Box1Abb[68],2.)*m_z12_2 - Box1Abb[1827]*Box1Abb[68]*m_z12_3 + Box1Abb[6018]*m_z12_4*m_z1k;

  Box1Abb[6020]=-55. + Box1Abb[4108]*m_z1k;

  Box1Abb[6021]=60. - Box1Abb[4112]*m_z12 + Box1Abb[6020]*m_z12_2 + 2.*Box1Abb[1172]*m_z12_3 + 4.*Box1Abb[1378]*m_z1k;

  Box1Abb[6022]=21. + 49.*m_z1k - 116.*m_z1k_3;

  Box1Abb[6023]=7. + 30.*m_z1k;

  Box1Abb[6024]=-4. + Box1Abb[6023]*m_z1k;

  Box1Abb[6025]=-2. + Box1Abb[4988]*m_z1k;

  Box1Abb[6026]=-3. + Box1Abb[6025]*m_z1k;

  Box1Abb[6027]=37. - 84.*m_z1k + 40.*m_z1k_2;

  Box1Abb[6028]=11. + 2.*Box1Abb[6027]*m_z1k;

  Box1Abb[6029]=9. + Box1Abb[6028]*m_z1k;

  Box1Abb[6030]=8.*Box1Abb[6026] + Box1Abb[6029]*m_z12 + Box1Abb[6022]*m_z12_2 - 3.*Box1Abb[6024]*m_z12_3 - m_z12_4*m_z1k;

  Box1Abb[6031]=13. - 5.*m_z1k + 80.*m_z1k_2;

  Box1Abb[6032]=-6. + Box1Abb[6031]*m_z1k;

  Box1Abb[6033]=-4.*Box1Abb[4114]*pow(Box1Abb[68],2.) - Box1Abb[4122]*Box1Abb[68]*m_z12 + Box1Abb[4118]*m_z12_2 + Box1Abb[6032]*m_z12_3 + Box1Abb[1843]*m_z12_4*m_z1k;

  Box1Abb[6034]=Box1Abb[6019]*m_x + Box1Abb[6033]*m_x_2 + Box1Abb[6030]*m_x_3 + Box1Abb[6021]*m_x_4 + Box1Abb[4100]*m_x_5 - Box1Abb[4096]*m_x_6 - Box1Abb[4]*Box1Abb[4087]*pow(Box1Abb[68],3.)*m_z12*m_z1k;

  Box1Abb[6035]=Box1Abb[4171]*m_cL + Box1Abb[6034]*m_cR*m_x;

  Box1Abb[6036]=2.*Box1Abb[4288]*m_cR - Box1Abb[4217]*Box1Abb[7]*m_cL*m_x;

  Box1Abb[6037]=Box1Abb[4462]*m_cR + Box1Abb[4341]*pow(Box1Abb[7],2.)*m_cL*m_x;

  Box1Abb[6038]=16. - Box1Abb[4553]*Box1Abb[68]*m_z1k;

  Box1Abb[6039]=-124. + 39.*Box1Abb[267]*m_z1k;

  Box1Abb[6040]=-96. + Box1Abb[6039]*m_z1k;

  Box1Abb[6041]=-145. + Box1Abb[6040]*m_z1k;

  Box1Abb[6042]=8.*Box1Abb[15]*pow(Box1Abb[68],2.) + 2.*Box1Abb[6038]*m_z12 + Box1Abb[6041]*m_z12_2 + Box1Abb[4557]*m_z12_3 + 4.*Box1Abb[4552]*m_z12_4 + 3.*m_z12_5*m_z1k;

  Box1Abb[6043]=-Box1Abb[4577]*m_x_3 + Box1Abb[6042]*m_x_4 - Box1Abb[4550]*m_x_5 + Box1Abb[4526]*m_x_6 - Box1Abb[4520]*m_x_7 + Box1Abb[4541]*m_x_2*m_z12 + Box1Abb[3318]*m_x_8*m_z12 - Box1Abb[4]*Box1Abb[4529]*pow(Box1Abb[68],2.)*m_x*m_z12_2 + 2.*pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*m_z12_3*m_z1k;

  Box1Abb[6044]=Box1Abb[4517]*m_cL + Box1Abb[6043]*m_cR;

  Box1Abb[6045]=Box1Abb[4517]*m_cL - Box1Abb[4629]*m_cR*m_x;

  Box1Abb[6046]=22. - Box1Abb[4632]*m_z12;

  Box1Abb[6047]=4. + Box1Abb[6046]*m_z12;

  Box1Abb[6048]=12. - 2.*Box1Abb[0]*Box1Abb[1621]*m_z12 + 8.*m_z1k + Box1Abb[6047]*m_z12*m_z1k + 2.*Box1Abb[4638]*m_z1k_2 + 5.*Box1Abb[1112]*m_z12*m_z1k_3;

  Box1Abb[6049]=1. + Box1Abb[1039]*m_z1k;

  Box1Abb[6050]=-Box1Abb[4645]*pow(Box1Abb[68],4.) - Box1Abb[4647]*pow(Box1Abb[68],2.)*m_z12 - 2.*Box1Abb[4649]*Box1Abb[68]*m_z12_2 + Box1Abb[6049]*m_z12_3*m_z1k;

  Box1Abb[6051]=-9. - Box1Abb[4651]*m_z12 + Box1Abb[234]*m_z12_2 - Box1Abb[1900]*m_z1k;

  Box1Abb[6052]=-8.*Box1Abb[15] + Box1Abb[6051]*m_z12;

  Box1Abb[6053]=1. - 7.*m_z1k;

  Box1Abb[6054]=2. - 7.*m_z1k;

  Box1Abb[6055]=2. + 3.*Box1Abb[6054]*m_z1k;

  Box1Abb[6056]=-2. + Box1Abb[6055]*m_z1k;

  Box1Abb[6057]=-46. + 33.*m_z1k;

  Box1Abb[6058]=-46. + Box1Abb[6057]*m_z1k;

  Box1Abb[6059]=108. + Box1Abb[6058]*m_z1k;

  Box1Abb[6060]=-51. + Box1Abb[6059]*m_z1k;

  Box1Abb[6061]=2. + Box1Abb[6060]*m_z1k_2;

  Box1Abb[6062]=Box1Abb[317]*Box1Abb[532]*pow(Box1Abb[68],4.) + Box1Abb[6061]*m_z12 + Box1Abb[15]*Box1Abb[4668]*m_z12_2 + 2.*Box1Abb[6056]*m_z12_3*m_z1k + Box1Abb[6053]*m_z12_4*m_z1k_2;

  Box1Abb[6063]=-Box1Abb[4688]*m_x_3 - Box1Abb[4662]*m_x_4 + Box1Abb[6048]*m_x_5 + Box1Abb[6052]*m_x_6 - Box1Abb[4643]*m_x_7 + Box1Abb[6062]*m_x_2*m_z12 + 3.*Box1Abb[79]*m_x_8*m_z12 + Box1Abb[6050]*m_x*m_z12_2*m_z1k + pow(Box1Abb[4],2.)*Box1Abb[457]*pow(Box1Abb[68],2.)*m_z12_3*m_z1k_2;

  Box1Abb[6064]=Box1Abb[6063]*m_cL - Box1Abb[4746]*m_cR*m_x;

  Box1Abb[6065]=Box1Abb[4782]*m_cL + Box1Abb[4763]*m_cR;

  Box1Abb[6066]=2.*Box1Abb[4854]*m_cL + Box1Abb[4936]*m_cR;

  Box1Abb[6067]=Box1Abb[5070]*m_cL - Box1Abb[4983]*Box1Abb[7]*m_cR*m_x;

  Box1Abb[6068]=-Box1Abb[5116]*Box1Abb[7]*m_cL + Box1Abb[5220]*m_cR*m_x;

  Box1Abb[6069]=Box1Abb[5311]*Box1Abb[7]*m_cL - 2.*Box1Abb[5421]*m_cR;

  Box1Abb[6070]=-4. + 45.*m_z1k;

  Box1Abb[6071]=-2. + 5.*Box1Abb[532]*m_z12 - 3.*m_z12_2 + Box1Abb[6070]*m_z1k;

  Box1Abb[6072]=9. + Box1Abb[6071]*m_z12 + 24.*m_z1k;

  Box1Abb[6073]=-3.*Box1Abb[174]*Box1Abb[3728] - 8.*Box1Abb[5512]*m_z12 - Box1Abb[5509]*m_z12_2 + 2.*Box1Abb[2780]*m_z12_3 + m_z12_4;

  Box1Abb[6074]=-9. + 22.*m_z1k;

  Box1Abb[6075]=-5. + 2.*Box1Abb[6074]*m_z1k;

  Box1Abb[6076]=-4. + Box1Abb[5516]*m_z1k;

  Box1Abb[6077]=-5. + Box1Abb[6076]*m_z1k;

  Box1Abb[6078]=3.*Box1Abb[6077] + 2.*Box1Abb[5524]*m_z12 + Box1Abb[5521]*m_z12_2 + Box1Abb[6075]*m_z12_3 - 5.*Box1Abb[15]*m_z12_4;

  Box1Abb[6079]=-3.*pow(Box1Abb[68],7.) - 2.*Box1Abb[5526]*pow(Box1Abb[68],4.)*m_z12 - Box1Abb[5537]*pow(Box1Abb[68],3.)*m_z12_2 - 2.*Box1Abb[5533]*pow(Box1Abb[68],2.)*m_z12_3 - 2.*Box1Abb[5529]*Box1Abb[68]*m_z12_4*m_z1k - 4.*Box1Abb[431]*m_z12_5*m_z1k_2;

  Box1Abb[6080]=5. + Box1Abb[1837]*m_z1k;

  Box1Abb[6081]=1. + m_z1k + 9.*m_z1k_2 - 14.*m_z1k_3;

  Box1Abb[6082]=9. - 71.*m_z1k;

  Box1Abb[6083]=9. + Box1Abb[6082]*m_z1k;

  Box1Abb[6084]=-10. + Box1Abb[6083]*m_z1k;

  Box1Abb[6085]=-4. + Box1Abb[5552]*m_z1k;

  Box1Abb[6086]=9. - 2.*Box1Abb[6085]*m_z12 - 3.*Box1Abb[5548]*m_z12_2 + 2.*Box1Abb[6084]*m_z12_3 + 2.*Box1Abb[6080]*m_z12_4 + 15.*Box1Abb[6081]*m_z1k;

  Box1Abb[6087]=-10. - 6.*m_z1k_2 + 51.*m_z1k_3;

  Box1Abb[6088]=-25. + 97.*m_z1k;

  Box1Abb[6089]=-15. + 2.*Box1Abb[6088]*m_z1k;

  Box1Abb[6090]=4. + Box1Abb[6089]*m_z1k;

  Box1Abb[6091]=35. + Box1Abb[6090]*m_z1k;

  Box1Abb[6092]=11. + 5.*Box1Abb[5568]*m_z1k;

  Box1Abb[6093]=-22. + Box1Abb[6092]*m_z1k;

  Box1Abb[6094]=-5. + Box1Abb[6093]*m_z1k;

  Box1Abb[6095]=5. + 2.*Box1Abb[6094]*m_z12 + Box1Abb[5566]*m_z12_2 + Box1Abb[6091]*m_z12_3 + Box1Abb[6087]*m_z12_4 - 14.*m_z1k + 3.*Box1Abb[5557]*m_z1k_3;

  Box1Abb[6096]=19. - 64.*m_z1k;

  Box1Abb[6097]=12. + Box1Abb[6096]*m_z1k;

  Box1Abb[6098]=-10. + Box1Abb[6097]*m_z1k;

  Box1Abb[6099]=5. + Box1Abb[6098]*m_z1k;

  Box1Abb[6100]=87. - 28.*m_z1k;

  Box1Abb[6101]=-95. + Box1Abb[6100]*m_z1k;

  Box1Abb[6102]=40. + Box1Abb[6101]*m_z1k;

  Box1Abb[6103]=-8. + 3.*Box1Abb[6102]*m_z1k;

  Box1Abb[6104]=217. - 5.*Box1Abb[5593]*m_z1k;

  Box1Abb[6105]=11. + Box1Abb[6104]*m_z1k;

  Box1Abb[6106]=-22. + Box1Abb[6105]*m_z1k;

  Box1Abb[6107]=26. + Box1Abb[6106]*m_z1k;

  Box1Abb[6108]=-1. - 4.*Box1Abb[5581]*m_z12 + Box1Abb[6107]*m_z12_2 - 2.*Box1Abb[5591]*m_z12_3 + Box1Abb[6099]*m_z12_4 - 3.*m_z1k + Box1Abb[6103]*m_z1k_2 - 6.*m_z12_5*m_z1k_3;

  Box1Abb[6109]=69. - 6.*m_z1k + 68.*m_z1k_2;

  Box1Abb[6110]=-112. + Box1Abb[6109]*m_z1k;

  Box1Abb[6111]=-6. + Box1Abb[6110]*m_z1k;

  Box1Abb[6112]=6. + Box1Abb[6111]*m_z1k;

  Box1Abb[6113]=5. + Box1Abb[6112]*m_z1k;

  Box1Abb[6114]=Box1Abb[5601]*pow(Box1Abb[68],3.) + Box1Abb[5622]*pow(Box1Abb[68],2.)*m_z12 + Box1Abb[5616]*Box1Abb[68]*m_z12_2 + Box1Abb[6113]*m_z12_3 + Box1Abb[5605]*m_z12_4 + 8.*Box1Abb[68]*m_z12_5*m_z1k_3;

  Box1Abb[6115]=Box1Abb[6114]*m_x_3 + Box1Abb[6108]*m_x_4 + Box1Abb[6095]*m_x_5 + Box1Abb[6086]*m_x_6 + Box1Abb[6078]*m_x_7 + Box1Abb[6073]*m_x_8 + Box1Abb[6072]*m_x_9 - Box1Abb[5501]*m_x_10 + m_x_11*m_z12 + Box1Abb[6079]*m_x_2*m_z1k + Box1Abb[4]*Box1Abb[5507]*pow(Box1Abb[68],2.)*m_x*m_z12*m_z1k_2 - pow(Box1Abb[4],2.)*pow(Box1Abb[68],4.)*Box1Abb[72]*m_z12_2*m_z1k_3;

  Box1Abb[6116]=Box1Abb[6115]*m_cL - Box1Abb[5500]*Box1Abb[7]*m_cR*m_x;

  Box1Abb[6117]=-48. + Box1Abb[5631]*m_z12;

  Box1Abb[6118]=261. - 46.*m_z12;

  Box1Abb[6119]=-140. + Box1Abb[6118]*m_z12;

  Box1Abb[6120]=95. + Box1Abb[6119]*m_z12;

  Box1Abb[6121]=-20. + Box1Abb[6120]*m_z12;

  Box1Abb[6122]=197. - 78.*m_z12;

  Box1Abb[6123]=-162. + Box1Abb[6122]*m_z12;

  Box1Abb[6124]=29. + Box1Abb[6123]*m_z12;

  Box1Abb[6125]=-3. + Box1Abb[5627]*m_z12 + 12.*m_z1k + Box1Abb[5630]*m_z12*m_z1k + Box1Abb[6117]*m_z12*m_z1k_2 + Box1Abb[6121]*m_z1k_3 + 5.*Box1Abb[6124]*m_z1k_4 - 14.*Box1Abb[5641]*m_z1k_5 - 210.*m_z12*m_z1k_6;

  Box1Abb[6126]=1. + Box1Abb[5645]*m_z12;

  Box1Abb[6127]=8. - Box1Abb[221]*m_z12;

  Box1Abb[6128]=14. + Box1Abb[6127]*m_z12;

  Box1Abb[6129]=3. + Box1Abb[6128]*m_z12;

  Box1Abb[6130]=2. + Box1Abb[6129]*m_z12;

  Box1Abb[6131]=-509. + 114.*m_z12;

  Box1Abb[6132]=850. + Box1Abb[6131]*m_z12;

  Box1Abb[6133]=-655. + Box1Abb[6132]*m_z12;

  Box1Abb[6134]=110. + Box1Abb[6133]*m_z12;

  Box1Abb[6135]=-1031. + 338.*m_z12;

  Box1Abb[6136]=996. + Box1Abb[6135]*m_z12;

  Box1Abb[6137]=-177. + Box1Abb[6136]*m_z12;

  Box1Abb[6138]=-pow(Box1Abb[0],2.)*Box1Abb[5643]*m_z12 - Box1Abb[6126]*m_z1k + 2.*Box1Abb[6130]*m_z1k_2 + 2.*Box1Abb[5654]*m_z1k_3 + Box1Abb[6134]*m_z1k_4 + Box1Abb[6137]*m_z1k_5 + 28.*Box1Abb[222]*Box1Abb[69]*m_z1k_6 + 120.*m_z12*m_z1k_7;

  Box1Abb[6139]=8. - 9.*m_z1k;

  Box1Abb[6140]=24. - Box1Abb[5665]*m_z12 + 7.*m_z12_2 + 5.*Box1Abb[6139]*m_z1k;

  Box1Abb[6141]=-13. + Box1Abb[6140]*m_z12 - 24.*m_z1k;

  Box1Abb[6142]=81. - 40.*m_z1k;

  Box1Abb[6143]=46. + Box1Abb[6142]*m_z1k;

  Box1Abb[6144]=61. - 420.*m_z1k;

  Box1Abb[6145]=16. + Box1Abb[6144]*m_z1k;

  Box1Abb[6146]=-50. + Box1Abb[6145]*m_z1k;

  Box1Abb[6147]=32. - 15.*m_z1k;

  Box1Abb[6148]=57. + 7.*Box1Abb[6147]*m_z1k;

  Box1Abb[6149]=-5. + 2.*Box1Abb[6148]*m_z1k;

  Box1Abb[6150]=-1. + m_z12 + Box1Abb[6146]*m_z12_2 + Box1Abb[6143]*m_z12_3 + Box1Abb[5678]*m_z12_4 - Box1Abb[5682]*m_z1k + Box1Abb[6149]*m_z12*m_z1k;

  Box1Abb[6151]=-27. - 44.*m_z1k + 30.*m_z1k_2;

  Box1Abb[6152]=-21. + 4.*Box1Abb[6151]*m_z1k;

  Box1Abb[6153]=18. + Box1Abb[6152]*m_z12 + Box1Abb[5671]*m_z12_2 - Box1Abb[5668]*m_z12_3 - 2.*m_z12_4 + Box1Abb[5669]*m_z1k;

  Box1Abb[6154]=-107. + 226.*m_z1k;

  Box1Abb[6155]=-98. + Box1Abb[6154]*m_z1k;

  Box1Abb[6156]=-10. + Box1Abb[6155]*m_z1k;

  Box1Abb[6157]=-485. + 588.*m_z1k;

  Box1Abb[6158]=-138. + Box1Abb[6157]*m_z1k;

  Box1Abb[6159]=68. + Box1Abb[6158]*m_z1k;

  Box1Abb[6160]=56. + Box1Abb[6159]*m_z1k;

  Box1Abb[6161]=-6. + Box1Abb[5705]*m_z12 + Box1Abb[6160]*m_z12_2 + Box1Abb[6156]*m_z12_3 - Box1Abb[5691]*m_z12_4 + 5.*Box1Abb[5696]*m_z1k - m_z12_5*m_z1k;

  Box1Abb[6162]=-6. + 7.*m_z1k;

  Box1Abb[6163]=-10. + Box1Abb[6162]*m_z1k;

  Box1Abb[6164]=10. + Box1Abb[6163]*m_z1k;

  Box1Abb[6165]=-5. + Box1Abb[6164]*m_z1k;

  Box1Abb[6166]=Box1Abb[4877]*pow(Box1Abb[68],6.) + Box1Abb[5715]*pow(Box1Abb[68],4.)*m_z12 + Box1Abb[5730]*pow(Box1Abb[68],3.)*m_z12_2 + Box1Abb[5720]*pow(Box1Abb[68],2.)*m_z12_3 + Box1Abb[5725]*Box1Abb[68]*m_z12_4 + Box1Abb[6165]*m_z12_5*m_z1k;

  Box1Abb[6167]=5. + Box1Abb[5021]*m_z1k_2;

  Box1Abb[6168]=113. - 84.*m_z1k;

  Box1Abb[6169]=4. + Box1Abb[6168]*m_z1k;

  Box1Abb[6170]=36. + Box1Abb[6169]*m_z1k;

  Box1Abb[6171]=-30. + Box1Abb[6170]*m_z1k;

  Box1Abb[6172]=7. + Box1Abb[6171]*m_z1k;

  Box1Abb[6173]=-24. - Box1Abb[5742]*Box1Abb[68]*m_z1k;

  Box1Abb[6174]=-Box1Abb[1163]*Box1Abb[5732]*pow(Box1Abb[68],4.) - Box1Abb[5749]*pow(Box1Abb[68],2.)*m_z12 - Box1Abb[5755]*Box1Abb[68]*m_z12_2 + Box1Abb[6173]*m_z12_3 + Box1Abb[6172]*m_z12_4 + 2.*Box1Abb[6167]*m_z12_5*m_z1k;

  Box1Abb[6175]=Box1Abb[6166]*m_x + Box1Abb[6174]*m_x_2 + Box1Abb[6138]*m_x_3 + Box1Abb[6125]*m_x_4 + Box1Abb[6161]*m_x_5 + Box1Abb[6150]*m_x_6 + Box1Abb[6153]*m_x_7 + Box1Abb[6141]*m_x_8 + Box1Abb[5664]*m_x_9 - m_x_10*m_z12 + Box1Abb[4]*Box1Abb[5677]*pow(Box1Abb[68],2.)*m_z12*m_z1k;

  Box1Abb[6176]=Box1Abb[5070]*Box1Abb[7]*m_cL + Box1Abb[6175]*m_cR*m_x;

  Box1Abb[6237]=(4.*m_cR*m_m*m_x)/(m_s_2*m_z1k);

  Box1Abb[6238]=(-4.*m_cR*m_m*m_x)/(m_s_2*m_z1k);

  Box1Abb[6239]=(4.*m_m*m_x)/m_s_2;

  Box1Abb[6240]=-m_m*m_s;

  return;
}

DivArrC Z_Decay_RV_Box_1::RV_Box_1(const int& ME, const int& LR,
				   const Vec4C& epsV, const Vec4C& epsP)
{
  // Box corrections for emission off leg 1
  // ME takes values 1..8, denoting m_standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return One
	*((Box1Abb[295]*Box1Abb[3138]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (epsP*m_p1)*((Box1Abb[3153]*Box1Abb[586]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[3163]*Box1Abb[5823]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)) + (Box1Abb[3143]*Box1Abb[5823]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3])) + (epsP*m_p2)*((Box1Abb[3143]*Box1Abb[5823]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[3150]*Box1Abb[5823]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3])))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((Box1Abb[422]*Box1Abb[5820]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[3039]*Box1Abb[5821]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[3041]*Box1Abb[93]*(epsP*epsV))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[296]*Box1Abb[3060]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[390]) + (epsP*m_p2)*((Box1Abb[3079]*Box1Abb[5822]*(epsV*m_p1))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[3094]*Box1Abb[5822]*(epsV*m_p2))/(pow(Box1Abb[3],2.)*Box1Abb[390])) + (epsP*m_p1)*((Box1Abb[296]*Box1Abb[3113]*(epsV*m_pP))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[3127]*Box1Abb[5822]*(epsV*m_p1))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[3137]*Box1Abb[5822]*(epsV*m_p2))/(pow(Box1Abb[3],2.)*Box1Abb[390])))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[1324]*Box1Abb[3138]*Box1Abb[92]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[295]*Box1Abb[3421]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (epsP*m_p1)*((Box1Abb[295]*Box1Abb[3648]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[295]*Box1Abb[3708]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[295]*Box1Abb[3762]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p2)*((Box1Abb[0]*Box1Abb[295]*Box1Abb[3485]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[3553]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58])))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[3763]*Box1Abb[93]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[296]*Box1Abb[3810]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p1)*((Box1Abb[3968]*Box1Abb[5830]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[4040]*Box1Abb[5830]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[296]*Box1Abb[4074]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p2)*((Box1Abb[296]*Box1Abb[3850]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[3889]*Box1Abb[5829]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[4086]*Box1Abb[92]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[295]*Box1Abb[4172]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (epsP*m_p1)*((Box1Abb[295]*Box1Abb[4289]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[4463]*Box1Abb[5823]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*pow(Box1Abb[7],2.)) + (Box1Abb[4579]*Box1Abb[5823]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57])) + (epsP*m_p2)*((Box1Abb[4630]*Box1Abb[5823]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[4747]*Box1Abb[5823]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57])))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[4783]*Box1Abb[5824]*Box1Abb[7]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],3.)) + (Box1Abb[4937]*Box1Abb[7]*Box1Abb[93]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[5422]*Box1Abb[93]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5625]*Box1Abb[5826]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5758]*Box1Abb[5826]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((Box1Abb[5071]*Box1Abb[5826]*Box1Abb[7]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5221]*Box1Abb[5826]*Box1Abb[7]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[3176]*Box1Abb[58]*Box1Abb[5824]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3191]*Box1Abb[58]*Box1Abb[92]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[3231]*Box1Abb[58]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3201]*Box1Abb[58]*Box1Abb[5825]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3209]*Box1Abb[58]*Box1Abb[5825]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[3201]*Box1Abb[58]*Box1Abb[5825]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[3209]*Box1Abb[58]*Box1Abb[5825]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[3176]*Box1Abb[5824]*Box1Abb[7]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3191]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[3231]*Box1Abb[7]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3201]*Box1Abb[5825]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3209]*Box1Abb[5825]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[3201]*Box1Abb[5825]*pow(Box1Abb[7],3.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[3209]*Box1Abb[5825]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[3241]*Box1Abb[7]*Box1Abb[94]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3263]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[3317]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3273]*Box1Abb[5826]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3327]*Box1Abb[5826]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[3273]*Box1Abb[5826]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[3283]*Box1Abb[5825]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[3241]*pow(Box1Abb[7],2.)*Box1Abb[95]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3263]*Box1Abb[585]*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[3317]*Box1Abb[585]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3273]*Box1Abb[5827]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[3327]*Box1Abb[5827]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[3273]*Box1Abb[5827]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[3283]*Box1Abb[5828]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)));
    }
    else if (LR == 1) {
      return One
	*((Box1Abb[296]*Box1Abb[5890]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (epsP*m_p1)*((Box1Abb[5898]*Box1Abb[6239]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[5822]*Box1Abb[5903]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)) + (Box1Abb[5822]*Box1Abb[5895]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3])) + (epsP*m_p2)*((Box1Abb[5822]*Box1Abb[5895]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[5822]*Box1Abb[5897]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3])))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((Box1Abb[422]*Box1Abb[6237]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[3039]*Box1Abb[6238]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[5833]*Box1Abb[93]*(epsP*epsV))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[295]*Box1Abb[5848]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[390]) + (epsP*m_p2)*((Box1Abb[5823]*Box1Abb[5857]*(epsV*m_p1))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[5823]*Box1Abb[5864]*(epsV*m_p2))/(pow(Box1Abb[3],2.)*Box1Abb[390])) + (epsP*m_p1)*((Box1Abb[295]*Box1Abb[5880]*(epsV*m_pP))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[5823]*Box1Abb[5887]*(epsV*m_p1))/(pow(Box1Abb[3],2.)*Box1Abb[390]) + (Box1Abb[5823]*Box1Abb[5889]*(epsV*m_p2))/(pow(Box1Abb[3],2.)*Box1Abb[390])))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[1324]*Box1Abb[5890]*Box1Abb[93]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[296]*Box1Abb[5972]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (epsP*m_p2)*((Box1Abb[0]*Box1Abb[296]*Box1Abb[5985]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[5986]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58])) + (epsP*m_p1)*((Box1Abb[296]*Box1Abb[5997]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[5998]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[5999]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[58])))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[6000]*Box1Abb[93]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[296]*Box1Abb[6001]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p2)*((Box1Abb[296]*Box1Abb[6002]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[5829]*Box1Abb[6003]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[5830]*Box1Abb[6004]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[5830]*Box1Abb[6005]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[296]*Box1Abb[6006]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[6016]*Box1Abb[93]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[296]*Box1Abb[6035]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (epsP*m_p1)*((Box1Abb[295]*Box1Abb[6036]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[5823]*Box1Abb[6037]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*pow(Box1Abb[7],2.)) + (Box1Abb[5822]*Box1Abb[6044]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57])) + (epsP*m_p2)*((Box1Abb[5822]*Box1Abb[6045]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[5822]*Box1Abb[6064]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57])))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[5824]*Box1Abb[6065]*Box1Abb[7]*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],3.)) + (Box1Abb[6066]*Box1Abb[7]*Box1Abb[93]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[6069]*Box1Abb[93]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5825]*Box1Abb[6116]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5825]*Box1Abb[6176]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((Box1Abb[5825]*Box1Abb[6067]*Box1Abb[7]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[5826]*Box1Abb[6068]*Box1Abb[7]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[58]*Box1Abb[5824]*Box1Abb[5904]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[58]*Box1Abb[5913]*Box1Abb[93]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[58]*Box1Abb[5919]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[58]*Box1Abb[5826]*Box1Abb[5917]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[58]*Box1Abb[5826]*Box1Abb[5918]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[58]*Box1Abb[5826]*Box1Abb[5917]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[58]*Box1Abb[5826]*Box1Abb[5918]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[5824]*Box1Abb[5904]*Box1Abb[7]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[5913]*Box1Abb[7]*Box1Abb[93]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[5919]*Box1Abb[7]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5826]*Box1Abb[5917]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5826]*Box1Abb[5918]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[5826]*Box1Abb[5917]*pow(Box1Abb[7],3.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[5826]*Box1Abb[5918]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[5824]*Box1Abb[5924]*Box1Abb[7]*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[5925]*Box1Abb[7]*Box1Abb[93]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[5935]*Box1Abb[93]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5826]*Box1Abb[5926]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5826]*Box1Abb[5927]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[5826]*Box1Abb[5926]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[5826]*Box1Abb[5927]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[5924]*Box1Abb[6240]*pow(Box1Abb[7],2.)*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[297]*Box1Abb[5925]*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[297]*Box1Abb[5935]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5827]*Box1Abb[5926]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[5827]*Box1Abb[5927]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[5827]*Box1Abb[5926]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[5827]*Box1Abb[5927]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)));
    }
  }
  // ubar1 \slashed{k} P_i v2
  if (ME == 2) {
    if (LR == 0) {
      return One
	*((Box1Abb[1270]*Box1Abb[2904]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (epsP*m_p2)*((Box1Abb[0]*Box1Abb[1277]*Box1Abb[2905]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1279]*Box1Abb[2904]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58])) + (epsP*m_p1)*((Box1Abb[1276]*Box1Abb[2904]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]*Box1Abb[7]) + (Box1Abb[0]*Box1Abb[2905]*Box1Abb[420]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1277]*Box1Abb[2905]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58])))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((Box1Abb[2902]*Box1Abb[418]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[2903]*Box1Abb[422]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[2906]*Box1Abb[384]*(epsP*epsV))/(Box1Abb[0]*Box1Abb[3]) + (Box1Abb[1287]*Box1Abb[2907]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p2)*((Box1Abb[1296]*Box1Abb[2908]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[1300]*Box1Abb[2908]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[1307]*Box1Abb[2907]*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1313]*Box1Abb[2908]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[3],2.)) + (Box1Abb[1322]*Box1Abb[2908]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[57])))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[1221]*Box1Abb[1270]*(epsP*epsV))/(Box1Abb[0]*Box1Abb[116]*Box1Abb[3]) + (Box1Abb[1536]*Box1Abb[2904]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p1)*((Box1Abb[1597]*Box1Abb[2904]*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1639]*Box1Abb[2904]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1690]*Box1Abb[2904]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p2)*((Box1Abb[1736]*Box1Abb[2904]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1787]*Box1Abb[2904]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1215]*Box1Abb[1789]*Box1Abb[206]*(epsP*epsV))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1813]*Box1Abb[2926]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (epsP*m_p2)*((Box1Abb[1855]*Box1Abb[2926]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[1884]*Box1Abb[2926]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[1925]*Box1Abb[2927]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*pow(Box1Abb[7],2.)) + (Box1Abb[1958]*Box1Abb[2927]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[2037]*Box1Abb[2926]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*Box1Abb[7])))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1216]*Box1Abb[192]*Box1Abb[2040]*(epsP*epsV))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)) + (Box1Abb[2077]*Box1Abb[2904]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (epsP*m_p2)*((Box1Abb[2108]*Box1Abb[2904]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[2132]*Box1Abb[2928]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[2185]*Box1Abb[2903]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*pow(Box1Abb[7],2.)) + (Box1Abb[2240]*Box1Abb[2903]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[2288]*Box1Abb[2908]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7])))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1223]*Box1Abb[2292]*Box1Abb[416]*pow(Box1Abb[7],2.)*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[2329]*Box1Abb[2929]*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[2597]*Box1Abb[2930]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[2698]*Box1Abb[2930]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[1216]*Box1Abb[2831]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((Box1Abb[1216]*Box1Abb[2417]*Box1Abb[7]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[2503]*Box1Abb[2929]*Box1Abb[7]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[2909]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],2.)*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[2910]*Box1Abb[384]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[0]*Box1Abb[1350]*Box1Abb[20]*Box1Abb[2912]*Box1Abb[57]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[1353]*Box1Abb[2913]*Box1Abb[57]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1216]*Box1Abb[1373]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1216]*Box1Abb[1341]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1347]*Box1Abb[2911]*Box1Abb[57]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[1219]*Box1Abb[1381]*(epsP*epsV))/(Box1Abb[1374]*pow(Box1Abb[3],2.)) + (Box1Abb[1390]*Box1Abb[2920]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[1428]*Box1Abb[2923]*(epsV*m_pP))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[117]*Box1Abb[1433]*Box1Abb[2923]*(epsV*m_p1))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[1450]*Box1Abb[2921]*(epsV*m_p2))/(Box1Abb[1374]*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((pow(Box1Abb[0],2.)*Box1Abb[1406]*Box1Abb[2921]*(epsV*m_p1))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[1416]*Box1Abb[2922]*(epsV*m_p2))/(Box1Abb[1374]*pow(Box1Abb[3],3.))))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1223]*Box1Abb[1455]*Box1Abb[384]*(epsP*epsV))/pow(Box1Abb[3],2.) + (Box1Abb[1462]*Box1Abb[2924]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[1480]*Box1Abb[2925]*(epsV*m_pP))/(pow(Box1Abb[3],3.)*Box1Abb[7]) + (Box1Abb[117]*Box1Abb[1433]*Box1Abb[2925]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1221]*Box1Abb[1450]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1221]*Box1Abb[1406]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1416]*Box1Abb[2924]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[2914]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],3.)*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[2915]*Box1Abb[384]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],3.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[0]*Box1Abb[1350]*Box1Abb[20]*Box1Abb[2918]*Box1Abb[57]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[1353]*Box1Abb[2919]*Box1Abb[57]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1373]*Box1Abb[2916]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1341]*Box1Abb[2916]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1347]*Box1Abb[2917]*Box1Abb[57]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)));
    }
    else if (LR == 1) {
      return One
	*((Box1Abb[1270]*Box1Abb[3003]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (epsP*m_p2)*((Box1Abb[0]*Box1Abb[1277]*Box1Abb[3004]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1279]*Box1Abb[3003]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58])) + (epsP*m_p1)*((Box1Abb[1276]*Box1Abb[3003]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]*Box1Abb[7]) + (Box1Abb[0]*Box1Abb[3004]*Box1Abb[420]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1277]*Box1Abb[3004]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58])))

	+B_0(0.,0.,m_m2,m_mu2)*(epsP*m_p1)*((Box1Abb[3001]*Box1Abb[418]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[3002]*Box1Abb[422]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[3005]*Box1Abb[384]*(epsP*epsV))/(Box1Abb[0]*Box1Abb[3]) + (Box1Abb[1287]*Box1Abb[3006]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p2)*((Box1Abb[1296]*Box1Abb[3007]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[1300]*Box1Abb[3007]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[1307]*Box1Abb[3006]*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1313]*Box1Abb[3007]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[3],2.)) + (Box1Abb[1322]*Box1Abb[3007]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[3],2.)*Box1Abb[57])))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[1263]*Box1Abb[1270]*(epsP*epsV))/(Box1Abb[0]*Box1Abb[116]*Box1Abb[3]) + (Box1Abb[1536]*Box1Abb[3003]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (epsP*m_p1)*((Box1Abb[1597]*Box1Abb[3003]*(epsV*m_pP))/(pow(Box1Abb[0],2.)*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1639]*Box1Abb[3003]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1690]*Box1Abb[3003]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p2)*((Box1Abb[1736]*Box1Abb[3003]*(epsV*m_p1))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[1787]*Box1Abb[3003]*(epsV*m_p2))/(Box1Abb[0]*pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1257]*Box1Abb[1789]*Box1Abb[206]*(epsP*epsV))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1813]*Box1Abb[3025]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (epsP*m_p2)*((Box1Abb[1855]*Box1Abb[3025]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]) + (Box1Abb[1884]*Box1Abb[3025]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[1925]*Box1Abb[3026]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*pow(Box1Abb[7],2.)) + (Box1Abb[1958]*Box1Abb[3026]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[2037]*Box1Abb[3025]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[57]*Box1Abb[7])))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1258]*Box1Abb[192]*Box1Abb[2040]*(epsP*epsV))/(Box1Abb[116]*Box1Abb[3]*pow(Box1Abb[7],2.)) + (Box1Abb[2077]*Box1Abb[3003]*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (epsP*m_p2)*((Box1Abb[2108]*Box1Abb[3003]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)) + (Box1Abb[2132]*Box1Abb[3027]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.))) + (epsP*m_p1)*((Box1Abb[2185]*Box1Abb[3002]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*pow(Box1Abb[7],2.)) + (Box1Abb[2240]*Box1Abb[3002]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7]) + (Box1Abb[2288]*Box1Abb[3007]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],2.)*Box1Abb[7])))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1265]*Box1Abb[2292]*Box1Abb[416]*pow(Box1Abb[7],2.)*(epsP*epsV))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[2329]*Box1Abb[3028]*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[2597]*Box1Abb[3029]*(epsV*m_pP))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[2698]*Box1Abb[3029]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[1258]*Box1Abb[2831]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((Box1Abb[1258]*Box1Abb[2417]*Box1Abb[7]*(epsV*m_p1))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.)) + (Box1Abb[2503]*Box1Abb[3028]*Box1Abb[7]*(epsV*m_p2))/(pow(Box1Abb[116],2.)*pow(Box1Abb[3],3.))))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1324]*Box1Abb[3008]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],2.)*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3009]*Box1Abb[384]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],2.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[0]*Box1Abb[1350]*Box1Abb[20]*Box1Abb[3011]*Box1Abb[57]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[1353]*Box1Abb[3012]*Box1Abb[57]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1258]*Box1Abb[1373]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1258]*Box1Abb[1341]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1347]*Box1Abb[3010]*Box1Abb[57]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[1261]*Box1Abb[1381]*(epsP*epsV))/(Box1Abb[1374]*pow(Box1Abb[3],2.)) + (Box1Abb[1390]*Box1Abb[3019]*(epsP*m_p2)*(epsV*m_pP))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (epsP*m_p1)*((Box1Abb[1428]*Box1Abb[3022]*(epsV*m_pP))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[117]*Box1Abb[1433]*Box1Abb[3022]*(epsV*m_p1))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[1450]*Box1Abb[3020]*(epsV*m_p2))/(Box1Abb[1374]*pow(Box1Abb[3],3.))) + (epsP*m_p2)*((pow(Box1Abb[0],2.)*Box1Abb[1406]*Box1Abb[3020]*(epsV*m_p1))/(Box1Abb[1374]*pow(Box1Abb[3],3.)) + (pow(Box1Abb[0],2.)*Box1Abb[1416]*Box1Abb[3021]*(epsV*m_p2))/(Box1Abb[1374]*pow(Box1Abb[3],3.))))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1265]*Box1Abb[1455]*Box1Abb[384]*(epsP*epsV))/pow(Box1Abb[3],2.) + (Box1Abb[1462]*Box1Abb[3023]*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[1480]*Box1Abb[3024]*(epsV*m_pP))/(pow(Box1Abb[3],3.)*Box1Abb[7]) + (Box1Abb[117]*Box1Abb[1433]*Box1Abb[3024]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1263]*Box1Abb[1450]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1263]*Box1Abb[1406]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1416]*Box1Abb[3023]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1324]*Box1Abb[3013]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],3.)*(epsP*epsV))/pow(Box1Abb[3],3.) + (Box1Abb[3014]*Box1Abb[384]*pow(Box1Abb[390],2.)*pow(Box1Abb[7],3.)*(epsP*m_p2)*(epsV*m_pP))/pow(Box1Abb[3],3.) + (epsP*m_p1)*((Box1Abb[0]*Box1Abb[1350]*Box1Abb[20]*Box1Abb[3017]*Box1Abb[57]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],3.) + (Box1Abb[117]*Box1Abb[1353]*Box1Abb[3018]*Box1Abb[57]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1373]*Box1Abb[3015]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],3.)) + (epsP*m_p2)*((Box1Abb[1341]*Box1Abb[3015]*pow(Box1Abb[7],2.)*(epsV*m_p1))/pow(Box1Abb[3],3.) + (Box1Abb[1347]*Box1Abb[3016]*Box1Abb[57]*pow(Box1Abb[7],2.)*(epsV*m_p2))/pow(Box1Abb[3],3.)));
    }
  }
  // ubar1 \slashed{epsP} P_i v2
  if (ME == 3) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[1214]*Box1Abb[390]*Box1Abb[978]*(epsV*m_pP))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1215]*Box1Abb[980]*(epsV*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1216]*Box1Abb[983]*(epsV*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[1026]*Box1Abb[1221]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1030]*Box1Abb[1221]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1037]*Box1Abb[1221]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1048]*Box1Abb[1049]*Box1Abb[1216]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1053]*Box1Abb[1216]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1059]*Box1Abb[1222]*Box1Abb[192]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1066]*Box1Abb[1215]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1075]*Box1Abb[1215]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[1085]*Box1Abb[1216]*Box1Abb[192]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1217]*Box1Abb[57]*Box1Abb[7]*Box1Abb[985]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[117]*Box1Abb[1218]*Box1Abb[7]*Box1Abb[990]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1218]*Box1Abb[7]*Box1Abb[996]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[1002]*Box1Abb[1218]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[1005]*Box1Abb[117]*Box1Abb[1218]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[1014]*Box1Abb[1218]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1047]*Box1Abb[1218]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[1005]*Box1Abb[117]*Box1Abb[1218]*Box1Abb[192]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1014]*Box1Abb[1218]*Box1Abb[192]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1108]*Box1Abb[1223]*Box1Abb[192]*(epsV*m_pP))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[1138]*Box1Abb[1223]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[1181]*Box1Abb[1223]*(epsV*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1219]*Box1Abb[57]*pow(Box1Abb[7],2.)*Box1Abb[985]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[117]*Box1Abb[1220]*pow(Box1Abb[7],2.)*Box1Abb[990]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1220]*pow(Box1Abb[7],2.)*Box1Abb[996]*(epsV*m_p2))/pow(Box1Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[1256]*Box1Abb[390]*Box1Abb[978]*(epsV*m_pP))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1257]*Box1Abb[980]*(epsV*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1258]*Box1Abb[983]*(epsV*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[1026]*Box1Abb[1263]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1030]*Box1Abb[1263]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[1037]*Box1Abb[1263]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[1048]*Box1Abb[1049]*Box1Abb[1258]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1053]*Box1Abb[1258]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1059]*Box1Abb[1264]*Box1Abb[192]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[1066]*Box1Abb[1257]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[1075]*Box1Abb[1257]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[1085]*Box1Abb[1258]*Box1Abb[192]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[1259]*Box1Abb[57]*Box1Abb[7]*Box1Abb[985]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[117]*Box1Abb[1260]*Box1Abb[7]*Box1Abb[990]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1260]*Box1Abb[7]*Box1Abb[996]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[1002]*Box1Abb[1260]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[1005]*Box1Abb[117]*Box1Abb[1260]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[1014]*Box1Abb[1260]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1047]*Box1Abb[1260]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[1005]*Box1Abb[117]*Box1Abb[1260]*Box1Abb[192]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1014]*Box1Abb[1260]*Box1Abb[192]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1108]*Box1Abb[1265]*Box1Abb[192]*(epsV*m_pP))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[1138]*Box1Abb[1265]*Box1Abb[192]*(epsV*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[1181]*Box1Abb[1265]*(epsV*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[1261]*Box1Abb[57]*pow(Box1Abb[7],2.)*Box1Abb[985]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[117]*Box1Abb[1262]*pow(Box1Abb[7],2.)*Box1Abb[990]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[1262]*pow(Box1Abb[7],2.)*Box1Abb[996]*(epsV*m_p2))/pow(Box1Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsV} P_i v2
  if (ME == 4) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[680]*Box1Abb[913]*(epsP*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[686]*Box1Abb[913]*(epsP*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[695]*Box1Abb[914]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[702]*Box1Abb[914]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[736]*Box1Abb[914]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[742]*Box1Abb[915]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[765]*Box1Abb[913]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[780]*Box1Abb[916]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[649]*Box1Abb[657]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[192]*Box1Abb[657]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[117]*Box1Abb[657]*Box1Abb[7]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[657]*pow(Box1Abb[7],2.)*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[721]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[674]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[846]*(epsP*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[889]*(epsP*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[117]*Box1Abb[674]*Box1Abb[912]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[674]*Box1Abb[7]*Box1Abb[912]*(epsP*m_p2))/pow(Box1Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[913]*Box1Abb[919]*(epsP*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[913]*Box1Abb[920]*(epsP*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[914]*Box1Abb[921]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[914]*Box1Abb[922]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[914]*Box1Abb[924]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[915]*Box1Abb[925]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[913]*Box1Abb[926]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[916]*Box1Abb[927]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[0]*Box1Abb[649]*Box1Abb[917]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[0]*Box1Abb[192]*Box1Abb[917]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[117]*Box1Abb[7]*Box1Abb[917]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (pow(Box1Abb[7],2.)*Box1Abb[917]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[923]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[918]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[952]*(epsP*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[953]*(epsP*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[117]*Box1Abb[912]*Box1Abb[918]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[7]*Box1Abb[912]*Box1Abb[918]*(epsP*m_p2))/pow(Box1Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsP} \slashed{k} P_i v2
  if (ME == 5) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[353]*Box1Abb[582]*(epsV*m_pP))/(Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[355]*(epsV*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[357]*(epsV*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[296]*Box1Abb[365]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[296]*Box1Abb[369]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[373]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[295]*Box1Abb[419]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[295]*Box1Abb[423]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[298]*Box1Abb[426]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[430]*Box1Abb[586]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[296]*Box1Abb[444]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[296]*Box1Abb[454]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[382]*Box1Abb[58]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[192]*Box1Abb[353]*Box1Abb[383]*Box1Abb[58]*Box1Abb[92]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[389]*Box1Abb[58]*Box1Abb[92]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[382]*Box1Abb[7]*Box1Abb[92]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[192]*Box1Abb[353]*Box1Abb[383]*Box1Abb[7]*Box1Abb[92]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[389]*Box1Abb[7]*Box1Abb[92]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[390]*Box1Abb[391]*Box1Abb[583]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[404]*Box1Abb[92]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[415]*Box1Abb[92]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[485]*Box1Abb[93]*(epsV*m_pP))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[520]*Box1Abb[93]*(epsV*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[552]*Box1Abb[93]*(epsV*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[390]*Box1Abb[391]*Box1Abb[584]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[404]*Box1Abb[585]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[415]*Box1Abb[585]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[582]*Box1Abb[587]*(epsV*m_pP))/(Box1Abb[3]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[588]*(epsV*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[589]*(epsV*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[296]*Box1Abb[599]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[296]*Box1Abb[600]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]) + (Box1Abb[0]*Box1Abb[296]*Box1Abb[601]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[295]*Box1Abb[605]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[295]*Box1Abb[606]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[298]*Box1Abb[607]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[586]*Box1Abb[608]*(epsV*m_pP))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[296]*Box1Abb[611]*(epsV*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]) + (Box1Abb[296]*Box1Abb[614]*(epsV*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[58]*Box1Abb[595]*Box1Abb[93]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[383]*Box1Abb[58]*Box1Abb[587]*Box1Abb[7]*Box1Abb[93]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[58]*Box1Abb[598]*Box1Abb[93]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[595]*Box1Abb[7]*Box1Abb[93]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[383]*Box1Abb[587]*pow(Box1Abb[7],2.)*Box1Abb[93]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[598]*Box1Abb[7]*Box1Abb[93]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[390]*Box1Abb[602]*Box1Abb[647]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[603]*Box1Abb[92]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[604]*Box1Abb[92]*(epsV*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[615]*Box1Abb[93]*(epsV*m_pP))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[616]*Box1Abb[93]*(epsV*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[617]*Box1Abb[93]*(epsV*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[390]*Box1Abb[602]*Box1Abb[648]*Box1Abb[7]*(epsV*m_pP))/pow(Box1Abb[3],2.) + (Box1Abb[585]*Box1Abb[603]*Box1Abb[7]*(epsV*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[585]*Box1Abb[604]*Box1Abb[7]*(epsV*m_p2))/pow(Box1Abb[3],2.));
    }
  }
  // ubar1 \slashed{k} \slashed{epsV} P_i v2
  if (ME == 6) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[111]*Box1Abb[295]*(epsP*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[115]*Box1Abb[295]*(epsP*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[122]*Box1Abb[296]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[128]*Box1Abb[296]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))
	
	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[191]*Box1Abb[295]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[194]*Box1Abb[298]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[205]*Box1Abb[295]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[213]*Box1Abb[295]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[137]*Box1Abb[58]*Box1Abb[92]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[148]*Box1Abb[58]*Box1Abb[92]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[137]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[148]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[164]*Box1Abb[93]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[182]*Box1Abb[93]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[248]*Box1Abb[92]*(epsP*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[274]*Box1Abb[92]*(epsP*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[164]*Box1Abb[297]*Box1Abb[7]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[182]*Box1Abb[297]*Box1Abb[7]*(epsP*m_p2))/pow(Box1Abb[3],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)*((Box1Abb[296]*Box1Abb[300]*(epsP*m_p1))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]) + (Box1Abb[296]*Box1Abb[305]*(epsP*m_p2))/(Box1Abb[3]*Box1Abb[57]*Box1Abb[58]))

	+B_0(m_s,m_m2,m_m2,m_mu2)*((Box1Abb[296]*Box1Abb[306]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]) + (Box1Abb[296]*Box1Abb[309]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[58]))

	+B_0(m_s1k,0.,m_m2,m_mu2)*((Box1Abb[295]*Box1Abb[314]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]) + (Box1Abb[192]*Box1Abb[298]*Box1Abb[315]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[7]))

	+B_0(m_m2,0.,m_m2,m_mu2)*((Box1Abb[296]*Box1Abb[325]*(epsP*m_p1))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]*Box1Abb[7]) + (Box1Abb[296]*Box1Abb[330]*(epsP*m_p2))/(Box1Abb[116]*Box1Abb[3]*Box1Abb[57]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)*((Box1Abb[310]*Box1Abb[58]*Box1Abb[92]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[311]*Box1Abb[58]*Box1Abb[92]*(epsP*m_p2))/pow(Box1Abb[3],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[310]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[311]*Box1Abb[7]*Box1Abb[92]*(epsP*m_p2))/pow(Box1Abb[3],2.))
	
	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)*((Box1Abb[312]*Box1Abb[93]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[313]*Box1Abb[93]*(epsP*m_p2))/pow(Box1Abb[3],2.))
	
	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)*((Box1Abb[331]*Box1Abb[93]*(epsP*m_p1))/(Box1Abb[116]*pow(Box1Abb[3],2.)) + (Box1Abb[332]*Box1Abb[93]*(epsP*m_p2))/(Box1Abb[116]*pow(Box1Abb[3],2.)))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)*((Box1Abb[297]*Box1Abb[312]*Box1Abb[7]*(epsP*m_p1))/pow(Box1Abb[3],2.) + (Box1Abb[297]*Box1Abb[313]*Box1Abb[7]*(epsP*m_p2))/pow(Box1Abb[3],2.));
    }
  }
  // ubar1 \slashed{epsP} \slashed{epsV} P_i v2
  if (ME == 7) {
    if (LR == 0) {
      return (Box1Abb[56]*Box1Abb[92]*B_0(m_m2,0.,m_m2,m_mu2))/Box1Abb[57]

	+(Box1Abb[56]*Box1Abb[93]*B_0(m_s12,m_m2,m_m2,m_mu2))/Box1Abb[57]

	+(Box1Abb[58]*Box1Abb[65]*Box1Abb[94]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[65]*Box1Abb[7]*Box1Abb[94]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[76]*Box1Abb[94]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[83]*Box1Abb[94]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[7]*Box1Abb[83]*Box1Abb[95]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[3];
    }
    else if (LR == 1) {
      return (Box1Abb[56]*Box1Abb[92]*B_0(m_m2,0.,m_m2,m_mu2))/Box1Abb[57]

	+(Box1Abb[56]*Box1Abb[93]*B_0(m_s12,m_m2,m_m2,m_mu2))/Box1Abb[57]

	+(Box1Abb[58]*Box1Abb[94]*Box1Abb[96]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[7]*Box1Abb[94]*Box1Abb[96]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[94]*Box1Abb[97]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[94]*Box1Abb[98]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[7]*Box1Abb[95]*Box1Abb[98]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[3];
    }
  }
  // ubar1 \slashed{epsP} \slashed{k} \slashed{epsV} P_i v2
  if (ME == 8) {
    if (LR == 0) {
      return (Box1Abb[0]*Box1Abb[6]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[11]*Box1Abb[7]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[19]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[23]*Box1Abb[40]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[28]*Box1Abb[41]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[3];
    }
    else if (LR == 1) {
      return (Box1Abb[0]*Box1Abb[43]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[44]*Box1Abb[7]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[45]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[3]

	+(Box1Abb[40]*Box1Abb[46]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[3]

	+(Box1Abb[41]*Box1Abb[47]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[3];
    }
  }
  else {
    msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
		<< "Values range from 1 to 8.";
  }
  return Zero;
}
