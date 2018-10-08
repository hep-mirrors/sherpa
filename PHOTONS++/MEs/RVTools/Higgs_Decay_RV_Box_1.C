#include "PHOTONS++/MEs/RVTools/Higgs_Decay_RV_Diagrams.H"
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

Higgs_Decay_RV_Box_1::Higgs_Decay_RV_Box_1
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Higgs_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Box1Abb = new Complex[364];
  Init_Coefficients();
}

Higgs_Decay_RV_Box_1::~Higgs_Decay_RV_Box_1()
{
  delete [] Box1Abb;
}

void Higgs_Decay_RV_Box_1::Init_Coefficients() 
{
  Init_Box_1_Coefficients();
  return;
}

// Set up coefficients multiplying scalar master integrals for box with emission off leg 1
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Box_1::Init_Box_1_Coefficients() 
{
  Box1Abb[0]=-1. + m_z12;

  Box1Abb[1]=-2.*m_x + m_z12;

  Box1Abb[2]=-1. + m_x;

  Box1Abb[3]=-1. - 2.*m_x + m_z12;

  Box1Abb[4]=m_x + Box1Abb[2]*m_x*m_z12 + Box1Abb[3]*m_z12*m_z1k + m_z12*m_z1k_2;

  Box1Abb[5]=m_x - m_z1k;

  Box1Abb[6]=2.*m_x - m_z12;

  Box1Abb[7]=1. + m_z12;

  Box1Abb[8]=1. + m_x - m_z12 + m_x*m_z12 - Box1Abb[7]*m_z1k;

  Box1Abb[9]=-1. - 2.*m_x + m_z12 + 2.*m_z1k;

  Box1Abb[10]=-2. + m_z12;

  Box1Abb[11]=Box1Abb[10]*m_x + m_z12*m_z1k;

  Box1Abb[20]=m_s*m_z12;

  Box1Abb[21]=m_s*m_x;

  Box1Abb[22]=m_s*m_z1k;

  Box1Abb[23]=-m_z12;

  Box1Abb[24]=-m_s;

  Box1Abb[33]=1. - m_z12;

  Box1Abb[34]=2.*m_z12 + m_z1k;

  Box1Abb[35]=1. + 2.*m_z1k;

  Box1Abb[36]=-Box1Abb[35]*m_x + m_x_2 + Box1Abb[34]*m_z1k;

  Box1Abb[37]=pow(Box1Abb[0],2.) + 2.*Box1Abb[36] - 2.*m_z1k;

  Box1Abb[38]=-1. + 3.*m_z1k;

  Box1Abb[39]=4. + 4.*m_z12 + 6.*m_z1k;

  Box1Abb[40]=-1. + m_z1k;

  Box1Abb[41]=1. + 2.*Box1Abb[40]*m_z1k;

  Box1Abb[42]=2. + 4.*Box1Abb[40]*m_z1k;

  Box1Abb[43]=-4. + m_z12;

  Box1Abb[44]=1. + Box1Abb[43]*m_z12 - 8.*m_z12*m_z1k - 6.*m_z1k_2;

  Box1Abb[45]=Box1Abb[40]*Box1Abb[41] + Box1Abb[44]*m_x + Box1Abb[39]*m_x_2 - 2.*m_x_3 + Box1Abb[42]*m_z12 + Box1Abb[38]*m_z12_2;

  Box1Abb[46]=pow(Box1Abb[0],2.) - 2.*m_z1k;

  Box1Abb[47]=-3. + m_z12 - 2.*m_z1k;

  Box1Abb[48]=-2. - 2.*Box1Abb[47]*m_z12 + 4.*m_z1k;

  Box1Abb[49]=Box1Abb[48]*m_x - 4.*Box1Abb[7]*m_x_2 + Box1Abb[46]*m_z12;

  Box1Abb[50]=1. - 3.*m_z12 + m_z12_2 - 2.*Box1Abb[7]*m_z1k;

  Box1Abb[51]=2.*Box1Abb[50]*m_x + 4.*Box1Abb[7]*m_x_2 - pow(Box1Abb[0],2.)*m_z12 + 2.*m_z12*m_z1k;

  Box1Abb[59]=m_m;

  Box1Abb[60]=-m_m*m_s;

  Box1Abb[68]=-1. - m_x + m_z12 + m_z1k;

  Box1Abb[69]=-m_x + m_z1k;

  Box1Abb[70]=1. + m_z1k;

  Box1Abb[71]=pow(Box1Abb[40],2.) - 2.*Box1Abb[70]*m_x + m_x_2;

  Box1Abb[72]=1. + 2.*m_z12;

  Box1Abb[73]=-1. + m_z12 + m_z1k;

  Box1Abb[74]=2. - 3.*m_z12 - 2.*m_z12*m_z1k;

  Box1Abb[75]=-Box1Abb[40]*Box1Abb[73] + Box1Abb[74]*m_x + Box1Abb[72]*m_x_2;

  Box1Abb[76]=-1. + m_x + m_z1k + 2.*m_z12*m_z1k;

  Box1Abb[77]=3. + 2.*m_x - 2.*m_z12;

  Box1Abb[78]=2. + m_z12;

  Box1Abb[79]=-1. + m_x + pow(Box1Abb[2],2.)*m_z12 + Box1Abb[77]*m_z1k - Box1Abb[78]*m_z1k_2;

  Box1Abb[80]=Box1Abb[78]*m_x + Box1Abb[40]*m_z12;

  Box1Abb[81]=1. + m_z12 + 3.*m_z12*m_z1k;

  Box1Abb[82]=1. + m_z12 + Box1Abb[10]*m_z1k + 3.*m_z1k_2;

  Box1Abb[83]=-1. + Box1Abb[82]*m_z12 + m_z1k;

  Box1Abb[84]=-Box1Abb[83]*m_x + Box1Abb[81]*m_x_2 + pow(Box1Abb[40],2.)*Box1Abb[73]*m_z12 - m_x_3*m_z12;

  Box1Abb[85]=-2. + m_z12 + m_z1k;

  Box1Abb[86]=pow(Box1Abb[0],2.) + 3.*m_z12*m_z1k;

  Box1Abb[87]=1. + Box1Abb[86]*m_z1k;

  Box1Abb[88]=-Box1Abb[87]*m_x + Box1Abb[81]*m_x_2 - m_x_3*m_z12 + Box1Abb[40]*Box1Abb[85]*m_z12*m_z1k;

  Box1Abb[89]=-13. + 8.*m_z12 + 5.*m_z1k;

  Box1Abb[90]=5. + Box1Abb[89]*m_z12;

  Box1Abb[91]=-4. + 3.*m_z12;

  Box1Abb[92]=2.*pow(Box1Abb[0],2.) + Box1Abb[91]*m_z1k + m_z1k_2;

  Box1Abb[93]=Box1Abb[92]*m_z12 + m_z1k;

  Box1Abb[94]=-Box1Abb[93]*m_x + Box1Abb[90]*m_x_2 - 3.*m_x_3*m_z12 - Box1Abb[73]*m_z12*m_z1k_2;

  Box1Abb[95]=-2. + 4.*m_z1k;

  Box1Abb[96]=pow(Box1Abb[40],2.) + Box1Abb[95]*m_z12 + m_z12_2;

  Box1Abb[97]=-4. + 3.*m_z12 + 11.*m_z1k;

  Box1Abb[98]=3. + Box1Abb[97]*m_z12;

  Box1Abb[99]=2. + 7.*m_z1k;

  Box1Abb[100]=-2. + m_z12 + Box1Abb[99]*m_z1k + 3.*m_z12*m_z1k;

  Box1Abb[101]=1. + Box1Abb[100]*m_z12 - m_z1k;

  Box1Abb[102]=-Box1Abb[101]*m_x + Box1Abb[98]*m_x_2 - 5.*m_x_3*m_z12 + Box1Abb[96]*m_z12*m_z1k;

  Box1Abb[103]=-14. + 15.*m_z12;

  Box1Abb[104]=-1. + m_z12_2;

  Box1Abb[105]=2.*pow(Box1Abb[0],2.)*m_z12 + 3.*Box1Abb[104]*m_z1k + Box1Abb[78]*m_z1k_2;

  Box1Abb[106]=-7. + 4.*m_z12 + 5.*m_z1k;

  Box1Abb[107]=-11. - 3.*Box1Abb[106]*m_z12 + 12.*m_z1k;

  Box1Abb[108]=2. + Box1Abb[107]*m_z12;

  Box1Abb[109]=Box1Abb[108]*m_x_2 + Box1Abb[105]*m_x*m_z12 + Box1Abb[103]*m_x_3*m_z12 - Box1Abb[73]*m_z12_2*m_z1k_2;

  Box1Abb[110]=-18. + 11.*m_z12;

  Box1Abb[111]=pow(Box1Abb[40],2.) + Box1Abb[10]*m_z12;

  Box1Abb[112]=2. + 3.*m_z12;

  Box1Abb[113]=pow(Box1Abb[0],2.) + m_z1k + 5.*Box1Abb[10]*m_z12*m_z1k - Box1Abb[112]*m_z1k_2;

  Box1Abb[114]=2. - 3.*m_z12 - 7.*m_z1k;

  Box1Abb[115]=5. + Box1Abb[114]*m_z12 + 20.*m_z1k;

  Box1Abb[116]=-2. + Box1Abb[115]*m_z12;

  Box1Abb[117]=Box1Abb[116]*m_x_2 + Box1Abb[113]*m_x*m_z12 + Box1Abb[110]*m_x_3*m_z12 - Box1Abb[111]*m_z12_2*m_z1k;

  Box1Abb[118]=2. + 4.*m_z12 - 7.*m_z1k;

  Box1Abb[119]=-5. + 2.*Box1Abb[118]*m_z12;

  Box1Abb[120]=-1. + 2.*m_z12 + m_z1k;

  Box1Abb[121]=28. - 25.*m_z1k;

  Box1Abb[122]=4. + 3.*m_z1k;

  Box1Abb[123]=-4. + 5.*Box1Abb[122]*m_z12 + 2.*m_z12_2 + Box1Abb[121]*m_z1k;

  Box1Abb[124]=15. - Box1Abb[123]*m_z12 + 16.*m_z1k;

  Box1Abb[125]=5. + m_z1k;

  Box1Abb[126]=1. + Box1Abb[125]*Box1Abb[40]*m_z1k;

  Box1Abb[127]=-15. + 23.*m_z1k;

  Box1Abb[128]=4. + Box1Abb[127]*m_z1k;

  Box1Abb[129]=-7. + m_z1k;

  Box1Abb[130]=3. + Box1Abb[129]*m_z1k;

  Box1Abb[131]=-1. + Box1Abb[130]*m_z1k;

  Box1Abb[132]=2.*Box1Abb[126]*pow(Box1Abb[40],2.)*m_z12 + Box1Abb[128]*Box1Abb[40]*m_z12_2 - 2.*Box1Abb[131]*m_z12_3 - pow(Box1Abb[40],3.)*m_z1k;

  Box1Abb[133]=1. + 6.*m_z1k;

  Box1Abb[134]=12. + 25.*m_z1k;

  Box1Abb[135]=-14. + 5.*m_z1k;

  Box1Abb[136]=5. + Box1Abb[135]*m_z1k;

  Box1Abb[137]=9. + 2.*Box1Abb[136]*m_z1k;

  Box1Abb[138]=1. - 2.*Box1Abb[137]*m_z12 + Box1Abb[134]*m_z12_2 + 6.*Box1Abb[70]*m_z12_3 - 3.*Box1Abb[133]*m_z1k;

  Box1Abb[139]=-7. + 8.*m_z1k;

  Box1Abb[140]=3. + Box1Abb[139]*m_z1k;

  Box1Abb[141]=-37. + 10.*m_z1k;

  Box1Abb[142]=-5. + Box1Abb[141]*m_z1k;

  Box1Abb[143]=4. + Box1Abb[142]*m_z1k;

  Box1Abb[144]=-8. + m_z1k;

  Box1Abb[145]=46. + 5.*Box1Abb[144]*m_z1k;

  Box1Abb[146]=-8. + Box1Abb[145]*m_z1k;

  Box1Abb[147]=5. + Box1Abb[146]*m_z1k;

  Box1Abb[148]=Box1Abb[140]*Box1Abb[40] + Box1Abb[147]*m_z12 + Box1Abb[143]*m_z12_2 - 6.*m_z12_3;

  Box1Abb[149]=-Box1Abb[132]*m_x - Box1Abb[148]*m_x_2 - Box1Abb[138]*m_x_3 - Box1Abb[124]*m_x_4 - Box1Abb[119]*m_x_5 - 3.*m_x_6*m_z12 + Box1Abb[120]*pow(Box1Abb[40],2.)*Box1Abb[73]*m_z12*m_z1k_2;

  Box1Abb[150]=-11. + 3.*m_z12 - 26.*m_z1k;

  Box1Abb[151]=-3. + Box1Abb[150]*m_z12;

  Box1Abb[152]=pow(Box1Abb[40],5) + Box1Abb[10]*pow(Box1Abb[40],3.)*m_z12;

  Box1Abb[153]=-5. + 4.*m_z12;

  Box1Abb[154]=-14. + m_z12;

  Box1Abb[155]=-3. + Box1Abb[154]*m_z12;

  Box1Abb[156]=-3. + m_z12 - 21.*m_z12_2;

  Box1Abb[157]=-6. + 3.*Box1Abb[153]*m_z12 + Box1Abb[155]*m_z12*m_z1k + 2.*Box1Abb[156]*m_z1k_2 - 60.*m_z12*m_z1k_3;

  Box1Abb[158]=-5. + 4.*m_z1k;

  Box1Abb[159]=23. + 55.*m_z1k;

  Box1Abb[160]=13. + 2.*Box1Abb[158]*m_z12 + Box1Abb[159]*m_z1k;

  Box1Abb[161]=10. + Box1Abb[160]*m_z12 + 8.*m_z1k;

  Box1Abb[162]=-3. + 7.*m_z1k;

  Box1Abb[163]=-3. + 5.*m_z1k + 24.*m_z1k_3;

  Box1Abb[164]=-4. + m_z1k - 34.*m_z1k_2 + 35.*m_z1k_3;

  Box1Abb[165]=10. + Box1Abb[164]*m_z1k;

  Box1Abb[166]=-2. + Box1Abb[165]*m_z12 + 2.*Box1Abb[163]*m_z12_2 + 8.*m_z1k + Box1Abb[162]*m_z12_3*m_z1k - 6.*m_z1k_2;

  Box1Abb[167]=2. + m_z1k - 5.*m_z1k_2 + 10.*m_z1k_3;

  Box1Abb[168]=2. + 9.*m_z1k;

  Box1Abb[169]=3. - Box1Abb[168]*m_z1k;

  Box1Abb[170]=-9. + 17.*m_z1k;

  Box1Abb[171]=-5. + Box1Abb[170]*m_z1k;

  Box1Abb[172]=1. + Box1Abb[171]*m_z1k;

  Box1Abb[173]=pow(Box1Abb[40],4.) - Box1Abb[167]*pow(Box1Abb[40],2.)*m_z12 - Box1Abb[172]*Box1Abb[40]*m_z12_2 + Box1Abb[169]*m_z12_3*m_z1k;

  Box1Abb[174]=-Box1Abb[173]*m_x - Box1Abb[166]*m_x_2 - Box1Abb[157]*m_x_3 - Box1Abb[161]*m_x_4 - Box1Abb[151]*m_x_5 - 5.*m_x_6*m_z12 - Box1Abb[152]*m_z12*m_z1k;

  Box1Abb[194]=(-8.*m_m*m_z12)/m_s_2;

  Box1Abb[195]=(8.*m_m*m_x)/m_s_2;

  Box1Abb[196]=(-8.*m_m*m_x)/m_s_2;

  Box1Abb[197]=(8.*m_m*m_z1k)/m_s_2;

  Box1Abb[198]=(-8.*m_m)/m_s_2;

  Box1Abb[199]=(2.*m_m)/m_s;

  Box1Abb[200]=(-2.*m_m)/m_s;

  Box1Abb[201]=-2.*m_m;

  Box1Abb[221]=-2. + 3.*m_z12 + 8.*m_z1k;

  Box1Abb[222]=-1. + m_z12 - 2.*m_z1k + 4.*m_z12*m_z1k + 4.*m_z1k_2;

  Box1Abb[223]=-Box1Abb[222]*m_x + Box1Abb[221]*m_x_2 - 4.*m_x_3 + Box1Abb[73]*m_z12*m_z1k;

  Box1Abb[224]=1. + m_x - m_z12 - m_z1k;

  Box1Abb[225]=3. - 2.*m_z12;

  Box1Abb[226]=-3. + m_z12;

  Box1Abb[227]=-16. + 5.*m_z12;

  Box1Abb[228]=2. + Box1Abb[226]*m_z12 + 8.*m_z1k + Box1Abb[227]*m_z12*m_z1k - 12.*m_z1k_2;

  Box1Abb[229]=3. + Box1Abb[43]*m_z12 - 3.*m_z12*m_z1k - 4.*m_z1k_2;

  Box1Abb[230]=2. + 3.*m_z12 + 12.*m_z1k;

  Box1Abb[231]=4. + Box1Abb[10]*Box1Abb[230]*m_z12;

  Box1Abb[232]=Box1Abb[231]*m_x_3 - Box1Abb[228]*m_x_2*m_z12 + 4.*Box1Abb[225]*m_x_4*m_z12 + Box1Abb[229]*m_x*m_z12_2*m_z1k + Box1Abb[73]*m_z12_3*m_z1k_2;

  Box1Abb[233]=4.*m_x - m_z12;

  Box1Abb[234]=2. - 3.*m_z12;

  Box1Abb[235]=Box1Abb[234]*m_x + Box1Abb[73]*m_z12;

  Box1Abb[236]=-Box1Abb[10]*m_x - m_z12*m_z1k;

  Box1Abb[237]=m_z12 - 2.*m_z1k;

  Box1Abb[238]=Box1Abb[40]*Box1Abb[73] + Box1Abb[237]*m_x + m_x_2;

  Box1Abb[239]=-Box1Abb[35]*m_x + m_x_2 + Box1Abb[120]*m_z1k;

  Box1Abb[240]=-1. + m_z12 - 2.*m_z1k;

  Box1Abb[241]=-4. + m_z12 + m_z12_2;

  Box1Abb[242]=m_z12 + Box1Abb[241]*m_z1k + 4.*Box1Abb[7]*m_z1k_2;

  Box1Abb[243]=-1. + 2.*m_z1k;

  Box1Abb[244]=-1. + 4.*m_z1k;

  Box1Abb[245]=Box1Abb[243]*pow(Box1Abb[40],2.) + Box1Abb[244]*Box1Abb[40]*Box1Abb[70]*m_z12 + 2.*Box1Abb[70]*m_z12_2*m_z1k;

  Box1Abb[246]=-Box1Abb[245]*m_x + Box1Abb[242]*m_x_2 + Box1Abb[240]*m_x_3 - m_x_4*m_z12 + pow(Box1Abb[40],2.)*Box1Abb[73]*m_z12*m_z1k;

  Box1Abb[247]=-1. + m_z12 - 4.*m_z1k + 7.*m_z12*m_z1k;

  Box1Abb[248]=-16. + 13.*m_z12;

  Box1Abb[249]=Box1Abb[248]*m_z12 + 24.*Box1Abb[0]*m_z1k;

  Box1Abb[250]=4. + Box1Abb[249]*m_z12;

  Box1Abb[251]=2. - 23.*m_z1k;

  Box1Abb[252]=8. - 3.*m_z1k;

  Box1Abb[253]=3. + 4.*Box1Abb[252]*m_z1k;

  Box1Abb[254]=-2. + Box1Abb[253]*m_z12 + Box1Abb[251]*m_z12_2 - 3.*m_z12_3 + 12.*Box1Abb[40]*m_z1k;

  Box1Abb[255]=Box1Abb[250]*m_x_3 + Box1Abb[254]*m_x_2*m_z12 - 12.*Box1Abb[0]*m_x_4*m_z12 + Box1Abb[247]*Box1Abb[73]*m_x*m_z12_2 - pow(Box1Abb[73],2.)*m_z12_3*m_z1k;

  Box1Abb[256]=2. + 5.*m_z12;

  Box1Abb[257]=-11. + 3.*m_z12;

  Box1Abb[258]=-1. + 2.*m_z12;

  Box1Abb[259]=4. + Box1Abb[257]*m_z12 + 2.*Box1Abb[0]*m_z12*m_z1k + 4.*Box1Abb[258]*m_z1k_2;

  Box1Abb[260]=-3. + m_z1k;

  Box1Abb[261]=8. + 3.*m_z12 + 20.*m_z1k;

  Box1Abb[262]=4.*Box1Abb[260] + Box1Abb[261]*m_z12;

  Box1Abb[263]=-7. + 6.*m_z1k;

  Box1Abb[264]=5. + 2.*Box1Abb[260]*m_z1k;

  Box1Abb[265]=-1. + m_z12 + 2.*Box1Abb[70]*m_z12_2*m_z1k + Box1Abb[264]*m_z1k_2 + Box1Abb[263]*m_z12*m_z1k_2;

  Box1Abb[266]=1. + 9.*m_z1k + 6.*m_z1k_2;

  Box1Abb[267]=-4. + m_z1k + m_z1k_2;

  Box1Abb[268]=1. + Box1Abb[267]*m_z1k;

  Box1Abb[269]=4.*pow(Box1Abb[40],3.) + 4.*Box1Abb[268]*m_z12 + Box1Abb[266]*m_z12_2 + m_z12_3*m_z1k;

  Box1Abb[270]=-Box1Abb[269]*m_x_2 - Box1Abb[259]*m_x_3 + Box1Abb[262]*m_x_4 - 2.*Box1Abb[256]*m_x_5 + Box1Abb[265]*m_x*m_z12 - pow(Box1Abb[40],2.)*Box1Abb[73]*m_z12_2*m_z1k;

  Box1Abb[271]=2. + m_z1k - 5.*m_z1k_2;

  Box1Abb[272]=-2.*pow(Box1Abb[40],2.) + Box1Abb[271]*m_z12 - m_z12_2*m_z1k;

  Box1Abb[273]=Box1Abb[272]*m_x + Box1Abb[78]*m_x_3 + 3.*Box1Abb[40]*m_x_2*m_z12 + Box1Abb[40]*Box1Abb[73]*m_z12*m_z1k;

  Box1Abb[274]=-6. + 7.*m_z12 - 24.*m_z1k;

  Box1Abb[275]=19. + 3.*m_z12;

  Box1Abb[276]=23. - Box1Abb[275]*m_z12 + 6.*m_z1k - 20.*m_z12*m_z1k + 60.*m_z1k_2;

  Box1Abb[277]=2. + 3.*m_z1k;

  Box1Abb[278]=3. + Box1Abb[277]*m_z1k;

  Box1Abb[279]=-4. + 13.*m_z1k;

  Box1Abb[280]=-1. + Box1Abb[279]*m_z1k;

  Box1Abb[281]=5. - 6.*m_z1k + 4.*m_z1k_2;

  Box1Abb[282]=1. + Box1Abb[281]*m_z1k;

  Box1Abb[283]=-Box1Abb[282]*pow(Box1Abb[40],2.) + Box1Abb[280]*Box1Abb[40]*m_z12 + Box1Abb[278]*m_z12_2*m_z1k;

  Box1Abb[284]=3. + m_z1k + 4.*m_z1k_2;

  Box1Abb[285]=-7. + 4.*m_z1k;

  Box1Abb[286]=5. + Box1Abb[285]*m_z1k;

  Box1Abb[287]=2. + 3.*Box1Abb[286]*m_z1k;

  Box1Abb[288]=-35. + 19.*m_z1k;

  Box1Abb[289]=8. + Box1Abb[288]*Box1Abb[40]*m_z1k;

  Box1Abb[290]=8. + Box1Abb[289]*m_z1k;

  Box1Abb[291]=2.*pow(Box1Abb[40],2.) - Box1Abb[290]*m_z12_2 + 2.*Box1Abb[284]*Box1Abb[70]*m_z12_3 - 2.*Box1Abb[287]*Box1Abb[40]*m_z12*m_z1k + 3.*Box1Abb[70]*m_z12_4*m_z1k;

  Box1Abb[292]=5. + 7.*m_z1k;

  Box1Abb[293]=29. + 11.*m_z1k;

  Box1Abb[294]=13. + Box1Abb[293]*m_z1k;

  Box1Abb[295]=-9. + 20.*m_z1k;

  Box1Abb[296]=12. + Box1Abb[295]*m_z1k;

  Box1Abb[297]=-34. + Box1Abb[294]*m_z12 + 2.*Box1Abb[292]*m_z12_2 - 4.*Box1Abb[296]*m_z1k;

  Box1Abb[298]=14. + Box1Abb[297]*m_z12;

  Box1Abb[299]=10. + 9.*m_z1k;

  Box1Abb[300]=6. + Box1Abb[299]*m_z1k;

  Box1Abb[301]=7. - 5.*m_z1k;

  Box1Abb[302]=-31. + 6.*Box1Abb[301]*m_z1k;

  Box1Abb[303]=2. + Box1Abb[302]*m_z1k;

  Box1Abb[304]=-17. + 8.*m_z1k;

  Box1Abb[305]=17. + 2.*Box1Abb[304]*m_z1k;

  Box1Abb[306]=5. + Box1Abb[305]*m_z1k;

  Box1Abb[307]=-14. - Box1Abb[306]*m_z12 + 2.*Box1Abb[300]*m_z12_2 + 2.*Box1Abb[303]*m_z1k + m_z12_3*m_z1k;

  Box1Abb[308]=8. + Box1Abb[307]*m_z12;

  Box1Abb[309]=Box1Abb[291]*m_x_2 - Box1Abb[308]*m_x_3 + Box1Abb[298]*m_x_4 - Box1Abb[283]*Box1Abb[73]*m_x*m_z12 + Box1Abb[276]*m_x_5*m_z12 + Box1Abb[274]*m_x_6*m_z12 + 4.*m_x_7*m_z12 + pow(Box1Abb[40],2.)*Box1Abb[70]*pow(Box1Abb[73],2.)*m_z12_2*m_z1k;

  Box1Abb[310]=-10. + 3.*m_z12 - 24.*m_z1k;

  Box1Abb[311]=17. - 10.*m_z12 + 26.*m_z1k + 60.*m_z1k_2;

  Box1Abb[312]=3. - 2.*m_z1k + 4.*m_z1k_2;

  Box1Abb[313]=4. - 6.*m_z1k + 8.*m_z1k_2;

  Box1Abb[314]=Box1Abb[312]*Box1Abb[40] + Box1Abb[313]*m_z12 + Box1Abb[244]*m_z12_2;

  Box1Abb[315]=12. - 29.*m_z1k_2;

  Box1Abb[316]=4. + m_z1k + 20.*m_z1k_2;

  Box1Abb[317]=-23. + Box1Abb[315]*m_z12 - 4.*Box1Abb[316]*m_z1k + m_z12_2*m_z1k;

  Box1Abb[318]=6. + Box1Abb[317]*m_z12;

  Box1Abb[319]=-3. + 4.*m_z1k;

  Box1Abb[320]=-6. + 8.*m_z1k + 56.*m_z1k_3;

  Box1Abb[321]=-1. + 2.*m_z1k - 44.*m_z1k_2 + 60.*m_z1k_3;

  Box1Abb[322]=15. + Box1Abb[320]*m_z12 + Box1Abb[321]*m_z1k + Box1Abb[319]*m_z12_2*m_z1k;

  Box1Abb[323]=8.*Box1Abb[40] + Box1Abb[322]*m_z12;

  Box1Abb[324]=-3. + m_z1k + 10.*m_z1k_2;

  Box1Abb[325]=2. - 36.*m_z1k + 39.*m_z1k_2;

  Box1Abb[326]=12. + Box1Abb[325]*m_z1k;

  Box1Abb[327]=-1. + Box1Abb[326]*m_z1k;

  Box1Abb[328]=-11. + 12.*m_z1k;

  Box1Abb[329]=1. + Box1Abb[328]*m_z1k;

  Box1Abb[330]=7. + 2.*Box1Abb[329]*m_z1k;

  Box1Abb[331]=-3. + Box1Abb[330]*m_z1k;

  Box1Abb[332]=-2.*pow(Box1Abb[40],2.) + Box1Abb[331]*Box1Abb[40]*m_z12 + Box1Abb[327]*m_z12_2 + Box1Abb[324]*m_z12_3*m_z1k;

  Box1Abb[333]=-Box1Abb[332]*m_x_2 + Box1Abb[323]*m_x_3 + Box1Abb[318]*m_x_4 + Box1Abb[311]*m_x_5*m_z12 + Box1Abb[310]*m_x_6*m_z12 + 4.*m_x_7*m_z12 + Box1Abb[314]*pow(Box1Abb[40],2.)*m_x*m_z12*m_z1k + pow(Box1Abb[40],3.)*Box1Abb[73]*m_z12_2*m_z1k_2;

  Box1Abb[355]=2.*m_z12;

  Box1Abb[356]=-2.*m_z12;

  Box1Abb[357]=2.*m_s;

  Box1Abb[358]=(8.*m_x)/m_s;

  Box1Abb[359]=(-8.*m_x)/m_s;

  Box1Abb[360]=4./m_s;

  Box1Abb[361]=(-8.*m_x*m_z1k)/m_s;

  Box1Abb[362]=-4./m_s;

  Box1Abb[363]=-1. - Box1Abb[7]*m_x + m_z12 + m_z1k + m_z12*m_z1k;

  return;
}

DivArrC Higgs_Decay_RV_Box_1::RV_Box_1(const int& ME, const int& LR, const Vec4C& epsP)
{
  // Box corrections for emission off leg 1
  // ME takes values 1..8, denoting standard MEs
  // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
  // ubar1 P_i v2
  if (ME == 1) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)
	*((Box1Abb[235]*Box1Abb[358]*(epsP*m_p1))/(Box1Abb[233]*Box1Abb[4]) 
	  + (Box1Abb[236]*Box1Abb[358]*(epsP*m_p2))/(Box1Abb[233]*Box1Abb[4]))

	+B_0(m_s,m_m2,m_m2,m_mu2)
	*((Box1Abb[238]*Box1Abb[359]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[71]) 
	  + (Box1Abb[239]*Box1Abb[359]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[71]))

	+B_0(m_s1k,0.,m_m2,m_mu2)
	*((Box1Abb[246]*Box1Abb[360]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71])
	  + (Box1Abb[361]*Box1Abb[8]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[71]))

	+B_0(m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[270]*Box1Abb[362]*(epsP*m_p1))/(Box1Abb[233]*Box1Abb[4]*Box1Abb[5]*Box1Abb[71])
	  + (Box1Abb[273]*Box1Abb[358]*(epsP*m_p2))/(Box1Abb[233]*Box1Abb[4]*Box1Abb[71]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	*((Box1Abb[0]*Box1Abb[223]*Box1Abb[355]*Box1Abb[68]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (Box1Abb[0]*Box1Abb[223]*Box1Abb[355]*Box1Abb[69]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[223]*Box1Abb[356]*Box1Abb[5]*Box1Abb[68]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[223]*Box1Abb[356]*Box1Abb[5]*Box1Abb[69]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	*((2.*Box1Abb[255]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (2.*Box1Abb[232]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((2.*Box1Abb[309]*(epsP*m_p1))/(pow(Box1Abb[4],2.)*Box1Abb[71]) 
	  + (2.*Box1Abb[333]*(epsP*m_p2))/(pow(Box1Abb[4],2.)*Box1Abb[71]))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[224]*Box1Abb[232]*Box1Abb[357]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[232]*Box1Abb[357]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)
	*((Box1Abb[235]*Box1Abb[358]*(epsP*m_p1))/(Box1Abb[233]*Box1Abb[4]) 
	  + (Box1Abb[236]*Box1Abb[358]*(epsP*m_p2))/(Box1Abb[233]*Box1Abb[4]))

	+B_0(m_s,m_m2,m_m2,m_mu2)
	*((Box1Abb[238]*Box1Abb[359]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[71]) 
	  + (Box1Abb[239]*Box1Abb[359]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[71]))

	+B_0(m_s1k,0.,m_m2,m_mu2)
	*((Box1Abb[246]*Box1Abb[360]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71])
	  + (Box1Abb[361]*Box1Abb[363]*Box1Abb[69]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71]))

	+B_0(m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[270]*Box1Abb[362]*(epsP*m_p1))/(Box1Abb[233]*Box1Abb[4]*Box1Abb[5]*Box1Abb[71]) 
	  + (Box1Abb[273]*Box1Abb[358]*(epsP*m_p2))/(Box1Abb[233]*Box1Abb[4]*Box1Abb[71]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	*((Box1Abb[0]*Box1Abb[223]*Box1Abb[355]*Box1Abb[68]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (Box1Abb[0]*Box1Abb[223]*Box1Abb[355]*Box1Abb[69]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[223]*Box1Abb[356]*Box1Abb[5]*Box1Abb[68]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[223]*Box1Abb[356]*Box1Abb[5]*Box1Abb[69]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	*((2.*Box1Abb[255]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (2.*Box1Abb[232]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((2.*Box1Abb[309]*(epsP*m_p1))/(pow(Box1Abb[4],2.)*Box1Abb[71]) 
	  + (2.*Box1Abb[333]*(epsP*m_p2))/(pow(Box1Abb[4],2.)*Box1Abb[71]))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[224]*Box1Abb[232]*Box1Abb[357]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[232]*Box1Abb[357]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.));
    }
  }
  // ubar1 \slashed{k} P_i v2
  else if (ME == 2) {
    if (LR == 0) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)
	*((Box1Abb[194]*Box1Abb[68]*(epsP*m_p1))/(Box1Abb[0]*Box1Abb[4]) 
	  + (Box1Abb[194]*Box1Abb[69]*(epsP*m_p2))/(Box1Abb[0]*Box1Abb[4]))

	+B_0(m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[195]*Box1Abb[75]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71]) 
	  + (Box1Abb[195]*Box1Abb[76]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[71]))

	+B_0(m_s1k,0.,m_m2,m_mu2)
	*((Box1Abb[196]*Box1Abb[79]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71])
	  + (Box1Abb[197]*Box1Abb[69]*Box1Abb[80]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71]))

	+B_0(m_s,m_m2,m_m2,m_mu2)
	*((Box1Abb[198]*Box1Abb[84]*(epsP*m_p1))/(Box1Abb[33]*Box1Abb[4]*Box1Abb[71]) 
	  + (Box1Abb[198]*Box1Abb[88]*(epsP*m_p2))/(Box1Abb[33]*Box1Abb[4]*Box1Abb[71]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	*((Box1Abb[199]*Box1Abb[33]*Box1Abb[94]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (Box1Abb[102]*Box1Abb[199]*Box1Abb[33]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[199]*Box1Abb[5]*Box1Abb[94]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (Box1Abb[102]*Box1Abb[199]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[109]*Box1Abb[200]*(epsP*m_p1))/pow(Box1Abb[4],2.) 
	  + (Box1Abb[117]*Box1Abb[200]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[149]*Box1Abb[200]*(epsP*m_p1))/(pow(Box1Abb[4],2.)*Box1Abb[71]) 
	  + (Box1Abb[174]*Box1Abb[200]*(epsP*m_p2))/(pow(Box1Abb[4],2.)*Box1Abb[71]))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[109]*Box1Abb[201]*Box1Abb[5]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[117]*Box1Abb[201]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.));
    }
    else if (LR == 1) {
      return B_0(m_s12,m_m2,m_m2,m_mu2)
	*((Box1Abb[194]*Box1Abb[68]*(epsP*m_p1))/(Box1Abb[0]*Box1Abb[4])
	  + (Box1Abb[194]*Box1Abb[69]*(epsP*m_p2))/(Box1Abb[0]*Box1Abb[4]))

	+B_0(m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[195]*Box1Abb[75]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71])
	  + (Box1Abb[195]*Box1Abb[76]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[71]))

	+B_0(m_s1k,0.,m_m2,m_mu2)
	*((Box1Abb[196]*Box1Abb[79]*(epsP*m_p1))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71]) 
	  + (Box1Abb[197]*Box1Abb[69]*Box1Abb[80]*(epsP*m_p2))/(Box1Abb[4]*Box1Abb[5]*Box1Abb[71]))

	+B_0(m_s,m_m2,m_m2,m_mu2)
	*((Box1Abb[198]*Box1Abb[84]*(epsP*m_p1))/(Box1Abb[33]*Box1Abb[4]*Box1Abb[71]) 
	  + (Box1Abb[198]*Box1Abb[88]*(epsP*m_p2))/(Box1Abb[33]*Box1Abb[4]*Box1Abb[71]))

	+C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	*((Box1Abb[199]*Box1Abb[33]*Box1Abb[94]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[102]*Box1Abb[199]*Box1Abb[33]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[199]*Box1Abb[5]*Box1Abb[94]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[102]*Box1Abb[199]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	*((Box1Abb[109]*Box1Abb[200]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[117]*Box1Abb[200]*(epsP*m_p2))/pow(Box1Abb[4],2.))

	+C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[149]*Box1Abb[200]*(epsP*m_p1))/(pow(Box1Abb[4],2.)*Box1Abb[71]) 
	  + (Box1Abb[174]*Box1Abb[200]*(epsP*m_p2))/(pow(Box1Abb[4],2.)*Box1Abb[71]))

	+D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2)
	*((Box1Abb[109]*Box1Abb[201]*Box1Abb[5]*(epsP*m_p1))/pow(Box1Abb[4],2.)
	  + (Box1Abb[117]*Box1Abb[201]*Box1Abb[5]*(epsP*m_p2))/pow(Box1Abb[4],2.));
    }
  }
  // ubar1 \slashed{m_epsP} P_i v2
  else if (ME == 3) {
    if (LR == 0) {
      return (Box1Abb[33]*Box1Abb[37]*Box1Abb[59]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[37]*Box1Abb[5]*Box1Abb[59]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[45]*Box1Abb[59]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[49]*Box1Abb[59]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[5]*Box1Abb[51]*Box1Abb[60]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[4];
    }
    else if (LR == 1) {
      return (Box1Abb[33]*Box1Abb[37]*Box1Abb[59]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[37]*Box1Abb[5]*Box1Abb[59]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[45]*Box1Abb[59]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[49]*Box1Abb[59]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[5]*Box1Abb[51]*Box1Abb[60]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[4];
    }
  }
  // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
  else if (ME == 4) {
    if (LR == 0) {
      return (pow(Box1Abb[0],2.)*Box1Abb[1]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[4]

	+(-Box1Abb[0]*Box1Abb[1]*Box1Abb[5]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[6]*Box1Abb[8]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[1]*Box1Abb[23]*Box1Abb[9]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[0]*Box1Abb[1]*Box1Abb[11]*Box1Abb[24]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[4];
    }
    else if (LR == 1) {
      return (pow(Box1Abb[0],2.)*Box1Abb[1]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box1Abb[4]

	+(-Box1Abb[0]*Box1Abb[1]*Box1Abb[5]*C_0(0.,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[6]*Box1Abb[8]*C_0(m_s,m_m2,m_s1k,m_m2,m_m2,0.,m_mu2))/Box1Abb[4]

	+(Box1Abb[1]*Box1Abb[23]*Box1Abb[9]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box1Abb[4]

	+(Box1Abb[0]*Box1Abb[1]*Box1Abb[11]*Box1Abb[24]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s1k,m_m2,m_m2,m_m2,0.,m_mu2))/Box1Abb[4];
    }
  }
  else {
    msg_Out() << "Standard matrix element not known. Value given: " << ME << "\n"
		<< "Values range from 1 to 4.";
  }
}

      
