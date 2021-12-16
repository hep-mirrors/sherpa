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
#include <algorithm>       
#include <assert.h> 

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)
#define D_0(A,B,C,D,E,F,G,H,I,J,M) Master_Box(A,B,C,D,E,F,G,H,I,J,M)


using namespace ATOOLS;
using namespace METOOLS;
using namespace PHOTONS;

Higgs_Decay_RV_Box_2::Higgs_Decay_RV_Box_2
(const double& m, const double& s, const Vec4D& p1, const Vec4D& p2, const Vec4D& pP,
 const Complex& cL, const Complex& cR, const double& mu2):
  Higgs_Decay_RV_Diagrams(m,s,p1,p2,pP,cL,cR,mu2)
{
  Box2Abb = new Complex[362];
  Init_Coefficients();
}

Higgs_Decay_RV_Box_2::~Higgs_Decay_RV_Box_2()
{
  delete [] Box2Abb;
}

void Higgs_Decay_RV_Box_2::Init_Coefficients() 
{
  Init_Box_2_Coefficients();
  return;
}

// Setup coefficients for box with emission of leg 2
// Calculated using FeynCalc 9.0.1
void Higgs_Decay_RV_Box_2::Init_Box_2_Coefficients() 
{
  Box2Abb[0]=-1. + m_z12;

  Box2Abb[1]=-2.*m_x + m_z12;

  Box2Abb[2]=-1. + m_x;

  Box2Abb[3]=-1. - 2.*m_x + m_z12;

  Box2Abb[4]=m_x + Box2Abb[2]*m_x*m_z12 + Box2Abb[3]*m_z12*m_z2k + m_z12*m_z2k_2;

  Box2Abb[5]=m_x - m_z2k;

  Box2Abb[6]=2.*m_x - m_z12;

  Box2Abb[7]=1. + m_z12;

  Box2Abb[8]=1. + m_x - m_z12 + m_x*m_z12 - Box2Abb[7]*m_z2k;

  Box2Abb[9]=-1. - 2.*m_x + m_z12 + 2.*m_z2k;

  Box2Abb[10]=-2. + m_z12;

  Box2Abb[11]=Box2Abb[10]*m_x + m_z12*m_z2k;

  Box2Abb[20]=m_s*m_z12;

  Box2Abb[21]=m_s*m_x;

  Box2Abb[22]=m_s*m_z2k;

  Box2Abb[23]=-m_z12;

  Box2Abb[24]=-m_s;

  Box2Abb[33]=1. - m_z12;

  Box2Abb[34]=2.*m_z12 + m_z2k;

  Box2Abb[35]=1. + 2.*m_z2k;

  Box2Abb[36]=-Box2Abb[35]*m_x + m_x_2 + Box2Abb[34]*m_z2k;

  Box2Abb[37]=pow(Box2Abb[0],2.) + 2.*Box2Abb[36] - 2.*m_z2k;

  Box2Abb[38]=-1. + 3.*m_z2k;

  Box2Abb[39]=4. + 4.*m_z12 + 6.*m_z2k;

  Box2Abb[40]=-1. + m_z2k;

  Box2Abb[41]=1. + 2.*Box2Abb[40]*m_z2k;

  Box2Abb[42]=2. + 4.*Box2Abb[40]*m_z2k;

  Box2Abb[43]=-4. + m_z12;

  Box2Abb[44]=1. + Box2Abb[43]*m_z12 - 8.*m_z12*m_z2k - 6.*m_z2k_2;

  Box2Abb[45]=Box2Abb[40]*Box2Abb[41] + Box2Abb[44]*m_x + Box2Abb[39]*m_x_2 - 2.*m_x_3 + Box2Abb[42]*m_z12 + Box2Abb[38]*m_z12_2;

  Box2Abb[46]=pow(Box2Abb[0],2.) - 2.*m_z2k;

  Box2Abb[47]=-3. + m_z12 - 2.*m_z2k;

  Box2Abb[48]=-2. - 2.*Box2Abb[47]*m_z12 + 4.*m_z2k;

  Box2Abb[49]=Box2Abb[48]*m_x - 4.*Box2Abb[7]*m_x_2 + Box2Abb[46]*m_z12;

  Box2Abb[50]=1. - 3.*m_z12 + m_z12_2 - 2.*Box2Abb[7]*m_z2k;

  Box2Abb[51]=2.*Box2Abb[50]*m_x + 4.*Box2Abb[7]*m_x_2 - pow(Box2Abb[0],2.)*m_z12 + 2.*m_z12*m_z2k;

  Box2Abb[59]=m_m;

  Box2Abb[60]=-m_m*m_s;

  Box2Abb[68]=-m_x + m_z2k;

  Box2Abb[69]=-1. - m_x + m_z12 + m_z2k;

  Box2Abb[70]=-1. + m_x + m_z2k + 2.*m_z12*m_z2k;

  Box2Abb[71]=1. + m_z2k;

  Box2Abb[72]=pow(Box2Abb[40],2.) - 2.*Box2Abb[71]*m_x + m_x_2;

  Box2Abb[73]=1. + 2.*m_z12;

  Box2Abb[74]=-1. + m_z12 + m_z2k;

  Box2Abb[75]=2. - 3.*m_z12 - 2.*m_z12*m_z2k;

  Box2Abb[76]=-Box2Abb[40]*Box2Abb[74] + Box2Abb[75]*m_x + Box2Abb[73]*m_x_2;

  Box2Abb[77]=2. + m_z12;

  Box2Abb[78]=Box2Abb[77]*m_x + Box2Abb[40]*m_z12;

  Box2Abb[79]=3. + 2.*m_x - 2.*m_z12;

  Box2Abb[80]=-1. + m_x + pow(Box2Abb[2],2.)*m_z12 + Box2Abb[79]*m_z2k - Box2Abb[77]*m_z2k_2;

  Box2Abb[81]=-2. + m_z12 + m_z2k;

  Box2Abb[82]=pow(Box2Abb[0],2.) + 3.*m_z12*m_z2k;

  Box2Abb[83]=1. + m_z12 + 3.*m_z12*m_z2k;

  Box2Abb[84]=m_x - Box2Abb[83]*m_x_2 + m_x_3*m_z12 + Box2Abb[82]*m_x*m_z2k - Box2Abb[40]*Box2Abb[81]*m_z12*m_z2k;

  Box2Abb[85]=1. + m_z12 + Box2Abb[10]*m_z2k + 3.*m_z2k_2;

  Box2Abb[86]=-1. + Box2Abb[85]*m_z12 + m_z2k;

  Box2Abb[87]=Box2Abb[86]*m_x - Box2Abb[83]*m_x_2 - pow(Box2Abb[40],2.)*Box2Abb[74]*m_z12 + m_x_3*m_z12;

  Box2Abb[88]=-2. + 4.*m_z2k;

  Box2Abb[89]=pow(Box2Abb[40],2.) + Box2Abb[88]*m_z12 + m_z12_2;

  Box2Abb[90]=-4. + 3.*m_z12 + 11.*m_z2k;

  Box2Abb[91]=3. + Box2Abb[90]*m_z12;

  Box2Abb[92]=2. + 7.*m_z2k;

  Box2Abb[93]=-2. + m_z12 + Box2Abb[92]*m_z2k + 3.*m_z12*m_z2k;

  Box2Abb[94]=1. + Box2Abb[93]*m_z12 - m_z2k;

  Box2Abb[95]=Box2Abb[94]*m_x - Box2Abb[91]*m_x_2 + 5.*m_x_3*m_z12 - Box2Abb[89]*m_z12*m_z2k;

  Box2Abb[96]=-13. + 8.*m_z12 + 5.*m_z2k;

  Box2Abb[97]=5. + Box2Abb[96]*m_z12;

  Box2Abb[98]=-4. + 3.*m_z12;

  Box2Abb[99]=2.*pow(Box2Abb[0],2.) + Box2Abb[98]*m_z2k + m_z2k_2;

  Box2Abb[100]=Box2Abb[99]*m_z12 + m_z2k;

  Box2Abb[101]=Box2Abb[100]*m_x - Box2Abb[97]*m_x_2 + 3.*m_x_3*m_z12 + Box2Abb[74]*m_z12*m_z2k_2;

  Box2Abb[102]=18. - 11.*m_z12;

  Box2Abb[103]=pow(Box2Abb[40],2.) + Box2Abb[10]*m_z12;

  Box2Abb[104]=2. + 3.*m_z12;

  Box2Abb[105]=pow(Box2Abb[0],2.) + m_z2k + 5.*Box2Abb[10]*m_z12*m_z2k - Box2Abb[104]*m_z2k_2;

  Box2Abb[106]=1. + 4.*m_z2k;

  Box2Abb[107]=-2. + 3.*m_z12 + 7.*m_z2k;

  Box2Abb[108]=-5.*Box2Abb[106] + Box2Abb[107]*m_z12;

  Box2Abb[109]=2. + Box2Abb[108]*m_z12;

  Box2Abb[110]=Box2Abb[109]*m_x_2 - Box2Abb[105]*m_x*m_z12 + Box2Abb[102]*m_x_3*m_z12 + Box2Abb[103]*m_z12_2*m_z2k;

  Box2Abb[111]=-14. + 15.*m_z12;

  Box2Abb[112]=-1. + m_z12_2;

  Box2Abb[113]=2.*pow(Box2Abb[0],2.)*m_z12 + 3.*Box2Abb[112]*m_z2k + Box2Abb[77]*m_z2k_2;

  Box2Abb[114]=-7. + 4.*m_z12 + 5.*m_z2k;

  Box2Abb[115]=-11. - 3.*Box2Abb[114]*m_z12 + 12.*m_z2k;

  Box2Abb[116]=2. + Box2Abb[115]*m_z12;

  Box2Abb[117]=-Box2Abb[116]*m_x_2 - Box2Abb[113]*m_x*m_z12 - Box2Abb[111]*m_x_3*m_z12 + Box2Abb[74]*m_z12_2*m_z2k_2;

  Box2Abb[118]=-11. + 3.*m_z12 - 26.*m_z2k;

  Box2Abb[119]=-3. + Box2Abb[118]*m_z12;

  Box2Abb[120]=pow(Box2Abb[40],5.) + Box2Abb[10]*pow(Box2Abb[40],3.)*m_z12;

  Box2Abb[121]=-5. + 4.*m_z12;

  Box2Abb[122]=-14. + m_z12;

  Box2Abb[123]=-3. + Box2Abb[122]*m_z12;

  Box2Abb[124]=-3. + m_z12 - 21.*m_z12_2;

  Box2Abb[125]=-6. + 3.*Box2Abb[121]*m_z12 + Box2Abb[123]*m_z12*m_z2k + 2.*Box2Abb[124]*m_z2k_2 - 60.*m_z12*m_z2k_3;

  Box2Abb[126]=-5. + 4.*m_z2k;

  Box2Abb[127]=23. + 55.*m_z2k;

  Box2Abb[128]=13. + 2.*Box2Abb[126]*m_z12 + Box2Abb[127]*m_z2k;

  Box2Abb[129]=10. + Box2Abb[128]*m_z12 + 8.*m_z2k;

  Box2Abb[130]=-3. + 7.*m_z2k;

  Box2Abb[131]=-3. + 5.*m_z2k + 24.*m_z2k_3;

  Box2Abb[132]=-4. + m_z2k - 34.*m_z2k_2 + 35.*m_z2k_3;

  Box2Abb[133]=10. + Box2Abb[132]*m_z2k;

  Box2Abb[134]=-2. + Box2Abb[133]*m_z12 + 2.*Box2Abb[131]*m_z12_2 + 8.*m_z2k + Box2Abb[130]*m_z12_3*m_z2k - 6.*m_z2k_2;

  Box2Abb[135]=2. + m_z2k - 5.*m_z2k_2 + 10.*m_z2k_3;

  Box2Abb[136]=2. + 9.*m_z2k;

  Box2Abb[137]=3. - Box2Abb[136]*m_z2k;

  Box2Abb[138]=-9. + 17.*m_z2k;

  Box2Abb[139]=-5. + Box2Abb[138]*m_z2k;

  Box2Abb[140]=1. + Box2Abb[139]*m_z2k;

  Box2Abb[141]=pow(Box2Abb[40],4.) - Box2Abb[135]*pow(Box2Abb[40],2.)*m_z12 - Box2Abb[140]*Box2Abb[40]*m_z12_2 + Box2Abb[137]*m_z12_3*m_z2k;

  Box2Abb[142]=Box2Abb[141]*m_x + Box2Abb[134]*m_x_2 + Box2Abb[125]*m_x_3 + Box2Abb[129]*m_x_4 + Box2Abb[119]*m_x_5 + 5.*m_x_6*m_z12 + Box2Abb[120]*m_z12*m_z2k;

  Box2Abb[143]=2. + 4.*m_z12 - 7.*m_z2k;

  Box2Abb[144]=-5. + 2.*Box2Abb[143]*m_z12;

  Box2Abb[145]=-1. + 2.*m_z12 + m_z2k;

  Box2Abb[146]=28. - 25.*m_z2k;

  Box2Abb[147]=4. + 3.*m_z2k;

  Box2Abb[148]=-4. + 5.*Box2Abb[147]*m_z12 + 2.*m_z12_2 + Box2Abb[146]*m_z2k;

  Box2Abb[149]=15. - Box2Abb[148]*m_z12 + 16.*m_z2k;

  Box2Abb[150]=5. + m_z2k;

  Box2Abb[151]=1. + Box2Abb[150]*Box2Abb[40]*m_z2k;

  Box2Abb[152]=-15. + 23.*m_z2k;

  Box2Abb[153]=4. + Box2Abb[152]*m_z2k;

  Box2Abb[154]=-7. + m_z2k;

  Box2Abb[155]=3. + Box2Abb[154]*m_z2k;

  Box2Abb[156]=-1. + Box2Abb[155]*m_z2k;

  Box2Abb[157]=2.*Box2Abb[151]*pow(Box2Abb[40],2.)*m_z12 + Box2Abb[153]*Box2Abb[40]*m_z12_2 - 2.*Box2Abb[156]*m_z12_3 - pow(Box2Abb[40],3.)*m_z2k;

  Box2Abb[158]=1. + 6.*m_z2k;

  Box2Abb[159]=12. + 25.*m_z2k;

  Box2Abb[160]=-14. + 5.*m_z2k;

  Box2Abb[161]=5. + Box2Abb[160]*m_z2k;

  Box2Abb[162]=9. + 2.*Box2Abb[161]*m_z2k;

  Box2Abb[163]=1. - 2.*Box2Abb[162]*m_z12 + Box2Abb[159]*m_z12_2 + 6.*Box2Abb[71]*m_z12_3 - 3.*Box2Abb[158]*m_z2k;

  Box2Abb[164]=-7. + 8.*m_z2k;

  Box2Abb[165]=3. + Box2Abb[164]*m_z2k;

  Box2Abb[166]=-37. + 10.*m_z2k;

  Box2Abb[167]=-5. + Box2Abb[166]*m_z2k;

  Box2Abb[168]=4. + Box2Abb[167]*m_z2k;

  Box2Abb[169]=-8. + m_z2k;

  Box2Abb[170]=46. + 5.*Box2Abb[169]*m_z2k;

  Box2Abb[171]=-8. + Box2Abb[170]*m_z2k;

  Box2Abb[172]=5. + Box2Abb[171]*m_z2k;

  Box2Abb[173]=Box2Abb[165]*Box2Abb[40] + Box2Abb[172]*m_z12 + Box2Abb[168]*m_z12_2 - 6.*m_z12_3;

  Box2Abb[174]=Box2Abb[157]*m_x + Box2Abb[173]*m_x_2 + Box2Abb[163]*m_x_3 + Box2Abb[149]*m_x_4 + Box2Abb[144]*m_x_5 + 3.*m_x_6*m_z12 - Box2Abb[145]*pow(Box2Abb[40],2.)*Box2Abb[74]*m_z12*m_z2k_2;

  Box2Abb[194]=(-8.*m_m*m_z12)/m_s_2;

  Box2Abb[195]=(8.*m_m*m_x)/m_s_2;

  Box2Abb[196]=(-8.*m_m*m_z2k)/m_s_2;

  Box2Abb[197]=(-8.*m_m*m_x)/m_s_2;

  Box2Abb[198]=(8.*m_m)/m_s_2;

  Box2Abb[199]=(-2.*m_m)/m_s;

  Box2Abb[200]=(2.*m_m)/m_s;

  Box2Abb[201]=2.*m_m;

  Box2Abb[221]=-2. + 3.*m_z12 + 8.*m_z2k;

  Box2Abb[222]=-1. + m_z12 - 2.*m_z2k + 4.*m_z12*m_z2k + 4.*m_z2k_2;

  Box2Abb[223]=-Box2Abb[222]*m_x + Box2Abb[221]*m_x_2 - 4.*m_x_3 + Box2Abb[74]*m_z12*m_z2k;

  Box2Abb[224]=3. - 2.*m_z12;

  Box2Abb[225]=-3. + m_z12;

  Box2Abb[226]=-16. + 5.*m_z12;

  Box2Abb[227]=2. + Box2Abb[225]*m_z12 + 8.*m_z2k + Box2Abb[226]*m_z12*m_z2k - 12.*m_z2k_2;

  Box2Abb[228]=3. + Box2Abb[43]*m_z12 - 3.*m_z12*m_z2k - 4.*m_z2k_2;

  Box2Abb[229]=2. + 3.*m_z12 + 12.*m_z2k;

  Box2Abb[230]=4. + Box2Abb[10]*Box2Abb[229]*m_z12;

  Box2Abb[231]=Box2Abb[230]*m_x_3 - Box2Abb[227]*m_x_2*m_z12 + 4.*Box2Abb[224]*m_x_4*m_z12 + Box2Abb[228]*m_x*m_z12_2*m_z2k + Box2Abb[74]*m_z12_3*m_z2k_2;

  Box2Abb[232]=1. + m_x - m_z12 - m_z2k;

  Box2Abb[233]=4.*m_x - m_z12;

  Box2Abb[234]=-2. + 3.*m_z12;

  Box2Abb[235]=m_z12 + m_z2k;

  Box2Abb[236]=Box2Abb[234]*m_x + m_z12 - Box2Abb[235]*m_z12;

  Box2Abb[237]=-Box2Abb[35]*m_x + m_x_2 + Box2Abb[145]*m_z2k;

  Box2Abb[238]=m_z12 - 2.*m_z2k;

  Box2Abb[239]=Box2Abb[40]*Box2Abb[74] + Box2Abb[238]*m_x + m_x_2;

  Box2Abb[240]=-1. + m_z12 - 4.*m_z2k + 7.*m_z12*m_z2k;

  Box2Abb[241]=-16. + 13.*m_z12;

  Box2Abb[242]=Box2Abb[241]*m_z12 + 24.*Box2Abb[0]*m_z2k;

  Box2Abb[243]=4. + Box2Abb[242]*m_z12;

  Box2Abb[244]=2. - 23.*m_z2k;

  Box2Abb[245]=8. - 3.*m_z2k;

  Box2Abb[246]=3. + 4.*Box2Abb[245]*m_z2k;

  Box2Abb[247]=-2. + Box2Abb[246]*m_z12 + Box2Abb[244]*m_z12_2 - 3.*m_z12_3 + 12.*Box2Abb[40]*m_z2k;

  Box2Abb[248]=Box2Abb[243]*m_x_3 + Box2Abb[247]*m_x_2*m_z12 - 12.*Box2Abb[0]*m_x_4*m_z12 + Box2Abb[240]*Box2Abb[74]*m_x*m_z12_2 - pow(Box2Abb[74],2.)*m_z12_3*m_z2k;

  Box2Abb[249]=1. - m_z12 + 2.*m_z2k;

  Box2Abb[250]=-4. + m_z12 + m_z12_2;

  Box2Abb[251]=m_z12 + Box2Abb[250]*m_z2k + 4.*Box2Abb[7]*m_z2k_2;

  Box2Abb[252]=-1. + 2.*m_z2k;

  Box2Abb[253]=-1. + 4.*m_z2k;

  Box2Abb[254]=Box2Abb[252]*pow(Box2Abb[40],2.) + Box2Abb[253]*Box2Abb[40]*Box2Abb[71]*m_z12 + 2.*Box2Abb[71]*m_z12_2*m_z2k;

  Box2Abb[255]=Box2Abb[254]*m_x - Box2Abb[251]*m_x_2 + Box2Abb[249]*m_x_3 + m_x_4*m_z12 - pow(Box2Abb[40],2.)*Box2Abb[74]*m_z12*m_z2k;

  Box2Abb[256]=2. + m_z2k - 5.*m_z2k_2;

  Box2Abb[257]=-2.*pow(Box2Abb[40],2.) + Box2Abb[256]*m_z12 - m_z12_2*m_z2k;

  Box2Abb[258]=Box2Abb[257]*m_x + Box2Abb[77]*m_x_3 + 3.*Box2Abb[40]*m_x_2*m_z12 + Box2Abb[40]*Box2Abb[74]*m_z12*m_z2k;

  Box2Abb[259]=2. + 5.*m_z12;

  Box2Abb[260]=-11. + 3.*m_z12;

  Box2Abb[261]=-1. + 2.*m_z12;

  Box2Abb[262]=4. + Box2Abb[260]*m_z12 + 2.*Box2Abb[0]*m_z12*m_z2k + 4.*Box2Abb[261]*m_z2k_2;

  Box2Abb[263]=-3. + m_z2k;

  Box2Abb[264]=8. + 3.*m_z12 + 20.*m_z2k;

  Box2Abb[265]=4.*Box2Abb[263] + Box2Abb[264]*m_z12;

  Box2Abb[266]=-7. + 6.*m_z2k;

  Box2Abb[267]=5. + 2.*Box2Abb[263]*m_z2k;

  Box2Abb[268]=-1. + m_z12 + 2.*Box2Abb[71]*m_z12_2*m_z2k + Box2Abb[267]*m_z2k_2 + Box2Abb[266]*m_z12*m_z2k_2;

  Box2Abb[269]=1. + 9.*m_z2k + 6.*m_z2k_2;

  Box2Abb[270]=-4. + m_z2k + m_z2k_2;

  Box2Abb[271]=1. + Box2Abb[270]*m_z2k;

  Box2Abb[272]=4.*pow(Box2Abb[40],3.) + 4.*Box2Abb[271]*m_z12 + Box2Abb[269]*m_z12_2 + m_z12_3*m_z2k;

  Box2Abb[273]=-Box2Abb[272]*m_x_2 - Box2Abb[262]*m_x_3 + Box2Abb[265]*m_x_4 - 2.*Box2Abb[259]*m_x_5 + Box2Abb[268]*m_x*m_z12 - pow(Box2Abb[40],2.)*Box2Abb[74]*m_z12_2*m_z2k;

  Box2Abb[274]=-10. + 3.*m_z12 - 24.*m_z2k;

  Box2Abb[275]=17. - 10.*m_z12 + 26.*m_z2k + 60.*m_z2k_2;

  Box2Abb[276]=3. - 2.*m_z2k + 4.*m_z2k_2;

  Box2Abb[277]=4. - 6.*m_z2k + 8.*m_z2k_2;

  Box2Abb[278]=Box2Abb[276]*Box2Abb[40] + Box2Abb[277]*m_z12 + Box2Abb[253]*m_z12_2;

  Box2Abb[279]=12. - 29.*m_z2k_2;

  Box2Abb[280]=4. + m_z2k + 20.*m_z2k_2;

  Box2Abb[281]=-23. + Box2Abb[279]*m_z12 - 4.*Box2Abb[280]*m_z2k + m_z12_2*m_z2k;

  Box2Abb[282]=6. + Box2Abb[281]*m_z12;

  Box2Abb[283]=-3. + 4.*m_z2k;

  Box2Abb[284]=-6. + 8.*m_z2k + 56.*m_z2k_3;

  Box2Abb[285]=-1. + 2.*m_z2k - 44.*m_z2k_2 + 60.*m_z2k_3;

  Box2Abb[286]=15. + Box2Abb[284]*m_z12 + Box2Abb[285]*m_z2k + Box2Abb[283]*m_z12_2*m_z2k;

  Box2Abb[287]=8.*Box2Abb[40] + Box2Abb[286]*m_z12;

  Box2Abb[288]=-3. + m_z2k + 10.*m_z2k_2;

  Box2Abb[289]=2. - 36.*m_z2k + 39.*m_z2k_2;

  Box2Abb[290]=12. + Box2Abb[289]*m_z2k;

  Box2Abb[291]=-1. + Box2Abb[290]*m_z2k;

  Box2Abb[292]=-11. + 12.*m_z2k;

  Box2Abb[293]=1. + Box2Abb[292]*m_z2k;

  Box2Abb[294]=7. + 2.*Box2Abb[293]*m_z2k;

  Box2Abb[295]=-3. + Box2Abb[294]*m_z2k;

  Box2Abb[296]=-2.*pow(Box2Abb[40],2.) + Box2Abb[295]*Box2Abb[40]*m_z12 + Box2Abb[291]*m_z12_2 + Box2Abb[288]*m_z12_3*m_z2k;

  Box2Abb[297]=-Box2Abb[296]*m_x_2 + Box2Abb[287]*m_x_3 + Box2Abb[282]*m_x_4 + Box2Abb[275]*m_x_5*m_z12 + Box2Abb[274]*m_x_6*m_z12 + 4.*m_x_7*m_z12 + Box2Abb[278]*pow(Box2Abb[40],2.)*m_x*m_z12*m_z2k + pow(Box2Abb[40],3.)*Box2Abb[74]*m_z12_2*m_z2k_2;

  Box2Abb[298]=-6. + 7.*m_z12 - 24.*m_z2k;

  Box2Abb[299]=19. + 3.*m_z12;

  Box2Abb[300]=23. - Box2Abb[299]*m_z12 + 6.*m_z2k - 20.*m_z12*m_z2k + 60.*m_z2k_2;

  Box2Abb[301]=2. + 3.*m_z2k;

  Box2Abb[302]=3. + Box2Abb[301]*m_z2k;

  Box2Abb[303]=-4. + 13.*m_z2k;

  Box2Abb[304]=-1. + Box2Abb[303]*m_z2k;

  Box2Abb[305]=5. - 6.*m_z2k + 4.*m_z2k_2;

  Box2Abb[306]=1. + Box2Abb[305]*m_z2k;

  Box2Abb[307]=-Box2Abb[306]*pow(Box2Abb[40],2.) + Box2Abb[304]*Box2Abb[40]*m_z12 + Box2Abb[302]*m_z12_2*m_z2k;

  Box2Abb[308]=3. + m_z2k + 4.*m_z2k_2;

  Box2Abb[309]=-7. + 4.*m_z2k;

  Box2Abb[310]=5. + Box2Abb[309]*m_z2k;

  Box2Abb[311]=2. + 3.*Box2Abb[310]*m_z2k;

  Box2Abb[312]=-35. + 19.*m_z2k;

  Box2Abb[313]=8. + Box2Abb[312]*Box2Abb[40]*m_z2k;

  Box2Abb[314]=8. + Box2Abb[313]*m_z2k;

  Box2Abb[315]=2.*pow(Box2Abb[40],2.) - Box2Abb[314]*m_z12_2 + 2.*Box2Abb[308]*Box2Abb[71]*m_z12_3 - 2.*Box2Abb[311]*Box2Abb[40]*m_z12*m_z2k + 3.*Box2Abb[71]*m_z12_4*m_z2k;

  Box2Abb[316]=5. + 7.*m_z2k;

  Box2Abb[317]=29. + 11.*m_z2k;

  Box2Abb[318]=13. + Box2Abb[317]*m_z2k;

  Box2Abb[319]=-9. + 20.*m_z2k;

  Box2Abb[320]=12. + Box2Abb[319]*m_z2k;

  Box2Abb[321]=-34. + Box2Abb[318]*m_z12 + 2.*Box2Abb[316]*m_z12_2 - 4.*Box2Abb[320]*m_z2k;

  Box2Abb[322]=14. + Box2Abb[321]*m_z12;

  Box2Abb[323]=10. + 9.*m_z2k;

  Box2Abb[324]=6. + Box2Abb[323]*m_z2k;

  Box2Abb[325]=7. - 5.*m_z2k;

  Box2Abb[326]=-31. + 6.*Box2Abb[325]*m_z2k;

  Box2Abb[327]=2. + Box2Abb[326]*m_z2k;

  Box2Abb[328]=-17. + 8.*m_z2k;

  Box2Abb[329]=17. + 2.*Box2Abb[328]*m_z2k;

  Box2Abb[330]=5. + Box2Abb[329]*m_z2k;

  Box2Abb[331]=-14. - Box2Abb[330]*m_z12 + 2.*Box2Abb[324]*m_z12_2 + 2.*Box2Abb[327]*m_z2k + m_z12_3*m_z2k;

  Box2Abb[332]=8. + Box2Abb[331]*m_z12;

  Box2Abb[333]=Box2Abb[315]*m_x_2 - Box2Abb[332]*m_x_3 + Box2Abb[322]*m_x_4 - Box2Abb[307]*Box2Abb[74]*m_x*m_z12 + Box2Abb[300]*m_x_5*m_z12 + Box2Abb[298]*m_x_6*m_z12 + 4.*m_x_7*m_z12 + pow(Box2Abb[40],2.)*Box2Abb[71]*pow(Box2Abb[74],2.)*m_z12_2*m_z2k;

  Box2Abb[355]=-2.*m_z12;

  Box2Abb[356]=2.*m_z12;

  Box2Abb[357]=-2.*m_s;

  Box2Abb[358]=(8.*m_x)/m_s;

  Box2Abb[359]=(8.*m_x*m_z2k)/m_s;

  Box2Abb[360]=4./m_s;

  Box2Abb[361]=(-8.*m_x)/m_s;

  return;
}


 DivArrC Higgs_Decay_RV_Box_2::RV_Box_2(const int& ME, const int& LR, const Vec4C& epsP)
 {
   // Box corrections for emission off leg 2
   // ME takes values 1..8, denoting m_standard MEs
   // LR takes values 0 for left-handed projection of ME, 1 for right-handed projection,
   // ubar1 P_i v2
   if (ME == 1) {
     if (LR == 0) {
       return B_0(m_s12,m_m2,m_m2,m_mu2)
	 *((Box2Abb[11]*Box2Abb[358]*(epsP*m_p1))/(Box2Abb[233]*Box2Abb[4]) 
	   + (Box2Abb[236]*Box2Abb[358]*(epsP*m_p2))/(Box2Abb[233]*Box2Abb[4]))

	 +B_0(m_s,m_m2,m_m2,m_mu2)
	 *((Box2Abb[237]*Box2Abb[358]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[239]*Box2Abb[358]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[72]))

	 +B_0(m_s2k,0.,m_m2,m_mu2)
	 *((Box2Abb[359]*Box2Abb[8]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72])
	   + (Box2Abb[255]*Box2Abb[360]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[258]*Box2Abb[361]*(epsP*m_p1))/(Box2Abb[233]*Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[273]*Box2Abb[360]*(epsP*m_p2))/(Box2Abb[233]*Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	 *((Box2Abb[0]*Box2Abb[223]*Box2Abb[355]*Box2Abb[68]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[0]*Box2Abb[223]*Box2Abb[355]*Box2Abb[69]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[223]*Box2Abb[356]*Box2Abb[5]*Box2Abb[68]*(epsP*m_p1))/pow(Box2Abb[4],2.)
	   + (Box2Abb[223]*Box2Abb[356]*Box2Abb[5]*Box2Abb[69]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	 *((-2.*Box2Abb[231]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   - (2.*Box2Abb[248]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((-2.*Box2Abb[297]*(epsP*m_p1))/(pow(Box2Abb[4],2.)*Box2Abb[72]) 
	   - (2.*Box2Abb[333]*(epsP*m_p2))/(pow(Box2Abb[4],2.)*Box2Abb[72]))

	 +D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[231]*Box2Abb[357]*Box2Abb[5]*(epsP*m_p1))/pow(Box2Abb[4],2.)
	   + (Box2Abb[231]*Box2Abb[232]*Box2Abb[357]*(epsP*m_p2))/pow(Box2Abb[4],2.));
     }
     else if (LR == 1) {
       return B_0(m_s12,m_m2,m_m2,m_mu2)
	 *((Box2Abb[11]*Box2Abb[358]*(epsP*m_p1))/(Box2Abb[233]*Box2Abb[4])
	   + (Box2Abb[236]*Box2Abb[358]*(epsP*m_p2))/(Box2Abb[233]*Box2Abb[4]))

	 +B_0(m_s,m_m2,m_m2,m_mu2)
	 *((Box2Abb[237]*Box2Abb[358]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72])
	   + (Box2Abb[239]*Box2Abb[358]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[72]))

	 +B_0(m_s2k,0.,m_m2,m_mu2)
	 *((Box2Abb[359]*Box2Abb[8]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72])
	   + (Box2Abb[255]*Box2Abb[360]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[258]*Box2Abb[361]*(epsP*m_p1))/(Box2Abb[233]*Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[273]*Box2Abb[360]*(epsP*m_p2))/(Box2Abb[233]*Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	 *((Box2Abb[0]*Box2Abb[223]*Box2Abb[355]*Box2Abb[68]*(epsP*m_p1))/pow(Box2Abb[4],2.)
	   + (Box2Abb[0]*Box2Abb[223]*Box2Abb[355]*Box2Abb[69]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[223]*Box2Abb[356]*Box2Abb[5]*Box2Abb[68]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[223]*Box2Abb[356]*Box2Abb[5]*Box2Abb[69]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	 *((-2.*Box2Abb[231]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   - (2.*Box2Abb[248]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((-2.*Box2Abb[297]*(epsP*m_p1))/(pow(Box2Abb[4],2.)*Box2Abb[72])
	   - (2.*Box2Abb[333]*(epsP*m_p2))/(pow(Box2Abb[4],2.)*Box2Abb[72]))

	 +D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[231]*Box2Abb[357]*Box2Abb[5]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[231]*Box2Abb[232]*Box2Abb[357]*(epsP*m_p2))/pow(Box2Abb[4],2.));
     }
   }
   // ubar1 \slashed{k} P_i v2
   else if (ME == 2) {
     if (LR == 0) {
       return B_0(m_s12,m_m2,m_m2,m_mu2)
	 *((Box2Abb[194]*Box2Abb[68]*(epsP*m_p1))/(Box2Abb[0]*Box2Abb[4]) 
	   + (Box2Abb[194]*Box2Abb[69]*(epsP*m_p2))/(Box2Abb[0]*Box2Abb[4]))

	 +B_0(m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[195]*Box2Abb[70]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72])
	   + (Box2Abb[195]*Box2Abb[76]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_s2k,0.,m_m2,m_mu2)
	 *((Box2Abb[196]*Box2Abb[78]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72])
	   + (Box2Abb[197]*Box2Abb[80]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_s,m_m2,m_m2,m_mu2)
	 *((Box2Abb[198]*Box2Abb[84]*(epsP*m_p1))/(Box2Abb[33]*Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[198]*Box2Abb[87]*(epsP*m_p2))/(Box2Abb[33]*Box2Abb[4]*Box2Abb[72]))
	 
	 +C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	 *((Box2Abb[199]*Box2Abb[33]*Box2Abb[95]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[101]*Box2Abb[199]*Box2Abb[33]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[199]*Box2Abb[5]*Box2Abb[95]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[101]*Box2Abb[199]*Box2Abb[5]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[110]*Box2Abb[200]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[117]*Box2Abb[200]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[142]*Box2Abb[200]*(epsP*m_p1))/(pow(Box2Abb[4],2.)*Box2Abb[72]) 
	   + (Box2Abb[174]*Box2Abb[200]*(epsP*m_p2))/(pow(Box2Abb[4],2.)*Box2Abb[72]))

	 +D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[110]*Box2Abb[201]*Box2Abb[5]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[117]*Box2Abb[201]*Box2Abb[5]*(epsP*m_p2))/pow(Box2Abb[4],2.));
     }
     else if (LR == 1) {
       return B_0(m_s12,m_m2,m_m2,m_mu2)
	 *((Box2Abb[194]*Box2Abb[68]*(epsP*m_p1))/(Box2Abb[0]*Box2Abb[4]) 
	   + (Box2Abb[194]*Box2Abb[69]*(epsP*m_p2))/(Box2Abb[0]*Box2Abb[4]))

	 +B_0(m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[195]*Box2Abb[70]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[195]*Box2Abb[76]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_s2k,0.,m_m2,m_mu2)
	 *((Box2Abb[196]*Box2Abb[78]*(epsP*m_p1))/(Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[197]*Box2Abb[80]*(epsP*m_p2))/(Box2Abb[4]*Box2Abb[5]*Box2Abb[72]))

	 +B_0(m_s,m_m2,m_m2,m_mu2)
	 *((Box2Abb[198]*Box2Abb[84]*(epsP*m_p1))/(Box2Abb[33]*Box2Abb[4]*Box2Abb[72]) 
	   + (Box2Abb[198]*Box2Abb[87]*(epsP*m_p2))/(Box2Abb[33]*Box2Abb[4]*Box2Abb[72]))

	 +C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2)
	 *((Box2Abb[199]*Box2Abb[33]*Box2Abb[95]*(epsP*m_p1))/pow(Box2Abb[4],2.)
	   + (Box2Abb[101]*Box2Abb[199]*Box2Abb[33]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[199]*Box2Abb[5]*Box2Abb[95]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[101]*Box2Abb[199]*Box2Abb[5]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2)
	 *((Box2Abb[110]*Box2Abb[200]*(epsP*m_p1))/pow(Box2Abb[4],2.)
	   + (Box2Abb[117]*Box2Abb[200]*(epsP*m_p2))/pow(Box2Abb[4],2.))

	 +C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[142]*Box2Abb[200]*(epsP*m_p1))/(pow(Box2Abb[4],2.)*Box2Abb[72])
	   + (Box2Abb[174]*Box2Abb[200]*(epsP*m_p2))/(pow(Box2Abb[4],2.)*Box2Abb[72]))

	 +D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2)
	 *((Box2Abb[110]*Box2Abb[201]*Box2Abb[5]*(epsP*m_p1))/pow(Box2Abb[4],2.) 
	   + (Box2Abb[117]*Box2Abb[201]*Box2Abb[5]*(epsP*m_p2))/pow(Box2Abb[4],2.));
     }
   }
   // ubar1 \slashed{m_epsP} P_i v2
   else if (ME == 3) {
     if (LR == 0) {
       return (Box2Abb[33]*Box2Abb[37]*Box2Abb[59]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[37]*Box2Abb[5]*Box2Abb[59]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[45]*Box2Abb[59]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[49]*Box2Abb[59]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[5]*Box2Abb[51]*Box2Abb[60]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[4];
     }
     else if (LR == 1) {
       return (Box2Abb[33]*Box2Abb[37]*Box2Abb[59]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[37]*Box2Abb[5]*Box2Abb[59]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[45]*Box2Abb[59]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[49]*Box2Abb[59]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[5]*Box2Abb[51]*Box2Abb[60]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[4];
     }
   }
   // ubar1 \slashed{m_epsP} \slashed{k} P_i v2
   else if (ME == 4) {
     if (LR == 0) {
       return (pow(Box2Abb[0],2.)*Box2Abb[1]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[4]

	 +(-Box2Abb[0]*Box2Abb[1]*Box2Abb[5]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[6]*Box2Abb[8]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[1]*Box2Abb[23]*Box2Abb[9]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[0]*Box2Abb[1]*Box2Abb[11]*Box2Abb[24]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[4];
     }
     else if (LR == 1) {
       return (pow(Box2Abb[0],2.)*Box2Abb[1]*C_0(0.,m_s,m_s12,m_m2,m_m2,m_m2,m_mu2))/Box2Abb[4]

	 +(-Box2Abb[0]*Box2Abb[1]*Box2Abb[5]*C_0(0.,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[6]*Box2Abb[8]*C_0(m_s,m_m2,m_s2k,m_m2,m_m2,0.,m_mu2))/Box2Abb[4]

	 +(Box2Abb[1]*Box2Abb[23]*Box2Abb[9]*C_0(m_m2,m_m2,m_s12,m_m2,0.,m_m2,m_mu2))/Box2Abb[4]

	 +(Box2Abb[0]*Box2Abb[1]*Box2Abb[11]*Box2Abb[24]*D_0(0.,m_s,m_m2,m_m2,m_s12,m_s2k,m_m2,m_m2,m_m2,0.,m_mu2))/Box2Abb[4];
     }
   }
   else {
     msg_Error() << "Standard matrix element not known. Value given: " << ME << "\n"
	       << "Values range from 1 to 8.";
   }
   return Zero;
 }


