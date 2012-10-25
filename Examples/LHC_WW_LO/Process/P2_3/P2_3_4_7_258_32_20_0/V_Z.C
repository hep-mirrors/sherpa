#include "V.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

void VP2_3_4_7_258_32_20_0::SetCouplFlav(vector<Complex>& coupl)
{
  f[0] = 24;
  f[1] = 2;
  f[2] = 11;
  f[3] = 1;

  for (int i=0;i<7;i++) c[i] = coupl[i];
  for (int i=0;i<32;i++)
    for (int j=0;j<4;j++) M[i][j] = Complex(0.,0.);

  Z[0] = Complex(0.,0.);
}

void VP2_3_4_7_258_32_20_0::Calculate()
{
  for(int i=0;i<32;i++) cl[i] = 0;

  Calculate_1();
}
