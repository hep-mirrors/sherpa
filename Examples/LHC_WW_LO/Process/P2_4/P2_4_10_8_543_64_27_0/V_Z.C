#include "V.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

void VP2_4_10_8_543_64_27_0::SetCouplFlav(vector<Complex>& coupl)
{
  f[0] = 24;
  f[1] = 2;
  f[2] = 1;
  f[3] = 11;

  for (int i=0;i<8;i++) c[i] = coupl[i];
  for (int i=0;i<64;i++)
    for (int j=0;j<10;j++) M[i][j] = Complex(0.,0.);

  Z[0] = Complex(0.,0.);
}

void VP2_4_10_8_543_64_27_0::Calculate()
{
  for(int i=0;i<64;i++) cl[i] = 0;

  Calculate_1();
  Calculate_2();
}
