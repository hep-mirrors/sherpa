#include "V.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

extern "C" Values* Getter_VP2_3_4_7_258_32_20_0(Basic_Sfuncs* bs) {
  return new VP2_3_4_7_258_32_20_0(bs);
}

VP2_3_4_7_258_32_20_0::VP2_3_4_7_258_32_20_0(Basic_Sfuncs* _BS) :
     Basic_Func(0,_BS),
     Basic_Zfunc(0,_BS),
     Basic_Xfunc(0,_BS),
     Basic_Vfunc(0,_BS),
     Basic_Pfunc(0,_BS),
     Basic_Epsilonfunc(0,_BS),
     Unitarityfunc(0,_BS)
{
  f = new int[4];
  c = new Complex[7];
  Z = new Complex[258];
  M = new Complex*[32];
  for(int i=0;i<32;i++) M[i] = new Complex[4];
  cl = new int[32];
}

VP2_3_4_7_258_32_20_0::~VP2_3_4_7_258_32_20_0()
{
  if (Z)  delete[] Z;
  if (f)  delete[] f;
  if (c)  delete[] c;
  if (cl) delete[] cl;
  if (M) {
    for(int i=0;i<32;i++) delete[] M[i];
    delete[] M;
  }
}

Complex VP2_3_4_7_258_32_20_0::Evaluate(int m,int n)
{
  if (cl[n]) return M[n][m];
  switch (n) {
    case 0: Calculate_M0(); break;
    case 1: Calculate_M1(); break;
    case 4: Calculate_M4(); break;
    case 5: Calculate_M5(); break;
  }
  cl[n]=1;
  return M[n][m];
}

void VP2_3_4_7_258_32_20_0::Calculate_M0()
{
  M[0][0] = (Z[9]*Z[8]-Z[4]*Z[3])*Z[12];
  M[0][1] = -(Z[67]*Z[66]+Z[92]*Z[91]+Z[117]*Z[116])*Z[142];
  M[0][2] = (Z[3]*Z[143]+Z[147]*Z[145])*Z[149];
  M[0][3] = Z[154]*Z[3]*Z[152];
}

void VP2_3_4_7_258_32_20_0::Calculate_M1()
{
  M[1][0] = (Z[155]*Z[8]-Z[4]*Z[6])*Z[12];
  M[1][1] = -(Z[67]*Z[161]+Z[92]*Z[162]+Z[117]*Z[163])*Z[142];
  M[1][2] = (Z[165]*Z[143]+Z[167]*Z[145])*Z[149];
  M[1][3] = Z[154]*Z[6]*Z[152];
}

void VP2_3_4_7_258_32_20_0::Calculate_M4()
{
  M[4][0] = -Z[12]*Z[171]*Z[3];
  M[4][1] = -(Z[67]*Z[200]+Z[92]*Z[214]+Z[117]*Z[232])*Z[142];
  M[4][2] = Z[149]*Z[3]*Z[249];
  M[4][3] = -(Z[151]*Z[250]-Z[3]*Z[251])*Z[154];
}

void VP2_3_4_7_258_32_20_0::Calculate_M5()
{
  M[5][0] = -Z[12]*Z[171]*Z[6];
  M[5][1] = -(Z[67]*Z[253]+Z[92]*Z[254]+Z[117]*Z[255])*Z[142];
  M[5][2] = Z[149]*Z[165]*Z[249];
  M[5][3] = -(Z[169]*Z[250]-Z[6]*Z[251])*Z[154];
}

