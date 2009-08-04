#include "METOOLS/Main/Four_Particle_Amplitudes.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// 0 -> VSSS
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

VSSS::VSSS(Flavour *fl,int* i,bool *out) :
  Partial_Amplitude_Base(fl,4,i,out)
{
  AssertIn(1);
  AssertSpins(2,0,0,0);
}

void VSSS::operator()(const Vec4D * moms,const bool anti)
{
  Vec4D pS1(moms[p_i[1]]);
  Vec4D pS2(moms[p_i[2]]);
  Vec4D pS3(moms[p_i[3]]);
  Vec4D pV(moms[p_i[0]]);
  Flavour flV(p_flavs[p_i[0]]);
  Polarization_Vector eps(pV,sqr(flV.HadMass()),flV.IsAnti()^anti,p_out[0]);
  int npol=IsZero(flV.HadMass())?2:3;
  for (int Vpol(0);Vpol<npol;Vpol++) {
    Insert(eps[Vpol]*cross(pS1,pS2,pS3),Vpol);
  }
}
