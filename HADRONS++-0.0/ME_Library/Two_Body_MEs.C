#include "Two_Body_MEs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

P_2Gamma::P_2Gamma(int _nout,Flavour * _flavs) :
  HD_ME_Base(_nout,_flavs)
{ }

double P_2Gamma::operator()(const Vec4D * moms)
{
  return sqr(p_masses2[0])/4.;
}


