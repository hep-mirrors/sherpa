#include "Four_Body_MEs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

P_2PFF::P_2PFF(int _nin,int _nout,Flavour * _flavs,string _met) :
  Hadron_Decay_ME(_nin,_nout,_flavs,_met)
{ }

double P_2PFF::operator()(const Vec4D * moms)
{

}


