#include "P2_3/2_3_8_4_101_0_1/P.H"
#include "Phase_Space_Generator.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace std;

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						   int chn,string& pID)
{
#ifdef P2_3_8_4_101_0_1_on
  if (pID==string("P2_3_8_4_101_0_1")) return (new P2_3_8_4_101_0_1(nin,nout,fl,chn));
#endif
  return 0;
}














































