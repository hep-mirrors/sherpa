#include "P2_2/2_2_2_6_42_2/P.H"
#include "P2_2/2_2_2_6_42_1/P.H"
#include "Phase_Space_Generator.H"

using namespace PHASIC;
using namespace std;

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,APHYTOOLS::Flavour* fl,
						   int chn,string& pID)
{
#ifdef P2_2_2_6_42_1_on
  if (pID==string("P2_2_2_6_42_1")) return (new P2_2_2_6_42_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_2_6_42_2_on
  if (pID==string("P2_2_2_6_42_2")) return (new P2_2_2_6_42_2(nin,nout,fl,chn));
#endif
  return 0;
}















































