#include "P2_4/2_4_16_10_378_2/P.H"
#include "P2_4/2_4_16_10_378_1/P.H"
#include "P2_4/2_4_32_7_623_2/P.H"
#include "P2_4/2_4_32_10_683_2/P.H"
#include "P2_4/2_4_11_14_521_2/P.H"
#include "P2_4/2_4_64_7_1051_2/P.H"
#include "P2_4/2_4_11_14_521_1/P.H"
#include "P2_4/2_4_32_10_683_1/P.H"
#include "P2_4/2_4_32_7_623_1/P.H"
#include "P2_4/2_4_43_15_973_1/P.H"
#include "P2_4/2_4_64_7_1051_1/P.H"
#include "P2_3/2_3_4_8_86_2/P.H"
#include "P2_3/2_3_4_8_86_1/P.H"
#include "P2_2/2_2_2_6_24_2/P.H"
#include "P2_2/2_2_2_6_24_1/P.H"
#include "Phase_Space_Generator.H"

using namespace PHASIC;
using namespace std;

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,APHYTOOLS::Flavour* fl,
						   int chn,string& pID)
{
#ifdef P2_2_2_6_24_1_on
  if (pID==string("P2_2_2_6_24_1")) return (new P2_2_2_6_24_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_2_6_24_2_on
  if (pID==string("P2_2_2_6_24_2")) return (new P2_2_2_6_24_2(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_86_1_on
  if (pID==string("P2_3_4_8_86_1")) return (new P2_3_4_8_86_1(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_86_2_on
  if (pID==string("P2_3_4_8_86_2")) return (new P2_3_4_8_86_2(nin,nout,fl,chn));
#endif
#ifdef P2_4_64_7_1051_1_on
  if (pID==string("P2_4_64_7_1051_1")) return (new P2_4_64_7_1051_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_43_15_973_1_on
  if (pID==string("P2_4_43_15_973_1")) return (new P2_4_43_15_973_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_32_7_623_1_on
  if (pID==string("P2_4_32_7_623_1")) return (new P2_4_32_7_623_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_32_10_683_1_on
  if (pID==string("P2_4_32_10_683_1")) return (new P2_4_32_10_683_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_11_14_521_1_on
  if (pID==string("P2_4_11_14_521_1")) return (new P2_4_11_14_521_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_64_7_1051_2_on
  if (pID==string("P2_4_64_7_1051_2")) return (new P2_4_64_7_1051_2(nin,nout,fl,chn));
#endif
#ifdef P2_4_11_14_521_2_on
  if (pID==string("P2_4_11_14_521_2")) return (new P2_4_11_14_521_2(nin,nout,fl,chn));
#endif
#ifdef P2_4_32_10_683_2_on
  if (pID==string("P2_4_32_10_683_2")) return (new P2_4_32_10_683_2(nin,nout,fl,chn));
#endif
#ifdef P2_4_32_7_623_2_on
  if (pID==string("P2_4_32_7_623_2")) return (new P2_4_32_7_623_2(nin,nout,fl,chn));
#endif
#ifdef P2_4_16_10_378_1_on
  if (pID==string("P2_4_16_10_378_1")) return (new P2_4_16_10_378_1(nin,nout,fl,chn));
#endif
#ifdef P2_4_16_10_378_2_on
  if (pID==string("P2_4_16_10_378_2")) return (new P2_4_16_10_378_2(nin,nout,fl,chn));
#endif
  return 0;
}




























































