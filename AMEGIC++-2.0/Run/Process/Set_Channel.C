#include "P2_3/2_3_4_8_55_6/P.H"
#include "P2_3/2_3_4_8_55_5/P.H"
#include "P2_3/2_3_4_8_55_4/P.H"
#include "P2_3/2_3_4_8_55_3/P.H"
#include "P2_3/2_3_4_8_55_2/P.H"
#include "P2_3/2_3_4_8_55_1/P.H"
#include "P2_2/2_2_2_6_14_2/P.H"
#include "P2_2/2_2_2_6_14_1/P.H"
#include "P2_2/2_2_3_4_37_1/P.H"
#include "P2_2/2_2_1_1_8_1/P.H"
#include "Phase_Space_Generator.H"

using namespace AMEGIC;
using namespace PHASIC;
using namespace std;

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,APHYTOOLS::Flavour* fl,
						   int chn,string& pID)
{
#ifdef P2_2_1_1_8_1_on
  if (pID==string("P2_2_1_1_8_1")) return (new P2_2_1_1_8_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_37_1_on
  if (pID==string("P2_2_3_4_37_1")) return (new P2_2_3_4_37_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_2_6_14_1_on
  if (pID==string("P2_2_2_6_14_1")) return (new P2_2_2_6_14_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_2_6_14_2_on
  if (pID==string("P2_2_2_6_14_2")) return (new P2_2_2_6_14_2(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_1_on
  if (pID==string("P2_3_4_8_55_1")) return (new P2_3_4_8_55_1(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_2_on
  if (pID==string("P2_3_4_8_55_2")) return (new P2_3_4_8_55_2(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_3_on
  if (pID==string("P2_3_4_8_55_3")) return (new P2_3_4_8_55_3(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_4_on
  if (pID==string("P2_3_4_8_55_4")) return (new P2_3_4_8_55_4(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_5_on
  if (pID==string("P2_3_4_8_55_5")) return (new P2_3_4_8_55_5(nin,nout,fl,chn));
#endif
#ifdef P2_3_4_8_55_6_on
  if (pID==string("P2_3_4_8_55_6")) return (new P2_3_4_8_55_6(nin,nout,fl,chn));
#endif
  return 0;
}























































