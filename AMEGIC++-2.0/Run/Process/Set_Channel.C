#include "P2_2/2_2_5_13_36_1/P.H"
#include "P2_2/2_2_4_11_31_1/P.H"
#include "P2_2/2_2_8_9_58_3/P.H"
#include "P2_2/2_2_8_9_58_2/P.H"
#include "P2_2/2_2_8_9_58_1/P.H"
#include "P2_2/2_2_6_8_48_2/P.H"
#include "P2_2/2_2_6_8_48_1/P.H"
#include "P2_2/2_2_4_9_30_3/P.H"
#include "P2_2/2_2_4_9_30_2/P.H"
#include "P2_2/2_2_3_4_25_8/P.H"
#include "P2_2/2_2_4_6_30_2/P.H"
#include "P2_2/2_2_5_11_36_2/P.H"
#include "P2_2/2_2_3_7_25_4/P.H"
#include "P2_2/2_2_1_2_7_5/P.H"
#include "P2_2/2_2_3_4_25_7/P.H"
#include "P2_2/2_2_3_4_25_6/P.H"
#include "P2_2/2_2_1_2_7_4/P.H"
#include "P2_2/2_2_4_9_31_4/P.H"
#include "P2_2/2_2_3_7_25_3/P.H"
#include "P2_2/2_2_1_2_7_3/P.H"
#include "P2_2/2_2_3_4_25_5/P.H"
#include "P2_2/2_2_1_2_7_2/P.H"
#include "P2_2/2_2_4_9_31_3/P.H"
#include "P2_2/2_2_3_4_35_1/P.H"
#include "P2_2/2_2_8_5_57_3/P.H"
#include "P2_2/2_2_4_9_29_1/P.H"
#include "P2_2/2_2_8_5_57_2/P.H"
#include "P2_2/2_2_4_6_29_1/P.H"
#include "P2_2/2_2_5_11_35_1/P.H"
#include "P2_2/2_2_8_5_57_1/P.H"
#include "P2_2/2_2_3_4_25_4/P.H"
#include "P2_2/2_2_6_4_48_4/P.H"
#include "P2_2/2_2_3_7_25_2/P.H"
#include "P2_2/2_2_3_4_25_3/P.H"
#include "P2_2/2_2_4_9_31_2/P.H"
#include "P2_2/2_2_6_4_48_3/P.H"
#include "P2_2/2_2_4_9_30_1/P.H"
#include "P2_2/2_2_4_6_30_1/P.H"
#include "P2_2/2_2_5_11_36_1/P.H"
#include "P2_2/2_2_3_4_25_2/P.H"
#include "P2_2/2_2_1_2_7_1/P.H"
#include "P2_2/2_2_3_7_25_1/P.H"
#include "P2_2/2_2_3_4_25_1/P.H"
#include "P2_2/2_2_4_9_31_1/P.H"
#include "P2_2/2_2_8_5_58_3/P.H"
#include "P2_2/2_2_8_5_58_2/P.H"
#include "P2_2/2_2_8_5_58_1/P.H"
#include "P2_2/2_2_6_4_48_2/P.H"
#include "P2_2/2_2_6_4_48_1/P.H"
#include "Phase_Space_Generator.H"

using namespace PHASIC;
using namespace std;

Single_Channel * Phase_Space_Generator::SetChannel(int nin,int nout,APHYTOOLS::Flavour* fl,
						   int chn,string& pID)
{
#ifdef P2_2_6_4_48_1_on
  if (pID==string("P2_2_6_4_48_1")) return (new P2_2_6_4_48_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_6_4_48_2_on
  if (pID==string("P2_2_6_4_48_2")) return (new P2_2_6_4_48_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_58_1_on
  if (pID==string("P2_2_8_5_58_1")) return (new P2_2_8_5_58_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_58_2_on
  if (pID==string("P2_2_8_5_58_2")) return (new P2_2_8_5_58_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_58_3_on
  if (pID==string("P2_2_8_5_58_3")) return (new P2_2_8_5_58_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_31_1_on
  if (pID==string("P2_2_4_9_31_1")) return (new P2_2_4_9_31_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_1_on
  if (pID==string("P2_2_3_4_25_1")) return (new P2_2_3_4_25_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_7_25_1_on
  if (pID==string("P2_2_3_7_25_1")) return (new P2_2_3_7_25_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_1_2_7_1_on
  if (pID==string("P2_2_1_2_7_1")) return (new P2_2_1_2_7_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_2_on
  if (pID==string("P2_2_3_4_25_2")) return (new P2_2_3_4_25_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_5_11_36_1_on
  if (pID==string("P2_2_5_11_36_1")) return (new P2_2_5_11_36_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_6_30_1_on
  if (pID==string("P2_2_4_6_30_1")) return (new P2_2_4_6_30_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_30_1_on
  if (pID==string("P2_2_4_9_30_1")) return (new P2_2_4_9_30_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_6_4_48_3_on
  if (pID==string("P2_2_6_4_48_3")) return (new P2_2_6_4_48_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_31_2_on
  if (pID==string("P2_2_4_9_31_2")) return (new P2_2_4_9_31_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_3_on
  if (pID==string("P2_2_3_4_25_3")) return (new P2_2_3_4_25_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_7_25_2_on
  if (pID==string("P2_2_3_7_25_2")) return (new P2_2_3_7_25_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_6_4_48_4_on
  if (pID==string("P2_2_6_4_48_4")) return (new P2_2_6_4_48_4(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_4_on
  if (pID==string("P2_2_3_4_25_4")) return (new P2_2_3_4_25_4(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_57_1_on
  if (pID==string("P2_2_8_5_57_1")) return (new P2_2_8_5_57_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_5_11_35_1_on
  if (pID==string("P2_2_5_11_35_1")) return (new P2_2_5_11_35_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_6_29_1_on
  if (pID==string("P2_2_4_6_29_1")) return (new P2_2_4_6_29_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_57_2_on
  if (pID==string("P2_2_8_5_57_2")) return (new P2_2_8_5_57_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_29_1_on
  if (pID==string("P2_2_4_9_29_1")) return (new P2_2_4_9_29_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_5_57_3_on
  if (pID==string("P2_2_8_5_57_3")) return (new P2_2_8_5_57_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_35_1_on
  if (pID==string("P2_2_3_4_35_1")) return (new P2_2_3_4_35_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_31_3_on
  if (pID==string("P2_2_4_9_31_3")) return (new P2_2_4_9_31_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_1_2_7_2_on
  if (pID==string("P2_2_1_2_7_2")) return (new P2_2_1_2_7_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_5_on
  if (pID==string("P2_2_3_4_25_5")) return (new P2_2_3_4_25_5(nin,nout,fl,chn));
#endif
#ifdef P2_2_1_2_7_3_on
  if (pID==string("P2_2_1_2_7_3")) return (new P2_2_1_2_7_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_7_25_3_on
  if (pID==string("P2_2_3_7_25_3")) return (new P2_2_3_7_25_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_31_4_on
  if (pID==string("P2_2_4_9_31_4")) return (new P2_2_4_9_31_4(nin,nout,fl,chn));
#endif
#ifdef P2_2_1_2_7_4_on
  if (pID==string("P2_2_1_2_7_4")) return (new P2_2_1_2_7_4(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_6_on
  if (pID==string("P2_2_3_4_25_6")) return (new P2_2_3_4_25_6(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_7_on
  if (pID==string("P2_2_3_4_25_7")) return (new P2_2_3_4_25_7(nin,nout,fl,chn));
#endif
#ifdef P2_2_1_2_7_5_on
  if (pID==string("P2_2_1_2_7_5")) return (new P2_2_1_2_7_5(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_7_25_4_on
  if (pID==string("P2_2_3_7_25_4")) return (new P2_2_3_7_25_4(nin,nout,fl,chn));
#endif
#ifdef P2_2_5_11_36_2_on
  if (pID==string("P2_2_5_11_36_2")) return (new P2_2_5_11_36_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_6_30_2_on
  if (pID==string("P2_2_4_6_30_2")) return (new P2_2_4_6_30_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_3_4_25_8_on
  if (pID==string("P2_2_3_4_25_8")) return (new P2_2_3_4_25_8(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_30_2_on
  if (pID==string("P2_2_4_9_30_2")) return (new P2_2_4_9_30_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_9_30_3_on
  if (pID==string("P2_2_4_9_30_3")) return (new P2_2_4_9_30_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_6_8_48_1_on
  if (pID==string("P2_2_6_8_48_1")) return (new P2_2_6_8_48_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_6_8_48_2_on
  if (pID==string("P2_2_6_8_48_2")) return (new P2_2_6_8_48_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_9_58_1_on
  if (pID==string("P2_2_8_9_58_1")) return (new P2_2_8_9_58_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_9_58_2_on
  if (pID==string("P2_2_8_9_58_2")) return (new P2_2_8_9_58_2(nin,nout,fl,chn));
#endif
#ifdef P2_2_8_9_58_3_on
  if (pID==string("P2_2_8_9_58_3")) return (new P2_2_8_9_58_3(nin,nout,fl,chn));
#endif
#ifdef P2_2_4_11_31_1_on
  if (pID==string("P2_2_4_11_31_1")) return (new P2_2_4_11_31_1(nin,nout,fl,chn));
#endif
#ifdef P2_2_5_13_36_1_on
  if (pID==string("P2_2_5_13_36_1")) return (new P2_2_5_13_36_1(nin,nout,fl,chn));
#endif
  return 0;
}






























































































