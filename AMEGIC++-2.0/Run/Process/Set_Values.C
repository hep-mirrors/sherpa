#include "P2_2/2_2_4_10_30_2/V.H"
#include "P2_2/2_2_4_10_30_1/V.H"
#include "P2_2/2_2_3_8_25_6/V.H"
#include "P2_2/2_2_3_8_25_5/V.H"
#include "P2_2/2_2_3_8_25_4/V.H"
#include "P2_2/2_2_1_3_7_2/V.H"
#include "P2_2/2_2_3_8_25_3/V.H"
#include "P2_2/2_2_1_3_7_1/V.H"
#include "P2_2/2_2_3_8_25_2/V.H"
#include "P2_2/2_2_3_8_25_1/V.H"
#include "P2_2/2_2_5_11_36_3/V.H"
#include "P2_2/2_2_5_13_36_1/V.H"
#include "P2_2/2_2_4_11_31_1/V.H"
#include "P2_2/2_2_8_9_58_3/V.H"
#include "P2_2/2_2_8_9_58_2/V.H"
#include "P2_2/2_2_8_9_58_1/V.H"
#include "P2_2/2_2_6_8_48_2/V.H"
#include "P2_2/2_2_6_8_48_1/V.H"
#include "P2_2/2_2_4_9_30_3/V.H"
#include "P2_2/2_2_4_9_30_2/V.H"
#include "P2_2/2_2_3_4_25_8/V.H"
#include "P2_2/2_2_4_6_30_2/V.H"
#include "P2_2/2_2_5_11_36_2/V.H"
#include "P2_2/2_2_3_7_25_4/V.H"
#include "P2_2/2_2_1_2_7_5/V.H"
#include "P2_2/2_2_3_4_25_7/V.H"
#include "P2_2/2_2_3_4_25_6/V.H"
#include "P2_2/2_2_1_2_7_4/V.H"
#include "P2_2/2_2_4_9_31_4/V.H"
#include "P2_2/2_2_3_7_25_3/V.H"
#include "P2_2/2_2_1_2_7_3/V.H"
#include "P2_2/2_2_3_4_25_5/V.H"
#include "P2_2/2_2_1_2_7_2/V.H"
#include "P2_2/2_2_4_9_31_3/V.H"
#include "P2_2/2_2_3_4_35_1/V.H"
#include "P2_2/2_2_8_5_57_3/V.H"
#include "P2_2/2_2_4_9_29_1/V.H"
#include "P2_2/2_2_8_5_57_2/V.H"
#include "P2_2/2_2_4_6_29_1/V.H"
#include "P2_2/2_2_5_11_35_1/V.H"
#include "P2_2/2_2_8_5_57_1/V.H"
#include "P2_2/2_2_3_4_25_4/V.H"
#include "P2_2/2_2_6_4_48_4/V.H"
#include "P2_2/2_2_3_7_25_2/V.H"
#include "P2_2/2_2_3_4_25_3/V.H"
#include "P2_2/2_2_4_9_31_2/V.H"
#include "P2_2/2_2_6_4_48_3/V.H"
#include "P2_2/2_2_4_9_30_1/V.H"
#include "P2_2/2_2_4_6_30_1/V.H"
#include "P2_2/2_2_5_11_36_1/V.H"
#include "P2_2/2_2_3_4_25_2/V.H"
#include "P2_2/2_2_1_2_7_1/V.H"
#include "P2_2/2_2_3_7_25_1/V.H"
#include "P2_2/2_2_3_4_25_1/V.H"
#include "P2_2/2_2_4_9_31_1/V.H"
#include "P2_2/2_2_8_5_58_3/V.H"
#include "P2_2/2_2_8_5_58_2/V.H"
#include "P2_2/2_2_8_5_58_1/V.H"
#include "P2_2/2_2_6_4_48_2/V.H"
#include "P2_2/2_2_6_4_48_1/V.H"
#include "String_Handler.H"

using namespace AMEGIC;
using namespace std;

Values* String_Handler::Set_Values(string& pID,Basic_Sfuncs* BS)
{
#ifdef V2_2_6_4_48_1_on
  if (pID==string("V2_2_6_4_48_1")) return (new V2_2_6_4_48_1(BS));
#endif
#ifdef V2_2_6_4_48_2_on
  if (pID==string("V2_2_6_4_48_2")) return (new V2_2_6_4_48_2(BS));
#endif
#ifdef V2_2_8_5_58_1_on
  if (pID==string("V2_2_8_5_58_1")) return (new V2_2_8_5_58_1(BS));
#endif
#ifdef V2_2_8_5_58_2_on
  if (pID==string("V2_2_8_5_58_2")) return (new V2_2_8_5_58_2(BS));
#endif
#ifdef V2_2_8_5_58_3_on
  if (pID==string("V2_2_8_5_58_3")) return (new V2_2_8_5_58_3(BS));
#endif
#ifdef V2_2_4_9_31_1_on
  if (pID==string("V2_2_4_9_31_1")) return (new V2_2_4_9_31_1(BS));
#endif
#ifdef V2_2_3_4_25_1_on
  if (pID==string("V2_2_3_4_25_1")) return (new V2_2_3_4_25_1(BS));
#endif
#ifdef V2_2_3_7_25_1_on
  if (pID==string("V2_2_3_7_25_1")) return (new V2_2_3_7_25_1(BS));
#endif
#ifdef V2_2_1_2_7_1_on
  if (pID==string("V2_2_1_2_7_1")) return (new V2_2_1_2_7_1(BS));
#endif
#ifdef V2_2_3_4_25_2_on
  if (pID==string("V2_2_3_4_25_2")) return (new V2_2_3_4_25_2(BS));
#endif
#ifdef V2_2_5_11_36_1_on
  if (pID==string("V2_2_5_11_36_1")) return (new V2_2_5_11_36_1(BS));
#endif
#ifdef V2_2_4_6_30_1_on
  if (pID==string("V2_2_4_6_30_1")) return (new V2_2_4_6_30_1(BS));
#endif
#ifdef V2_2_4_9_30_1_on
  if (pID==string("V2_2_4_9_30_1")) return (new V2_2_4_9_30_1(BS));
#endif
#ifdef V2_2_6_4_48_3_on
  if (pID==string("V2_2_6_4_48_3")) return (new V2_2_6_4_48_3(BS));
#endif
#ifdef V2_2_4_9_31_2_on
  if (pID==string("V2_2_4_9_31_2")) return (new V2_2_4_9_31_2(BS));
#endif
#ifdef V2_2_3_4_25_3_on
  if (pID==string("V2_2_3_4_25_3")) return (new V2_2_3_4_25_3(BS));
#endif
#ifdef V2_2_3_7_25_2_on
  if (pID==string("V2_2_3_7_25_2")) return (new V2_2_3_7_25_2(BS));
#endif
#ifdef V2_2_6_4_48_4_on
  if (pID==string("V2_2_6_4_48_4")) return (new V2_2_6_4_48_4(BS));
#endif
#ifdef V2_2_3_4_25_4_on
  if (pID==string("V2_2_3_4_25_4")) return (new V2_2_3_4_25_4(BS));
#endif
#ifdef V2_2_8_5_57_1_on
  if (pID==string("V2_2_8_5_57_1")) return (new V2_2_8_5_57_1(BS));
#endif
#ifdef V2_2_5_11_35_1_on
  if (pID==string("V2_2_5_11_35_1")) return (new V2_2_5_11_35_1(BS));
#endif
#ifdef V2_2_4_6_29_1_on
  if (pID==string("V2_2_4_6_29_1")) return (new V2_2_4_6_29_1(BS));
#endif
#ifdef V2_2_8_5_57_2_on
  if (pID==string("V2_2_8_5_57_2")) return (new V2_2_8_5_57_2(BS));
#endif
#ifdef V2_2_4_9_29_1_on
  if (pID==string("V2_2_4_9_29_1")) return (new V2_2_4_9_29_1(BS));
#endif
#ifdef V2_2_8_5_57_3_on
  if (pID==string("V2_2_8_5_57_3")) return (new V2_2_8_5_57_3(BS));
#endif
#ifdef V2_2_3_4_35_1_on
  if (pID==string("V2_2_3_4_35_1")) return (new V2_2_3_4_35_1(BS));
#endif
#ifdef V2_2_4_9_31_3_on
  if (pID==string("V2_2_4_9_31_3")) return (new V2_2_4_9_31_3(BS));
#endif
#ifdef V2_2_1_2_7_2_on
  if (pID==string("V2_2_1_2_7_2")) return (new V2_2_1_2_7_2(BS));
#endif
#ifdef V2_2_3_4_25_5_on
  if (pID==string("V2_2_3_4_25_5")) return (new V2_2_3_4_25_5(BS));
#endif
#ifdef V2_2_1_2_7_3_on
  if (pID==string("V2_2_1_2_7_3")) return (new V2_2_1_2_7_3(BS));
#endif
#ifdef V2_2_3_7_25_3_on
  if (pID==string("V2_2_3_7_25_3")) return (new V2_2_3_7_25_3(BS));
#endif
#ifdef V2_2_4_9_31_4_on
  if (pID==string("V2_2_4_9_31_4")) return (new V2_2_4_9_31_4(BS));
#endif
#ifdef V2_2_1_2_7_4_on
  if (pID==string("V2_2_1_2_7_4")) return (new V2_2_1_2_7_4(BS));
#endif
#ifdef V2_2_3_4_25_6_on
  if (pID==string("V2_2_3_4_25_6")) return (new V2_2_3_4_25_6(BS));
#endif
#ifdef V2_2_3_4_25_7_on
  if (pID==string("V2_2_3_4_25_7")) return (new V2_2_3_4_25_7(BS));
#endif
#ifdef V2_2_1_2_7_5_on
  if (pID==string("V2_2_1_2_7_5")) return (new V2_2_1_2_7_5(BS));
#endif
#ifdef V2_2_3_7_25_4_on
  if (pID==string("V2_2_3_7_25_4")) return (new V2_2_3_7_25_4(BS));
#endif
#ifdef V2_2_5_11_36_2_on
  if (pID==string("V2_2_5_11_36_2")) return (new V2_2_5_11_36_2(BS));
#endif
#ifdef V2_2_4_6_30_2_on
  if (pID==string("V2_2_4_6_30_2")) return (new V2_2_4_6_30_2(BS));
#endif
#ifdef V2_2_3_4_25_8_on
  if (pID==string("V2_2_3_4_25_8")) return (new V2_2_3_4_25_8(BS));
#endif
#ifdef V2_2_4_9_30_2_on
  if (pID==string("V2_2_4_9_30_2")) return (new V2_2_4_9_30_2(BS));
#endif
#ifdef V2_2_4_9_30_3_on
  if (pID==string("V2_2_4_9_30_3")) return (new V2_2_4_9_30_3(BS));
#endif
#ifdef V2_2_6_8_48_1_on
  if (pID==string("V2_2_6_8_48_1")) return (new V2_2_6_8_48_1(BS));
#endif
#ifdef V2_2_6_8_48_2_on
  if (pID==string("V2_2_6_8_48_2")) return (new V2_2_6_8_48_2(BS));
#endif
#ifdef V2_2_8_9_58_1_on
  if (pID==string("V2_2_8_9_58_1")) return (new V2_2_8_9_58_1(BS));
#endif
#ifdef V2_2_8_9_58_2_on
  if (pID==string("V2_2_8_9_58_2")) return (new V2_2_8_9_58_2(BS));
#endif
#ifdef V2_2_8_9_58_3_on
  if (pID==string("V2_2_8_9_58_3")) return (new V2_2_8_9_58_3(BS));
#endif
#ifdef V2_2_4_11_31_1_on
  if (pID==string("V2_2_4_11_31_1")) return (new V2_2_4_11_31_1(BS));
#endif
#ifdef V2_2_5_13_36_1_on
  if (pID==string("V2_2_5_13_36_1")) return (new V2_2_5_13_36_1(BS));
#endif
#ifdef V2_2_5_11_36_3_on
  if (pID==string("V2_2_5_11_36_3")) return (new V2_2_5_11_36_3(BS));
#endif
#ifdef V2_2_3_8_25_1_on
  if (pID==string("V2_2_3_8_25_1")) return (new V2_2_3_8_25_1(BS));
#endif
#ifdef V2_2_3_8_25_2_on
  if (pID==string("V2_2_3_8_25_2")) return (new V2_2_3_8_25_2(BS));
#endif
#ifdef V2_2_1_3_7_1_on
  if (pID==string("V2_2_1_3_7_1")) return (new V2_2_1_3_7_1(BS));
#endif
#ifdef V2_2_3_8_25_3_on
  if (pID==string("V2_2_3_8_25_3")) return (new V2_2_3_8_25_3(BS));
#endif
#ifdef V2_2_1_3_7_2_on
  if (pID==string("V2_2_1_3_7_2")) return (new V2_2_1_3_7_2(BS));
#endif
#ifdef V2_2_3_8_25_4_on
  if (pID==string("V2_2_3_8_25_4")) return (new V2_2_3_8_25_4(BS));
#endif
#ifdef V2_2_3_8_25_5_on
  if (pID==string("V2_2_3_8_25_5")) return (new V2_2_3_8_25_5(BS));
#endif
#ifdef V2_2_3_8_25_6_on
  if (pID==string("V2_2_3_8_25_6")) return (new V2_2_3_8_25_6(BS));
#endif
#ifdef V2_2_4_10_30_1_on
  if (pID==string("V2_2_4_10_30_1")) return (new V2_2_4_10_30_1(BS));
#endif
#ifdef V2_2_4_10_30_2_on
  if (pID==string("V2_2_4_10_30_2")) return (new V2_2_4_10_30_2(BS));
#endif
  return 0;
}




















































































































































































































