#include "P2_3/2_3_4_8_55_6/V.H"
#include "P2_3/2_3_4_8_55_5/V.H"
#include "P2_3/2_3_4_8_55_4/V.H"
#include "P2_3/2_3_4_8_55_3/V.H"
#include "P2_3/2_3_4_8_55_2/V.H"
#include "P2_3/2_3_4_8_55_1/V.H"
#include "P2_2/2_2_2_6_14_2/V.H"
#include "P2_2/2_2_2_6_14_1/V.H"
#include "P2_2/2_2_3_4_37_1/V.H"
#include "P2_2/2_2_1_1_8_1/V.H"
#include "String_Handler.H"

using namespace AMEGIC;
using namespace std;

Values* String_Handler::Set_Values(string& pID,Basic_Sfuncs* BS)
{
#ifdef V2_2_1_1_8_1_on
  if (pID==string("V2_2_1_1_8_1")) return (new V2_2_1_1_8_1(BS));
#endif
#ifdef V2_2_3_4_37_1_on
  if (pID==string("V2_2_3_4_37_1")) return (new V2_2_3_4_37_1(BS));
#endif
#ifdef V2_2_2_6_14_1_on
  if (pID==string("V2_2_2_6_14_1")) return (new V2_2_2_6_14_1(BS));
#endif
#ifdef V2_2_2_6_14_2_on
  if (pID==string("V2_2_2_6_14_2")) return (new V2_2_2_6_14_2(BS));
#endif
#ifdef V2_3_4_8_55_1_on
  if (pID==string("V2_3_4_8_55_1")) return (new V2_3_4_8_55_1(BS));
#endif
#ifdef V2_3_4_8_55_2_on
  if (pID==string("V2_3_4_8_55_2")) return (new V2_3_4_8_55_2(BS));
#endif
#ifdef V2_3_4_8_55_3_on
  if (pID==string("V2_3_4_8_55_3")) return (new V2_3_4_8_55_3(BS));
#endif
#ifdef V2_3_4_8_55_4_on
  if (pID==string("V2_3_4_8_55_4")) return (new V2_3_4_8_55_4(BS));
#endif
#ifdef V2_3_4_8_55_5_on
  if (pID==string("V2_3_4_8_55_5")) return (new V2_3_4_8_55_5(BS));
#endif
#ifdef V2_3_4_8_55_6_on
  if (pID==string("V2_3_4_8_55_6")) return (new V2_3_4_8_55_6(BS));
#endif
  return 0;
}


































































































































































