#include "P2_4/2_4_16_10_378_2/V.H"
#include "P2_4/2_4_16_10_378_1/V.H"
#include "P2_4/2_4_32_7_623_2/V.H"
#include "P2_4/2_4_32_10_683_2/V.H"
#include "P2_4/2_4_11_14_521_2/V.H"
#include "P2_4/2_4_64_7_1051_2/V.H"
#include "P2_4/2_4_11_14_521_1/V.H"
#include "P2_4/2_4_32_10_683_1/V.H"
#include "P2_4/2_4_32_7_623_1/V.H"
#include "P2_4/2_4_43_15_973_1/V.H"
#include "P2_4/2_4_64_7_1051_1/V.H"
#include "P2_3/2_3_4_8_86_2/V.H"
#include "P2_3/2_3_4_8_86_1/V.H"
#include "P2_2/2_2_2_6_24_2/V.H"
#include "P2_2/2_2_2_6_24_1/V.H"
#include "String_Handler.H"

using namespace AMEGIC;
using namespace std;

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
#ifdef V2_2_2_6_24_1_on
  if (pID==string("V2_2_2_6_24_1")) return (new V2_2_2_6_24_1(BS));
#endif
#ifdef V2_2_2_6_24_2_on
  if (pID==string("V2_2_2_6_24_2")) return (new V2_2_2_6_24_2(BS));
#endif
#ifdef V2_3_4_8_86_1_on
  if (pID==string("V2_3_4_8_86_1")) return (new V2_3_4_8_86_1(BS));
#endif
#ifdef V2_3_4_8_86_2_on
  if (pID==string("V2_3_4_8_86_2")) return (new V2_3_4_8_86_2(BS));
#endif
#ifdef V2_4_64_7_1051_1_on
  if (pID==string("V2_4_64_7_1051_1")) return (new V2_4_64_7_1051_1(BS));
#endif
#ifdef V2_4_43_15_973_1_on
  if (pID==string("V2_4_43_15_973_1")) return (new V2_4_43_15_973_1(BS));
#endif
#ifdef V2_4_32_7_623_1_on
  if (pID==string("V2_4_32_7_623_1")) return (new V2_4_32_7_623_1(BS));
#endif
#ifdef V2_4_32_10_683_1_on
  if (pID==string("V2_4_32_10_683_1")) return (new V2_4_32_10_683_1(BS));
#endif
#ifdef V2_4_11_14_521_1_on
  if (pID==string("V2_4_11_14_521_1")) return (new V2_4_11_14_521_1(BS));
#endif
#ifdef V2_4_64_7_1051_2_on
  if (pID==string("V2_4_64_7_1051_2")) return (new V2_4_64_7_1051_2(BS));
#endif
#ifdef V2_4_11_14_521_2_on
  if (pID==string("V2_4_11_14_521_2")) return (new V2_4_11_14_521_2(BS));
#endif
#ifdef V2_4_32_10_683_2_on
  if (pID==string("V2_4_32_10_683_2")) return (new V2_4_32_10_683_2(BS));
#endif
#ifdef V2_4_32_7_623_2_on
  if (pID==string("V2_4_32_7_623_2")) return (new V2_4_32_7_623_2(BS));
#endif
#ifdef V2_4_16_10_378_1_on
  if (pID==string("V2_4_16_10_378_1")) return (new V2_4_16_10_378_1(BS));
#endif
#ifdef V2_4_16_10_378_2_on
  if (pID==string("V2_4_16_10_378_2")) return (new V2_4_16_10_378_2(BS));
#endif
  return 0;
}







































































































































































