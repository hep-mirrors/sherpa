#include "P2_2/2_2_2_6_14_1_4/V.H"
#include "P2_2/2_2_2_6_14_1_3/V.H"
#include "P2_2/2_2_2_6_14_1_2/V.H"
#include "P2_2/2_2_2_6_14_1_1/V.H"
#include "String_Handler.H"

using namespace AMEGIC;
using namespace std;

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
#ifdef V2_2_2_6_14_1_1_on
  if (pID==string("V2_2_2_6_14_1_1")) return (new V2_2_2_6_14_1_1(BS));
#endif
#ifdef V2_2_2_6_14_1_2_on
  if (pID==string("V2_2_2_6_14_1_2")) return (new V2_2_2_6_14_1_2(BS));
#endif
#ifdef V2_2_2_6_14_1_3_on
  if (pID==string("V2_2_2_6_14_1_3")) return (new V2_2_2_6_14_1_3(BS));
#endif
#ifdef V2_2_2_6_14_1_4_on
  if (pID==string("V2_2_2_6_14_1_4")) return (new V2_2_2_6_14_1_4(BS));
#endif
  return 0;
}




























































































































































