#include "P2_3/2_3_4_8_56_1_2/V.H"
#include "P2_3/2_3_4_8_56_1_1/V.H"
#include "P2_2/2_2_1_3_5_1_1/V.H"
#include "String_Handler.H"

using namespace AMEGIC;
using namespace std;

Values* String_Handler::Set_Values(std::string& pID,Basic_Sfuncs* BS)
{
#ifdef V2_2_1_3_5_1_1_on
  if (pID==string("V2_2_1_3_5_1_1")) return (new V2_2_1_3_5_1_1(BS));
#endif
#ifdef V2_3_4_8_56_1_1_on
  if (pID==string("V2_3_4_8_56_1_1")) return (new V2_3_4_8_56_1_1(BS));
#endif
#ifdef V2_3_4_8_56_1_2_on
  if (pID==string("V2_3_4_8_56_1_2")) return (new V2_3_4_8_56_1_2(BS));
#endif
  return 0;
}



























































































































































