//Header file for process VP2_3_4_7_258_32_20_0

#ifndef VP2_3_4_7_258_32_20_0_on
#define VP2_3_4_7_258_32_20_0_on
#include "AMEGIC++/String/Values.H"

extern "C" AMEGIC::Values* Getter_VP2_3_4_7_258_32_20_0(AMEGIC::Basic_Sfuncs* bs);

namespace AMEGIC {

class VP2_3_4_7_258_32_20_0 : public Values,
  public Basic_Zfunc,
  public Basic_Xfunc,
  public Basic_Vfunc,
  public Basic_Pfunc,
  public Basic_Epsilonfunc,
  public Unitarityfunc {
  Complex*  Z;
  int*      f;
  Complex*  c;
  Complex** M;
  int*      cl;
public:
  VP2_3_4_7_258_32_20_0(Basic_Sfuncs* _BS);
  ~VP2_3_4_7_258_32_20_0();
  void SetCouplFlav(std::vector<Complex>&);
  int NumberOfCouplings() { return 7; }
  Complex Evaluate(int,int);
  void    Calculate();
  void Calculate_1();
  void Calculate_M0();
  void Calculate_M1();
  void Calculate_M4();
  void Calculate_M5();
};
}

#endif
