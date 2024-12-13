#ifndef YFS_NLO_Amp_H
#define YFS_NLO_Amp_H

#include "YFS/Main/YFS_Base.H"


namespace YFS {

class NLO_Amp : YFS_Base {
public:
  NLO_Amp();
  ~NLO_Amp();

  double Calculate(const Vec4D_Vector & p);

  };

}
#endif