#ifndef Coulomb_h
#define Coulomb_h

#include "Single_Amplitude.H"

namespace AMEGIC {
class Coulomb {
  double** COULM;
  Single_Amplitude* first_ampl;
  double fcoul;
  int icoulomb;

  void Build_Matrix();
public:
  Coulomb(Single_Amplitude*);
  void Calculate(int*,ATOOLS::Vec4D*);
  inline double Factor(int i,int j) {return 1.+fcoul*COULM[i][j];}
};
}
#endif
