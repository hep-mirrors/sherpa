// modul for parameters

#ifndef Couplings_LED_h
#define Couplings_LED_h

#include "MyComplex.H"

namespace AMEGIC {

  class Couplings_LED {
    double gn,ms,rad;
    int ed,kk_sum_mode;
  public:

    void   Init();
    double Gn();
    double Ms();
    double Msq();
    int Ned();
    double Kappa();
    double Omega();
    double R();
    inline bool DoSum(){return kk_sum_mode!=0;}
    double Ifunc(double);
    double IEfunc(double);
    Complex KKProp(double);
 };
}
#endif








