// modul for parameters

#ifndef Couplings_EW_h
#define Couplings_EW_h

//#include "Flavour.H"
#include "MyComplex.H"

namespace AMEGIC {

  class Couplings_EW {
    double  sinthetaW_MZ,vev;
    Complex** _CKM;
  public:
    Couplings_EW();
    ~Couplings_EW();

    void   Init();
    double VEV();
    double SinThetaW();
    double CosThetaW();
    Complex CKM(short int,short int);
 };
}
#endif








