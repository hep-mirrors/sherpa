#ifndef AMEGIC_String_Values_H
#define AMEGIC_String_Values_H

#include "ATOOLS/Math/MyComplex.H"
#include <vector>
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"

namespace AMEGIC {
  class Values {
  public:
    virtual ~Values() {}
    virtual Complex Evaluate(int,int) {return Complex(0.,0.);}
    virtual void    Calculate() 
    { std::cerr<<"Wrong Method in Values::Calculate()"<<std::endl;}
    virtual void    SetCouplFlav(std::vector<Complex>&) 
    { std::cerr<<"Wrong Method in Values::SetCouplFlav()"<<std::endl;}
    virtual int    NumberOfCouplings() { return 0; }
  };
}
#endif









