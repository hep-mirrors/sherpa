#ifndef Values_h
#define Values_h

#include "MyComplex.H"
#include <vector>
#include "Basic_Func.H"

namespace AMEGIC {
  class Values {
  public:
    virtual ~Values() {}
    virtual Complex Evaluate(int&,int&) {return Complex(0.,0.);}
    virtual void    Calculate() 
    { std::cerr<<"Wrong Method in Values::Calculate()"<<std::endl;}
    virtual void    SetCouplFlav(std::vector<Complex>&) 
    { std::cerr<<"Wrong Method in Values::SetCouplFlav()"<<std::endl;}
  };
}
#endif









