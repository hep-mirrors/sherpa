#ifndef Function_Base_H
#define Function_Base_H

namespace AMATOOLS {
  class Function_Base {
    // Interval for the values ?
  public:
    //  virtual double Get(double)=0; -> GetValue;
    virtual void SetParameters(double *) { };
    // obsolete, use operator() instead!!!
    virtual double GetValue(double xx) {return (*this)(xx);}   
    virtual double operator() (double)=0;  // {return GetValue(xx);};
  };
}

#endif
