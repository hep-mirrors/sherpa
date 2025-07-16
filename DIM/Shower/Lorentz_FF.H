#ifndef DIM__Shower__Lorentz_FF_H
#define DIM__Shower__Lorentz_FF_H

#include "DIM/Shower/Lorentz.H"

namespace DIM {

  class Lorentz_FF: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_FF

}// end of namepsace DIM

#endif
