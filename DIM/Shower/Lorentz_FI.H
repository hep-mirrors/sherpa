#ifndef DIM__Shower__Lorentz_FI_H
#define DIM__Shower__Lorentz_FI_H

#include "DIM/Shower/Lorentz.H"

namespace DIM {

  class Lorentz_FI: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_FI(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_FI

}// end of namepsace DIM

#endif
