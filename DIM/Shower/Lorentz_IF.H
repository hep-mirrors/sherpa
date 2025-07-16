#ifndef DIM__Shower__Lorentz_IF_H
#define DIM__Shower__Lorentz_IF_H

#include "DIM/Shower/Lorentz.H"

namespace DIM {

  class Lorentz_IF: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_IF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_IF

}// end of namepsace DIM

#endif
