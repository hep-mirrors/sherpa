#ifndef DIM__Shower__Lorentz_II_H
#define DIM__Shower__Lorentz_II_H

#include "DIM/Shower/Lorentz.H"

namespace DIM {

  class Lorentz_II: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_II(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

    double MEPSWeight(const Splitting &s) const;

  };// end of class Lorentz_II

}// end of namepsace DIM

#endif
