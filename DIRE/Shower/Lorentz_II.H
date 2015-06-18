#ifndef DIRE__Shower__Lorentz_II_H
#define DIRE__Shower__Lorentz_II_H

#include "DIRE/Shower/Lorentz.H"

namespace DIRE {

  class Lorentz_II: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_II(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

  };// end of class Lorentz_II

}// end of namepsace DIRE

#endif
