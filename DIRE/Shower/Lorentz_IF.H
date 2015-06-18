#ifndef DIRE__Shower__Lorentz_IF_H
#define DIRE__Shower__Lorentz_IF_H

#include "DIRE/Shower/Lorentz.H"

namespace DIRE {

  class Lorentz_IF: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

    double PDFEstimate(const Splitting &s) const;

  public:

    Lorentz_IF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

  };// end of class Lorentz_IF

}// end of namepsace DIRE

#endif
