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

    bool Compute(Splitting &s,const int mode=1) const;

  };// end of class Lorentz_II

  class Lorentz_II_123: public Lorentz_II {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_II_123(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

  };// end of class Lorentz_II_123

}// end of namepsace DIRE

#endif
