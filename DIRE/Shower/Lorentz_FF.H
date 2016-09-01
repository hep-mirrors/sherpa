#ifndef DIRE__Shower__Lorentz_FF_H
#define DIRE__Shower__Lorentz_FF_H

#include "DIRE/Shower/Lorentz.H"

namespace DIRE {

  class Lorentz_FF: public Lorentz {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

    bool Compute(Splitting &s,const int mode=1) const;

  };// end of class Lorentz_FF

  class Lorentz_FF_123: public Lorentz_FF {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF_123(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

  };// end of class Lorentz_FF_123

}// end of namepsace DIRE

#endif
