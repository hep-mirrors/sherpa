#ifndef DIRE__Shower__Lorentz_FF_H
#define DIRE__Shower__Lorentz_FF_H

#include "DIRE/Shower/Lorentz.H"
#include "MCATNLO/Showers/Kinematics_Base.H"

namespace DIRE {

  class Lorentz_FF: public Lorentz {
  private:

    MCATNLO::Kinematics_FF m_ff;

  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

    double GetVimax(const Splitting &s) const;
    double GetVi(const Splitting &s, const double &vitilde_max) const;
    double GetViab(const Splitting &s) const;
    double KT2(const Splitting &s, const double &vi) const;
    double JacobianResAware(const Splitting &s) const;

  };// end of class Lorentz_FF

  class Lorentz_FF_123: public Lorentz_FF {
  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF_123(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Cluster(Splitting &s,const int mode) const;

    void SetMS(ATOOLS::Mass_Selector *const ms);

    bool Allowed(const Splitting &s) const;

  };// end of class Lorentz_FF_123

}// end of namepsace DIRE

#endif
