#ifndef DIM__Shower__Lorentz_FF_H
#define DIM__Shower__Lorentz_FF_H

#include "DIM/Shower/Lorentz.H"
#include "MCATNLO/Showers/Kinematics_Base.H"

namespace DIM {

  class Lorentz_FF: public Lorentz {
  private:

    MCATNLO::Kinematics_FF m_ff;

  protected:

    double Jacobian(const Splitting &s) const;

  public:

    Lorentz_FF(const Kernel_Key &k);

    int Construct(Splitting &s,const int mode) const;

    bool Compute(Splitting &s) const;

    double MEPSWeight(const Splitting &s) const;

    double GetVimax(const Splitting &s) const;

    double GetVi(const Splitting &s, const double &vitilde_max) const;

    double GetViab(const Splitting &s) const;

    double KT2(const Splitting &s, const double &vi) const;

    double JacobianResAware(const Splitting &s) const;

  };// end of class Lorentz_FF

}// end of namepsace DIM

#endif
