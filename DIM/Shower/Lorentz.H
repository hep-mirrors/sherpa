#ifndef DIM__Shower__Lorentz_H
#define DIM__Shower__Lorentz_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "DIM/Tools/Kernel_Key.H"
#include "DIM/Tools/Splitting.H"

namespace DIM {

  class Shower;

  class Lorentz {
  protected:

    Kernel *p_sk;
    int     m_type;

    ATOOLS::Flavour m_fl[3];

    ATOOLS::Mass_Selector *p_ms;

    inline double Lam(const double &a,
		      const double &b,
		      const double &c) const
    {  return (a-b-c)*(a-b-c)-4.0*b*c; }

    void SetParams(Splitting &s,const PHASIC::Kin_Args &ff) const;

    int Update(Splitting &s,const PHASIC::Kin_Args &ff,
	       const int mode) const;

  public:

    Lorentz(const Kernel_Key &k,const int type);

    ~Lorentz();

    virtual double Value(const Splitting &s) const = 0;
    virtual double Jacobian(const Splitting &s) const = 0;

    virtual double Estimate(const Splitting &s) const = 0;
    virtual double Integral(const Splitting &s) const = 0;

    virtual bool GeneratePoint(Splitting &s) const = 0;
    virtual int Construct(Splitting &s,const int mode) const = 0;

    virtual bool Compute(Splitting &s) const = 0;

    virtual double MEPSWeight(const Splitting &s) const = 0;

    virtual double AsymmetryFactor(const Splitting &s) const;

    bool Allowed(const Splitting &s) const;

    bool SetLimits(Splitting &s) const;

    inline const ATOOLS::Flavour &Flav(const int i) const { return m_fl[i]; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

  };// end of class Lorentz

  typedef ATOOLS::Getter_Function
  <Lorentz,Kernel_Key,std::less<std::string> > Lorentz_Getter;

}// end of namepsace DIM

#endif
