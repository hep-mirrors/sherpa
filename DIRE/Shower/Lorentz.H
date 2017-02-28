#ifndef DIRE__Shower__Lorentz_H
#define DIRE__Shower__Lorentz_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "DIRE/Tools/Kernel_Key.H"
#include "DIRE/Tools/Splitting.H"

namespace DIRE {

  class Shower;

  class Lorentz {
  protected:

    Kernel *p_sk;
    int     m_type;

    ATOOLS::Flavour_Vector m_fl;

    ATOOLS::Mass_Selector *p_ms;

    inline double Lam(const double &a,
		      const double &b,
		      const double &c) const
    {  return (a-b-c)*(a-b-c)-4.0*b*c; }

    void SetParams(Splitting &s,const PHASIC::Kin_Args &ff) const;

    int Update(Splitting &s,const PHASIC::Kin_Args &ff,const int mode,
	       const ATOOLS::Vec4D &pl=ATOOLS::Vec4D()) const;

  public:

    Lorentz(const Kernel_Key &k,const int type);

    virtual ~Lorentz();

    virtual double Value(const Splitting &s) const = 0;
    virtual double Jacobian(const Splitting &s) const = 0;

    virtual double Estimate(const Splitting &s) const = 0;
    virtual double Integral(const Splitting &s) const = 0;

    virtual bool GeneratePoint(Splitting &s) const = 0;
    virtual int Construct(Splitting &s,const int mode) const = 0;

    virtual bool Compute(Splitting &s,const int mode=1) const = 0;
    virtual bool Cluster(Splitting &s,const int mode) const = 0;

    bool Allowed(const Splitting &s) const;

    bool SetLimits(Splitting &s) const;

    inline const ATOOLS::Flavour &Flav(const int i) const { return m_fl[i]; }
    inline const ATOOLS::Flavour_Vector &Flavs() const { return m_fl; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

  };// end of class Lorentz

  typedef ATOOLS::Getter_Function
  <Lorentz,Kernel_Key,std::less<std::string> > Lorentz_Getter;

}// end of namepsace DIRE

#endif
