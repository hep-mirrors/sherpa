#ifndef ALARIC__Shower__Lorentz_H
#define ALARIC__Shower__Lorentz_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "PHASIC++/Channels/Antenna_Kinematics.H"
#include "ALARIC++/Tools/Kernel_Key.H"
#include "ALARIC++/Tools/Splitting.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace ALARIC {

  class Shower;

  class Lorentz {
  protected:

    Kernel *p_sk;
    int     m_type;

    ATOOLS::Flavour_Vector m_fl;

    const ATOOLS::Mass_Selector *p_ms;

    inline double Lam(const double &a,
		      const double &b,
		      const double &c) const
    {  return (a-b-c)*(a-b-c)-4.0*b*c; }

    void SetParams(Splitting &s) const;

    int Update(Splitting &s,const int mode) const;

  public:

    Lorentz(const Kernel_Key &k,const int type);

    virtual ~Lorentz();

    virtual double Value(const Splitting &s) const = 0;
    virtual double Jacobian(const Splitting &s) const = 0;

    virtual double Estimate(const Splitting &s) const = 0;
    virtual double Integral(const Splitting &s) const = 0;

    virtual bool GeneratePoint(Splitting &s) const = 0;
    virtual int Construct(Splitting &s,const int mode) const = 0;
    virtual bool Cluster(Splitting &s,const int mode,
			 const int i, const int j,const int k,
			 const ATOOLS::Cluster_Amplitude *a) const = 0;

    virtual void SetMS(const ATOOLS::Mass_Selector *ms);

    virtual bool Allowed(const Splitting &s) const;

    bool SetLimits(Splitting &s) const;

    inline const ATOOLS::Flavour &Flav(const int i) const { return m_fl[i]; }
    inline const ATOOLS::Flavour_Vector &Flavs() const { return m_fl; }

  };// end of class Lorentz

  typedef ATOOLS::Getter_Function
  <Lorentz,Kernel_Key,std::less<std::string> > Lorentz_Getter;

}// end of namepsace ALARIC

#endif
