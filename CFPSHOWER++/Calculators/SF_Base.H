#ifndef CFPSHOWER_Calculators_SF_Base_H
#define CFPSHOWER_Calculators_SF_Base_H

#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace CFPSHOWER {
  class SF_Base {
  protected:
    ATOOLS::Flavour_Vector m_flavs;
    double m_CF, m_CA, m_TR, m_zeta3;
    int    m_orderA, m_orderB;
    bool   m_swap;

    std::string m_name;
    virtual double S1(const double & q2) const;
    virtual double S2(const double & q2) const;
    virtual double Lambda(const double & a,const double & b,const double & c) const;
    virtual double Lambda2(const double & a,const double & b,const double & c) const;
  public:
    SF_Base(const Kernel_Info & info);
    virtual double operator()(const Splitting & split)   const=0;
    virtual double Jacobean(const Splitting & split)     const=0;
    virtual double Integral(const Splitting & split)     const=0;
    virtual double OverEstimate(const Splitting & split) const=0;
    virtual void   GeneratePoint(Splitting & split)      const=0;
    virtual bool   InitKinematics(Splitting & split)     const=0;
    virtual int    Construct(Splitting & split)          const=0;

    virtual inline const std::string & Name()                const { return m_name; }
    virtual inline void SetName(const std::string & name)          { m_name = name; }
    virtual inline const ATOOLS::Flavour_Vector & GetFlavs() const { return m_flavs; }
    virtual inline const bool IsSwapped()                    const { return m_swap; }
  };

  typedef ATOOLS::Getter_Function<SF_Base,Kernel_Info,std::less<std::string> > SF_Getter;
}

#endif
