#ifndef CFPSHOWER_Calculators_SF_Base_H
#define CFPSHOWER_Calculators_SF_Base_H

#include "CFPSHOWER++/Shower/Splitting.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace CFPSHOWER {
  struct ffff_mode {
    enum code {
      same = 0,
      diff = 1
    };
  };
  std::ostream &operator<<(std::ostream &s,const ffff_mode::code & mode);

  struct subtract {
    enum code {
      none = 0,
      coll = 1,
      soft = 2,
      both = 3
    };
  };
  std::ostream &operator<<(std::ostream &s,const subtract::code & sub);
  
  class SF_Base {
  protected:
    static double Nc, CF, CA, TR;
    ATOOLS::Flavour        m_split;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<size_t>    m_tagsequence;
    std::string            m_name;
    int                    m_CMW;
  
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
    
    virtual inline const ATOOLS::Flavour_Vector & GetFlavs() const { return m_flavs; }
    virtual inline const std::vector<size_t> & TagSequence() const { return m_tagsequence; }
    virtual inline const std::string & Name()                const { return m_name; }
    virtual inline void SetName(const std::string & name)          { m_name = name; }
  };
  
  typedef ATOOLS::Getter_Function<SF_Base,Kernel_Info,std::less<std::string> > SF_Getter;
}

#endif
