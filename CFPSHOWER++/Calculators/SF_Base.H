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
    ATOOLS::Flavour            m_split;
    ATOOLS::Flavour_Vector     m_flavs;
    std::vector<ATOOLS::Vec4D> m_moms;
    ATOOLS::Vec4D              m_specmom;
    std::vector<size_t>        m_tags, m_invtags;
    std::string                m_name;
    int                        m_CMW;
    subtract::code             m_subtract;

    virtual double Lambda(const double & a,const double & b,const double & c) const;
    virtual double Lambda2(const double & a,const double & b,const double & c) const;
  public:
    SF_Base(const Kernel_Info & info);
    virtual double operator()(const Splitting & split)   const=0;
    virtual double Jacobean(const Splitting & split)     const=0;
    virtual double Integral(const Splitting & split)     const=0;
    virtual double OverEstimate(const Splitting & split) const=0;
    virtual void   GeneratePoint(Splitting & split)      const=0;
    virtual bool   Construct(Splitting & split)               =0;
    
    virtual inline const ATOOLS::Flavour & GetSplit()           const { return m_split; }
    virtual inline const ATOOLS::Flavour_Vector & GetFlavs()    const { return m_flavs; }
    virtual inline const std::vector<ATOOLS::Vec4D> & GetMoms() const { return m_moms;  }
    virtual inline const ATOOLS::Vec4D & SpecMom()              const { return m_specmom; }
    virtual inline const std::vector<size_t> & Tags()           const { return m_tags; }
    virtual inline const size_t & Tags(const size_t & pos)      const { return m_tags[pos]; }
    virtual inline const std::vector<size_t> & InvTags()        const { return m_invtags; }
    virtual inline const size_t & InvTags(const size_t & pos)   const { return m_invtags[pos]; }
    virtual inline const std::string & Name()                   const { return m_name; }
    virtual inline void SetName(const std::string & name)             { m_name = name; }

    virtual inline const bool IsSubtract(subtract::code sub) const {
      return (int(sub & m_subtract)>0);
    }
  };
  
  typedef ATOOLS::Getter_Function<SF_Base,Kernel_Info,std::less<std::string> > SF_Getter;
}

#endif
