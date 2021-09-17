#ifndef CFPSHOWER_SplittingFunctions_SF_Base_H
#define CFPSHOWER_SplittingFunctions_SF_Base_H

#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace CFPSHOWER {
  static double s_CF = 4./3.;
  static double s_CA = 3.;
  static double s_Nc = 3.;
  static double s_TR = 1./2.;
  
  struct subtract {
    enum code {
      undefined = -1,
      none      = 0,
      coll      = 1,
      soft      = 2,
      both      = 3
    };
  };
  inline subtract::code operator|(const subtract::code& c1,const subtract::code& c2)
  { return (subtract::code)((int)c1|(int)c2); }
  inline subtract::code operator&(const subtract::code& c1,const subtract::code& c2)
  { return (subtract::code)((int)c1&(int)c2); }
  std::ostream &operator<<(std::ostream &s,const subtract::code & sub);

  
  class SF_Base {
  protected:
    ATOOLS::Flavour        m_split;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<size_t>    m_tags;
    size_t                 m_nout;
    std::string            m_name;
    kernel_type::code      m_type;
    double                 m_weight;
    int                    m_CMW, m_softcorr, m_endpoint;
    subtract::code         m_subtract;
  public:
    SF_Base(const Kernel_Info & info);
    ~SF_Base() {}
    virtual double operator()(const Splitting & split) = 0;
    virtual double Integral(const Splitting & split)      const = 0;
    virtual double OverEstimate(const Splitting & split)  const = 0;
    virtual void   GeneratePoint(Splitting & split)       const = 0;    
    virtual inline const size_t size()                    const { return m_nout; } 
    virtual inline const ATOOLS::Flavour & GetSplit()     const { return m_split; }
    virtual inline const ATOOLS::Flavour_Vector & Flavs() const { return m_flavs; }
    virtual inline const std::vector<size_t> & Tags()     const { return m_tags; }
    virtual inline const kernel_type::code & Type()       const { return m_type; }
    virtual inline const std::string & Name()             const { return m_name; }
  };
  
  typedef ATOOLS::Getter_Function<SF_Base,Kernel_Info,
				  std::less<std::string> > SF_Getter;
}

#endif
