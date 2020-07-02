#ifndef CFPSHOWER_Calculators_SF_Base_H
#define CFPSHOWER_Calculators_SF_Base_H

#include "CFPSHOWER++/Calculators/Kinematics_Base.H"
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
    ATOOLS::Flavour        m_split, m_spect;
    double                 m_msplit, m_msplit2, m_mspect, m_mspect2, m_msum;
    ATOOLS::Flavour_Vector m_flavs;
    std::vector<size_t>    m_tags;
    std::vector<double>    m_m, m_m2;
    std::vector<ATOOLS::Vec4D>        m_moms;
    std::vector<std::vector<double> > m_pp;
    size_t                 m_nout;
    bool                   m_ismassive;
    double                 m_shat, m_Ehat;
    ATOOLS::Vec4D          m_psplit, m_pspect, m_pboth;

    
    std::string            m_name;
    log_type::code         m_logtype;
    kernel_type::code      m_type;
    double                 m_weight;
    int                    m_CMW, m_softcorr, m_endpoint;
    subtract::code         m_subtract;

    static ATOOLS::Vec4D SumMomenta(const Configuration & config); 
    static ATOOLS::Vec4D LT(const ATOOLS::Vec4D &a,
			    const ATOOLS::Vec4D &b,
			    const ATOOLS::Vec4D &c);
    static double Lambda(const double & a,const double & b,const double & c);
    static double Lambda2(const double & a,const double & b,const double & c);
  public:
    SF_Base(const Kernel_Info & info);
    ~SF_Base() {}
    virtual bool   Init(const Splitting & split,const ATOOLS::Mass_Selector * msel);
    virtual double operator()(const Splitting & split)                    = 0;
    virtual double Integral(const Splitting & split)                const = 0;
    virtual double OverEstimate(const Splitting & split)            const = 0;
    virtual void   GeneratePoint(Splitting & split)                 const = 0;
    virtual bool   Construct(Splitting & split,Configuration & config)    = 0;
    virtual void   CalculateJacobean(Splitting & split)                   = 0;
    virtual void   CalculateInvariants(Splitting & split)                 = 0;
    virtual bool   UpdateSystem(Splitting & split,Configuration & config) = 0;

    virtual inline const size_t size()                              const { return m_flavs.size(); }
    virtual inline const ATOOLS::Flavour & GetSplit()               const { return m_split; }
    virtual inline const ATOOLS::Flavour_Vector & GetFlavs()        const { return m_flavs; }
    virtual inline const ATOOLS::Vec4D & Mom(const size_t & pos)    const { return m_moms[pos]; }
    virtual inline const ATOOLS::Flavour & Flav(const size_t & pos) const { return m_flavs[pos]; }
    virtual inline const std::vector<size_t> & Tags()               const { return m_tags; }
    virtual inline const size_t & Tags(const size_t & pos)          const { return m_tags[pos]; }
    virtual inline const kernel_type::code & Type()                 const { return m_type; }
    virtual inline const log_type::code & LogType()                 const { return m_logtype; }
    virtual inline const std::string & Name()                       const { return m_name; }
    virtual inline const double & Jacobean()                        const { return m_weight; }
    virtual inline void SetName(const std::string & name)                 { m_name = name; }

    virtual inline const bool IsSubtract(subtract::code sub=subtract::undefined) const {
      return (sub==subtract::undefined ? m_subtract>0 : int(sub & m_subtract)>0);
    }
  };
  
  typedef ATOOLS::Getter_Function<SF_Base,Kernel_Info,std::less<std::string> > SF_Getter;
}

#endif
