#ifndef CFPSHOWER_Kinematics_Kinematics_Base_H
#define CFPSHOWER_Kinematics_Kinematics_Base_H

#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Tools/Configuration.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace CFPSHOWER {
  class Kernel;
  
  class Kinematics_Base {
  protected:
    std::string    m_name;
    kin_type::code m_scheme;
    Kernel *       p_kernel;
    ATOOLS::Vec4D  m_psplit, m_pspect, m_pboth;
    double         m_msplit, m_mspect, m_msplit2, m_mspect2, m_m[2], m_m2[2];
    double         m_Q2, m_Q, m_z[2], m_y, m_weight;
    bool           m_ismassive;

    std::map<Parton *, ATOOLS::Vec4D> m_allmomenta;
    std::map<Parton *, double>        m_allmasses2;
    
    static double Lambda(const double & a,const double & b,const double & c);
    static double Lambda2(const double & a,const double & b,const double & c);
    virtual void InitSystem(Splitting & split, const ATOOLS::Mass_Selector * msel);
  public:
    Kinematics_Base(const Kernel_Info & info);
    ~Kinematics_Base();

    inline void  SetKernel(Kernel * kernel) { p_kernel = kernel; }
    virtual bool Init(Splitting & split, Configuration & config,
		      const ATOOLS::Mass_Selector * msel)                       = 0;
    virtual bool operator()(Splitting & split, Configuration & config)          = 0;
    virtual void CalculateJacobean(Splitting & split,Configuration & config)    = 0;
    virtual bool UpdateSystem(Splitting & split, Configuration & config)        = 0;

    
    inline virtual const double            & Weight()        const { return m_weight; }
    inline virtual const std::string       & Name()          const { return m_name; }
    inline virtual const kin_type::code    & Scheme()        const { return m_scheme; }
    inline std::map<Parton*,ATOOLS::Vec4D> * GetMomentaMap()       { return &m_allmomenta; }
    inline std::map<Parton*,double>        * GetMasses2Map()       { return &m_allmasses2; }
  };

  typedef ATOOLS::Getter_Function<Kinematics_Base,Kernel_Info,
				  std::less<std::string> > Kinematics_Getter;

}

#endif
