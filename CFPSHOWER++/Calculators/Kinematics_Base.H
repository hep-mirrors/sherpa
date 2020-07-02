#ifndef CFPSHOWER_Calculators_Kinematics_Base_H
#define CFPSHOWER_Calculators_Kinematics_Base_H

#include "CFPSHOWER++/Tools/Splitting.H"
#include "CFPSHOWER++/Tools/Configuration.H"
#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/Getter_Function.H"

namespace CFPSHOWER {
  class Kinematics_Base {
  protected:
    log_type::code    m_logtype;
    kernel_type::code m_type;
    std::string       m_name;
    double            m_weight;
        
    virtual void Init(const Splitting & split);

    static ATOOLS::Vec4D SumMomenta(const Configuration & config); 
    static ATOOLS::Vec4D LT(const ATOOLS::Vec4D &a,
			    const ATOOLS::Vec4D &b,
			    const ATOOLS::Vec4D &c);
    static double Lambda(const double & a,const double & b,const double & c);
    static double Lambda2(const double & a,const double & b,const double & c);
  public:
    Kinematics_Base(log_type::code logtype=log_type::none,
		    kernel_type::code type=kernel_type::none);
    ~Kinematics_Base();

    virtual bool operator()(Splitting & split,Configuration & config)   = 0;
    virtual void CalculateWeight(Splitting & split)                     = 0;
    virtual bool UpdateSystem(Splitting & split,Configuration & config) = 0;
    inline virtual const double & Weight()    { return m_weight; }
    inline virtual const std::string & Name() { return m_name; }
  };

  typedef ATOOLS::Getter_Function<Kinematics_Base,Kernel_Info,
				  std::less<std::string> > Kinematics_Getter;

}

#endif
