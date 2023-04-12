#ifndef NEUTRINOS_Tools_Form_Factor_Library_H
#define NEUTRINOS_Tools_Form_Factor_Library_H

#include "ATOOLS/Org/Getter_Function.H"
#include <iostream>

namespace NEUTRINOS {
  struct ff_type {
    enum code {
      unknown     = 0,
      none        = 1,
      dipole      = 2,
      exponential = 3,
      Gaussian    = 4
    };
  };
  std::ostream & operator<<(std::ostream & s,const ff_type::code & type);

  struct cpl_info {
    enum code {
      unknown      = 0,
      scalar       = 1,
      pseudoscalar = 2,
      vector       = 3,
      axialvector  = 4,
      tensor       = 5
    };
  };
  std::ostream & operator<<(std::ostream & s,const cpl_info::code & cpl);
  
  struct ff_info {
    cpl_info::code      m_cpl;
    ff_type::code       m_type;
    std::vector<double> m_params;
    ff_info(const cpl_info::code & cpl, const ff_type::code & type, 
	    const size_t & nparams=2) :
      m_cpl(cpl), m_type(type) { m_params.resize(nparams); }
  };
  std::ostream & operator<<(std::ostream & s,const ff_info & info);

  struct ModelParameters: public std::map<std::string, double> {
  };
  
  
  class Form_Factor_Base {
  protected:
    std::string    m_name;
    ff_type::code  m_type;
    cpl_info::code m_cpl;
  public:
    Form_Factor_Base(const std::string & name,const ff_info & info) :
      m_name(name), m_type(info.m_type), m_cpl(info.m_cpl) {}
    ~Form_Factor_Base() {}
    virtual double Calc(const double & q2) = 0;
    const std::string    & Name() const { return m_name; }
    const ff_type::code  & Type() const { return m_type; } 
    const cpl_info::code & Cpl()  const { return m_cpl; } 
  };
}

//#define DECLARE_NN_FORMFACTOR_GETTER(NAME,TAG)		\
//  DECLARE_GETTER(NAME,TAG,Form_Factor_Base,ff_info);		\
//  void ATOOLS::Getter<Form_Factor_Base,ff_info,NAME>::		\
//  PrintInfo(std::ostream &str,const size_t width) const {       \
//    str<<#TAG;                                                  \
//  }

#endif
