#ifndef NEUTRINOS_Tools_Form_Factor_Library_H
#define NEUTRINOS_Tools_Form_Factor_Library_H

#include "ATOOLS/Org/Getter_Function.H"
#include <iostream>

namespace NEUTRINOS {
  struct ff_type {
    enum code {
      none             = 0,
      constant         = 1,
      dipole           = 2,
      neutron_electric = 3,
      exponential      = 4,
      Gaussian         = 5,
      Kelly            = 6,
      BBBA             = 7,
      ArringtonHill    = 8,
      Helm             = 9,
      Lovato           = 10,
      unknown          = 99
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
      tensor       = 5,
      GE           = 10,
      GM           = 11,
      F1           = 20,
      F2           = 21
    };
  };
  std::ostream & operator<<(std::ostream & s,const cpl_info::code & cpl);
  
  struct ff_info {
    cpl_info::code      m_cpl;
    ff_type::code       m_type;
    std::vector<double> m_params;
    ff_info(const cpl_info::code & cpl, const ff_type::code & type=ff_type::none, 
	    const size_t & nparams=0) :
      m_cpl(cpl), m_type(type) { m_params.resize(nparams); }
  };
  std::ostream & operator<<(std::ostream & s,const ff_info & info);

  struct ModelParameters: public std::map<std::string, double> {
  };
  //std::ostream & operator<<(std::ostream & s,const ff_type::code & type);
  
  
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

  class Zero_Form_Factor : public Form_Factor_Base {
  public:
    Zero_Form_Factor(cpl_info::code & cpl) : Form_Factor_Base("Zero_FF", ff_info(cpl)) {}
    ~Zero_Form_Factor() {}
    double Calc(const double & q2) { return 0.; } 
  };

  class Dummy_Form_Factor : public Form_Factor_Base {
  public:
    Dummy_Form_Factor(cpl_info::code & cpl) : Form_Factor_Base("No_FF", ff_info(cpl)) {}
    ~Dummy_Form_Factor() {}
    double Calc(const double & q2) { return 1.; } 
  };

  class Constant_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm;
  public:
    Constant_Form_Factor(const ff_info & info);
    ~Constant_Form_Factor() {}
    double Calc(const double & q2);
  };

  class Dipole_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_invlambda2;
  public:
    Dipole_Form_Factor(const ff_info & info);
    ~Dipole_Form_Factor() {}
    double Calc(const double & q2);
  };

  class Neutron_Electric_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_invlambda2, m_mass2, m_pref;
  public:
    Neutron_Electric_Form_Factor(const ff_info & info);
    ~Neutron_Electric_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class Exponential_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_invlambda2;
  public:
    Exponential_Form_Factor(const ff_info & info);
    ~Exponential_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class Gaussian_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_invlambda2;
  public:
    Gaussian_Form_Factor(const ff_info & info);
    ~Gaussian_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class Kelly_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_invlambda2, m_mass2, m_A[4], m_En[2];
    const double ModifiedDipole(const double & q2) const;
    const double Polynomial(const double & q2) const;
  public:
    Kelly_Form_Factor(const ff_info & info);
    ~Kelly_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class BBBA_Form_Factor : public Form_Factor_Base {
  private:
    double m_norm, m_mass2, m_A[4], m_B[4];
    const double Numerator(const double & tau) const;
    const double Denominator(const double & tau) const;
  public:
    BBBA_Form_Factor(const ff_info & info);
    ~BBBA_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class ArringtonHill_Form_Factor : public Form_Factor_Base {
  private:
    double m_mass2, m_tcut, m_t0, m_A[13];
    const double ZExpand(const double & z) const;
  public:
    ArringtonHill_Form_Factor(const ff_info & info);
    ~ArringtonHill_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class Helm_Form_Factor : public Form_Factor_Base {
  private:
    double m_s, m_r, m_kappa, m_hbarc;
  public:
    Helm_Form_Factor(const ff_info & info);
    ~Helm_Form_Factor() {}
    double Calc(const double & q2);
  };
  
  class Lovato_Form_Factor : public Form_Factor_Base {
  private:
    double m_b, m_c[5], m_hbarc;
    const double Polynomial(const double & bx) const;
  public:
    Lovato_Form_Factor(const ff_info & info);
    ~Lovato_Form_Factor() {}
    double Calc(const double & q2);
  };

}

#endif
