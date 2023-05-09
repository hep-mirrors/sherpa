#ifndef NEUTRINOS_Tools_Propagator_Library_H
#define NEUTRINOS_Tools_Propagator_Library_H

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Getter_Function.H"
#include <iostream>

using namespace std;
using namespace ATOOLS;

namespace ATOOLS {
  class Scoped_Settings;
}

namespace NEUTRINOS {
  struct prop_type {
    enum code {
      none             = 0,
      massless         = 1,
      massive          = 2,
      unstable         = 3, 
      unknown          = 99
    };
  };
  std::ostream & operator<<(std::ostream & s,const prop_type::code & type);

  struct prop_info {
    prop_type::code       m_type;
    double                m_mass;
    double                m_width;
    prop_info(const prop_type::code & type=prop_type::none, const double & mass=0.0, const double & width=0.0) :
      m_type(type), m_mass(mass) , m_width(width) {}
  };
  std::ostream & operator<<(std::ostream & s,const prop_info & info);

  class Propagator_Base {
  protected:
    std::string    m_name;
    double         m_mass;
    double         m_width;
    prop_type::code  m_type;
  public:
    Propagator_Base(const std::string & name,const prop_info & info) :
      m_name(name), m_type(info.m_type){}
    ~Propagator_Base() {}
    virtual Complex Calc(const double & p2) = 0;
    const std::string    & Name() const { return m_name; }
    const prop_type::code  & Type() const { return m_type; } 
  };

  class Dummy_Prop : public Propagator_Base {
  public:
    Dummy_Prop(const prop_info & info);
    ~Dummy_Prop() {}
    Complex Calc(const double & p2);
  };

  class Massless_Prop : public Propagator_Base {
  public:
    Massless_Prop(const prop_info & info);
    ~Massless_Prop() {}
    Complex Calc(const double & p2);
  };

  class Massive_Prop : public Propagator_Base {
  private:
    double mass;
  public:
    Massive_Prop(const prop_info & info);
    ~Massive_Prop() {}
    Complex Calc(const double & p2);
  };

  class Unstable_Prop : public Propagator_Base {
  private:
    double mass, width;
  public:
    Unstable_Prop(const prop_info & info);
    ~Unstable_Prop() {}
    Complex Calc(const double & p2);
  };

}

#endif
