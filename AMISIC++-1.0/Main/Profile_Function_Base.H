#ifndef Profile_Function_Base_H
#define Profile_Function_Base_H

#include "Message.H"
#include <vector>

namespace AMISIC {

  struct pft {

    enum code {
      none            = 0,
      flat            = 1,
      exponential     = 2,
      gaussian        = 3,
      double_gaussian = 4
    };

  };// end of struct pft

  std::ostream &operator<<(std::ostream &ostr,const pft::code code);

  class Profile_Function_Base {
  protected:

    pft::code m_type;

    // constructor
    Profile_Function_Base(const pft::code code);

    template <class Profile_Type>
    static Profile_Function_Base *CreateProfile(const std::string &type,
						const std::vector<double> &parameters);

  public:
  
    static Profile_Function_Base *SelectProfile(const std::string &type,
						const std::vector<double> &parameters);

    // member functions
    virtual double Value(const double b) const = 0;
    virtual double Integral(const double b) const = 0;

    // inline functions
    inline pft::code Type() const { return m_type; }

    inline double operator()(const double b) const { return Value(b);    }
    inline double operator[](const double b) const { return Integral(b); }

  };// end of class Profile_Function_Base

}// end of namespace AMISIC

#endif
