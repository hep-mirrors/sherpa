#ifndef Regulator_Base_H
#define Regulator_Base_H

namespace EXTRAXS { class Regulator_Base; }

#include "Single_XS.H"

namespace EXTRAXS {

  struct rf {
    
    enum code {
      none               = 1,
      massive_propagator = 2
    };
    
  };// end of struct rf

  std::ostream &operator<<(std::ostream &ostr,const rf::code code);

  class Regulator_Base {
  protected:

    const rf::code m_type;

    std::vector<double> m_parameters;

    EXTRAXS::Single_XS *const p_xs;

    // constructor
    Regulator_Base(Single_XS *const xs,const std::vector<double> &parameters);

    template <class Type>
    static Regulator_Base *SelectRegulator(Single_XS *const xs,const std::string &regulator,
					   const std::vector<double> &parameters);

  public:

    // destructor
    virtual ~Regulator_Base();

    // member functions
    virtual double operator()(const double dsigma) const;
    virtual double operator[](const double scale) const;
    
    static Regulator_Base *GetRegulator(Single_XS *const xs,const std::string &regulator,
					const std::vector<double> &parameters);

    // inline functions
    inline rf::code Type() const { return m_type; }
    
  };// end of class Regulator_Base

}// end of namespace EXTRAXS

#endif
