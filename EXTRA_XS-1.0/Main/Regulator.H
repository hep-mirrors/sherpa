#ifndef Regulator_H
#define Regulator_H

#include "Regulator_Base.H"

namespace EXTRAXS {

  class Identity: public Regulator_Base {
  protected:

    // constructor
    Identity(XS_Base *const xs,const std::vector<double> &parameters);
    
    // member functions
    double operator()(const double dsigma) const;
    double operator[](const double scale) const;
    
    friend class Regulator_Base;

  };// end of class Identity

  class Massive_Propagator: public Regulator_Base {
  protected:

    // constructor
    Massive_Propagator(XS_Base *const xs,
		       const std::vector<double> &parameters);

    // member functions
    double operator()(const double dsigma) const;
    double operator[](const double scale) const;
    
    friend class Regulator_Base;

  };// end of class Massive_Propagator

}// end of namespace EXTRAXS

#endif
