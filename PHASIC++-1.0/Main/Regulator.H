#ifndef Regulator_H
#define Regulator_H

#include "Regulator_Base.H"

namespace PHASIC {

  class Identity: public Regulator_Base {
  protected:

    // constructor
    Identity(Integrable_Base *const xs,const std::vector<double> &parameters);
    
    // member functions
    double operator()(const double dsigma) const;
    double operator[](const double scale) const;
    
    friend class Regulator_Base;

  };// end of class Identity

  class LO_QCD_Regulator: public Regulator_Base {
  protected:

    // constructor
    LO_QCD_Regulator(Integrable_Base *const xs,
		     const std::vector<double> &parameters);

    // member functions
    double operator()(const double dsigma) const;
    double operator[](const double scale) const;
    
    friend class Regulator_Base;

  };// end of class LO_QCD_Regulator

}// end of namespace EXTRAXS

#endif
