#ifndef PHASIC_EWSudakov_CoefficientChecker_H
#define PHASIC_EWSudakov_CoefficientChecker_H

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"
#include "PHASIC++/EWSudakov/EWSudakov_Calculator.H"

#include <string>
#include <map>

namespace PHASIC {

  struct Mandelstam_Variables {
    Mandelstam_Variables(double _s, double _t, double _u):
      s {_s}, t {_t}, u {_u}
    {};
    double s {0.0};
    double t {0.0};
    double u {0.0};
  };

  class Coefficient_Checker {

  public:

    Coefficient_Checker(const std::string& _procname,
                        const std::set<EWSudakov_Log_Type>& _activecoeffs):
      procname{ _procname },
      activecoeffs{ _activecoeffs }
    { }

    // check coefficients against a ref value when available
    bool CheckCoeffs(const Coeff_Map& coeffs,
                     const METOOLS::Spin_Amplitudes&,
                     const Mandelstam_Variables&);

  private:

    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    std::map<Coeff_Map_Key, HelicityCoeffMap> ReferenceCoeffs(
      const Mandelstam_Variables&);

    const std::string procname;
    const std::set<EWSudakov_Log_Type> activecoeffs;

    bool CheckCoeff(const Coeff_Value&,
                    Complex ref,
                    const std::vector<int>& helicities,
		    const EWSudakov_Log_Type ewlt) const;

  };

}

#endif
