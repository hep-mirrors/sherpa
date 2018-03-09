#ifndef PHASIC_EWSudakov_CoefficientChecker_H
#define PHASIC_EWSudakov_CoefficientChecker_H

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/EWSudakov/Sudakov.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include <string>
#include <map>

namespace PHASIC {

  class Coefficient_Checker {

  public:

    Coefficient_Checker(const std::string& _procname):
      procname{ _procname }
    { }

    // check coefficients against a ref value when available
    bool CheckCoeffs(const LogType_Coeff_Map& coeffs,
                     const METOOLS::Spin_Amplitudes&);
    bool CheckAngularCoeffs(const std::vector<LegIndizes_Coeff_Map>& coeffs,
                            const METOOLS::Spin_Amplitudes&);

  private:

    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    const std::map<std::string, HelicityCoeffMap>& ReferenceCoeffs();
    const std::map<Two_Leg_Indizes, HelicityCoeffMap>& AngularReferenceCoeffs();

    const std::string procname;

  };

}

#endif
