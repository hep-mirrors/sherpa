#ifndef PHASIC_EWSudakov_CoefficientChecker_H
#define PHASIC_EWSudakov_CoefficientChecker_H

#include "METOOLS/Main/Spin_Structure.H"

#include <string>
#include <map>

namespace PHASIC {

  class Coefficient_Checker {

  public:

    Coefficient_Checker(const std::string& _procname):
      procname{ _procname }
    { }

    // check coefficients against a ref value when available
    bool CheckCoeffs(const std::map<std::string, std::vector<Complex>>& coeffs,
                     const METOOLS::Spin_Amplitudes&);

  private:

    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    const std::map<std::string, HelicityCoeffMap>& ReferenceCoeffs();

    const std::string procname;

  };

}

#endif
