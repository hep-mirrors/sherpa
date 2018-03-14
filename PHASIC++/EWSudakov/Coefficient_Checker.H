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

    Coefficient_Checker(const std::string& _procname,
                        const std::set<EWSudakov_Log_Type>& _activecoeffs):
      procname{ _procname },
      activecoeffs{ _activecoeffs }
    { }

    // check coefficients against a ref value when available
    bool CheckCoeffs(const Coeff_Map& coeffs,
                     const METOOLS::Spin_Amplitudes&);

  private:

    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    const std::map<Coeff_Map_Key, HelicityCoeffMap>& ReferenceCoeffs();

    const std::string procname;
    const std::set<EWSudakov_Log_Type> activecoeffs;

  };

}

#endif
