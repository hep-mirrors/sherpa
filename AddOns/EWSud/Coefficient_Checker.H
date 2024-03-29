#ifndef EWSudakov_CoefficientChecker_H
#define EWSudakov_CoefficientChecker_H

#include "AddOns/EWSud/Amplitudes.H"
#include "AddOns/EWSud/Calculator.H"

#include <string>
#include <map>

namespace EWSud {

  class Coefficient_Checker {

  public:

    Coefficient_Checker(
        const std::string& _procname,
        const std::set<EWSudakov_Log_Type>& _activecoeffs)
        : procname {_procname}, activecoeffs {_activecoeffs}
    { }

    // set a log file name for writing out checked results and references for
    // later analysis
    void SetLogFileName(const std::string& f) { logfilename = f; }

    // check coefficients against a ref value when available
    bool CheckCoeffs(const Coeff_Map& coeffs,
                     const METOOLS::Spin_Amplitudes&,
                     const Mandelstam_Variables&,
                     const EWGroupConstants&);

  private:

    using HelicityCoeffMap = std::map<std::vector<int>, double>;
    std::map<Coeff_Map_Key, HelicityCoeffMap>
    ReferenceCoeffs(const Mandelstam_Variables&, const EWGroupConstants&);

    const std::string procname;
    const std::set<EWSudakov_Log_Type> activecoeffs;

    std::string logfilename;

    bool CheckCoeff(const Coeff_Value&,
                    Complex ref,
                    const std::vector<int>& helicities,
		    const EWSudakov_Log_Type ewlt) const;

  };

}

#endif
