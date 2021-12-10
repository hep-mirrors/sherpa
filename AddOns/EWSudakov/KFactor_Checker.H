#ifndef PHASIC_EWSudakov_KFactorChecker_H
#define PHASIC_EWSudakov_KFactorChecker_H

#include "AddOns/EWSudakov/EWSudakov_Amplitudes.H"
#include "AddOns/EWSudakov/EWSudakov_Calculator.H"

#include <string>
#include <map>

namespace PHASIC {

  class KFactor_Checker {

  public:
    KFactor_Checker(const std::string& _procname) : procname {_procname} {}

    // set a log file name for writing out checked results and references for
    // later analysis
    void SetLogFileName(const std::string& f) { logfilename = f; }

    // check coefficients against a ref value when available
    bool CheckKFactor(double kfac,
                      const ATOOLS::Mandelstam_Variables&,
                      const EWGroupConstants&);

  private:
    double ReferenceKFactor(const ATOOLS::Mandelstam_Variables&,
                            const EWGroupConstants&);

    const std::string procname;
    const std::set<ATOOLS::EWSudakov_Log_Type> activecoeffs;

    std::string logfilename;

    bool CheckKFactor(double kfac, double ref) const;

  };

}

#endif
