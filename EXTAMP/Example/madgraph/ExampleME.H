#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "ATOOLS/Phys/Flavour.H"

#include "CPPProcess.h"

namespace PHASIC {
  class Process_Info;
}


class ExampleME: public PHASIC::Tree_ME2_Base {

public:
  
  ExampleME(const PHASIC::Process_Info &);
    
    ~ExampleME() {}
  
  double Calc(const ATOOLS::Vec4D_Vector& momenta);

  int OrderQCD(const int &id=-1) const { return m_oqcd; }
  int OrderEW (const int &id=-1) const { return m_oew;  }

  double RunningCouplingFactor() const;
  
private:

  size_t m_nin, m_nout;
  size_t m_oew,m_oqcd;

  PHASIC::Process_Info m_pinfo;

  double m_symfac;

  CPPProcess m_mg_proc;
};
