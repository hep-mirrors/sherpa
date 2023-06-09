#ifndef YFS_REAL_H
#define YFS_REAL_H


#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "YFS/Main/YFS_Base.H"


namespace PHASIC {
class Tree_ME2_Base;
}

namespace ATOOLS {
class NLO_subevtlist;
}
namespace EXTAMP {
class CS_Dipole;
class Dipole_Wrapper_Process;
}
namespace YFS {

class Real : public YFS_Base {

public:

  Real(const PHASIC::Process_Info& pi);

  ~Real();
  // void InitME(MODEL::Model_Base *const model);
  double Calc(const ATOOLS::Vec4D_Vector momenta);
   PHASIC::Tree_ME2_Base     *p_real_me;
  double Calc_R(const ATOOLS::Vec4D_Vector& p);
  MODEL::Coupling_Map m_cpls;
  PHASIC::Scale_Setter_Base *p_scale;
  PHASIC::Tree_ME2_Base* p_born_me;
  double m_sym, m_factor;
};
}

#endif