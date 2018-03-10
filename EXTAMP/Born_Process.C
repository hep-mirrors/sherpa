#include "EXTAMP/External_ME_Interface.H"
#include "EXTAMP/Born_Process.H"

#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace EXTAMP {

  Born_Process::Born_Process(const PHASIC::Process_Info& pi) : Process(pi)
  {
    p_born_me = External_ME_Interface::GetExternalBornME(pi);
    p_born_me->SetCouplings(m_cpls);
  }

  double Born_Process::Partonic(const ATOOLS::Vec4D_Vector &p,
				const int mode)
  {

    /* Maybe move to PHASIC::Single_Process */
    ScaleSetter()->CalculateScale(p);

    double dxs = p_born_me->Calc(p)/NormFac();

    /* Single_Process derivatives are responsible for storing the
       return value in m_lastdxs and for filling the m_mewgtinfo
       inherited from PHASIC::Process_Base */
    m_mewgtinfo.m_K = 1.0;
    m_mewgtinfo.m_B = dxs;
    m_lastxs        = dxs;

    return dxs;
  }

}
