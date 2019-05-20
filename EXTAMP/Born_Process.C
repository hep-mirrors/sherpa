#include "EXTAMP/External_ME_Interface.H"
#include "EXTAMP/Born_Process.H"

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Selectors/Combined_Selector.H"

namespace EXTAMP {

  Born_Process::Born_Process(const PHASIC::Process_Info& pi) : Process(pi)
  {
    PHASIC::External_ME_Args args (pi.m_ii.GetExternal(),
				   pi.m_fi.GetExternal(),
				   pi.m_maxcpl);
    p_born_me = 
      External_ME_Interface::GetExternalBornME(args);
    p_born_me->SetCouplings(m_cpls);
  }

  double Born_Process::Partonic(const ATOOLS::Vec4D_Vector &p,
				const int mode)
  {
    if (!Selector()->Result()) return m_mewgtinfo.m_B=m_lastbxs=m_lastxs=0.0;

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
