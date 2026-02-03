#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Virtual.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"


using namespace MODEL;
using namespace YFS;

Virtual::Virtual(const PHASIC::Process_Info& pi)
  {
    /* Load loop ME */
    PHASIC::Process_Info loop_pi(pi);
    loop_pi.m_fi.m_nlotype=ATOOLS::nlo_type::loop;
    loop_pi.m_mincpl[0] = pi.m_mincpl[0];
    loop_pi.m_maxcpl[0] = pi.m_maxcpl[0];
    loop_pi.m_mincpl[1] = pi.m_mincpl[1];
    loop_pi.m_maxcpl[1] = pi.m_maxcpl[1]+1;
    p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
    if (!p_loop_me)  {
      msg_Error()<<loop_pi<<std::endl;
      THROW(not_implemented, "Couldn't find virtual ME for this process.");
    }
    MODEL::s_model->GetCouplings(m_cpls);
    p_loop_me->SetSubType(ATOOLS::sbt::qed);
    p_loop_me->SetCouplings(m_cpls);
    m_factor = p_loop_me->AlphaQED()/2.0/M_PI;
  }

Virtual::~Virtual()
{
 if(p_loop_me) delete p_loop_me;
}


double Virtual::Calc(const ATOOLS::Vec4D_Vector momenta, double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()));
}




double Virtual::Calc_V(const ATOOLS::Vec4D_Vector& p,
           const double B,
           const double mur)
  {
    p_loop_me->SetRenScale(100);
    m_failcut = false;
    // if(m_nlocuts && !p_vproc->Trigger(p)) {
    //   m_failcut = true;
    //   return 0;
    // }
    double V(0.0), run_corr(0.0), scale(0.0);
    if(aqed->m_mode!=vpmode::off) {
     if(m_tchannel==2) scale = -(p[0]-p[2]).Abs2();  
     else scale = (p[0]+p[1]).Abs2();
     double dalpha = ((*aqed)(scale) - aqed->AqedThomson());
    //  run_corr = 4.*dalpha*B;
    }
    p_loop_me->Calc(p,B);
    switch(p_loop_me->Mode())
      {
      case 0:
        V =  m_factor *  p_loop_me->ME_Finite() * B ; break;
      case 1:
        // For Griffin
        V =  p_loop_me->ME_Finite();
        break;
      default:
        THROW(not_implemented, "Loop ME mode not implemented: "+ATOOLS::ToString(p_loop_me->Mode()));
      }
    return V-run_corr*m_factor;
  }