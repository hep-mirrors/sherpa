#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/RealVirtual.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace MODEL;
using namespace YFS;

RealVirtual::RealVirtual(const PHASIC::Process_Info& pi)
  {

    /* Load loop ME */
    PHASIC::Process_Info loop_pi(pi);
    loop_pi.m_fi.m_nlotype=ATOOLS::nlo_type::loop;
    loop_pi.m_mincpl[0] = pi.m_mincpl[0];
    loop_pi.m_maxcpl[0] = pi.m_maxcpl[0];
    loop_pi.m_mincpl[1] = pi.m_mincpl[1]+2;
    loop_pi.m_maxcpl[1] = pi.m_maxcpl[1]+2;
    p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
    if (!p_loop_me)  THROW(not_implemented, "Couldn't find RealVirtual ME for this process.");
    MODEL::s_model->GetCouplings(m_cpls);
    p_loop_me->SetSubType(ATOOLS::sbt::qed);
    /* Load color-correlated ME. TODO: orders */
    PHASIC::External_ME_Args args(loop_pi.m_ii.GetExternal(),
          loop_pi.m_fi.GetExternal(),
          loop_pi.m_maxcpl);
    p_corr_me = PHASIC::Color_Correlated_ME2::GetME2(args);  
    p_loop_me->SetCouplings(m_cpls);
    m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
    m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
    m_factor = p_loop_me->AlphaQED()/2.0/M_PI;

  }

RealVirtual::~RealVirtual()
{
 // if(p_loop_me) delete p_loop_me;
 // if(p_scale)   delete p_scale;
}


double RealVirtual::Calc(const ATOOLS::Vec4D_Vector momenta, double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()));
}




double RealVirtual::Calc_V(const ATOOLS::Vec4D_Vector& p,
           const double B,
           const double mur)
  {
    double V(0.0);
    // p_loop_me->SetRenScale(mur);
    p_loop_me->Calc(p,B);
    V = p_loop_me->ME_Finite();
    return V*m_rescale_alpha/m_sym*m_factor*m_factor;
  }
