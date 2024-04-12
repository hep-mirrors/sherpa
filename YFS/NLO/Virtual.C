#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Virtual.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
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
    loop_pi.m_mincpl[1] = pi.m_mincpl[1]+1;
    loop_pi.m_maxcpl[1] = pi.m_maxcpl[1]+1;
    p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
    if (!p_loop_me)  THROW(not_implemented, "Couldn't find virtual ME for this process.");
    MODEL::s_model->GetCouplings(m_cpls);
    p_loop_me->SetSubType(ATOOLS::sbt::qed);
    PHASIC::External_ME_Args args(loop_pi.m_ii.GetExternal(),
          loop_pi.m_fi.GetExternal(),
          loop_pi.m_maxcpl);

    double sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
    sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
    /* Load color-correlated ME. TODO: orders */
    p_corr_me = PHASIC::Color_Correlated_ME2::GetME2(args);  
    p_loop_me->SetCouplings(m_cpls);
    m_factor = p_loop_me->AlphaQED()/2.0/M_PI;
    // m_factor = 1./2.0/M_PI/137.036;
  }

Virtual::~Virtual()
{
 if(p_loop_me) delete p_loop_me;
 // if(p_scale)   delete p_scale;
}


double Virtual::Calc(const ATOOLS::Vec4D_Vector momenta, double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()));
}




double Virtual::Calc_V(const ATOOLS::Vec4D_Vector& p,
           const double B,
           const double mur)
  {
    double V(0.0), run_corr(0.0), scale(0.0);
    if(s_model->IsQEDRunning()) {
     if(m_tchannel) scale = -(p[0]-p[2]).Abs2();  
     else scale = (p[0]+p[1]).Abs2();
     double dalpha = ((*aqed)(scale) - aqed->AqedThomson());
     // PRINT_VAR(dalpha);
     run_corr = 4.*dalpha*B;
    }
    p_loop_me->Calc(p,B);
    V = p_loop_me->ME_Finite()*B-run_corr;
    return V*m_rescale_alpha*m_factor;
  }