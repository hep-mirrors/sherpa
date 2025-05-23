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
    PHASIC::Process_Info rr_pi(pi);
    rr_pi.m_fi.m_nlotype=ATOOLS::nlo_type::rvirt;
    rr_pi.m_mincpl[0] = pi.m_mincpl[0];
    rr_pi.m_maxcpl[0] = pi.m_maxcpl[0];
    rr_pi.m_mincpl[1] = pi.m_mincpl[1]+1;
    rr_pi.m_maxcpl[1] = pi.m_maxcpl[1]+1;
    PRINT_VAR(rr_pi.m_maxcpl);
    // MODEL::Coupling_Data* aqcd=m_cpls.Get("Alpha_QCD");
    // MODEL::Coupling_Data* aqed=m_cpls.Get("Alpha_QED");
    p_loop_me = PHASIC::Virtual_ME2_Base::GetME2(rr_pi);
    if (!p_loop_me)  THROW(not_implemented, "Couldn't find RealVirtual ME for this process.");
    MODEL::s_model->GetCouplings(m_cpls);
    p_loop_me->SetSubType(ATOOLS::sbt::qed);
    /* Load color-correlated ME. TODO: orders */
    PHASIC::External_ME_Args args(rr_pi.m_ii.GetExternal(),
          rr_pi.m_fi.GetExternal(),
          rr_pi.m_maxcpl);
    p_corr_me = PHASIC::Color_Correlated_ME2::GetME2(args);  
    p_loop_me->SetCouplings(m_cpls);
    m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
    m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
    m_factor  = p_loop_me->AlphaQED()/2.0/M_PI;
    // m_factor  = 1.0/2.0/M_PI;
    double cplfac(1.0);
    // cplfac *= pow(p_loop_me->AlphaQCD(),rr_pi.m_mincpl[0]);
    cplfac *= pow(p_loop_me->AlphaQED(),rr_pi.m_mincpl[1]);
    // m_factor = cplfac/2.0/M_PI;
    m_factor = p_loop_me->AlphaQED()/2.0/M_PI;
    // m_factor =1.;
    // PRINT_VAR(m_factor);
    // PRINT_VAR(cplfac/2.0/M_PI);

    // m_factor *= p_loop_me->AlphaQED()/2.0/M_PI;

  }

RealVirtual::~RealVirtual()
{
 if(p_loop_me) delete p_loop_me;
 // if(p_scale)   delete p_scale;
}


double RealVirtual::Calc(const ATOOLS::Vec4D_Vector momenta, double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()));
}



  
double RealVirtual::Calc_V(const ATOOLS::Vec4D_Vector& p,
           const double B,
           const double mur)
  {
    double V(0.0), run_corr(0.0), scale(0.0);
    m_failcut = false;
    if(!p_rvproc->Trigger(p)) {
      m_failcut = true;
      return 0;
    }
    // p_loop_me->SetRenScale(mur);
    if(aqed->m_mode!=vpmode::off) {
     if(m_tchannel) scale = -(p[0]-p[2]).Abs2();  
     else scale = (p[0]+p[1]).Abs2();
     double dalpha = ((*aqed)(scale) - aqed->AqedThomson());
     run_corr = 4.*dalpha*B;
    }
    p_loop_me->Calc(p,B);
    double gammaborn = p_loop_me->ME_Born();
    switch(p_loop_me->Mode())
      {
      case 0:
        V =  m_factor * p_loop_me->ME_Finite()*gammaborn; break;

      case 1:
        V =  m_factor *  p_loop_me->ME_Finite(); break;
      case 2:
        // For Griffin
        THROW(not_implemented, "No Real-Virtuals implemented for this mode")
        break;
      default:
        THROW(not_implemented, "Loop ME mode not implemented: "+ATOOLS::ToString(p_loop_me->Mode()));
      }
    return V-run_corr;
  }
