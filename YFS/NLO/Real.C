#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Real.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "EXTAMP/External_ME_Interface.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace YFS;
using namespace MODEL;


Real::Real(const PHASIC::Process_Info& pi)  {
   /* Load Real ME */
   PRINT_VAR(pi);
   PHASIC::Process_Info real_pi(pi);
   real_pi.m_fi.m_nlotype = ATOOLS::nlo_type::lo;
   real_pi.m_mincpl[0] = pi.m_mincpl[0];
   real_pi.m_maxcpl[0] = pi.m_maxcpl[0];
   real_pi.m_mincpl[1] = pi.m_mincpl[1] + 1;
   real_pi.m_maxcpl[1] = pi.m_maxcpl[1] + 1;
   p_real_me = PHASIC::Tree_ME2_Base::GetME2(real_pi);
   if (!p_real_me)  THROW(not_implemented, "Couldn't find real ME for this process.");
   MODEL::s_model->GetCouplings(m_cpls);
   /* Load color-correlated ME. TODO: orders */
   PHASIC::External_ME_Args args(real_pi.m_ii.GetExternal(),
                                 real_pi.m_fi.GetExternal(),
                                 real_pi.m_maxcpl);
   p_real_me->SetCouplings(m_cpls);
   m_sym  = ATOOLS::Flavour::FSSymmetryFactor(args.m_outflavs);
   m_sym *= ATOOLS::Flavour::ISSymmetryFactor(args.m_inflavs);
   ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
  m_factor = p_real_me->AlphaQED()/2.0/M_PI;

}

Real::~Real() {
}

double Real::Calc_R(const ATOOLS::Vec4D_Vector& p)
  {
    double R = p_real_me->Calc(p);
    return R*m_rescale_alpha/m_sym;
  }