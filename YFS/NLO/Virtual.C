#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/Virtual.H"

#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"
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
    // massification registers its massless species before the provider
    // registers the process
    if (Massification::Enabled())
      p_massif.reset(new Massification(loop_pi.ExtractFlavours(), "virtual"));
    p_loop_me.reset(PHASIC::Virtual_ME2_Base::GetME2(loop_pi));
    if (!p_loop_me)  {
      msg_Error()<<loop_pi<<std::endl;
      THROW(not_implemented, "Couldn't find virtual ME for this process.");
    }
    if (p_massif) p_massif->Attach(p_loop_me.get());
    MODEL::s_model->GetCouplings(m_cpls);
    p_loop_me->SetSubType(ATOOLS::sbt::qed);
    p_loop_me->SetCouplings(m_cpls);
    m_factor = p_loop_me->AlphaQED()/2.0/M_PI;
    /*
      YFS: EW_HIGHER_ORDER griffin: GRIFFIN's NNLO - NLO on top of the loop
      provider's full NLO. Not GRIFFIN as the provider: its IR-subtracted
      matrix elements carry no photonic QED virtual, so Loop_Generator Griffin
      loses the O(alpha) photonic virtual (Z-pole mu mu: -5.5% on sigma_fid,
      A_FB changes sign; 2026-09-30).
    */
    const ewhigherorder::code ewho(ATOOLS::Settings::GetMainSettings()["YFS"]
        ["EW_HIGHER_ORDER"].SetDefault(ewhigherorder::none).Get<ewhigherorder::code>());
    if (ewho == ewhigherorder::griffin) {
      PHASIC::Process_Info ho_pi(loop_pi);
      ho_pi.m_loopgenerator = "GriffinHO";
      p_ewho_me.reset(PHASIC::Virtual_ME2_Base::GetME2(ho_pi));
      if (!p_ewho_me)
        THROW(fatal_error, "YFS: EW_HIGHER_ORDER griffin needs Griffin in "
                           "ME_GENERATORS and a 2->2 s-channel process");
      p_ewho_me->SetCouplings(m_cpls);
      msg_Info() << "YFS: EW_HIGHER_ORDER griffin - adding GRIFFIN NNLO - NLO to "
                    "the Born-level virtual\n";
    }
  }

double Virtual::EWHigherOrder(const ATOOLS::Vec4D_Vector &momenta)
{
  if (!p_ewho_me) return 0.;
  p_ewho_me->Calc(momenta);
  return p_ewho_me->ME_Finite();
}

Virtual::~Virtual()
{
}


double Virtual::Calc(const ATOOLS::Vec4D_Vector momenta, const double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()));
}




double Virtual::CalcInFrame(const ATOOLS::Vec4D_Vector &momenta,
                            const ATOOLS::Vec4D_Vector &peval, double born){
  return Calc_V(momenta,born,sqr(rpa->gen.Ecms()),&peval);
}

double Virtual::Calc_V(const ATOOLS::Vec4D_Vector& p,
           const double B,
           const double mur,
           const ATOOLS::Vec4D_Vector *peval)
  {
    // Renormalisation scale as specified by SCALES in the runcard, for an
    // arbitrary (momentum-dependent) expression. The Born process's scale setter
    // has already evaluated it for THIS event: Single_Process::Differential runs
    // Partonic() - which calls ScaleSetter()->CalculateScale(p) - before the YFS
    // weight is built, so the value is simply read here. Do NOT call
    // CalculateScale again: it would overwrite the Born's cached scales while
    // Differential is still using them, and the virtual multiplies that same
    // Born so it must use the same mu. Already squared; SetRenScale takes mu^2.
    // ScaleSetter(1) is the mapped-process setter, which is the one
    // Single_Process actually evaluates.
    p_loop_me->SetRenScale(p_vproc ? p_vproc->ScaleSetter(1)->Scale(ATOOLS::stp::ren)
                                   : mur);
    m_failcut = false;
    if(m_nlocuts && !p_vproc->Trigger(p)) {
      m_failcut = true;
      return 0;
    }
    double V(0.0);
  /*
    No running-coupling subtraction. The old run_corr = 4 (alpha(s) - alpha(0)) B
    (Alpha_QED: VPMODE != off) removed a term the virtual never double counted:
    the Born does NOT run with VPMODE (YFS.LO and YFS.BR are bit-identical with
    VPMODE off and Full, 2026-09-30), and the loop provider's virtual already
    contains the full photon vacuum polarisation. It also used the absolute
    alpha difference, ~70x below a relative 2 Delta alpha. Removed; see
    NOTES-hvp-design-2026-09-30.md.
  */
    const ATOOLS::Vec4D_Vector &pe(peval ? *peval : p);
    if (p_massif) {
      // the provider's V/B at the massless image, rebuilt into the massive
      // one: a ratio, so it multiplies Sherpa's massive B
      p_massif->Count(pe);
      p_loop_me->Calc(p_massif->MasslessImage(pe),B);
      p_massif->CheckDoublePole(p_loop_me->ME_E2());
    }
    else p_loop_me->Calc(pe,B);
    switch(p_loop_me->Mode())
      {
      case 0:
        V =  m_factor * (p_loop_me->ME_Finite()
                         + (p_massif ? p_massif->Z(sqr(p_loop_me->IRscale())) : 0.)) * B;
        break;
      case 1:
        // For Griffin
        if (p_massif) THROW(not_implemented, "YFS: VIRTUAL_MASSIFICATION needs a"
                            " provider returning V/B (Mode 0)");
        V =  p_loop_me->ME_Finite();
        break;
      default:
        THROW(not_implemented, "Loop ME mode not implemented: "+ATOOLS::ToString(p_loop_me->Mode()));
      }
    return V;
  }
double Virtual::E2() const
{
  return p_loop_me->ME_E2() + (p_massif ? p_massif->ZPole2() : 0.);
}

double Virtual::E1() const
{
  return p_loop_me->ME_E1()
    + (p_massif ? p_massif->ZPole1(sqr(p_loop_me->IRscale())) : 0.);
}
