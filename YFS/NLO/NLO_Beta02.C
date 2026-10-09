/*!
  \file NLO_Beta02.C

  YFS: BETA0_2 / ACRAIC: BETA0_2 - the photonic two-loop beta_0^(2) of one
  event (formulas and derivation in YFS/Main/Beta0_2.H). This file reads
  the settings, builds the analytic photonic factors at the Born point of
  the loop and gives each column the delta^(2) that goes with its own
  one-loop delta.
*/

#include "YFS/NLO/NLO_Base.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"
#include "MODEL/Main/Model_Base.H"

#include <iomanip>
#include <iostream>
#include <sstream>

using namespace YFS;
using namespace ATOOLS;

namespace {

  Settings &Main() { return Settings::GetMainSettings(); }

  //! the coupling of delta_f: alpha(0) as YFS reads it, or the model's
  double Beta02Alpha(beta02alpha::code which)
  {
    if (which == beta02alpha::loop)
      return MODEL::s_model->ScalarConstant("alpha_QED");
    return 1./Main()["YFS"]["1/ALPHAQED(0)"].SetDefault(137.03599976).Get<double>();
  }

}


void NLO_Base::InitBeta02()
{
  Beta02_Config shared;
  shared.pairs = Main()["YFS"]["BETA0_2_PAIRS"].Get<beta02pairs::code>();
  shared.cross = Main()["YFS"]["BETA0_2_CROSS"].Get<bool>();
  shared.weak  = Main()["YFS"]["BETA0_2_WEAK_PRODUCT"].Get<bool>();
  shared.alpha = Beta02Alpha(Main()["YFS"]["BETA0_2_ALPHA"].Get<beta02alpha::code>());
  shared.imc   = Main()["ACRAIC"]["VIRT_HEL_IM_C"].SetDefault(0.5).Get<double>();
  m_b02acraic = m_b02nlo = shared;
  m_b02acraic.mode = Main()["ACRAIC"]["BETA0_2"].Get<beta02::code>();
  m_b02nlo.mode  = Main()["YFS"]["BETA0_2"].Get<beta02::code>();
  const dysonvp::code vp(Main()["YFS"]["DYSON_VP"].Get<dysonvp::code>());
  if (vp != dysonvp::none) {
    p_dysonvp = std::make_unique<Dyson_VP>(vp == dysonvp::complete);
    msg_Info()<<"YFS: DYSON_VP: "<<p_dysonvp->Describe()<<std::endl;
  }
  m_ewhohel = Main()["ACRAIC"]["EW_HIGHER_ORDER_HELICITY"].Get<ewhohel::code>();
  m_m1sq = Main()["YFS"]["M1_SQUARE"].Get<m1square::code>();
  if (!(m_b02acraic.On() || m_b02nlo.On())) return;
  /*
    delta_f is the photonic remainder in the B convention of the
    subtraction actually used (Beta0_2.H). With Dim_Reg false the s-channel
    dipoles are subtracted with BVirtGeneral, the one-in-one-out formula of
    YFS_Form_Factor::A2 applied to two incoming or two outgoing legs: it is
    (alpha/pi)(L^2/2 + L/2 + 1 + pi^2/2) below KKMC's SBvirt per s-channel
    pair (310 alpha/pi for the Z-pole beams), so the one-loop delta itself is
    wrong there and no delta^(2) can be built on it.
  */
  if (!m_dim_reg)
    THROW(fatal_error, "YFS/ACRAIC: BETA0_2 needs YFS: Dim_Reg: true. The"
          " photon-mass subtraction (BVirtGeneral) of an s-channel dipole is"
          " not KKMC's B, see NLO_Beta02.C.");
  msg_Info()<<"YFS: BETA0_2 (YFS.NLO): "<<Beta02Describe(m_b02nlo)<<"\n"
            <<"ACRAIC: BETA0_2 (ACRAIC): "<<Beta02Describe(m_b02acraic)<<std::endl;
}


void NLO_Base::BuildBeta02()
{
  m_b02 = Photonic_Factors();
  if (!(m_b02acraic.On() || m_b02nlo.On()) || m_born == 0. || IsBad(m_oneloop))
    return;
  m_b02 = PhotonicFactors(m_flavs, m_plab, m_b02acraic, 2);
  static const bool boxes(Main()["YFS"]["BETA0_2_BOXES"].Get<bool>());
  if (boxes && m_b02.ok && (m_b02acraic.mode == beta02::two_loop ||
                            m_b02nlo.mode == beta02::two_loop))
    AddPhotonicTwoLoopBoxes(m_b02, m_flavs, m_plab, m_b02acraic.alpha);
  // e+e- -> gamma gamma: the AGT two-loop hard function on top of two_loop
  static const bool diphoton(Main()["YFS"]["BETA0_2_DIPHOTON"].Get<bool>());
  double daa(0.);
  static long aawarn(0);
  if (diphoton && m_b02.ok && (m_b02acraic.mode == beta02::two_loop ||
                               m_b02nlo.mode == beta02::two_loop)) {
    // the hard function comes from the double-virtual provider
    // (NLO_Part with M, VV_Generator: Internal)
    PHASIC::Virtual_ME2_Base *vv(p_vv != NULL ? p_vv->p_loop_me.get() : NULL);
    static const double powercut(Main()["YFS"]["BETA0_2_DIPHOTON_POWER_CUT"]
                                 .SetDefault(1e-3).Get<double>());
    bool provided(false), skipped(false);
    daa = AddPhotonicTwoLoopDiphoton(m_b02, vv, m_plab, m_b02acraic.alpha,
                                     m_flavs[0].Mass(), powercut, provided, skipped);
    if (provided) { ++m_b02aa_n; if (skipped) ++m_b02aa_skip; }
    if (daa != 0. && !(m_b02acraic.weak && m_b02nlo.weak) && aawarn++ == 0)
      msg_Error()<<METHOD<<"(): YFS: BETA0_2_DIPHOTON assumes"
                 <<" BETA0_2_WEAK_PRODUCT true; the photonic delta_f w of"
                 <<" e+e- -> gamma gamma is missing.\n";
    static long novv(0);
    if (!provided && m_flavs.size() == 4 && m_flavs[2].IsPhoton()
        && m_flavs[3].IsPhoton() && novv++ == 0)
      msg_Error()<<METHOD<<"(): e+e- -> gamma gamma with BETA0_2 two_loop but"
                 <<" no VV_Generator: Internal double virtual (NLO_Part with M):"
                 <<" the two-loop hard function (BETA0_2_DIPHOTON) is not added.\n";
  }
  if (!m_b02.ok) {
    if (m_b02_fail++ == 0)
      msg_Error()<<METHOD<<"(): a charged leg is massless (or nothing is"
                 <<" charged); beta_0^(2) needs Massive: true on every charged"
                 <<" leg. It is 0 for such events.\n";
    return;
  }
  ++m_b02_n;
  static const int probe(Main()["YFS"]["BETA0_2_PROBE"].Get<int>());
  static long nprobe(0);
  if (probe > 0 && nprobe < probe) {
    ++nprobe;
    const double v1(VOneLoop());
    std::ostringstream o;
    o<<std::setprecision(12)<<"@@@ B02 v1="<<v1
     <<" deltaf="<<m_b02.deltaf.real()<<","<<m_b02.deltaf.imag()
     <<" nll="<<m_b02.nll.real()<<","<<m_b02.nll.imag()
     <<" d2ext="<<Beta02External()<<" v2nlo="<<Beta02SquaredLevel()
     <<" npairs="<<m_b02.f.size()<<" daa="<<daa;
    for (size_t d(0); d < m_b02.f.size(); ++d)
      o<<" f"<<d<<"="<<m_b02.f[d].real()<<","<<m_b02.f[d].imag();
    o<<"\n";
    std::cerr<<o.str();
  }
}


void NLO_Base::BuildDysonVP()
{
  m_vp = VP_Factors();
  if (!p_dysonvp || m_born == 0. || IsBad(m_oneloop)) return;
  m_vp = p_dysonvp->Factors(m_flavs, m_plab);
  if (!m_vp.ok) {
    if (m_vp_fail++ == 0)
      msg_Error()<<METHOD<<"(): YFS: DYSON_VP dresses one neutral-current"
                 <<" e+e- -> f fbar Born only; this event is left undressed.\n";
    return;
  }
  ++m_vp_n;
  static const int probe(Main()["YFS"]["DYSON_VP_PROBE"].Get<int>());
  static long nprobe(0);
  if (probe > 0 && nprobe < probe) {
    ++nprobe;
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ DYSONVP s="<<m_vp.s<<" t="<<m_vp.t
     <<" Ds="<<m_vp.Ds.real()<<","<<m_vp.Ds.imag()
     <<" Dt="<<m_vp.Dt.real()<<","<<m_vp.Dt.imag()
     <<" v2="<<m_vp.v2<<" rbar="<<m_vp.rbar.real()<<","<<m_vp.rbar.imag();
    for (const VP_Helicity &h : m_vp.hel)
      o<<" ["<<h.l[0]<<h.l[1]<<h.l[2]<<h.l[3]<<" w="<<h.w
       <<" r="<<h.r.real()<<","<<h.r.imag()<<"]";
    o<<"\n";
    std::cerr<<o.str();
  }
}


double NLO_Base::VOneLoop() const
{
  if (m_born == 0.) return 0.;
  return m_oneloop/m_born - m_vewho;
}


double NLO_Base::Beta02External() const
{
  if (!m_b02.ok || !m_b02acraic.On()) return 0.;
  return Beta02DeltaReal(m_b02, 0.5*VOneLoop(), m_b02acraic);
}


double NLO_Base::Beta02SquaredLevel() const
{
  if (!m_b02.ok || !m_b02nlo.On()) return 0.;
  const double v1(VOneLoop());
  /*
    YFS: M1_SQUARE helicity: the Born-weighted sum_h |1 + delta_h + d2_h|^2
    - 1 - v1 over the provider's helicity factors (complex, the same phase
    convention as delta_f, so the sum is c independent at O(alpha^2)); the
    weighted 2 Re delta_h is v1 by construction (BuildLoopHelicityFactors).
  */
  if (m_m1sq == m1square::helicity && m_lhelok && m_la0.size() == m_ldelta.size()) {
    double s0(0.), s2(0.);
    for (size_t c(0); c < m_la0.size(); ++c) {
      const double w(std::norm(m_la0[c]));
      Complex d2(Beta02Delta(m_b02, m_ldelta[c], m_b02nlo));
      if (m_b02nlo.mode == beta02::two_loop && c < m_lhel.size())
        d2 += Beta02BoxFor(m_b02, m_lhel[c]);
      s0 += w;
      s2 += w*std::norm(1. + m_ldelta[c] + d2);
    }
    if (s0 > 0.) { ++m_m1sq_hel; return s2/s0 - 1. - v1; }
  }
  ++m_m1sq_avg;
  const double d2(Beta02DeltaReal(m_b02, 0.5*v1, m_b02nlo));
  return sqr(1. + 0.5*v1 + d2) - 1. - v1;
}


void NLO_Base::BuildEWHigherOrderHelicity()
{
  m_ghelok = false;
  m_ghel.clear(); m_gdelta.clear();
  if (m_ewhohel != ewhohel::helicity || !p_virt || !p_virt->HasEWHigherOrder()
      || m_born == 0.) return;
  std::vector<Complex> a0, a1;
  if (!p_virt->p_ewho_me->HelicityAmplitudes(a0, a1, m_ghel) || a0.empty() ||
      a0.size() != a1.size()) {
    if (m_ghel_fail++ == 0)
      msg_Error()<<METHOD<<"(): the EW_HIGHER_ORDER provider returned no"
                 <<" helicity factors; v_HO/2 on every helicity instead.\n";
    return;
  }
  // delta_G(h) = a1/a0, rescaled so that sum_h |a0|^2 2 Re delta_G / sum |a0|^2
  // is v_HO, the term YFS.NLO and the external column carry
  double s0(0.), s1(0.);
  for (size_t c(0); c < a0.size(); ++c) {
    s0 += std::norm(a0[c]);
    s1 += 2.*std::real(std::conj(a0[c])*a1[c]);
  }
  if (!(s0 > 0.) || s1 == 0.) { ++m_ghel_fail; return; }
  const double kappa(m_vewho/(s1/s0));
  m_gdelta.resize(a0.size());
  for (size_t c(0); c < a0.size(); ++c)
    m_gdelta[c] = std::abs(a0[c]) > 0. ? kappa*a1[c]/a0[c] : Complex(0.5*m_vewho, 0.);
  m_ghelok = true;
  ++m_ghel_n;
  static const int probe(Main()["YFS"]["DYSON_VP_PROBE"].Get<int>());
  static long nprobe(0);
  if (probe > 0 && nprobe < probe) {
    ++nprobe;
    std::ostringstream o;
    // GRIFFIN's NLO helicity sum in Griffin_Virtual's normalisation (s^2/pi)
    // against the Born the virtual is divided by
    const double sh((m_plab[0] + m_plab[1]).Abs2());
    o<<std::setprecision(10)<<"@@@ GRIFFINHEL vho="<<m_vewho<<" kappa="<<kappa
     <<" nloG_over_born="<<0.25*s0*sh*sh/M_PI/m_rescale_alpha/m_born;
    for (size_t c(0); c < a0.size(); ++c)
      o<<" ["<<m_ghel[c][0]<<m_ghel[c][2]<<" w="<<std::norm(a0[c])/s0
       <<" dG="<<m_gdelta[c].real()<<"]";
    o<<"\n";
    std::cerr<<o.str();
  }
}


void NLO_Base::ReportBeta02()
{
  if (p_dysonvp) {
    double vb[2] = {(double)m_vp_n, (double)m_vp_fail};
#ifdef USING__MPI
    if (mpi->Size() > 1) mpi->Allreduce(vb, 2, MPI_DOUBLE, MPI_SUM);
#endif
    msg_Out()<<"DYSON_VP: dressed "<<(long)vb[0]<<" events, "<<(long)vb[1]
             <<" not a neutral-current 2 -> 2 Born.\n";
  }
  if (m_ewhohel == ewhohel::helicity || m_m1sq == m1square::helicity) {
    double hb[4] = {(double)m_ghel_n, (double)m_ghel_fail, (double)m_m1sq_hel,
                    (double)m_m1sq_avg};
#ifdef USING__MPI
    if (mpi->Size() > 1) mpi->Allreduce(hb, 4, MPI_DOUBLE, MPI_SUM);
#endif
    msg_Out()<<"EW_HIGHER_ORDER_HELICITY: "<<(long)hb[0]<<" events per helicity, "
             <<(long)hb[1]<<" without; M1_SQUARE: "<<(long)hb[2]<<" helicity, "
             <<(long)hb[3]<<" averaged.\n";
  }
  if (!(m_b02acraic.On() || m_b02nlo.On())) return;
  double buf[2] = {(double)m_b02_n, (double)m_b02_fail};
#ifdef USING__MPI
  if (mpi->Size() > 1) mpi->Allreduce(buf, 2, MPI_DOUBLE, MPI_SUM);
#endif
  msg_Out()<<"BETA0_2: photonic factors built on "<<(long)buf[0]
           <<" events, "<<(long)buf[1]<<" without (massless charged leg).\n";
}
