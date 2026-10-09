/*!
  \file Acraic_Resonance.C

  ACRAIC: RESONANCE_EXPONENTIATION - the resonant initial-final log of a
  stage's resonance, exponentiated into that resonance's propagator in every
  amplitude of the partition sum (Greco, Pancheri, Srivastava), with

    a_g = alpha/pi sum_{i outside g, j in g} w_i w_j ln(2 |p_i.p_j|).

  Pairs between two decay stages take the log with each stage's own
  off-shellness; the mixed non-factorisable terms are not reproduced (an
  ansatz). Charged resonance stages (ACRAIC: W_STAGES) are refused.
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"

#include <iomanip>

using namespace YFS;
using namespace ATOOLS;

double Acraic_Base::StageResonanceCoefficient(int g, const Vec4D_Vector &p) const
{
  if (g < 0 || g >= (int)m_stagelegs.size()) return 0.;
  double sum(0.);
  for (size_t h(0); h < m_stagelegs.size(); ++h) {
    if ((int)h == g) continue;
    for (const StageLeg &i : m_stagelegs[h]) {
      if (i.w == 0. || IsResonanceLeg(i.leg) || i.leg >= (int)p.size()) continue;
      for (const StageLeg &j : m_stagelegs[g]) {
        if (j.w == 0. || IsResonanceLeg(j.leg) || j.leg >= (int)p.size()) continue;
        sum += i.w * j.w * log(2. * std::abs(p[i.leg] * p[j.leg]));
      }
    }
  }
  return m_alpha / M_PI * sum;
}

/*
  With greco, the resonant propagator spanning each final stage g is
  multiplied by exp(I_g) in every amplitude of the partition sum, with
    I_g(q^2) = a_g ln((M~^2 - q^2)/M~^2).
*/
void Acraic_Base::BuildResonanceExponents()
{
  m_resexp.clear();
  if (!ResonanceExponentiation()) return;
  // RESONANCE_EXPONENTIATION auto: no exponent on this event
  if (WStagesActive() && m_resexpmode == acraicresexp::automatic) return;
  // ISR only: the resonant log is initial-final interference
  if (m_isronly) {
    if (m_resexpmode == acraicresexp::automatic) return;
    THROW(fatal_error, "ACRAIC: RESONANCE_EXPONENTIATION: greco exponentiates "
          "the resonant initial-final log, which YFS: MODE ISR (ISR only) "
          "leaves out. Set it to none or auto.");
  }
  if (WStagesActive())
    THROW(fatal_error, "ACRAIC: RESONANCE_EXPONENTIATION: greco is not defined "
          "for charged resonance stages, and this event activated ACRAIC: "
          "W_STAGES (auto): the W radiates itself and its terms in the "
          "resonant log are not derived. Set ACRAIC: W_STAGES off or "
          "RESONANCE_EXPONENTIATION none.");
  for (int g(0); g < m_nstages; ++g) {
    if (g == m_initstage || m_stagelegs[g].empty()) continue;
    const size_t legs(StageShiftMask(g));
    if (legs == 0) continue;
    m_resexp.push_back({g, legs, StageResonanceCoefficient(g, m_pacraic)});
  }
  if (m_redborn && !m_resexp.empty()) {
    static bool warned(false);
    if (!warned) { warned = true;
      msg_Error()<<METHOD<<"(): ACRAIC: TCHANNEL_REDUCED_BORN evaluates the "
                 <<"partition Born without propagator shifts, so it does not "
                 <<"carry the resonance factor that M_1 carries. Reported "
                 <<"once."<<std::endl; }
  }
}

void Acraic_Base::AddResonanceEntries(PropShifts &sh) const
{
  for (const ResonanceStage &r : m_resexp)
    sh.push_back(PHASIC::Process_Base::ResonanceEntry(r.legs, r.a));
}

/*
  The resonant log's share of the provider's virtual at the Born point:
  delta_I(h) = (A_B(h; 1 + I) - A_B(h))/A_B(h) per helicity, with the
  resonant propagators times the linear factor 1 + I_g, and its
  Born-weighted 2 Re.
*/
void Acraic_Base::BuildResonanceBornDelta()
{
  const int nh(Amplitude::NHel());
  m_resexpdv = 0.;
  m_resexpdh.assign(nh, Complex(0., 0.));
  if (m_resexp.empty() || p_bornproc == NULL) return;
  // The provider's Born point is m_prefsr, not m_bornmomenta. The amplitude
  // is linear in each resonant propagator, so one stage at a time the
  // difference is exactly I_g times the resonant part of the Born.
  const Vec4D_Vector &pb(m_prefsr);
  if (pb.size() != m_flavs.size()) return;
  Amplitude B;
  if (!ComixBornShifted(pb, B, PropShifts())) return;
  Amplitude dI;
  for (const ResonanceStage &r : m_resexp) {
    PropShifts lin;
    lin.push_back(PHASIC::Process_Base::ResonanceEntry
                  (r.legs, StageResonanceCoefficient(r.stage, pb), true));
    Amplitude L;
    if (!ComixBornShifted(pb, L, lin)) return;
    for (int c = 0; c < nh; ++c) dI.m_A[c] += L.m_A[c] - B.m_A[c];
  }
  double num(0.), den(0.);
  for (int c = 0; c < nh; ++c) {
    num += 2. * std::real(std::conj(B.m_A[c]) * dI.m_A[c]);
    den += std::norm(B.m_A[c]);
  }
  if (!(den > 0.)) return;
  m_resexpdv = num / den;
  // per ACRAIC helicity, in the index the provider's delta^h is mapped to
  const int flip(m_comixflip >= 0 ? ComixLegFlips(Amplitude::s_nlegs) : 0);
  for (int f = 0; f < nh; ++f) {
    const int c(FlipLegs(f, flip));
    if (std::abs(B.m_A[c]) > 1e-12 * sqrt(den))
      m_resexpdh[f] = dI.m_A[c] / B.m_A[c];
  }
}

void Acraic_Base::ResonanceProbe(int partition) const
{
  static const int np(Settings::GetMainSettings()["ACRAIC"]
                      ["RESONANCE_EXPONENTIATION_PROBE"].Get<int>());
  static int nev(0), last(-1);
  if (np <= 0 || m_resexp.empty()) return;
  if (partition == 1) { if (nev >= np) return; ++nev; last = nev; }
  if (last != nev || nev > np) return;
  std::ostringstream o;
  o<<std::setprecision(12)<<"@@@ RESEXP part="<<partition
   <<" nphot="<<m_allphotons.size()<<" dIB="<<m_resexpdv;
  for (const ResonanceStage &r : m_resexp) {
    const Vec4D q(StageSystemMomentum(r.stage) + StagePhotonSum(r.stage, -1));
    double sB(-1.);
    if (m_prefsr.size() == m_flavs.size()) {
      Vec4D qb;
      for (int l : StageExternalLegs(r.stage)) qb += m_prefsr[l];
      sB = qb.Abs2();
    }
    o<<" stage="<<r.stage<<" a="<<r.a
     <<" aB="<<StageResonanceCoefficient(r.stage, m_prefsr)
     <<" q2="<<q.Abs2()<<" sB="<<sB;
  }
  if (m_pacraic.size() >= 4)
    o<<" t="<<(m_pacraic[0] - m_pacraic[m_if1]).Abs2()
     <<" u="<<(m_pacraic[0] - m_pacraic[m_if2]).Abs2()
     <<" s="<<(m_pacraic[0] + m_pacraic[1]).Abs2();
  o<<"\n";
  std::cerr<<o.str();
}
