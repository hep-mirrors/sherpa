#include "YFS/Main/YFS_Weight_Builder.H"
#include "YFS/Main/YFS_Handler.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/Beta0_2.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;

// YFS: DIFF_BUDGET_STREAMS (diagnostic): the named columns BuildNamedWeights
// adds, in this order (see the comment there)
static const char *const s_diffbudgetnames[] = {
  "DB_CRUDE", "DB_A1", "DB_A0", "DB_V", "DB_X", "DB_R", "DB_RIF",
  "DB_X1", "DB_X2", "DB_NPH", "DB_SXI", "DB_SXF",
  "DB_VSQ", "DB_REST", "DB_NPV"};

static bool DiffBudgetStreamsOn()
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["DIFF_BUDGET_STREAMS"].Get<bool>());
  return on;
}

namespace {
  //! ACRAIC: WEIGHT_O0 - the O(alpha^0) ACRAIC column (rho_0/rho_crude) as a
  //! named weight
  bool AcraicO0Stream()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["WEIGHT_O0"].SetDefault(false).Get<bool>());
    return on;
  }
  //! ACRAIC: FIXED_ORDER_STREAM - the strict O(alpha) truncation of ACRAIC as
  //! the named weight ACRAIC_FixedOrder (see GenerateWeight)
  bool AcraicStrictO1Stream()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["STRICT_O1_STREAM"].SetDefault(false).Get<bool>());
    return on;
  }
  bool AcraicFixedOrderStream()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["FIXED_ORDER_STREAM"].SetDefault(false).Get<bool>());
    return on;
  }
  //! ACRAIC: TNP_BAND_STREAMS - the TNP band as shifted named weights
  bool AcraicTnpBandStreams()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["TNP_BAND_STREAMS"].SetDefault(false).Get<bool>());
    return on;
  }
  //! ACRAIC: TNP_BAND_INCOMPLETE - the incomplete-second-order parameters
  //! RV_ISR and B0_2L (the real-virtual beyond REAL_VIRTUAL factorisable, the
  //! non-logarithmic two-loop vertex) in the O(alpha^2) band
  bool AcraicTnpBandIncomplete()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["TNP_BAND_INCOMPLETE"].SetDefault(true).Get<bool>());
    return on;
  }
  /*
    ACRAIC: TNP_BAND_STREAMS. One theory nuisance parameter: its norm N and
    the template classes it scales (the ACRAIC_TNP1_/ACRAIC_TNP2_ columns of
    TNP_STREAMS), for the initial state. The named weights
    ACRAIC_TNPB_<name>_UP/_DN are nominal +- N x the summed templates, so
    rivet-mkhtml's BandComponentEnv (one component per parameter, added in
    quadrature) draws sigma +- sqrt(sum (N T)^2).
  */
  struct TnpBandParam {
    std::string name;
    double norm;
    std::vector<int> classes;
  };
  /*
    The initial-state parameters of a run truncated at O(alpha^k), with
    g = (alpha(0)/pi)(ln(s/m_e^2) - 1). Halved: T_c responds to an amplitude
    rescaling, which moves |B_c|^2 twice as fast. beta02: ACRAIC: BETA0_2 is
    in the prediction, so beta_0 first misses order k+1 at k = 2 too.
  */
  std::vector<TnpBandParam> TnpBandParams(double s, double me, int k,
                                          bool beta02, bool incomplete)
  {
    typedef YFS::Acraic_Base A;
    const double a(1./137.035999/M_PI), g(a*(log(s/(me*me)) - 1.));
    const int mb(k == 1 || beta02 ? k + 1 : k);
    auto fact = [](int n) { double f(1.); for (int i(2); i <= n; ++i) f *= i; return f; };
    std::vector<TnpBandParam> b;
    b.push_back({"B0_ISR", 0.5*pow(g, mb)/fact(mb), {A::tnp_b0}});
    b.push_back({"ISR",    0.5*pow(g, k)/fact(k),   {A::tnp_isr}});
    b.push_back({"MULTI",  0.5*pow(g, k)/fact(k),   {A::tnp_multi}});
    if (k == 2) {
      // the two-photon classes share one parameter
      b.push_back({"B2", 0.5*g, {A::tnp_b2i, A::tnp_b2f, A::tnp_b2m}});
      if (incomplete) {
        b.push_back({"RV_ISR", 0.5*a,                {A::tnp_isr}});
        b.push_back({"B0_2L",  0.5*a*a*M_PI*M_PI,    {A::tnp_b0}});
      }
    }
    return b;
  }
}

void YFS_Weight_Builder::CalculateBeta() {
  m_h.m_ev.m_nlo_current = false;
  m_h.m_ev.m_nlo_real = m_h.m_ev.m_nlo_virtual = m_h.m_ev.m_nlo_rv = m_h.m_ev.m_nlo_rr = 0.;
  m_h.m_ev.m_nlo_vxr = 0.;
  m_h.m_ev.m_nlo_real_ifisoft = 0.;
  if(!m_h.m_rmode && !m_h.m_int_nlo) return;
  double realISR(0), realFSR(0);
  if (m_h.m_betaorder > 0) {
    if(m_h.m_real_only) {
      if(!m_h.m_no_born) m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEX()+1;
      else m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEX();
    }
    else if(m_h.m_virtual_only) {
      if(!m_h.m_no_born) m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEXVirtual();
      else m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEXVirtual()-1;
    }
    else {
      if(!m_h.m_no_born) m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEX()+m_h.p_dipoles->CalculateEEXVirtual();
      else m_h.m_ev.m_real = m_h.p_dipoles->CalculateEEX()/m_h.m_born+m_h.p_dipoles->CalculateEEXVirtual()/m_h.m_born;
    }
    m_h.m_ev.m_eex = m_h.m_ev.m_real;
    if(IsNan(m_h.m_ev.m_eex)) m_h.m_ev.m_eex=0;
  }
  double acraicfac(1.);
  bool   haveacraic(false);
  if (m_h.m_useacraic) {
    m_h.MakeACRAIC();
    /*
      YFS: REAL_BORN_PHOTON_MULTICHANNEL, the ACRAIC half: rho_crude becomes the
      sum of the crudes of every assignment of photons to Born pair and ISR
      that the generator could have produced this final state with
      (NLO_Base::BornPhotonChannelSum, all Born-pair choices). on applies
      it, report only computes it. G = 1 exactly without
      Born photons or when no ISR photon passes the Born cuts.
    */
    m_h.m_bpmc_evtG = 1.; m_h.m_bpmc_evtn = 0;
    static const bornphotonmc::code bpmc(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(bornphotonmc::on)
      .Get<bornphotonmc::code>());
    if (bpmc != bornphotonmc::off && m_h.p_acraic && m_h.p_nlo && m_h.p_dipoles->HasDipoleII()
        && !m_h.m_ev.m_ISRPhotons.empty()) {
      InitNLO();
      const std::vector<Vec4D> isr(m_h.m_ev.m_ISRPhotons.begin(), m_h.m_ev.m_ISRPhotons.end());
      m_h.m_bpmc_evtG = m_h.p_nlo->BornPhotonChannelSum(m_h.m_ev.m_reallab, isr,
                                                m_h.p_dipoles->GetDipoleII(), -1,
                                                &m_h.m_bpmc_evtn);
      if (bpmc == bornphotonmc::on && m_h.m_bpmc_evtG != 1.)
        m_h.p_acraic->ScaleRhoCrude(m_h.m_bpmc_evtG);
    }
    if (m_h.p_acraic) {
      /*
        Denominator = the INCOHERENT partition sum.

        rho0 is the coherent sum, so the partitions can cancel and it can come
        arbitrarily close to zero, while the crude weight this factor
        multiplies is a factorised eikonal times Born and does not.
      */
      const double r0(m_h.p_acraic->GetRhoCrude()), r1(m_h.p_acraic->GetResult());
      if (r0 > 0. && !IsBad(r1/r0)) { acraicfac = r1/r0; haveacraic = true; }
      else ++m_h.m_acraicstats.m_bad;
      m_h.m_ev.m_acraicfactor = haveacraic ? acraicfac : 0.;
      if (m_h.p_nlo) m_h.p_nlo->SetAcraicVirtual(m_h.p_acraic->VirtualFactor());
    }
  }

  if(m_h.m_nlotype!=nlo_type::born && m_h.p_nlo) {
    if(m_h.m_no_born) m_h.m_ev.m_real=CalculateNLO()/m_h.m_born;
    else m_h.m_ev.m_real=(m_h.m_born+CalculateNLO())/m_h.m_born;
    m_h.m_ev.m_nlo_current = true;
    if (m_h.m_acraic_compare && haveacraic) AcraicCompare();
  }

  if (haveacraic) {
    m_h.m_acraicstats.AddOverEex(m_h.m_ev.m_real != 0. ? acraicfac/m_h.m_ev.m_real : 0.);
  }
}
/*!
  ACRAIC against the EXTERNAL providers, event by event.

  Both sides are O(alpha) corrections relative to the Born, so they are
  directly comparable without any normalisation being matched by hand.

      external virtual / Born   <->   rho(Born+virtual)/rho(Born) - 1
      external real    / Born   <->   rho(Born+real)   /rho(Born) - 1

  The total row is the meaningful one: both sides are the O(alpha) correction
  factor to the same resummed Born, and it is what drives the weight.

*/
void YFS_Weight_Builder::AcraicCompare() {
  if (!m_h.p_acraic || !m_h.p_nlo || m_h.m_born == 0.) return;
  const double rc(m_h.p_acraic->RealFactor() - 1.);
  const double tc(m_h.m_ev.m_acraicfactor - 1.);
  const double vc(tc - rc);
  const double ve(m_h.m_ev.m_nlo_virtual/m_h.m_born);
  const double re(m_h.m_ev.m_nlo_real   /m_h.m_born);
  const double te(ve + re);
  auto rel = [](double a, double b) {
    const double s(std::abs(a)+std::abs(b));
    return s > 0. ? std::abs(a-b)/s : 0.;
  };
  const double dv(rel(vc,ve)), dr(rel(rc,re)), dt(rel(tc,te));
  m_h.m_acraicstats.AddCompare(dv, dr, dt);
  const bool dumped(m_h.p_acraic->JustDumped());
  if (m_h.m_acraicstats.m_cmp_n <= (long)m_h.m_acraic_compare || dumped) {
    size_t ng(m_h.m_ev.m_ISRPhotons.size() + m_h.m_ev.m_FSRPhotons.size());
    msg_Out()<<std::setprecision(8)
             <<"@@@ ACRAICCMP"<<(dumped?"-DUMPED":"")<<" n="<<m_h.m_acraicstats.m_cmp_n<<" ngam="<<ng
             <<" born="<<m_h.m_born<<"\n"
             <<"    virtual  acraic="<<vc<<"  ext="<<ve<<"  reldiff="<<dv<<"\n"
             <<"    real*    acraic="<<rc<<"  ext="<<re<<"  reldiff="<<dr<<"\n"
             <<"    total    acraic="<<tc<<"  ext="<<te<<"  reldiff="<<dt<<"\n"
             <<"      virt pieces: raw/born="<<(m_h.m_born!=0.?m_h.p_nlo->m_virt_raw/m_h.m_born:0.)
             <<"  sub/born="<<(m_h.m_born!=0.?m_h.p_nlo->m_virt_subval/m_h.m_born:0.)
             <<"  raw-sub="<<(m_h.m_born!=0.?(m_h.p_nlo->m_virt_raw-m_h.p_nlo->m_virt_subval)/m_h.m_born:0.)
             <<"\n      formfactor="<<m_h.m_ev.m_formfactor<<"  log(FF)="<<(m_h.m_ev.m_formfactor>0.?log(m_h.m_ev.m_formfactor):0.)
             <<"\n";
    // The Born point of the virtual: cos theta of the first final-state leg
    // in the beam CM frame, and the virtual subtraction split by dipole type
    // (taken in NLO_Base::CalculateVirtual; the default, non-eps, s-channel
    // branch only).
    if (m_h.p_nlo->m_plab.size() > 2) {
      const Vec4D P(m_h.p_nlo->m_plab[0] + m_h.p_nlo->m_plab[1]);
      Poincare cms(P);
      Vec4D pa(m_h.p_nlo->m_plab[0]), pf(m_h.p_nlo->m_plab[2]);
      cms.Boost(pa); cms.Boost(pf);
      double xmax(0.);
      // hardest photon: 1 - |cos| to the beam axis and 1 - cos to the nearest
      // final-state charged leg, in the beam CM
      double th1c(2.), thf(2.);
      auto hardest = [&](const Vec4D &k) {
        if (!(2.*k.E()/P.Mass() > xmax)) return;
        xmax = 2.*k.E()/P.Mass();
        Vec4D kc(k); cms.Boost(kc); th1c = 1. - dabs(kc.CosTheta());
        thf = 2.;
        for (size_t i(2); i < m_h.m_ev.m_plab.size(); ++i) {
          if (m_h.m_flavs.size() > i && m_h.m_flavs[i].Charge() == 0.) continue;
          Vec4D pc(m_h.m_ev.m_plab[i]); cms.Boost(pc);
          thf = Min(thf, 1. - kc.CosTheta(pc));
        }
      };
      for (const Vec4D &k : m_h.m_ev.m_ISRPhotons) hardest(k);
      for (const Vec4D &k : m_h.m_ev.m_FSRPhotons) hardest(k);
      msg_Out()<<"      @@@ VPT cos="<<pf.CosTheta(pa)
               <<" sqrts="<<P.Mass()<<" raw/born="<<m_h.p_nlo->m_virt_raw/m_h.m_born
               <<" sub/born="<<m_h.p_nlo->m_virt_subval/m_h.m_born
               <<" subII="<<m_h.p_nlo->m_virt_sub_ii/m_h.m_born<<" subFF="<<m_h.p_nlo->m_virt_sub_ff/m_h.m_born
               <<" subIF="<<m_h.p_nlo->m_virt_sub_if/m_h.m_born
               <<" v_used="<<ve<<" dIB="<<m_h.p_acraic->ResonanceBornDeltaV()<<" ngam="<<ng<<" xmax="<<xmax<<" th1c="<<th1c<<" thf="<<thf
               <<" nISR="<<m_h.m_ev.m_ISRPhotons.size()<<" nFSR="<<m_h.m_ev.m_FSRPhotons.size()
               <<" xisr="<<[&]{ double x(0.); for (const Vec4D &k : m_h.m_ev.m_ISRPhotons) x = Max(x, 2.*k.E()/P.Mass()); return x; }()
               <<" xfsr="<<[&]{ double x(0.); for (const Vec4D &k : m_h.m_ev.m_FSRPhotons) x = Max(x, 2.*k.E()/P.Mass()); return x; }()
               <<" mff="<<(m_h.m_ev.m_plab.size() > 3 ? (m_h.m_ev.m_plab[2] + m_h.m_ev.m_plab[3]).Mass() : -1.)
               // the pair plus one photon (any label), the mass closest to M_Z: an
               // ISR-labelled photon that would be resonant FSR off the Z return
               <<" mffk="<<[&]{ double best(-1.), mz(Flavour(kf_Z).Mass());
                   if (m_h.m_ev.m_plab.size() < 4) return best;
                   const Vec4D q(m_h.m_ev.m_plab[2] + m_h.m_ev.m_plab[3]);
                   for (const auto *ks : {&m_h.m_ev.m_ISRPhotons, &m_h.m_ev.m_FSRPhotons})
                     for (const Vec4D &k : *ks) {
                       const double mk((q + k).Mass());
                       if (best < 0. || std::abs(mk - mz) < std::abs(best - mz)) best = mk;
                     }
                   return best; }()
               <<" acraicfac="<<m_h.m_ev.m_acraicfactor<<" eexfac="<<m_h.m_ev.m_real
               <<" re="<<re<<" rho0crude="<<(m_h.p_acraic->GetRhoCrude() > 0. ?
                                            m_h.p_acraic->GetResult0()/m_h.p_acraic->GetRhoCrude() : 0.)
               <<" rho01crude="<<(m_h.p_acraic->GetRhoCrude() > 0. ?
                                  m_h.p_acraic->GetResult01()/m_h.p_acraic->GetRhoCrude() : 0.)
               <<" rho1="<<m_h.p_acraic->GetResult()<<" rho0="<<m_h.p_acraic->GetResult0()
               <<" rhocrude="<<m_h.p_acraic->GetRhoCrude()<<" born="<<m_h.m_born
               <<" y_r="<<m_h.p_nlo->m_dbg_r<<" y_flux="<<m_h.p_nlo->m_dbg_flux
               <<" y_subloc="<<m_h.p_nlo->m_dbg_subloc<<" y_subb="<<m_h.p_nlo->m_dbg_subb
               <<" y_tot="<<m_h.p_nlo->m_dbg_tot<<" y_dk="<<m_h.p_nlo->m_dbg_dk
               <<" y_dp="<<m_h.p_nlo->m_dbg_dp
               <<" x2pt="<<m_h.p_nlo->m_dbg_x2pt<<" x2gen="<<m_h.p_nlo->m_dbg_x2gen
               <<" x2isr="<<m_h.p_nlo->m_dbg_x2isr<<" facisr="<<m_h.p_nlo->m_dbg_facisr
               <<" facfsr="<<m_h.p_nlo->m_dbg_facfsr<<"\n";
      m_h.p_nlo->m_dbg_x2pt = m_h.p_nlo->m_dbg_x2gen = m_h.p_nlo->m_dbg_x2isr = 0.;
      m_h.p_nlo->m_dbg_facisr = m_h.p_nlo->m_dbg_facfsr = 0.;
      // one-photon events: the momenta, for an offline soft-factor check
      if (ng == 1) {
        const Vec4D k(m_h.m_ev.m_ISRPhotons.size() ? m_h.m_ev.m_ISRPhotons[0] : m_h.m_ev.m_FSRPhotons[0]);
        msg_Out()<<std::setprecision(16)<<"      @@@ VPTMOM k="<<k[0]<<","<<k[1]<<","<<k[2]<<","<<k[3];
        for (size_t i(0); i < m_h.m_ev.m_plab.size(); ++i)
          msg_Out()<<" lab"<<i<<"="<<m_h.m_ev.m_plab[i][0]<<","<<m_h.m_ev.m_plab[i][1]<<","<<m_h.m_ev.m_plab[i][2]<<","<<m_h.m_ev.m_plab[i][3];
        for (size_t i(0); i < m_h.p_nlo->m_plab.size(); ++i)
          msg_Out()<<" nlo"<<i<<"="<<m_h.p_nlo->m_plab[i][0]<<","<<m_h.p_nlo->m_plab[i][1]<<","<<m_h.p_nlo->m_plab[i][2]<<","<<m_h.p_nlo->m_plab[i][3];
        msg_Out()<<std::setprecision(8)<<"\n";
      }
    }
    /*
      Photon by photon. The external stores its own per-photon real in
      YFS::Photon::m_beta10 (NLO_Base::CalculateReal), ACRAIC in m_realphot;
      matched by MOMENTUM rather than by index, because the two lists are
      built independently and an ordering assumption would silently pair the
      wrong emissions.

      Neither column sums to the total: rho is |B + sum_j R_j|^2, so the cross
      terms between photons belong to no single photon.
    */
    const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
    const std::vector<std::pair<Vec4D,double> > &eph(m_h.m_ev.m_extrealphot);
    /*
      The two sides are NOT handed the same photons. ACRAIC gets m_h.m_ev.m_ISRPhotons +
      m_h.m_ev.m_FSRPhotons; m_h.p_nlo->m_photons is m_h.m_ev.m_ISRPhotons + m_h.m_ev.m_me_photons, and the
      final-state halves come from different places - Dipole::GetPhotons() for
      one, GetMEPhotons() for the other. Report both counts, because if they
      differ there is no per-photon comparison to be made and the real
      correction is not even over the same emissions.
    */
    msg_Out()<<"      photons: acraic="<<cph.size()<<"  ext="<<eph.size()
             <<(cph.size()!=eph.size() ? "   <-- DIFFERENT LISTS" : "")<<"\n";
    // Both lists in full, so a photon present on one side and absent on the
    // other is visible directly rather than inferred from a failed match.
    msg_Out()<<"        acraic list:";
    for (size_t a(0); a < cph.size(); ++a) msg_Out()<<" "<<cph[a].E();
    msg_Out()<<"\n        ext  list:";
    for (size_t b(0); b < eph.size(); ++b) msg_Out()<<" "<<eph[b].first.E();
    msg_Out()<<"\n        nISR="<<m_h.m_ev.m_ISRPhotons.size()
             <<" nFSR="<<m_h.m_ev.m_FSRPhotons.size()
             <<" nFSRforME="<<m_h.m_ev.m_fsrphotonsforME.size()<<"\n";
    /*
      Which list belongs to THIS event? Momentum conservation decides it with
      no reference to either matrix element:  p_a + p_b - q_c - q_d - sum k
      must vanish. A list captured in the wrong frame, or at the wrong stage
      of the dipole's boost, cannot balance.
    */
    auto bal = [&](const Vec4D &sum) {
      Vec4D b(m_h.m_ev.m_plab[0] + m_h.m_ev.m_plab[1] - sum);
      for (size_t i(2); i < m_h.m_ev.m_plab.size(); ++i) b -= m_h.m_ev.m_plab[i];
      return Max(Max(dabs(b[0]),dabs(b[1])),Max(dabs(b[2]),dabs(b[3])));
    };
    Vec4D sc, se;
    for (size_t a(0); a < cph.size(); ++a) sc += cph[a];
    for (size_t b(0); b < eph.size(); ++b) se += eph[b].first;
    msg_Out()<<"        momentum balance:  acraic list "<<bal(sc)
             <<"   ext list "<<bal(se)<<"\n";
    for (size_t a(0); a < cph.size(); ++a) {
      long match(-1); double best(1e30);
      for (size_t b(0); b < eph.size(); ++b) {
        const Vec4D d(cph[a]-eph[b].first);
        const double m(Max(Max(dabs(d[0]),dabs(d[1])),Max(dabs(d[2]),dabs(d[3]))));
        if (m < best) { best = m; match = (long)b; }
      }
      const bool ok(match >= 0 && best < 1e-6*Max(cph[a].E(),1e-30));
      const double rce(m_h.p_acraic->RealFactorPhoton(a));
      const double rex(ok && m_h.m_born != 0. ? eph[match].second/m_h.m_born : 0.);
      msg_Out()<<"      gam["<<a<<"] E="<<cph[a].E()
               <<(a < m_h.m_ev.m_ISRPhotons.size() ? " ISR" : " FSR")
               <<"  acraic="<<rce
               <<(ok ? "  ext=" : "  ext=(no match) ")<<rex
               <<"  reldiff="<<(ok && std::abs(rce)+std::abs(rex) > 0. ?
                                std::abs(rce-rex)/(std::abs(rce)+std::abs(rex)) : -1.)
               <<"\n";
    }
  }
}
void YFS_Weight_Builder::InitNLO(){
  m_h.p_nlo->Init(m_h.m_flavs,m_h.m_ev.m_reallab,m_h.m_ev.m_bornMomenta);
  m_h.p_nlo->m_bpmc_vmax = m_h.m_vmax; m_h.p_nlo->m_bpmc_s = m_h.m_s;
  m_h.p_nlo->SetPostEmissionMomenta(m_h.m_ev.m_plab);
  m_h.p_nlo->p_dipoles = m_h.p_dipoles.get();
  m_h.p_nlo->SetBorn(m_h.m_born);
  m_h.p_nlo->SetFSR(m_h.p_fsr.get());
  m_h.p_nlo->m_ISRPhotons = m_h.m_ev.m_ISRPhotons;
  if (m_h.m_nlo_fsr_photons)
    m_h.p_nlo->m_FSRPhotons = m_h.m_nlo_fsr_from_event ? m_h.m_ev.m_FSRPhotons : m_h.m_ev.m_fsrphotonsforME;
  else
    m_h.p_nlo->m_FSRPhotons.clear();
  // Mirror the two lines above: same photons, now carrying their dipole.
  m_h.p_nlo->m_photons.clear();
  if (m_h.p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_h.m_ev.m_ISRPhotons)
      m_h.p_nlo->m_photons.push_back(YFS::Photon(k, &m_h.p_dipoles->GetDipoleII()));
  if (m_h.m_nlo_fsr_photons) {
    // m_ev.m_photons holds the event-record photons WITH their dipoles: ISR first,
    // then the FF ones, so skip the ISR head to take only the final-state tail.
    if (m_h.m_nlo_fsr_from_event) {
      const size_t skip(m_h.p_dipoles->HasDipoleII() ? m_h.m_ev.m_ISRPhotons.size() : 0);
      for (size_t i(skip); i < m_h.m_ev.m_photons.size(); ++i)
        m_h.p_nlo->m_photons.push_back(m_h.m_ev.m_photons[i]);
    } else {
      for (const YFS::Photon &k : m_h.m_ev.m_me_photons) m_h.p_nlo->m_photons.push_back(k);
    }
  }
}
/*
  The ACRAIC photon list (AllPhotonsLab, initial-state photons first) and
  YFS.NLO's are built independently, so each ACRAIC photon takes the dv of
  the NLO photon with the same lab momentum (to 1e-6 of its energy). An
  initial-state photon's v_B is the provider's at the post-ISR Born point,
  while ACRAIC's v has its initial vertex moved to s (ISR_VERTEX_AT_S), so
  its dv drops that shift. The shift goes to the beams of the photon's own
  (n+1)-point, s_pt (m_rvdv_spt), which lies below s when other
  initial-state photons were taken out of the beams: the provider's v_{n+1}
  has its initial vertex there. Returns the number of photons matched; entry
  (if given) gets each ACRAIC photon's m_rvdv index, -1 without a match.
*/
int YFS_Weight_Builder::MatchRealVirtualDV(std::vector<double> &dv, std::vector<long> *entry)
{
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  const std::vector<double> &spt(m_h.p_nlo->m_rvdv_spt);
  const size_t nisr(m_h.p_acraic->NISRPhotons());
  if (entry) entry->assign(cph.size(), -1);
  int nmatch(0);
  for (size_t b(0); b < nd.size(); ++b) {
    long a(-1); double best(1e30);
    for (size_t c(0); c < cph.size(); ++c) {
      const Vec4D d(cph[c] - nd[b].first);
      const double m(Max(Max(dabs(d[0]),dabs(d[1])),Max(dabs(d[2]),dabs(d[3]))));
      if (m < best) { best = m; a = (long)c; }
    }
    if (a < 0 || !(best < 1e-6*Max(nd[b].first.E(), 1e-30))) { ++m_h.m_acrv_nomatch; continue; }
    if (!m_h.p_acraic->PhotonHasM1(a)) { ++m_h.m_acrv_nom1; continue; }
    const double dvisr((size_t)a < nisr && b < spt.size() ? IsrVertexShift(spt[b]) : 0.);
    dv[a] = nd[b].second - dvisr;
    if (entry) (*entry)[a] = (long)b;
    ++nmatch;
  }
  return nmatch;
}

double YFS_Weight_Builder::IsrVertexShift() const
{
  return IsrVertexShift(m_h.m_s);
}

double YFS_Weight_Builder::IsrVertexShift(double spt) const
{
  static const bool isrvs(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["ISR_VERTEX_AT_S"].Get<bool>());
  if (!isrvs || m_h.m_ev.m_plab.size() < 2) return 0.;
  return m_h.p_acraic->InitialVertexShift(spt, (m_h.m_ev.m_plab[0] + m_h.m_ev.m_plab[1]).Abs2());
}

/*
  ACRAIC: REAL_VIRTUAL 2. The ACRAIC photon list (m_h.m_ev.m_ISRPhotons +
  m_h.m_ev.m_FSRPhotons, AllPhotonsLab) and YFS.NLO's (m_ISRPhotons +
  m_me_photons) are built independently, so each ACRAIC photon takes the dv of
  the NLO photon with the same lab momentum (to 1e-6 of its energy); photons
  without a match or below RV_SOFT_CUT keep dv = 0 (their remainder -> 0 in
  the soft limit anyway), and are counted. YFS: RV_PROBE 1 prints per event
  "@@@ ACRV": ACRAIC's addition over rho_crude next to YFS.NLO's rv/B.
*/
double YFS_Weight_Builder::AcraicRealVirtualRemainder()
{
  if (!m_h.p_acraic || !m_h.p_nlo || !m_h.p_nlo->RVRemainder()) return 0.;
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  std::vector<double> dv(cph.size(), 0.);
  const int nmatch(MatchRealVirtualDV(dv));
  ++m_h.m_acrv_n;
  const double r(m_h.p_acraic->RealVirtualRemainderRho(dv));
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  /*
    YFS: RV_PROBE_MIN (default 0): print only events with |acrv| above it, and
    then per photon its energy fraction, its smallest lab angle to a beam and
    to an outgoing charged leg, and its dv (heavy real-virtual weights, Bhabha
    O(alpha^2), parsed by yfs-helpers/pion-virtual/isr-rv-vs-v2/compare.py)
  */
  static const double probemin(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE_MIN"]
                               .SetDefault(0.).Get<double>());
  const double rcr0(m_h.p_acraic->GetRhoCrude());
  if (probe && probemin > 0. && rcr0 > 0. && std::abs(r/rcr0) > probemin) {
    const Vec4D_Vector pl(m_h.GetMomenta());
    auto ang = [](const Vec4D &a, const Vec4D &b) {
      const double c(Vec3D(a)*Vec3D(b)/(Vec3D(a).Abs()*Vec3D(b).Abs()));
      return acos(Max(-1., Min(1., c)));
    };
    std::ostringstream o;
    o<<std::setprecision(6)<<"@@@ ACRVBIG acrv="<<r/rcr0<<" nph="<<cph.size()
     <<" nnlo="<<nd.size()<<" matched="<<nmatch;
    for (size_t c(0); c < cph.size(); ++c) {
      double tb(1e9), tf(1e9);
      for (size_t i(0); i < pl.size(); ++i) {
        if (i < 2) tb = Min(tb, ang(cph[c], pl[i]));
        else if (m_h.m_flavs.size() > i && m_h.m_flavs[i].Charge() != 0.) tf = Min(tf, ang(cph[c], pl[i]));
      }
      o<<" | x="<<2.*cph[c].E()/sqrt(m_h.m_s)<<" thb="<<tb<<" thf="<<tf<<" dv="<<dv[c];
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  if (probe && !(probemin > 0.)) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ ACRV nph="<<cph.size()<<" nnlo="<<nd.size()
     <<" matched="<<nmatch<<" acrv="<<(rcr > 0. ? r/rcr : 0.)
     <<" nlorv="<<(m_h.m_born != 0. ? m_h.m_ev.m_nlo_rv/m_h.m_born : 0.)
     <<" rho1/rcr="<<(rcr > 0. ? m_h.p_acraic->GetResult()/rcr : 0.)
     <<" nloreal="<<(m_h.m_born != 0. ? 1. + m_h.m_ev.m_nlo_real/m_h.m_born : 0.)
     <<" dv="<<(nd.size() == 1 ? nd[0].second : 0.)
     <<" x="<<(nd.size() == 1 ? 2.*nd[0].first.E()/sqrt(m_h.m_s) : 0.)
     <<" nisr="<<m_h.m_ev.m_ISRPhotons.size()<<"\n";
    std::cerr<<o.str();
  }
  return r;
}
/*
  ACRAIC: REAL_VIRTUAL helicity. As AcraicRealVirtualRemainder, with photon
  j's remainder resolved by helicity: r_j^h from NLO_Base::m_rvhel (filled
  next to m_rvdv, same index), mapped to ACRAIC's index at the photon's drawn
  helicity. A photon without per-helicity factors keeps dv_j/2 on every h,
  i.e. the averaged remainder. ACRAIC: RV_HEL_FORCE_AVERAGE 1 sets dv_j/2 on
  every photon: the regression against REAL_VIRTUAL averaged through this
  code path.
*/
double YFS_Weight_Builder::AcraicRealVirtualHelicityRemainder()
{
  if (!m_h.p_acraic || !m_h.p_nlo || !m_h.p_nlo->RVRemainder()) return 0.;
  static const bool force(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["RV_HEL_FORCE_AVERAGE"].SetDefault(false).Get<bool>());
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  const std::vector<RVHelicityFactors> &nh(m_h.p_nlo->m_rvhel);
  const int nhel(Amplitude::NHel());
  std::vector<double> dv(cph.size(), 0.);
  std::vector<std::vector<Complex> > r(cph.size(),
                                       std::vector<Complex>(nhel, Complex(0., 0.)));
  bool anyhel(false);
  std::vector<long> entry;
  MatchRealVirtualDV(dv, &entry);
  for (size_t a(0); a < cph.size(); ++a) {
    if (entry[a] < 0) continue;
    const size_t b(entry[a]);
    // the ISR vertex shift MatchRealVirtualDV took off dv is helicity-blind:
    // half of it off every amplitude-level remainder
    const Complex shift(0.5*(nd[b].second - dv[a]), 0.);
    const Complex half(0.5*dv[a], 0.);
    if (force || b >= nh.size() ||
        !m_h.p_acraic->RealVirtualHelicityToAcraic(nh[b], m_h.p_acraic->PhotonHelicity(a),
                                               half, r[a]))
      r[a].assign(nhel, half);
    else {
      for (Complex &rh : r[a]) rh -= shift;
      anyhel = true;
    }
  }
  ++m_h.m_acrv_n;
  // ACRAIC: LOOP_BETA_MAX 2: the two-photon one-loop residual, built on the
  // same per-photon remainders r (Acraic_LoopN.C); 0 otherwise
  const double rloop(m_h.p_acraic->LoopNRemainderRho(m_h.p_nlo->m_loopnhel, r));
  const double rhel(m_h.p_acraic->HelicityRealVirtualRemainderRho(r) + rloop);
  const double ravg(m_h.p_acraic->RealVirtualRemainderRho(dv));
  const double rcr(m_h.p_acraic->GetRhoCrude());
  if (rcr > 0.) m_h.p_acraic->AccumulateRealVirtualHelicity(anyhel, rhel/rcr, ravg/rcr);
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  if (probe) {
    std::ostringstream o;
    o<<std::setprecision(12)<<"@@@ ACRVH nph="<<cph.size()<<" nnlo="<<nd.size()
     <<" anyhel="<<anyhel<<" hel="<<(rcr > 0. ? rhel/rcr : 0.)
     <<" avg="<<(rcr > 0. ? ravg/rcr : 0.)
     <<" nlorv="<<(m_h.m_born != 0. ? m_h.m_ev.m_nlo_rv/m_h.m_born : 0.)
     <<" rho1/rcr="<<(rcr > 0. ? m_h.p_acraic->GetResult()/rcr : 0.)<<"\n";
    std::cerr<<o.str();
  }
  return rhel;
}
double YFS_Weight_Builder::CalculateNLO(){
  InitNLO();
  // one-shot fixed-point dump for the KKMC CEEX comparison (YFS: ACRAIC_Compare)
  m_h.p_nlo->ACRAICComparePoint();
  InitNLO();
  m_h.m_ev.m_nlo_real = m_h.p_nlo->CalculateReal();
  if (m_h.m_acraic_compare) {
    m_h.m_ev.m_extrealphot.clear();
    for (const YFS::Photon &g : m_h.p_nlo->m_photons)
      m_h.m_ev.m_extrealphot.push_back(std::make_pair(g.K(), g.beta10()));
  }
  // Hardest-photon-only contributions are captured as a side effect of the
  // nominal sums above (see NLO_Base::CalculateReal/CalculateRealVirtual/
  // CalculateRealReal) - no extra ME evaluation needed here.
  m_h.m_ev.m_nlo_real_hardest = m_h.p_nlo->m_real_hard1;
  m_h.m_ev.m_nlo_real_2hardest = m_h.p_nlo->m_real_hard2;
  // read now: the double real calls CalculateReal(k) again, moving m_ifi_prod
  m_h.m_ev.m_nlo_real_ifisoft = m_h.p_nlo->m_real_ifi_soft;
  InitNLO();
  m_h.m_ev.m_nlo_virtual = m_h.p_nlo->CalculateVirtual();
  // YFS: VV_POLE_CHECK with a double-virtual provider (diagnostic)
  if (m_h.p_nlo->HasVV()) m_h.p_nlo->CheckVVPoles();
  InitNLO();
  m_h.m_ev.m_nlo_rv = m_h.p_nlo->CalculateRealVirtual();
  m_h.m_ev.m_nlo_rv_hardest = m_h.p_nlo->m_rv_hard1;
  m_h.m_ev.m_nlo_rv_2hardest = m_h.p_nlo->m_rv_hard2;
  InitNLO();
  m_h.m_ev.m_nlo_rr = m_h.p_nlo->CalculateRealReal();
  m_h.m_ev.m_nlo_rr_2hardest = m_h.p_nlo->m_rr_hard2;
  // Everything above the double real. This has to live HERE, not in
  // NLO_Base::CalculateNLO(): that function is never called - this is the
  // driver - so anything added to it silently does nothing.
  // CalculateRealMultiplicity returns 0 for any n without a provider, so a run
  // that leaves YFS: NLO_MAX_PHOTONS at its default of 2 pays nothing.
  m_h.m_ev.m_nlo_rn = 0.;
  m_h.m_ev.m_nlo_vxr = 0.;
  // NLO_MAX_PHOTONS bounds it: providers above it exist only for ACRAIC's
  // tree residuals (ACRAIC: TREE_BETA_MAX)
  const size_t nrmax(Min(m_h.p_nlo->MaxRealPhotons(),
                         YFS::NLO_Base::RequestedMaxRealPhotons()));
  for (size_t n(3); n <= nrmax; ++n) {
    InitNLO();
    m_h.m_ev.m_nlo_rn += m_h.p_nlo->CalculateRealMultiplicity(n);
  }
  /*
    YFS: VIRTUAL_COMBINE. Where the virtual v = V/B enters the squared-level
    weight (name or old integer). sum (0): 1 + v + sum_j delta_j, v times the
    sampled Born on every event. product (1, default): the product
    (1 + v)(1 + sum_j delta_j), v times
    the event's own real correction - the placement ACRAIC has, where the
    virtual sits inside sum_h |A_1 + (v/2) A_0|^2. The two agree at O(alpha);
    the product adds v * sum_j delta_j, nothing ACRAIC lacks.

    Why (e+e- -> gamma gamma, Z pole, BVR, 2026-09-26): with the sum,
    YFS.NLO/ACRAIC is 1.22 at hardest-ISR x 0.6-0.8 and 1.33 at x 0.8-1.
    70-76% of that excess comes from events with rho_1/rho_crude < 0.01
    (0.3-0.5% of ACRAIC): Borns on the generator's reduced-frame t-channel
    pole, which the exact real empties but the sum still pays v * B.
    Stands aside when a real-virtual is requested: the legacy RV (YFS:
    RV_MODE 0) supplies the v x real term itself. The RV_MODE 1 remainder
    does not - it is defined as the real-virtual beyond (1 + v) x real - so
    there the term is kept, whatever VIRTUAL_COMBINE says.
    v x RR and v x R_n are O(alpha^3). With RR_MODE singles_product they are
    left out: they put the double real, and its soft cut, into the O(alpha)
    columns YFS.NLO and YFS.Virtual (Z-pole mu mu: 0.4% per decade of
    RR_SOFT_CUT, NOTES-yfsnlo-rr-cutdep-2026-09-30.md). The other modes keep
    them, unchanged.
  */
  static const virtualcombine::code vcomb(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["VIRTUAL_COMBINE"].SetDefault(virtualcombine::product).Get<virtualcombine::code>());
  const bool vxr_real_only(m_h.p_nlo->RRSinglesProduct());
  const bool addvxr((vcomb == virtualcombine::product && !m_h.p_nlo->HasRealVirtual())
                    || m_h.p_nlo->RVRemainder());
  if (addvxr && m_h.m_born != 0.) {
    const double beyond_real(vxr_real_only ? 0. : m_h.m_ev.m_nlo_rr + m_h.m_ev.m_nlo_rn);
    const double vxr((m_h.m_ev.m_nlo_virtual/m_h.m_born)*(m_h.m_ev.m_nlo_real + beyond_real));
    if (!IsBad(vxr)) m_h.m_ev.m_nlo_vxr = vxr;
  }
  /*
    YFS: BETA0_2 - the squared-level beta_0^(2), v2 = |1 + delta + d2|^2
    - 1 - v with delta = v1/2 (NLO_Base::Beta02SquaredLevel): it carries
    the |m_1|^2 = (v1/2)^2 YFS.NLO lacks and 2 Re delta^(2). It goes where
    v goes: on the Born, and on the real where v multiplies the real, so
    that 1 + v + v2 is the Born-level virtual factor of the product.
  */
  m_h.m_ev.m_nlo_beta02 = 0.;
  if (m_h.p_nlo->Beta02Nlo().On() && m_h.p_nlo->Beta02Ok() && m_h.m_born != 0.) {
    const double v2(m_h.p_nlo->Beta02SquaredLevel());
    const double beyond_real(vxr_real_only ? 0. : m_h.m_ev.m_nlo_rr + m_h.m_ev.m_nlo_rn);
    const double b02(v2*(m_h.m_born + (addvxr ? m_h.m_ev.m_nlo_real + beyond_real : 0.)));
    if (!IsBad(b02)) m_h.m_ev.m_nlo_beta02 = b02;
  }
  // YFS: DYSON_VP - sum_h w_h (|1 + r_h|^2 - 1), placed like v2
  m_h.m_ev.m_nlo_vp = 0.;
  if (m_h.p_nlo->DysonVPOn() && m_h.p_nlo->DysonVPFactors().ok && m_h.m_born != 0.) {
    const double beyond_real(vxr_real_only ? 0. : m_h.m_ev.m_nlo_rr + m_h.m_ev.m_nlo_rn);
    const double vp(m_h.p_nlo->DysonVPFactors().v2
                    *(m_h.m_born + (addvxr ? m_h.m_ev.m_nlo_real + beyond_real : 0.)));
    if (!IsBad(vp)) m_h.m_ev.m_nlo_vp = vp;
  }
  const double total(m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_rv
                     + m_h.m_ev.m_nlo_rr + m_h.m_ev.m_nlo_rn + m_h.m_ev.m_nlo_vxr);
  const double with02(m_h.p_nlo->Beta02Nlo().On() ? total + m_h.m_ev.m_nlo_beta02 : total);
  return m_h.p_nlo->DysonVPOn() ? with02 + m_h.m_ev.m_nlo_vp : with02;
}
/*
  delta^(2) of the ACRAIC column's external virtual (beta_0^(2), ACRAIC: BETA0_2),
  of v' = v - dIB with ACRAIC: RESONANCE_EXPONENTIATION greco (its delta_f w
  must not see the resonant log, which the propagators carry). With dIB = 0
  this is NLO_Base::Beta02External() bit for bit.
*/
static double AcraicBeta02External(const YFS::NLO_Base &nlo, double dIB)
{
  return Beta02DeltaReal(nlo.Beta02Factors(), 0.5*(nlo.VOneLoop() - dIB),
                         nlo.Beta02Acraic());
}

double YFS_Weight_Builder::AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const
{
  return 2.*d2*(r0x + 0.5*v*r0v) + d2*d2*m_h.p_acraic->GetResult0();
}
void YFS_Weight_Builder::GenerateWeight() {
  if (m_h.m_dump_dipoles) m_h.p_dipoles->DumpDipoles();
  m_h.AddFormFactor();
  if (m_h.m_mode == yfsmode::isrfsr) m_h.m_ev.m_yfsweight = m_h.m_isrWeight * m_h.m_fsrWeight;
  else if (m_h.m_mode == yfsmode::fsr) m_h.m_ev.m_yfsweight = m_h.m_fsrWeight;
  else m_h.m_ev.m_yfsweight = m_h.m_isrWeight;
  if (m_h.m_coulomb) m_h.m_ev.m_yfsweight *= m_h.p_coulomb->GetWeight();
  if (m_h.m_coulombho != coulombho::none) m_h.m_ev.m_yfsweight *= m_h.m_ev.m_coulho;
  if (m_h.m_formWW) m_h.m_ev.m_yfsweight *= m_h.m_ev.m_ww_formfact; //*exp(m_ev.m_coulSub);
  CalculateBeta();

  double wif = 1.;
  if (m_h.m_ifireal && m_h.m_mode == yfsmode::isrfsr && m_h.m_nlotype == nlo_type::born &&
      m_h.p_nlo && m_h.p_nlo->HasReal()) {
    Vec4D_Vector allphotons(m_h.m_ev.m_ISRPhotons);
    allphotons.insert(allphotons.end(), m_h.m_ev.m_FSRPhotons.begin(), m_h.m_ev.m_FSRPhotons.end());
    // Same eikonal ratios the real correction subtracted with (see
    // NLO_Base::m_wifterms); the event-dipole ratio only when the fixed-order
    // side did not see every photon.
    const std::vector<double> &wt(m_h.p_nlo->WIFTerms());
    static const bool wifnlo(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["IFI_REAL_FROM_NLO"].SetDefault(true).Get<bool>());
    if (wifnlo && wt.size() == allphotons.size() && !wt.empty()) {
      wif = 1.; for (double t : wt) wif *= t;
    } else
      wif = m_h.p_dipoles->RealIFWeight(allphotons);
  }
  // The Born-level YFS weight: ISR x FSR crude (plus Coulomb/WW if on) times
  // the form factor, with NO NLO correction applied. This is what YFS.LO has
  // to reproduce -- built here directly rather than recovered downstream as
  // 1/m_ev.m_real, so the LO column cannot inherit anything m_ev.m_real does.
  double w_lo = m_h.m_ev.m_yfsweight * m_h.m_ev.m_formfactor * (1.-m_h.m_v);
  /*
    YFS.LO without IFI_Real. With IFI_Real on, the IF exponent stops at
    IFIOmega() (the generation cutoff) and the interference above it is
    restored only by wif, which the LO column does not carry. The LO weight is
    then the "exponent lowered, nothing restoring it" configuration of
    Define_Dipoles::IFIOmega(), whose A_FB runs with log(IR_CUTOFF): at the Z
    pole it made dsigma/dcos(theta) of YFS.LO strongly forward-peaked. So LO
    takes the IF form factor at sqrt(s)/2, the whole soft integral, which is
    exactly the LO of an IFI_Real 0 run (the event generation does not depend
    on IFI_Real). Same conditions as the ACRAIC IF cut below.
    YFS: LO_WITHOUT_IFI: false restores the old column.
  */
  static const bool lonoifi(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["LO_WITHOUT_IFI"].SetDefault(true).Get<bool>());
  if (lonoifi && m_h.m_ifireal && m_h.m_ifisub && m_h.m_fullform >= 1 &&
      m_h.p_dipoles) {
    const double dy(m_h.p_dipoles->FormFactorSumIF(0.5*sqrt(m_h.m_s))
                    - m_h.p_dipoles->FormFactorSumIF());
    if (!IsBad(dy)) w_lo *= exp(dy);
  }
  /*
    IFI_Real is an EEX-only correction and must NOT reach the ACRAIC weight.

    RealIFWeight() supplies the real initial-final interference that the EEX
    beta expansion does not have. ACRAIC's partition sum IS that interference:
    the 2^n assignments of each photon to the initial or the final line are
    summed COHERENTLY, and the cross terms between them are real IFI by
    construction. Adding wif on top counts it twice.

    So the two corrections are built separately - the nominal one with wif, the
    ACRAIC one without - and whichever drives the weight is chosen after.
  */
  const double corr_eex (m_h.m_ev.m_real + (wif - 1.));

  AcraicColumns c;
  c.corr_acraic = m_h.m_ev.m_acraicfactor;            // 0 if ACRAIC produced nothing
  /*
    ACRAIC: ORDER 2 keeps the O(alpha^1) column of the same event alongside
    (named weight ACRAIC_O1), composed with the external virtual exactly as an
    ORDER 1 run composes it, so ORDER2/ORDER1 is a same-event ratio.
    GetResultVV() is |A_v|^2 with A_v the amplitude (v/2) multiplies: A_0
    (= GetResult0(), bit for bit) unless ORDER2_VIRTUAL_ON_BETA1.
  */
  c.acraico2 = (m_h.p_acraic && m_h.p_acraic->Order() == 2);
  c.corr_acraic_o0 = AcraicO0Stream() && m_h.p_acraic && c.corr_acraic != 0.
                     && m_h.p_acraic->GetRhoCrude() > 0.
                     ? m_h.p_acraic->GetResult0()/m_h.p_acraic->GetRhoCrude() : 0.;
  c.corr_acraic_o0_ifcut = c.corr_acraic_o0;
  /*
    ACRAIC: BETA0_2 - delta^(2) on A_0 (YFS/Main/Beta0_2.H), the increments of
    the nominal and the two ORDER 2 compositions, added after the virtual and
    before the IF cut; corr_acraic_nob02 keeps the column without it
    (ACRAIC_NO_B02). The O(alpha^1) column ACRAIC_O1 never carries it.
  */
  c.acraicb02 = (m_h.p_acraic && m_h.p_nlo && m_h.p_nlo->Beta02Acraic().On());
  /*
    YFS: DYSON_VP in ACRAIC: r(h) on A_0 alone, with delta^(2) when that is
    on too (one amplitude-level insertion e(h) = d2 + r(h)); vp_* are the
    parts beyond beta_0^(2) (ACRAIC_NO_VP keeps the column without them).
  */
  c.acraicvp = (m_h.p_acraic && m_h.p_nlo && m_h.p_nlo->DysonVPOn() &&
               m_h.p_nlo->DysonVPFactors().ok);
  if (c.acraico2 && c.corr_acraic != 0.) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    c.corr_acraic_o1 = rcr > 0. ? m_h.p_acraic->GetResult1()/rcr : 0.;
    c.corr_acraic_o1v1 = c.corr_acraic_o1;
    c.corr_acraic_v0 = c.corr_acraic_v1 = c.corr_acraic;
  }
  /*
    ACRAIC_Virtual: external - the provider's virtual in the ACRAIC column.
    See ApplyExternalVirtual for the derivation.
  */
  if (c.corr_acraic != 0. && (m_h.m_acraicvirtsrc == acraicvirt::external ||
                          m_h.m_acraicvirtsrc == acraicvirt::automatic) && m_h.p_acraic &&
      m_h.p_nlo && m_h.p_nlo->HasVirtual() && m_h.m_ev.m_nlo_current && m_h.m_born != 0.)
    ApplyExternalVirtual(c);
  // With ACRAIC: ORDER 2 the helicity virtual goes into A_2, on the amplitude
  // the external one multiplies (Acraic_Base::HelicityVirtualRho).
  if (m_h.p_acraic) for (double &x : m_h.p_acraic->m_lhcol) x = 0.;
  /*
    ACRAIC_Virtual: helicity (ACRAIC) - the same insertion with a complex,
    helicity-resolved factor from the loop provider's helicity amplitudes.
    See ApplyHelicityVirtual for the derivation.
  */
  if (c.corr_acraic != 0. && m_h.m_acraicvirtsrc == acraicvirt::helicity && m_h.p_acraic &&
      m_h.p_nlo && m_h.p_nlo->HasVirtual() && m_h.m_ev.m_nlo_current && m_h.m_born != 0.)
    ApplyHelicityVirtual(c);
  /*
    The IF part of the form factor, cut where ACRAIC needs it.

    With IFI_Real off, IFIOmega() is sqrt(s)/2: the exponent holds the WHOLE
    soft initial-final integral, because no explicit photon of the EEX/NLO
    weights carries interference. ACRAIC's coherent partition sum does carry it,
    for every generated photon above the generation cutoff, so for the ACRAIC
    weight the region between that cutoff and sqrt(s)/2 was counted twice -
    once in exp(Y_IF) and once in |sum_g s_g|^2. The symptom is an A_FB that
    runs with log(IR_CUTOFF): e+e- -> mu mu at 0.7 GeV (CMD cuts) gave -3.4,
    -4.9, -6.1, -7.0% at IR_CUTOFF 1e-4..1e-7 against KKMC's -0.08 +- 0.14%;
    with the exponent cut at the generation cutoff it is +1.03/+1.06% at 1e-5
    and 1e-7. So the ACRAIC column swaps exp(Y_IF(IFIOmega)) for
    exp(Y_IF(omega_gen)); with IFI_Real on the two coincide and this is 1.
    ACRAIC: IF_FORMFACTOR_CUT: false restores the old behaviour.
  */
  FinalizeAcraicBeta02VP(c);
  ApplyIFFormFactorCut(c);
  /*
    ACRAIC: FIXED_ORDER_STREAM (2026-10-11) - ACRAIC truncated to O(alpha), on
    the resummed events, as the named weight ACRAIC_FixedOrder. See
    ComputeFixedOrderColumn for the derivation.
  */
  ComputeFixedOrderColumn(c);

  // NLO_Part: B is beta_0 only: neither ACRAIC nor the EEX beta_1 may correct
  // the nominal; their columns stay as named weights
  const bool   bornonly   (m_h.m_nlotype == nlo_type::born);
  if (bornonly && m_h.m_useacraic) {
    static bool shown(false);
    if (!shown) { shown = true;
      msg_Info()<<"YFS: NLO_Part: B - ACRAIC is not used for the nominal "
                <<"weight (beta_0 only); its correction is kept in the named "
                <<"weights only. Add R (and V) to NLO_Part to use it."<<std::endl; }
  }
  const bool   acraic_nom (!bornonly && m_h.m_acraic_weight && c.corr_acraic != 0.);
  const double corr_nom   (bornonly ? 1. : acraic_nom ? c.corr_acraic : corr_eex);
  m_h.m_ev.m_yfsweight *= corr_nom;
  m_h.m_ev.m_yfsweight *= m_h.m_ev.m_formfactor*(1.-m_h.m_v);
  // What the named ACRAIC column has to divide by to become a ratio.
  m_h.m_ev.m_corr_nominal = corr_nom;
  m_h.m_ev.m_corr_acraic    = c.corr_acraic;
  m_h.m_ev.m_corr_eex     = corr_eex;
  m_h.m_ev.m_corr_acraic_o1 = c.corr_acraic_o1;
  m_h.m_ev.m_corr_acraic_o0 = c.corr_acraic_o0_ifcut;
  m_h.m_ev.m_corr_acraic_fo = c.corr_acraic_fo;
  m_h.m_ev.m_corr_acraic_o1s = c.corr_acraic_o1s;
  m_h.m_ev.m_corr_acraic_o1v1 = c.corr_acraic_o1v1;
  m_h.m_ev.m_corr_acraic_v0 = c.corr_acraic_v0;
  m_h.m_ev.m_corr_acraic_v1 = c.corr_acraic_v1;
  m_h.m_ev.m_corr_acraic_v2 = c.corr_acraic_v2;
  m_h.m_ev.m_corr_acraic_nob02 = c.corr_acraic_nob02;
  m_h.m_ev.m_corr_acraic_novp = c.corr_acraic_novp;
  m_h.m_ev.m_tnp1.swap(c.tnp1);
  m_h.m_ev.m_tnp2.swap(c.tnp2);
  // Captured before the IsBad/negative-weight clamps below, since the named
  // weights are ratios against the weight the event actually carries.
  const double w_full = m_h.m_ev.m_yfsweight;
  if(m_h.m_isr_debug) {
    Vec4D ele;
    for (int i = 2; i < m_h.m_flavs.size(); ++i)
    {
      if(IsEqual(m_h.m_flavs[i],kf_e)) {
        ele = m_h.m_ev.m_plab[m_h.p_dipoles->m_flav_label[m_h.m_flavs[i]]];
        m_h.p_beams->BoostBackLab(ele);
        m_h.p_debug->FillHist("Form_Factor_FS_Angle", ele.Theta()*1000,m_h.m_ev.m_formfactor,1);
      }
    }
  }
  DEBUG_FUNC("\nISR Weight = " << m_h.m_isrWeight << "\n" <<
             "  FSR Weight = " << m_h.m_fsrWeight << "\n" <<
             "  WW form Weight = " << m_h.m_ev.m_ww_formfact << "\n" <<
             "  Total form Weight = " << m_h.m_ev.m_formfactor << "\n" <<
             "  Coulomb Weight = " << m_h.p_coulomb->GetWeight() << "\n" <<
             " Coulomb Subtraction Weight = " << exp(m_h.m_ev.m_coulSub) << "\n" <<
             "Total Weight = " << m_h.m_ev.m_yfsweight << "\n");
  if(IsBad(m_h.m_ev.m_yfsweight)){
    msg_Error()<<"\nISR Weight = " << m_h.m_isrWeight << "\n" <<
             "  FSR Weight = " << m_h.m_fsrWeight << "\n" <<
             "  Form Factor = " << m_h.m_ev.m_formfactor << "\n" <<
             "  NLO  Correction = " << m_h.m_ev.m_real << "\n" <<
             "Total Weight = " << m_h.m_ev.m_yfsweight << "\n";
    m_h.m_ev.m_yfsweight = 0;
    ++m_h.m_nonfinite;
  }
  if(m_h.m_ev.m_yfsweight < 0 && m_h.m_skipNegWeights){
    msg_Debugging()<<"Skipping negative Weight in YFS"<<std::endl;
    m_h.m_ev.m_yfsweight=0;
    m_h.m_negskip++;
  }

  ProbeSign(w_full);
  BuildNamedWeights(w_lo, w_full);
}

//! ACRAIC_Virtual: external/automatic, split out of GenerateWeight.
//!
//! The squared-level provider's virtual v on the Born-level amplitude A_0
//! inside the coherent sum:
//!     rho = sum_h |A_0 (1 + delta) + (A_1 - A_0)|^2
//!         = rho_1 + 2 Re(delta) Re<A_0, A_1> + |delta|^2 rho_0,
//! with delta = v/2. Also folds in beta_0^(2), Dyson VP and the ACRAIC:
//! REAL_VIRTUAL non-factorisable remainder on the same v.
void YFS_Weight_Builder::ApplyExternalVirtual(AcraicColumns &c)
{
  const double rcr(m_h.p_acraic->GetRhoCrude());
  // ACRAIC: ISR_VERTEX_AT_S - the initial vertex from the post-ISR Born point
  // (m_plab, where the provider's virtual was taken) to the full s
  const double dvisr(IsrVertexShift());
  /*
    ACRAIC: RESONANCE_EXPONENTIATION greco: the resonant initial-final log is
    in the resonance propagators of the ACRAIC amplitudes now, so it leaves
    the provider's virtual: v' = v - dI_B, dI_B its exact Born-weighted
    share at the Born point (Acraic_Base::BuildResonanceBornDelta; 0, and v
    bit for bit, with none). Everything below - beta_0^(2), Dyson VP, the
    ORDER 2 streams - takes v'.
  */
  const double v(m_h.m_ev.m_nlo_virtual/m_h.m_born - m_h.p_acraic->ResonanceBornDeltaV() + dvisr);
  c.tnp_vhalf = 0.5*v;
  c.vext_fo = v;
  if (AcraicStrictO1Stream() && rcr > 0.) {
    const double s1((m_h.p_acraic->GetResult1() + v*m_h.p_acraic->GetResult01Order1())/rcr);
    if (!IsBad(s1)) c.corr_acraic_o1s = s1;
  }
  const double add(rcr > 0. ? (v*m_h.p_acraic->GetResult01()
                               + 0.25*v*v*m_h.p_acraic->GetResultVV())/rcr : 0.);
  if (!IsBad(add)) c.corr_acraic += add;
  else ++m_h.m_acraicstats.m_bad;
  if (c.acraicb02 && m_h.p_nlo->Beta02Ok() && rcr > 0.) {
    // delta^(2) of v' (its delta_f w must not see the resonant log either)
    const double d2(AcraicBeta02External(*m_h.p_nlo, m_h.p_acraic->ResonanceBornDeltaV()));
    c.b02_acraic = AcraicBeta02Rho(d2, v, m_h.p_acraic->GetResult0X(), m_h.p_acraic->GetResult0V())/rcr;
    if (c.acraico2) {
      c.b02_v0 = AcraicBeta02Rho(d2, v, m_h.p_acraic->GetResult02(), m_h.p_acraic->GetResult0())/rcr;
      c.b02_v1 = AcraicBeta02Rho(d2, v, m_h.p_acraic->GetResult02(), m_h.p_acraic->GetResult01Order1())/rcr;
    }
  }
  if (c.acraicvp && rcr > 0.) {
    const double d2(c.acraicb02 && m_h.p_nlo->Beta02Ok()
                    ? AcraicBeta02External(*m_h.p_nlo, m_h.p_acraic->ResonanceBornDeltaV()) : 0.);
    std::vector<Complex> e(m_h.p_acraic->HelicityShares(m_h.p_nlo->DysonVPFactors()));
    for (Complex &x : e) x += d2;
    const std::vector<Complex> e2(e.size(), Complex(d2, 0.));
    const bool vb1(m_h.p_acraic->VirtualOnBeta1());
    c.vp_acraic = (m_h.p_acraic->A0InsertionRho(e, 0.5*v, vb1)
               - m_h.p_acraic->A0InsertionRho(e2, 0.5*v, vb1))/rcr;
    if (c.acraico2) {
      c.vp_v0 = (m_h.p_acraic->A0InsertionRho(e, 0.5*v, false)
               - m_h.p_acraic->A0InsertionRho(e2, 0.5*v, false))/rcr;
      c.vp_v1 = (m_h.p_acraic->A0InsertionRho(e, 0.5*v, true)
               - m_h.p_acraic->A0InsertionRho(e2, 0.5*v, true))/rcr;
    }
  }
  // ACRAIC: REAL_VIRTUAL 2 - the non-factorisable real-virtual, YFS.NLO's
  // per-photon v_{n+1} - v_B on each photon's M_1 (Acraic_Base.C)
  const acraicrv::code rvmode(YFS::Acraic_Base::RealVirtualMode());
  if ((rvmode == acraicrv::averaged || rvmode == acraicrv::helicity) && rcr > 0.) {
    // the provider's per-photon loop has both vertices and the
    // initial-final boxes, which a one-sided MODE leaves out
    if (m_h.p_acraic->SingleSided())
      THROW(fatal_error, "ACRAIC: REAL_VIRTUAL averaged takes the loop "
            "provider's real-virtual, which contains both vertices and the "
            "initial-final loops; YFS: MODE ISR or FSR (one side only) needs ACRAIC: "
            "{REAL_VIRTUAL: factorisable} with ACRAIC's own virtual "
            "(YFS: ACRAIC_Virtual: acraic).");
    const double arv((rvmode == acraicrv::helicity ? AcraicRealVirtualHelicityRemainder()
                                                  : AcraicRealVirtualRemainder())/rcr);
    if (!IsBad(arv)) c.corr_acraic += arv;
  }
  if (c.acraico2) {
    const double add1(rcr > 0. ? (v*m_h.p_acraic->GetResult01Order1()
                                  + 0.25*v*v*m_h.p_acraic->GetResult0())/rcr : 0.);
    if (!IsBad(add1)) c.corr_acraic_o1 += add1;
    const double add1v1(rcr > 0. ? (v + 0.25*v*v)*m_h.p_acraic->GetResult1()/rcr : 0.);
    if (!IsBad(add1v1)) c.corr_acraic_o1v1 += add1v1;
    const double a0(rcr > 0. ? (v*m_h.p_acraic->GetResult02()
                                + 0.25*v*v*m_h.p_acraic->GetResult0())/rcr : 0.);
    const double a1(rcr > 0. ? (v*m_h.p_acraic->GetResult12()
                                + 0.25*v*v*m_h.p_acraic->GetResult1())/rcr : 0.);
    if (!IsBad(a0)) c.corr_acraic_v0 += a0;
    if (!IsBad(a1)) c.corr_acraic_v1 += a1;
    const double rho2(m_h.p_acraic->GetResult());   // |A_2|^2 at ORDER 2
    const double d2(rcr > 0. ? ((v + 0.25*v*v)*rho2
                                - v*m_h.p_acraic->GetResult01()
                                - 0.25*v*v*m_h.p_acraic->GetResultVV())/rcr : 0.);
    if (!IsBad(d2)) c.dv2_acraic = d2;
  }
}

//! ACRAIC_Virtual: helicity, split out of GenerateWeight: the same insertion
//! with a complex, helicity-resolved factor from the loop provider
//! (NLO_Base::BuildLoopHelicityFactors), delta^h in place of v/2:
//!     rho_V = 1/4 sum_h |A_1^h + delta^h A_0^h|^2 ,
//! added as (rho_V - rho_1)/rho_crude. Also fills ACRAIC's m_lhcol (the
//! ACRAIC_EXTV/L/R named-weight columns).
void YFS_Weight_Builder::ApplyHelicityVirtual(AcraicColumns &c)
{
  static const bool forcehalf(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                              ["VIRT_HEL_FORCE_HALF_V"].SetDefault(0).Get<int>() != 0);
  const double rcr(m_h.p_acraic->GetRhoCrude());
  const double vfull(m_h.m_ev.m_nlo_virtual/m_h.m_born);
  std::vector<Complex> d;
  bool ok(!forcehalf && m_h.p_nlo->m_lhelok &&
          m_h.p_acraic->LoopHelicityToAcraic(m_h.p_nlo->m_lhel, m_h.p_nlo->m_la0,
                                     m_h.p_nlo->m_ldelta, m_h.p_nlo->m_lhelmom,
                                     Complex(0.5*vfull, 0.), d));
  if (!ok) d.assign(Amplitude::NHel(), Complex(0.5*vfull, 0.));
  /*
    ACRAIC: RESONANCE_EXPONENTIATION greco: the resonant log leaves delta^h
    helicity by helicity, exactly (Acraic_Base::BuildResonanceBornDelta),
    and v' = v - dI_B is its Born-weighted 2 Re; with none both are zero.
  */
  const double v(vfull - m_h.p_acraic->ResonanceBornDeltaV());
  c.tnp_vhalf = 0.5*v;
  if (m_h.p_acraic->ResonanceExponentiation()) {
    const std::vector<Complex> &dI(m_h.p_acraic->ResonanceBornDeltaHelicity());
    for (size_t f(0); f < d.size() && f < dI.size(); ++f) d[f] -= dI[f];
  }
  const double rho1(m_h.p_acraic->GetResult());
  const double add(rcr > 0. ? (m_h.p_acraic->HelicityVirtualRho(d) - rho1)/rcr : 0.);
  // the external form on the same event, and rho_V(v/2) against it
  const double ext(v*m_h.p_acraic->GetResult01() + 0.25*v*v*m_h.p_acraic->GetResultVV());
  const std::vector<Complex> dh(Amplitude::NHel(), Complex(0.5*v, 0.));
  /*
    ACRAIC: BETA0_2 per helicity: delta^h -> delta^h + delta^(2)(delta^h),
    complex, with delta_f in the same phase convention (Beta0_2.H), so that
    the O(alpha^2) sum does not depend on VIRT_HEL_IM_C. Events without
    helicity factors use the external, real form on v/2.
  */
  std::vector<Complex> db(d), dhb(dh);
  const bool hb02(c.acraicb02 && m_h.p_nlo->Beta02Ok() && rcr > 0.);
  if (hb02) {
    const double d2ext(AcraicBeta02External(*m_h.p_nlo, m_h.p_acraic->ResonanceBornDeltaV()));
    for (size_t f(0); f < db.size(); ++f)
      db[f] += ok ? Beta02Delta(m_h.p_nlo->Beta02Factors(), d[f], m_h.p_nlo->Beta02Acraic())
                  : Complex(d2ext, 0.);
    // BETA0_2 two_loop: the photon-exchange two-loop boxes per helicity
    // (the external fallback above has them Born-weighted already)
    if (ok && m_h.p_nlo->Beta02Acraic().mode == beta02::two_loop &&
        !m_h.p_nlo->Beta02Factors().boxhel.empty()) {
      const std::vector<Complex> bx(m_h.p_acraic->PhysicalLabelsToAcraic
        (m_h.p_nlo->Beta02Factors().boxhel, m_h.p_nlo->Beta02Factors().boxd2, Complex(0., 0.)));
      for (size_t f(0); f < db.size() && f < bx.size(); ++f) db[f] += bx[f];
    }
    for (Complex &x : dhb) x += Complex(d2ext, 0.);
    c.b02_acraic = (m_h.p_acraic->HelicityVirtualRho(db) - m_h.p_acraic->HelicityVirtualRho(d))/rcr;
  }
  /*
    ACRAIC: EW_HIGHER_ORDER_HELICITY - EW_HIGHER_ORDER (in v, so in the v/2
    fallback) is otherwise missing from delta^h: averaged puts v_HO/2 on
    every helicity, helicity GRIFFIN's own per-helicity factor. After
    delta^(2), whose delta_f w must not see this O(alpha^2) term.
  */
  double gho_acraic(0.);
  if (ok && rcr > 0. && m_h.p_nlo->EWHOHelicityMode() != ewhohel::none &&
      m_h.p_nlo->VEWHigherOrder() != 0.) {
    const Complex half(0.5*m_h.p_nlo->VEWHigherOrder(), 0.);
    std::vector<Complex> g(Amplitude::NHel(), half);
    if (m_h.p_nlo->EWHOHelicityMode() == ewhohel::helicity && m_h.p_nlo->m_ghelok)
      g = m_h.p_acraic->PhysicalLabelsToAcraic(m_h.p_nlo->m_ghel, m_h.p_nlo->m_gdelta, half);
    std::vector<Complex> dg(db);
    for (size_t f(0); f < dg.size() && f < g.size(); ++f) dg[f] += g[f];
    gho_acraic = (m_h.p_acraic->HelicityVirtualRho(dg) - m_h.p_acraic->HelicityVirtualRho(db))/rcr;
    db.swap(dg);
  }
  // YFS: DYSON_VP per helicity, on top of delta^h (+ delta^(2))
  const bool hvp(c.acraicvp && rcr > 0.);
  if (hvp) {
    const std::vector<Complex> r(m_h.p_acraic->HelicityShares(m_h.p_nlo->DysonVPFactors()));
    std::vector<Complex> dv(db);
    for (size_t f(0); f < dv.size() && f < r.size(); ++f) {
      dv[f] += r[f];
      dhb[f] += r[f];
    }
    c.vp_acraic = (m_h.p_acraic->HelicityVirtualRho(dv) - m_h.p_acraic->HelicityVirtualRho(db))/rcr;
    db.swap(dv);
  }
  const double reg(std::abs(m_h.p_acraic->HelicityVirtualRho(dh) - rho1 - ext)
                   /Max(std::abs(rho1), 1e-300));
  double split[3] = {0., 0., 0.};
  if (ok && rcr > 0.) {
    m_h.p_acraic->HelicityVirtualSplit(d, 0.5*v, split);
    for (double &x : split) x /= rcr;
  }
  m_h.p_acraic->AccumulateLoopHelicity(ok, add, rcr > 0. ? ext/rcr : 0., reg, split);
  /*
    The same event for the extra named weights: ACRAIC_EXTV (external
    virtual), ACRAIC_L/R and ACRAIC_EXTV_L/R (beam 1 fully left/right
    polarised, beam 2 unpolarised; the unpolarised column is their mean).
    Stored as absolute factors before the IF cut, which multiplies all.
  */
  const bool hextra(hb02 || hvp || gho_acraic != 0.);
  if (rcr > 0. && !IsBad(add) && !hextra) {
    m_h.p_acraic->m_lhcol[0] = c.corr_acraic + add;
    m_h.p_acraic->m_lhcol[1] = c.corr_acraic + ext/rcr;
    m_h.p_acraic->m_lhcol[2] = m_h.p_acraic->HelicityVirtualRho(d,  -1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[3] = m_h.p_acraic->HelicityVirtualRho(d,  +1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[4] = m_h.p_acraic->HelicityVirtualRho(dh, -1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[5] = m_h.p_acraic->HelicityVirtualRho(dh, +1., 0.)/rcr;
  }
  // the same columns with beta_0^(2) in each (the external one on v/2)
  if (rcr > 0. && !IsBad(add) && hextra) {
    m_h.p_acraic->m_lhcol[0] = c.corr_acraic + add + c.b02_acraic + c.vp_acraic + gho_acraic;
    m_h.p_acraic->m_lhcol[1] = c.corr_acraic
      + (m_h.p_acraic->HelicityVirtualRho(dhb) - rho1)/rcr;
    m_h.p_acraic->m_lhcol[2] = m_h.p_acraic->HelicityVirtualRho(db,  -1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[3] = m_h.p_acraic->HelicityVirtualRho(db,  +1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[4] = m_h.p_acraic->HelicityVirtualRho(dhb, -1., 0.)/rcr;
    m_h.p_acraic->m_lhcol[5] = m_h.p_acraic->HelicityVirtualRho(dhb, +1., 0.)/rcr;
  }
  if (!IsBad(add)) c.corr_acraic += add;
  else ++m_h.m_acraicstats.m_bad;
  if (!IsBad(gho_acraic)) c.corr_acraic += gho_acraic;
  /*
    ORDER 2: the helicity analogues of the external compositions of the same
    event, delta^h in place of v/2. ACRAIC_O1 |A_1 + d A_0|^2, ACRAIC_O1V1
    |A_1 + d A_1|^2, ACRAIC_O2V0/V1 |A_2 + d A_0|^2 and |A_2 + d A_1|^2, and
    ACRAIC_O2V2 the nominal with |A_2 + d A_2|^2 in place of |A_2 + d A_v|^2.
  */
  if (c.acraico2 && rcr > 0.) {
    YFS::Acraic_Base &a(*m_h.p_acraic);
    const double rho2(a.GetResult()), r1(a.GetResult1());
    const double o1((a.HelicityInsertionRho(d, 1, 0) - r1)/rcr);
    const double o1v1((a.HelicityInsertionRho(d, 1, 1) - r1)/rcr);
    const double v0((a.HelicityInsertionRho(d, 2, 0) - rho2)/rcr);
    const double v1((a.HelicityInsertionRho(d, 2, 1) - rho2)/rcr);
    const double d2((a.HelicityInsertionRho(d, 2, 2) - a.HelicityVirtualRho(d))/rcr);
    if (!IsBad(o1)) c.corr_acraic_o1 += o1;
    if (!IsBad(o1v1)) c.corr_acraic_o1v1 += o1v1;
    if (!IsBad(v0)) c.corr_acraic_v0 += v0;
    if (!IsBad(v1)) c.corr_acraic_v1 += v1;
    if (!IsBad(d2)) c.dv2_acraic = d2;
  }
}

void YFS_Weight_Builder::FinalizeAcraicBeta02VP(AcraicColumns &c)
{
  c.corr_acraic_nob02 = c.corr_acraic;
  if (c.acraicb02 && c.corr_acraic != 0.) {
    if (!IsBad(c.b02_acraic)) c.corr_acraic += c.b02_acraic;
    if (c.acraico2 && !IsBad(c.b02_v0)) c.corr_acraic_v0 += c.b02_v0;
    if (c.acraico2 && !IsBad(c.b02_v1)) c.corr_acraic_v1 += c.b02_v1;
  }
  c.corr_acraic_novp = c.corr_acraic;
  if (c.acraicvp && c.corr_acraic != 0.) {
    if (!IsBad(c.vp_acraic)) { c.corr_acraic += c.vp_acraic; c.corr_acraic_nob02 += c.vp_acraic; }
    if (c.acraico2 && !IsBad(c.vp_v0)) c.corr_acraic_v0 += c.vp_v0;
    if (c.acraico2 && !IsBad(c.vp_v1)) c.corr_acraic_v1 += c.vp_v1;
  }
  c.corr_acraic_v2 = c.acraico2 && c.corr_acraic != 0. ? c.corr_acraic + c.dv2_acraic : 0.;
  // ACRAIC: TNP_STREAMS - templates, 0 on events without an ACRAIC weight
  if (m_h.p_acraic && m_h.p_acraic->TnpOn() && c.corr_acraic != 0. && m_h.p_acraic->GetRhoCrude() > 0.) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    if (m_h.p_acraic->TnpTemplates(1, c.tnp_vhalf, c.tnp1)) for (double &x : c.tnp1) x /= rcr;
    else c.tnp1.clear();
    if (c.acraico2 && m_h.p_acraic->TnpTemplates(2, c.tnp_vhalf, c.tnp2)) for (double &x : c.tnp2) x /= rcr;
    else c.tnp2.clear();
  }
}

void YFS_Weight_Builder::ApplyIFFormFactorCut(AcraicColumns &c)
{
  static const bool ifcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["IF_FORMFACTOR_CUT"].SetDefault(true).Get<bool>());
  if (ifcut && c.corr_acraic != 0. && m_h.m_ifisub && m_h.m_fullform >= 1 &&
      m_h.p_dipoles) {
    const double wgen(0.5*sqrt(m_h.m_s)*m_h.m_isrcut);
    const double dy(m_h.p_dipoles->FormFactorSumIF(wgen)
                    - m_h.p_dipoles->FormFactorSumIF());
    if (!IsBad(dy)) { c.dyif_fo = dy; c.corr_acraic_o1s *= exp(dy);
                      c.corr_acraic *= exp(dy); c.corr_acraic_o1 *= exp(dy); c.corr_acraic_o1v1 *= exp(dy);
                      c.corr_acraic_o0_ifcut *= exp(dy);
                      c.corr_acraic_v0 *= exp(dy); c.corr_acraic_v1 *= exp(dy); c.corr_acraic_v2 *= exp(dy);
                      c.corr_acraic_nob02 *= exp(dy); c.corr_acraic_novp *= exp(dy);
                      for (double &x : c.tnp1) x *= exp(dy);
                      for (double &x : c.tnp2) x *= exp(dy); }
  }
}

void YFS_Weight_Builder::ComputeFixedOrderColumn(AcraicColumns &c)
{
  if (AcraicFixedOrderStream() && m_h.p_acraic && c.corr_acraic != 0. && m_h.p_dipoles &&
      m_h.p_acraic->GetRhoCrude() > 0.) {
    if (m_h.m_fullform < 1 || m_h.m_fullform > 2 || m_h.m_coulomb ||
        m_h.m_acraicvirtsrc == acraicvirt::helicity || m_h.p_dipoles->FormFactorSumDecay() != 0.)
      THROW(not_implemented, "ACRAIC: FIXED_ORDER_STREAM needs FULL_FORM 1 or 2, "
            "no Coulomb factor, no W decay stages and ACRAIC_Virtual other than "
            "helicity.");
    const size_t nphot(m_h.p_acraic->NPhot());
    const double yfo(m_h.m_ev.m_formfactor_sum + c.dyif_fo
                     + m_h.p_dipoles->FormFactorShiftToRegion(0.5*sqrt(m_h.m_s)*m_h.m_isrcut));
    double rhofo(0.);
    if (nphot == 0)
      rhofo = m_h.p_acraic->GetResult1() + (c.vext_fo + yfo)*m_h.p_acraic->GetResult0();
    else if (nphot == 1)
      rhofo = m_h.p_acraic->GetResult1();
    c.corr_acraic_fo = exp(c.dyif_fo - yfo)*rhofo/m_h.p_acraic->GetRhoCrude();
    if (IsBad(c.corr_acraic_fo)) c.corr_acraic_fo = 0.;
  }
}
// The YFS.* named weights: each is a ratio that turns the nominal weight into
// one of the truncated or reordered matchings. Split out of GenerateWeight,
// which was 216 lines with more than half of them here.
//
// w_lo   : the Born-level YFS weight, before any NLO correction
// w_full : the weight the event actually carries, before the clamps below
void YFS_Weight_Builder::BuildNamedWeights(double w_lo, double w_full) {
  m_h.m_ev.m_nlo_weightsmap = Weights_Map{1.0};
  Weights wyfs{1.0};
  bool any(false);

  std::vector<std::string> names, tnpnames;
  RegisterColumnNames(names, tnpnames, wyfs);
  if (!names.empty()) any = true;

  FillAcraicColumns(wyfs);
  FillTnpColumns(wyfs, tnpnames);
  FillDiffBudgetColumns(wyfs);
  FillWithoutColumns(wyfs);
  if (FillNLOBreakdownColumns(wyfs, w_lo, w_full)) any = true;

  if (BuildLadderWeights(wyfs)) any = true;

  if (any) m_h.m_ev.m_nlo_weightsmap["YFS"] = wyfs;
}

//! The YFS.* column names this configuration emits. Fixed once per event
//! from configuration alone, and the only names emit() will write -- a
//! per-event code path must not be able to conjure a column that other
//! events lack, since HepMC binds the set from event one.
void YFS_Weight_Builder::RegisterColumnNames(std::vector<std::string> &names,
                                             std::vector<std::string> &tnpnames,
                                             ATOOLS::Weights &wyfs)
{
  const bool hasnlo (m_h.p_nlo && m_h.p_nlo->HasNLO());
  const bool hasnnlo(m_h.p_nlo && m_h.p_nlo->HasNNLO());
  if (m_h.m_nlotype != nlo_type::born && (hasnlo || hasnnlo)) {
    if (hasnlo) {
      names.push_back("LO");
      if (m_h.m_nlo_columns) names.push_back("NLO");
    }
    if (hasnnlo) names.push_back("NNLO");
    if (m_h.m_nlo_weight_breakdown) {
      if (hasnlo)
        for (const char *n : {"Real","Virtual","BR","BV","NLO_1g",
                              "NLO_2g","NLO_FixedOrder"})
          names.push_back(n);
      // EEX only when the beta expansion ran: m_ev.m_eex is assigned solely inside
      // CalculateBeta's `if (m_betaorder > 0)`, so at BETA:0 the column would
      // be a no-op 1.0 masquerading as a measurement. m_betaorder is
      // configuration, so gating the name on it keeps the set constant.
      if (hasnlo && m_h.m_betaorder > 0) names.push_back("EEX");
      if (hasnnlo)
        for (const char *n : {"RealVirtual","RealReal","NLO+RR","NLO+RV",
                              "NNLO_1g","NNLO_2g","NNLO_FixedOrder"})
          names.push_back(n);
      // the EEX estimate of the double virtual, retired by YFS: BETA0_2
      // (its LL part, without the NLL and the 100% envelope); never both
      if (hasnnlo && !m_h.p_nlo->Beta02Nlo().On())
        for (const char *n : {"VV_EEX","NNLO_VV","NNLO_VV_up","NNLO_VV_down"})
          names.push_back(n);
    }
  }
  /*
    ACRAIC gets its own column whenever ACRAIC is on, independently of the NLO
    names above: it is defined at Born level too, where none of those exist.
  */
  if (m_h.m_useacraic) names.push_back("ACRAIC");
  // ACRAIC nominal (YFS: ACRAIC_WEIGHT, default with ACRAIC on): the squared-level
  // YFS.NLO weight the event would otherwise carry, as a named column
  if (m_h.m_useacraic && m_h.m_acraic_weight && m_h.m_nlo_columns)
    names.push_back("NLO_EW");
  // YFS: DIFF_BUDGET_STREAMS: the per-event ingredients (filled further down)
  if (DiffBudgetStreamsOn())
    for (const char *n : s_diffbudgetnames) names.push_back(n);
  // BETA0_2: the column without beta_0^(2) on the same event, per column
  if (m_h.m_useacraic && m_h.p_nlo && m_h.p_nlo->Beta02Acraic().On() &&
      m_h.m_acraicvirtsrc != acraicvirt::acraic)
    names.push_back("ACRAIC_NO_B02");
  // COULOMB_HIGHER_ORDER: every column without the beyond-one-loop W-pair
  // Coulomb factor, on the same event
  if (m_h.m_coulombho != coulombho::none) names.push_back("NoCoulombHO");
  if (m_h.p_nlo && m_h.p_nlo->Beta02Nlo().On()) names.push_back("NLO_EW_NO_B02");
  if (m_h.p_nlo && m_h.p_nlo->DysonVPOn()) {
    if (m_h.m_useacraic && m_h.m_acraicvirtsrc != acraicvirt::acraic) names.push_back("ACRAIC_NO_VP");
    names.push_back("NLO_EW_NO_VP");
  }
  // ACRAIC_Virtual: helicity - the external-virtual ACRAIC column on the same
  // events, and the beam-1-polarised columns, see GenerateWeight
  if (m_h.m_useacraic && m_h.m_acraicvirtsrc == acraicvirt::helicity)
    for (const char *n : {"ACRAIC_EXTV","ACRAIC_L","ACRAIC_R","ACRAIC_EXTV_L","ACRAIC_EXTV_R"})
      names.push_back(n);
  if (m_h.m_useacraic && AcraicO0Stream()) names.push_back("ACRAIC_O0");
  if (m_h.m_useacraic && AcraicFixedOrderStream()) names.push_back("ACRAIC_FixedOrder");
  if (m_h.m_useacraic && AcraicStrictO1Stream()) names.push_back("ACRAIC_O1S");
  // ACRAIC: ORDER 2 also carries the O(alpha^1) ACRAIC column of each event
  if (m_h.m_useacraic && ATOOLS::Settings::GetMainSettings()["ACRAIC"]["ORDER"]
      .SetDefault(1).Get<int>() == 2) {
    names.push_back("ACRAIC_O1");
    names.push_back("ACRAIC_O1V1");
    names.push_back("ACRAIC_O2V0");
    names.push_back("ACRAIC_O2V1");
    names.push_back("ACRAIC_O2V2");
  }
  // ACRAIC: TNP_STREAMS - templates, 0 on events without an ACRAIC weight
  if (m_h.m_useacraic && m_h.p_acraic && m_h.p_acraic->TnpOn()) {
    for (int c(0); c <= YFS::Acraic_Base::tnp_multi; ++c)
      tnpnames.push_back(std::string("ACRAIC_TNP1_") + YFS::Acraic_Base::TnpClassName(c));
    if (m_h.p_acraic->Order() == 2)
      for (int c(0); c < YFS::Acraic_Base::tnp_nclass; ++c)
        tnpnames.push_back(std::string("ACRAIC_TNP2_") + YFS::Acraic_Base::TnpClassName(c));
    names.insert(names.end(), tnpnames.begin(), tnpnames.end());
  }
  // ACRAIC: TNP_BAND_STREAMS - nominal +- N T per parameter, 1 (no shift) on
  // events without an ACRAIC weight
  if (AcraicTnpBandStreams() && m_h.m_useacraic && m_h.p_acraic) {
    if (!m_h.p_acraic->TnpOn())
      THROW(fatal_error, "ACRAIC: TNP_BAND_STREAMS needs ACRAIC: TNP_STREAMS true.");
    if (m_h.HasFSR())
      THROW(not_implemented, "ACRAIC: TNP_BAND_STREAMS has the initial-state "
            "parameters only (YFS: MODE ISR); with FSR combine the "
            "ACRAIC_TNP1_/ACRAIC_TNP2_ template weights instead.");
    const int k(m_h.p_acraic->Order());
    static bool shown(false);
    if (!shown) {
      shown = true;
      if (k == 1)
        msg_Info()<<"ACRAIC: TNP_BAND_STREAMS at O(alpha): the two-photon classes "
                  <<"missing entirely at this order need an ORDER 2 run "
                  <<"and are not in the band."<<std::endl;
    }
    const double me(m_h.m_flavs[0].Mass() > 0. ? m_h.m_flavs[0].Mass() : 0.000510999);
    for (const TnpBandParam &b :
           TnpBandParams(sqr(rpa->gen.Ecms()), me, k,
                         m_h.p_nlo && m_h.p_nlo->Beta02Acraic().On(),
                         AcraicTnpBandIncomplete())) {
      names.push_back("ACRAIC_TNPB_" + b.name + "_UP");
      names.push_back("ACRAIC_TNPB_" + b.name + "_DN");
    }
  }
  if (m_h.m_ladder_weights) {
    if (m_h.m_coulomb && m_h.p_coulomb) names.push_back("NoCoulomb");
    if (m_h.m_ifisub && m_h.m_fullform >= 1 &&
        m_h.FixedOrder() != fixed_order::nlo && m_h.p_dipoles)
      names.push_back("NoIFI");
  }
  for (const std::string &n : names) wyfs[n] = 1.;
  for (const std::string &n : tnpnames) wyfs[n] = 0.;
  m_h.m_wnames.clear();
  m_h.m_wnames.insert(names.begin(), names.end());
}

void YFS_Weight_Builder::FillAcraicColumns(ATOOLS::Weights &wyfs)
{
  if (m_h.m_useacraic && m_h.m_wnames.count("ACRAIC") && m_h.m_ev.m_corr_acraic != 0. &&
      !IsZero(m_h.m_ev.m_corr_nominal)) {
    const double r(m_h.m_ev.m_corr_acraic/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC"] = r;
    if (m_h.m_wnames.count("ACRAIC_EXTV") && m_h.p_acraic) {
      const char *nm[5] = {"ACRAIC_EXTV","ACRAIC_L","ACRAIC_R","ACRAIC_EXTV_L","ACRAIC_EXTV_R"};
      const double c0(m_h.p_acraic->m_lhcol[0]);
      for (int i(0); i < 5; ++i) {
        const double ri(c0 != 0. ? r*m_h.p_acraic->m_lhcol[i+1]/c0 : r);
        if (!IsBad(ri)) wyfs[nm[i]] = ri;
      }
    }
    if (m_h.m_wnames.count("ACRAIC_O1S") && m_h.m_ev.m_corr_acraic_o1s != 0.) {
      const double ro1s(m_h.m_ev.m_corr_acraic_o1s/m_h.m_ev.m_corr_nominal);
      if (!IsBad(ro1s)) wyfs["ACRAIC_O1S"] = ro1s;
    }
    if (m_h.m_wnames.count("ACRAIC_FixedOrder")) {
      const double rfo(m_h.m_ev.m_corr_acraic_fo/m_h.m_ev.m_corr_nominal);
      if (!IsBad(rfo)) wyfs["ACRAIC_FixedOrder"] = rfo;
    }
    if (m_h.m_wnames.count("ACRAIC_O0")) {
      const double r0(m_h.m_ev.m_corr_acraic_o0/m_h.m_ev.m_corr_nominal);
      if (!IsBad(r0)) wyfs["ACRAIC_O0"] = r0;
    }
    if (m_h.m_wnames.count("ACRAIC_O1")) {
      const double r1(m_h.m_ev.m_corr_acraic_o1/m_h.m_ev.m_corr_nominal);
      if (!IsBad(r1)) wyfs["ACRAIC_O1"] = r1;
      const double r1v1(m_h.m_ev.m_corr_acraic_o1v1/m_h.m_ev.m_corr_nominal);
      if (!IsBad(r1v1)) wyfs["ACRAIC_O1V1"] = r1v1;
      const double rv0(m_h.m_ev.m_corr_acraic_v0/m_h.m_ev.m_corr_nominal);
      const double rv1(m_h.m_ev.m_corr_acraic_v1/m_h.m_ev.m_corr_nominal);
      if (!IsBad(rv0)) wyfs["ACRAIC_O2V0"] = rv0;
      if (!IsBad(rv1)) wyfs["ACRAIC_O2V1"] = rv1;
      const double rv2(m_h.m_ev.m_corr_acraic_v2/m_h.m_ev.m_corr_nominal);
      if (!IsBad(rv2) && m_h.m_ev.m_corr_acraic_v2 != 0.) wyfs["ACRAIC_O2V2"] = rv2;
    }
  }
  if (m_h.m_wnames.count("NLO_EW") && !IsZero(m_h.m_ev.m_corr_nominal)) {
    const double r(m_h.m_ev.m_corr_eex/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW"] = r;
  }
}

void YFS_Weight_Builder::FillTnpColumns(ATOOLS::Weights &wyfs,
                                        const std::vector<std::string> &tnpnames)
{
  if (!tnpnames.empty() && !IsZero(m_h.m_ev.m_corr_nominal)) {
    // the template contribution T_c of this event: nominal weight x T_c/corr
    const size_t n1(YFS::Acraic_Base::tnp_multi + 1);
    for (size_t c(0); c < m_h.m_ev.m_tnp1.size() && c < n1; ++c)
      if (!IsBad(m_h.m_ev.m_tnp1[c])) wyfs[tnpnames[c]] = m_h.m_ev.m_tnp1[c]/m_h.m_ev.m_corr_nominal;
    for (size_t c(0); c < m_h.m_ev.m_tnp2.size() && n1 + c < tnpnames.size(); ++c)
      if (!IsBad(m_h.m_ev.m_tnp2[c])) wyfs[tnpnames[n1 + c]] = m_h.m_ev.m_tnp2[c]/m_h.m_ev.m_corr_nominal;
  }
  // ACRAIC: TNP_BAND_STREAMS - 1 +- N sum_c T_c/corr, from the templates of
  // the run's own order (registered, and so checked, in RegisterColumnNames)
  if (AcraicTnpBandStreams() && m_h.p_acraic && m_h.p_acraic->TnpOn() &&
      !IsZero(m_h.m_ev.m_corr_nominal)) {
    const int k(m_h.p_acraic->Order());
    const std::vector<double> &t(k == 2 ? m_h.m_ev.m_tnp2 : m_h.m_ev.m_tnp1);
    if (t.empty()) return;
    const double me(m_h.m_flavs[0].Mass() > 0. ? m_h.m_flavs[0].Mass() : 0.000510999);
    for (const TnpBandParam &b :
           TnpBandParams(sqr(rpa->gen.Ecms()), me, k,
                         m_h.p_nlo && m_h.p_nlo->Beta02Acraic().On(),
                         AcraicTnpBandIncomplete())) {
      double r(0.);
      for (int c : b.classes) if (c < (int)t.size()) r += t[c];
      r *= b.norm/m_h.m_ev.m_corr_nominal;
      if (IsBad(r)) continue;
      wyfs["ACRAIC_TNPB_" + b.name + "_UP"] = 1. + r;
      wyfs["ACRAIC_TNPB_" + b.name + "_DN"] = 1. - r;
    }
  }
}

/*
  YFS: DIFF_BUDGET_STREAMS (diagnostic, default false). The ingredients of
  the two O(alpha) correction factors of this event,
      ACRAIC      C = a1 (1 + v'/2)^2   (ACRAIC: REAL_VIRTUAL factorisable, ORDER 1)
      YFS.NLO_EW  N = (1 + v)(1 + r)    (YFS: VIRTUAL_COMBINE product)
  with v' = v - x (x = dI_B, 0 without RESONANCE_EXPONENTIATION greco), as
  columns DB_<q> = w_crude q, i.e. the ratio q/corr_nominal. A per-event
  dump (paper/rivet CRAIC_DUMP) divides by DB_CRUDE and splits C - N into
  named pieces offline. a1 = rho_1/rho_crude and a0 = rho_0/rho_crude carry
  the IF form-factor cut that multiplies the ACRAIC column (1 with IFI_Real
  true). r is m_nlo_real/B (with the IF restoration), RIF = 1 + that
  restoration (prod_j r_IF,j). X1, X2 (hardest and second-hardest visible
  photon, x = 2E/sqrt(s), + ISR / - FSR), NPH (their number) and SXI, SXF
  (the sum of x over the ISR and over the FSR photons) are payloads, not
  weights. A column is 0 when its event has no ACRAIC factor.
  Histogram-ready budget columns (2026-10-05, agree1005; true weights, in
  the same w_crude units, so a CRAIC_EE histogram of each is the piece per
  bin with its own error). With N = corr_eex (YFS.NLO_EW):
      VSQ  = a1 v^2/4          the amplitude-square piece V
      REST = a1 (1 + v) - N    C - N - VSQ - (greco part)
      NPV  = N + a1 v^2/4      YFS.NLO_EW + V, the "N + V" prediction
  In a RESONANCE_EXPONENTIATION none run a1 = a1n, so these are exactly V,
  the remainder R = C - N - V and N + V. In a greco run they use a1g; there
  G = C_g - C_n comes from the none twin on the same events.
*/
void YFS_Weight_Builder::FillDiffBudgetColumns(ATOOLS::Weights &wyfs)
{
  if (DiffBudgetStreamsOn()) {
    std::map<std::string, double> q;
    for (const char *n : s_diffbudgetnames) q[n] = 0.;
    const double cnom(m_h.m_ev.m_corr_nominal);
    if (m_h.p_acraic && m_h.m_ev.m_corr_acraic != 0. && !IsZero(cnom) && m_h.m_born != 0. &&
        m_h.p_acraic->GetRhoCrude() > 0.) {
      double ifcut(1.);
      static const bool ifcuton(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["IF_FORMFACTOR_CUT"].Get<bool>());
      if (ifcuton && m_h.m_ifisub && m_h.m_fullform >= 1 && m_h.p_dipoles) {
        const double dy(m_h.p_dipoles->FormFactorSumIF(0.5*sqrt(m_h.m_s)*m_h.m_isrcut)
                        - m_h.p_dipoles->FormFactorSumIF());
        if (!IsBad(dy)) ifcut = exp(dy);
      }
      const double rcr(m_h.p_acraic->GetRhoCrude());
      std::vector<double> xs;   // signed photon x, + ISR / - FSR
      double sxi(0.), sxf(0.);
      for (const Vec4D &k : m_h.m_ev.m_ISRPhotons) { xs.push_back( 2.*k[0]/sqrt(m_h.m_s)); sxi += xs.back(); }
      for (const Vec4D &k : m_h.m_ev.m_FSRPhotons) { xs.push_back(-2.*k[0]/sqrt(m_h.m_s)); sxf -= xs.back(); }
      std::sort(xs.begin(), xs.end(),
                [](double a, double b) { return std::abs(a) > std::abs(b); });
      q["DB_CRUDE"] = 1.;
      q["DB_A1"]  = m_h.p_acraic->GetResult()/rcr*ifcut;
      q["DB_A0"]  = m_h.p_acraic->GetResult0()/rcr*ifcut;
      q["DB_V"]   = m_h.m_ev.m_nlo_virtual/m_h.m_born;
      q["DB_X"]   = m_h.p_acraic->ResonanceBornDeltaV();
      q["DB_R"]   = m_h.m_ev.m_nlo_real/m_h.m_born;
      q["DB_RIF"] = 1. + m_h.m_ev.m_nlo_real_ifisoft/m_h.m_born;
      q["DB_X1"]  = xs.size() > 0 ? xs[0] : 0.;
      q["DB_X2"]  = xs.size() > 1 ? xs[1] : 0.;
      q["DB_NPH"] = xs.size();
      q["DB_SXI"] = sxi;
      q["DB_SXF"] = sxf;
      q["DB_VSQ"] = q["DB_A1"]*q["DB_V"]*q["DB_V"]/4.;
      q["DB_REST"] = q["DB_A1"]*(1. + q["DB_V"]) - m_h.m_ev.m_corr_eex;
      q["DB_NPV"] = m_h.m_ev.m_corr_eex + q["DB_VSQ"];
      for (auto &kv : q) kv.second /= cnom;
    }
    for (const auto &kv : q) wyfs[kv.first] = IsBad(kv.second) ? 0. : kv.second;
  }
}

void YFS_Weight_Builder::FillWithoutColumns(ATOOLS::Weights &wyfs)
{
  if (m_h.m_wnames.count("NoCoulombHO") && !IsZero(m_h.m_ev.m_coulho) && !IsBad(m_h.m_ev.m_coulho))
    wyfs["NoCoulombHO"] = 1./m_h.m_ev.m_coulho;
  if (m_h.m_wnames.count("ACRAIC_NO_B02") && m_h.m_ev.m_corr_acraic != 0. &&
      !IsZero(m_h.m_ev.m_corr_nominal)) {
    const double r(m_h.m_ev.m_corr_acraic_nob02/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC_NO_B02"] = r;
  }
  if (m_h.m_wnames.count("ACRAIC_NO_VP") && m_h.m_ev.m_corr_acraic != 0. &&
      !IsZero(m_h.m_ev.m_corr_nominal)) {
    const double r(m_h.m_ev.m_corr_acraic_novp/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC_NO_VP"] = r;
  }
  if (m_h.m_wnames.count("NLO_EW_NO_VP") && !IsZero(m_h.m_ev.m_corr_nominal) &&
      m_h.m_born != 0.) {
    const double r((m_h.m_ev.m_corr_eex - m_h.m_ev.m_nlo_vp/m_h.m_born)/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW_NO_VP"] = r;
  }
  if (m_h.m_wnames.count("NLO_EW_NO_B02") && !IsZero(m_h.m_ev.m_corr_nominal) &&
      m_h.m_born != 0.) {
    const double r((m_h.m_ev.m_corr_eex - m_h.m_ev.m_nlo_beta02/m_h.m_born)/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW_NO_B02"] = r;
  }
}

bool YFS_Weight_Builder::FillNLOBreakdownColumns(ATOOLS::Weights &wyfs, double w_lo, double w_full)
{
  if (!(m_h.m_ev.m_nlo_current && m_h.m_nlotype != nlo_type::born && !IsZero(m_h.m_ev.m_real) &&
      (m_h.p_nlo->HasNLO() || m_h.p_nlo->HasNNLO())))
    return false;
  auto ratio = [this](double term, double denom) -> double {
    return term / denom;
  };
  /*
    The YFS.NLO breakdown columns are factors x against m_real, the
    squared-level correction. They were written as x/m_real because the
    nominal WAS that correction; with ACRAIC as the nominal the event carries
    corr_acraic instead, so the denominator is rescaled by corr_nominal/corr_eex
    and every column keeps exactly the absolute weight it had (identical to
    m_real when YFS.NLO is the nominal).
  */
  const double nomden(IsZero(m_h.m_ev.m_corr_eex) ? m_h.m_ev.m_real
                      : m_h.m_ev.m_real*(m_h.m_ev.m_corr_nominal/m_h.m_ev.m_corr_eex));

  // Values only. Writing through operator[] CREATES a column, which is how
  // YFS.EEX kept appearing at BETA:0 after it was dropped from the
  // registration list -- and a column that is zero on every event comes back
  // as NaN once Rivet divides by its own sumOfWeights. So refuse any name the
  // configuration did not register, and leave a non-finite value at the
  // registered 1.0 rather than poisoning the column.
  auto emit = [this, &wyfs](const std::string &name, double value) {
    if (!m_h.m_wnames.count(name)) return;
    if (!IsBad(value)) wyfs[name] = value;
  };

  // NB (2026-10-11): NLO_FixedOrder / NNLO_FixedOrder are not strict
  // truncations: every event keeps its photons, the form factor becomes
  // 1+Y, and multi-photon events still contribute with the hardest photon's
  // real. Never compared with an independent fixed-order result; the strict
  // O(alpha) truncation is ACRAIC_FixedOrder (ACRAIC: FIXED_ORDER_STREAM).
  const bool have_fixed_order_ff =
      m_h.m_fullform >= 1 && !IsZero(m_h.m_ev.m_formfactor);
  // The cached exponent, not a fresh FormFactorSum(): the two must be the
  // same Y or this is not a truncation of anything.
  const double ff_fixedorder_ratio =
      have_fixed_order_ff ? (1. + m_h.m_ev.m_formfactor_sum) / m_h.m_ev.m_formfactor : 1.;

  /*
    YFS: IFI_Real true cuts the IF part of the exponent at IFIOmega()
    instead of sqrt(s)/2 and gives every photon above it the IF eikonal
    ratio r_j; the soft (beta_0) level paired with that exponent is
    B prod_j r_j. m_nlo_real carries B (prod_j r_j - 1); every column NOT
    built from m_nlo_real adds it back, or its exponent has lost the IF
    interference above omega (A_FB +0.15 at the Z pole). 0 with IFI_Real
    false, where all columns are unchanged.
  */
  const double ifisoft = m_h.m_ev.m_nlo_real_ifisoft;

  // NLO: Real + Virtual
  if (m_h.p_nlo->HasNLO()) {
    // m_nlo_vxr: YFS: VIRTUAL_COMBINE 1 (zero otherwise), so the named
    // column equals the nominal weight it is compared against
    const double nlo_sum   = (m_h.m_born + m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual
                              + m_h.m_ev.m_nlo_vxr)/m_h.m_born;
    const double real_sum  = (m_h.m_born + m_h.m_ev.m_nlo_real)/m_h.m_born;
    const double virt_sum  = (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_virtual)/m_h.m_born;
    if (!IsZero(m_h.m_ev.m_real)) {
      emit("Real", ratio((m_h.m_ev.m_nlo_real)/m_h.m_born, nomden));
      emit("Virtual", ratio((m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_vxr)/m_h.m_born, nomden));
      emit("NLO", ratio(nlo_sum, nomden));
      emit("BR", ratio(real_sum, nomden));
      emit("BV", ratio(virt_sum, nomden));
      // LO = (Born-level YFS weight) / (full weight), so that
      // nominal * YFS.LO == w_lo identically. Algebraically 1/m_ev.m_real, but
      // built from the two weights themselves.
      if (!IsZero(w_full)) emit("LO", w_lo/w_full);
      emit("EEX", ratio(m_h.m_ev.m_eex + ifisoft/m_h.m_born, nomden));
      // Matching truncated to a fixed real-photon multiplicity, to see the
      // result "as if" only the 1 or 2 hardest photons were used in the
      // matching (full "NLO" above keeps all generated photons). Real is
      // summed over the 1 / 2 hardest photons; Virtual is always full, and so
      // is the IF restoration, which belongs to the exponent, not to beta_1.
      // With one photon NLO_1g is NLO without m_nlo_vxr (either REAL_COMBINE).
      const double nlo_1g = (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_hardest
                             + m_h.m_ev.m_nlo_virtual)/m_h.m_born;
      const double nlo_2g = (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_2hardest
                             + m_h.m_ev.m_nlo_virtual)/m_h.m_born;
      emit("NLO_1g", ratio(nlo_1g, nomden));
      emit("NLO_2g", ratio(nlo_2g, nomden));
      // Fixed-order comparison point: 1-photon NLO correction with the
      // resummed exp(form) form factor undone in favour of its 1+form
      // fixed-order truncation - matches a plain (non-YFS-resummed) NLO EW
      // calculation, which only ever has at most one real photon.
      if (have_fixed_order_ff)
        emit("NLO_FixedOrder", ratio(nlo_1g, nomden) * ff_fixedorder_ratio);
    }
  }

  // NNLO: RealVirtual + RealReal
  if (m_h.p_nlo->HasNNLO()) {
    // YFS: BETA0_2's v2 belongs to the O(alpha^2) column (0 when off)
    const double nnlo_total = (m_h.m_born + m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_rv + m_h.m_ev.m_nlo_rr
                               + m_h.m_ev.m_nlo_vxr + m_h.m_ev.m_nlo_beta02 + m_h.m_ev.m_nlo_vp)/m_h.m_born;
    const double RR_total = (m_h.m_born + m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_rr
                             + m_h.m_ev.m_nlo_vxr)/m_h.m_born;
    /*
      NLO+RV is YFS.NLO plus the real-virtual, so it carries m_nlo_vxr like
      YFS.NLO and NLO+RR. Without it the column was YFS.NLO + rv - v x real
      whenever vxr is set (YFS: RV_MODE remainder, or VIRTUAL_COMBINE product
      without an RV): Z-pole mu mu NLO+RV/NLO 1.003 with IFI_Real false and
      0.991 with a 2% cos theta tilt with true, where the real carries the IF
      restoration B (prod r_IF - 1) (2026-10-04). With the legacy RV vxr is 0,
      so nothing changes there.
    */
    const double RV_total = (m_h.m_born + m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_rv
                             + m_h.m_ev.m_nlo_vxr)/m_h.m_born;
    // No separate denominator: every column here is x/m_ev.m_real. NNLO used to be
    // make_ratio(nnlo_total, m_born + nnlo_total), which is
    // (m_born + nnlo_total)/(m_born + nnlo_total) -- identically 1, whatever
    // the physics did. It also added a dimensionful m_born to a ratio.
    if (!IsZero(m_h.m_ev.m_real)) {
      emit("RealVirtual", ratio((m_h.m_ev.m_nlo_rv)/m_h.m_born, nomden));
      emit("RealReal", ratio((m_h.m_ev.m_nlo_rr)/m_h.m_born, nomden));
      emit("NLO+RR", ratio(RR_total, nomden));
      emit("NLO+RV", ratio(RV_total, nomden));
      emit("NNLO", ratio(nnlo_total, nomden));
      // Matching truncated to a fixed real-photon multiplicity, the NNLO
      // analogue of NLO_1g/NLO_2g. 1 photon: Real + RealVirtual on the
      // single hardest, RealReal = 0 (a pair needs two photons). 2 photons:
      // Real + RealVirtual summed over the two hardest, plus the RealReal
      // pair formed by them. Virtual is always full.
      const double nnlo_1g =
          (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_hardest + m_h.m_ev.m_nlo_virtual
           + m_h.m_ev.m_nlo_rv_hardest)/m_h.m_born;
      const double nnlo_2g =
          (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_2hardest + m_h.m_ev.m_nlo_virtual
           + m_h.m_ev.m_nlo_rv_2hardest + m_h.m_ev.m_nlo_rr_2hardest)/m_h.m_born;
      emit("NNLO_1g", ratio(nnlo_1g, nomden));
      emit("NNLO_2g", ratio(nnlo_2g, nomden));
      // Fixed-order NNLO comparison: the 2-photon truncation (fixed-order
      // NNLO EW allows up to two real photons) with the resummed form factor
      // undone to its 1+form truncation.
      if (have_fixed_order_ff)
        emit("NNLO_FixedOrder", ratio(nnlo_2g, nomden) * ff_fixedorder_ratio);

      // ---- approximate double-virtual (VV) ----
      // The NNLO weights above are RV + RR only: there is no exact
      // double-virtual provider, so at O(alpha^2) the VV is simply MISSING and
      // "NNLO" is incomplete. Estimate it from the EEX virtual series, whose
      // O(alpha^2) term is the 0.125*gamma^2 in Dipole::VirtualEEX. Taking the
      // difference of the EEX virtual evaluated at order 2 and order 1 isolates
      // exactly that term:
      //     VV_EEX = prod(1 + 0.5g + 0.125g^2) - prod(1 + 0.5g)
      // (products over the II and FF dipoles, so the cross terms between
      // dipoles are kept). It is a leading-log/eikonal estimate of a term
      // whose exact form is unknown here, NOT a calculation of it.
      //
      // Because it is an estimate, it is shipped with an explicit envelope
      // rather than silently folded into the nominal: VV_up/VV_down scale it
      // by 1 +/- VV_Approx_Uncertainty (default 1, i.e. the band runs from
      // "twice the estimate" down to "no VV at all"). That is the honest
      // statement of ignorance - the term is known to be of this size, but its
      // coefficient is not - and it lets the VV uncertainty be propagated as a
      // normal variation weight instead of quoted by hand.
      //
      // The nominal "NNLO" weight is deliberately left VV-free so it keeps
      // meaning what it meant before; NNLO_VV is the one including the estimate.
      // p_nlo's flag: the handler's own m_vvtool is a YFS_Base member that
      // nothing on the handler side ever sets (see YFS_Base.H)
      if (!m_h.p_nlo->HasVV()) {
        const double vv = m_h.p_dipoles->CalculateEEXVirtual(2)
                        - m_h.p_dipoles->CalculateEEXVirtual(1);
        if (!IsBad(vv)) {
          // Flat envelope: the band is (1 +/- d) times the estimate, with
          // d = VV_Approx_Uncertainty. At the default d = 1 it runs from "the
          // term is absent" to "twice the estimate", which is the conventional
          // missing-higher-order convention and the most ignorance an envelope
          // can honestly express - d > 1 would push the down variation to a
          // NEGATIVE VV, i.e. assert the opposite sign rather than absence.
          //
          // A calibrated alternative was tried (scaling d by EEX's measured
          // error on the O(alpha) virtual, where the exact result is known)
          // and dropped: on this process that ratio came out ~2.3, i.e. EEX
          // does not predict the exact virtual at all, so it carries no
          // information about the order above and only degenerated to the flat
          // band once capped.
          //
          // Keep the size of this band in perspective: VV_EEX is ~0.3% of the
          // cross section while RR alone is ~5% and the NLO->NNLO shift ~7.8%,
          // so this is NOT the dominant NNLO uncertainty.
          const double d = m_h.m_vv_approx_unc;
          emit("VV_EEX", ratio(vv, nomden));
          emit("NNLO_VV", ratio(nnlo_total + vv, nomden));
          emit("NNLO_VV_up", ratio(nnlo_total + (1.+d)*vv, nomden));
          emit("NNLO_VV_down", ratio(nnlo_total + (1.-d)*vv, nomden));
        } else {
          msg_Error() << METHOD << ": EEX double-virtual estimate is "
                      << vv << ", skipping the VV weights\n";
        }
      }
    }
  }

  return true;
}
bool YFS_Weight_Builder::BuildLadderWeights(ATOOLS::Weights &w) {
  if (!m_h.m_ladder_weights) return false;
  bool any(false);

  if (m_h.m_coulomb && m_h.p_coulomb) {
    const double wc(m_h.p_coulomb->GetWeight());
    if (!IsZero(wc) && !IsBad(wc)) { w["NoCoulomb"] = 1./wc; any = true; }
  }

  if (m_h.m_ifisub && m_h.m_fullform >= 1 &&
      m_h.FixedOrder() != fixed_order::nlo && m_h.p_dipoles) {
    const double fif(m_h.p_dipoles->FormFactorSumIF());
    if (!IsBad(fif)) { w["NoIFI"] = exp(-fif); any = true; }
  }

  return any;
}
