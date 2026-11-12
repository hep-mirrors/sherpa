#include "YFS/Main/YFS_Weight_Builder.H"
#include "YFS/Main/YFS_Handler.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/Beta0_2.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"
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
    ACRAIC: TNP_BAND_STREAMS. One nuisance parameter: norm N and the template
    classes it scales. ACRAIC_TNPB_<name>_UP/_DN = nominal +- N x summed templates.
  */
  struct TnpBandParam {
    std::string name;
    double norm;
    std::vector<int> classes;
  };

  std::vector<TnpBandParam> TnpBandParams(double s, double me, int k,
                                          bool beta02, bool incomplete)
  {
    typedef YFS::Acraic_Base A;
    const double a(1./137.035999/M_PI), g(a*(log(s/(me*me)) - 1.));
    const int mb(k == 1 || beta02 ? k + 1 : k);
    std::vector<TnpBandParam> b;
    b.push_back({"B0_ISR", 0.5*pow(g, mb)/Factorial(mb), {A::tnp_b0}});
    b.push_back({"ISR",    0.5*pow(g, k)/Factorial(k),   {A::tnp_isr}});
    b.push_back({"MULTI",  0.5*pow(g, k)/Factorial(k),   {A::tnp_multi}});
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
    m_h.m_bpmc_evtG = 1.; m_h.m_bpmc_evtn = 0;
    static const bornphotonmc::code bpmc(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(bornphotonmc::on)
      .Get<bornphotonmc::code>());
    // YFS: BORN_PHOTON_CHANNELS stage: the FSR photons are alternatives too, on
    // the post-emission legs (NLO_Base::BornPhotonChannelSum)
    static const bornphotonchannels::code bpch(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["BORN_PHOTON_CHANNELS"].SetDefault(bornphotonchannels::isr)
      .Get<bornphotonchannels::code>());
    const bool allchannels(bpch != bornphotonchannels::isr);
    std::vector<Vec4D> photons(m_h.m_ev.m_ISRPhotons.begin(), m_h.m_ev.m_ISRPhotons.end());
    if (allchannels)
      photons.insert(photons.end(), m_h.m_ev.m_FSRPhotons.begin(), m_h.m_ev.m_FSRPhotons.end());
    if (bpmc != bornphotonmc::off && m_h.p_acraic && m_h.p_nlo && m_h.p_dipoles->HasDipoleII()
        && !photons.empty()) {
      InitNLO();
      // each assignment weighed with ACRAIC's crude of the swapped event (the
      // photons above are in ACRAIC's order, ISR then FSR, as m_allphotons)
      const NLO_Base::ChannelDensity density(
          [this](const std::vector<size_t> &slot, const std::vector<size_t> &pick)
          { return m_h.p_acraic->CrudeOfPhotonAssignment(slot, pick); });
      m_h.m_bpmc_evtG = m_h.p_nlo->BornPhotonChannelSum(allchannels ? m_h.m_ev.m_plab
                                                                    : m_h.m_ev.m_reallab,
                                                photons, m_h.p_dipoles->GetDipoleII(), -1,
                                                &m_h.m_bpmc_evtn,
                                                bpch == bornphotonchannels::stage
                                                ? &density : nullptr);
      if (bpmc == bornphotonmc::on && m_h.m_bpmc_evtG != 1.)
        m_h.p_acraic->ScaleRhoCrude(m_h.m_bpmc_evtG);
    }
    if (m_h.p_acraic) {
      // Denominator is the crude (incoherent) weight: the coherent rho can cancel to ~0.
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
  ACRAIC against the external providers, event by event (diagnostic):
      external virtual / Born   <->   rho(Born+virtual)/rho(Born) - 1
      external real    / Born   <->   rho(Born+real)   /rho(Born) - 1
  The total row is the meaningful one.
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
    // Born point of the virtual: cos theta of the first final-state leg in the
    // beam CM, and the virtual subtraction split by dipole type.
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
    // Photon by photon, matched by momentum (the two lists are built independently).
    // The columns do not sum to the total: cross terms belong to no single photon.
    const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
    const std::vector<std::pair<Vec4D,double> > &eph(m_h.m_ev.m_extrealphot);
    // The sides can get different photons: ACRAIC takes ISR + FSR photons, the
    // external ISR + ME photons (final-state halves from different dipole getters).
    msg_Out()<<"      photons: acraic="<<cph.size()<<"  ext="<<eph.size()
             <<(cph.size()!=eph.size() ? "   <-- DIFFERENT LISTS" : "")<<"\n";
    msg_Out()<<"        acraic list:";
    for (size_t a(0); a < cph.size(); ++a) msg_Out()<<" "<<cph[a].E();
    msg_Out()<<"\n        ext  list:";
    for (size_t b(0); b < eph.size(); ++b) msg_Out()<<" "<<eph[b].first.E();
    msg_Out()<<"\n        nISR="<<m_h.m_ev.m_ISRPhotons.size()
             <<" nFSR="<<m_h.m_ev.m_FSRPhotons.size()
             <<" nFSRforME="<<m_h.m_ev.m_fsrphotonsforME.size()<<"\n";
    // Momentum balance p_a + p_b - q_c - q_d - sum k must vanish for the right list.
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
  m_h.p_nlo->m_photons.clear();
  if (m_h.p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_h.m_ev.m_ISRPhotons)
      m_h.p_nlo->m_photons.push_back(YFS::Photon(k, &m_h.p_dipoles->GetDipoleII()));
  if (m_h.m_nlo_fsr_photons) {
    // m_ev.m_photons: ISR first, then FF; skip the ISR head.
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
  ACRAIC: REAL_VIRTUAL isr_nll. Non-soft part V_2 of the one-loop ISR photon
  emission, LL + NLL, small-mass limit (KKMC MakeVini, Phys. Rev. D65 (2002)
  073030 eq. 2.18), r_i = p_i.k/p_1.p_2.
  Returns dv = 2 Re V_2 on the photon's M_1, on top of the factorisable (1 + v/2).
*/
static double IsrNllV2(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, double alpha)
{
  const double p1p2(p1*p2), r1((p1*k)/p1p2), r2((p2*k)/p1p2);
  if (!(r1 > 0. && r2 > 0. && r1 < 1. && r2 < 1.)) return 0.;
  return alpha/M_PI*0.5*(log(r1)*log(1.-r2) + log(r2)*log(1.-r1)
                         + DiLog(r1) + DiLog(r2)
                         - 0.5*sqr(log(1.-r1)) - 0.5*sqr(log(1.-r2))
                         + 1.5*log(1.-r1) + 1.5*log(1.-r2)
                         + 0.5*r1*(1.-r1)/(1.+sqr(1.-r1))
                         + 0.5*r2*(1.-r2)/(1.+sqr(1.-r2)));
}

double YFS_Weight_Builder::AcraicIsrNllRemainder()
{
  if (!m_h.p_acraic) return 0.;
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const Vec4D_Vector beams(m_h.BornMomenta());
  std::vector<double> dv(cph.size(), 0.);
  for (size_t j(0); j < m_h.p_acraic->NISRPhotons() && j < cph.size(); ++j)
    if (m_h.p_acraic->PhotonHasM1(j)) dv[j] = 2.*IsrNllV2(beams[0], beams[1], cph[j], m_h.m_alpha);
  return m_h.p_acraic->RealVirtualRemainderRho(dv);
}

// dv of the NLO photon matching each ACRAIC photon; ISR photons drop the vertex shift at s_pt
int YFS_Weight_Builder::MatchRealVirtualDV(std::vector<double> &dv, size_t first,
                                           std::vector<long> *entry)
{
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  const size_t nisr(m_h.p_acraic->NISRPhotons());
  const std::vector<double> &spt(m_h.p_nlo->m_rvdv_spt);
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
    if ((size_t)a < first) { ++m_h.m_acrv_isr; continue; }
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

double YFS_Weight_Builder::TnpVirtualIF() const
{
  // e+e- -> f fbar: ACRAIC's own boxes; the mirrored loop is exact only for a
  // charge-symmetric Born, which the gamma-Z interference is not
  static const tnpifivirtual::code mode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                       ["TNP_IFI_VIRTUAL"].Get<tnpifivirtual::code>());
  if (mode == tnpifivirtual::none) return 0.;
  double vown(0.);
  const bool own(mode != tnpifivirtual::mirrored && m_h.p_acraic->OwnBoxVirtualIF(vown));
  if (mode == tnpifivirtual::own_boxes && !own)
    THROW(fatal_error, "ACRAIC: TNP_IFI_VIRTUAL own_boxes needs e+e- -> f fbar (no Bhabha).");
  static const int nprobe(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["TNP_IFI_CHECK"].Get<int>());
  static int nprinted(0);
  if (nprinted < (nprobe > 0 ? 20 : 0)) {
    ++nprinted;
    std::cerr<<std::setprecision(10)<<"@@@ VIFCMP own="<<(own ? vown : 0.)
             <<" mirrored="<<m_h.p_nlo->VirtualIFOdd()<<" dI="<<m_h.p_acraic->ResonanceBornDeltaV()<<"\n";
  }
  const double dI(m_h.p_acraic->ResonanceBornDeltaV());
  if (own) return vown - dI;
  if (!m_h.p_nlo->VirtualIFOddOk()) return 0.;
  return m_h.p_nlo->VirtualIFOdd() - dI;
}

double YFS_Weight_Builder::IsrVertexShift(double spt) const
{
  static const bool isrvs(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["ISR_VERTEX_AT_S"].Get<bool>());
  if (!isrvs || m_h.m_ev.m_plab.size() < 2) return 0.;
  return m_h.p_acraic->InitialVertexShift(spt, (m_h.m_ev.m_plab[0] + m_h.m_ev.m_plab[1]).Abs2());
}

/*
  ACRAIC: REAL_VIRTUAL hybrid. The analytic LL + NLL remainder V_2 (as isr_nll)
  on the initial-state photons, YFS.NLO's exact remainder (as averaged) on the
  final-state photons. HYBRID_ISR_RMIN r_c >= 0: V_2 only for min(r_1, r_2) < r_c.
*/
double YFS_Weight_Builder::AcraicHybridRemainder()
{
  if (!m_h.p_acraic || !m_h.p_nlo || !m_h.p_nlo->RVRemainder()) return 0.;
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const Vec4D_Vector beams(m_h.BornMomenta());
  const size_t nisr(Min(m_h.p_acraic->NISRPhotons(), cph.size()));
  static const double rc(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["HYBRID_ISR_RMIN"]
                         .SetDefault(-1.).Get<double>());
  std::vector<double> dv(cph.size(), 0.);
  for (size_t j(0); j < nisr; ++j) {
    if (!m_h.p_acraic->PhotonHasM1(j)) continue;
    // as NLO_Base::InHybridIsrCone
    if (rc < 0. || Min(cph[j][0] - cph[j][3], cph[j][0] + cph[j][3])/sqrt(m_h.m_s) < rc)
      dv[j] = 2.*IsrNllV2(beams[0], beams[1], cph[j], m_h.m_alpha);
  }
  // photons above r_c take the exact dv
  MatchRealVirtualDV(dv, rc < 0. ? nisr : 0);
  ++m_h.m_acrv_n;
  const double r(m_h.p_acraic->RealVirtualRemainderRho(dv));
  // YFS: RV_PROBE 1: per event the addition over rho_crude and per photon
  // I/F, x, the smallest angle to a beam and the dv it carries
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  if (probe) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    std::ostringstream o;
    // rho: event weight, sh: relative effect of one unit of dv
    const double rho(m_h.p_acraic->GetResult());
    o<<std::setprecision(8)<<"@@@ ACRVH acrv="<<(rcr > 0. ? r/rcr : 0.)<<" nph="<<cph.size()
     <<" nisr="<<nisr<<" rho="<<(rcr > 0. ? rho/rcr : 0.);
    for (size_t c(0); c < cph.size(); ++c) {
      const double tb(Min(cph[c].Theta(beams[0]), cph[c].Theta(beams[1])));
      std::vector<double> e(cph.size(), 0.);
      e[c] = 1.;
      o<<" | "<<(c < nisr ? "I" : "F")<<" x="<<2.*cph[c].E()/sqrt(m_h.m_s)<<" thb="<<tb
       <<" dv="<<dv[c]<<" sh="<<(rho > 0. ? m_h.p_acraic->RealVirtualRemainderRho(e)/rho : 0.);
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  return r;
}

/*
  ACRAIC: REAL_VIRTUAL 2. Each ACRAIC photon takes the dv of the NLO photon with
  the same lab momentum (1e-6 of its energy); unmatched or sub-RV_SOFT_CUT photons
  keep dv = 0 and are counted. YFS: RV_PROBE 1 prints "@@@ ACRV" per event.
*/
double YFS_Weight_Builder::AcraicRealVirtualRemainder()
{
  if (!m_h.p_acraic || !m_h.p_nlo || !m_h.p_nlo->RVRemainder()) return 0.;
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  std::vector<double> dv(cph.size(), 0.);
  const int nmatch(MatchRealVirtualDV(dv, 0));
  ++m_h.m_acrv_n;
  const double r(m_h.p_acraic->RealVirtualRemainderRho(dv));
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  // YFS: RV_PROBE_MIN: print only events with |acrv| above it, per photon x, angles, dv
  // (parsed by yfs-helpers/pion-virtual/isr-rv-vs-v2/compare.py)
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
  ACRAIC: REAL_VIRTUAL helicity. As AcraicRealVirtualRemainder with photon j's
  remainder resolved by helicity (NLO_Base::m_rvhel); photons without per-helicity
  factors keep dv_j/2 on every h.
*/
double YFS_Weight_Builder::AcraicRealVirtualHelicityRemainder()
{
  if (!m_h.p_acraic || !m_h.p_nlo || !m_h.p_nlo->RVRemainder()) return 0.;
  const Vec4D_Vector &cph(m_h.p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(m_h.p_nlo->m_rvdv);
  const std::vector<RVHelicityFactors> &nh(m_h.p_nlo->m_rvhel);
  const int nhel(Amplitude::NHel());
  std::vector<double> dv(cph.size(), 0.);
  std::vector<std::vector<Complex> > r(cph.size(),
                                       std::vector<Complex>(nhel, Complex(0., 0.)));
  bool anyhel(false);
  std::vector<long> entry;
  MatchRealVirtualDV(dv, 0, &entry);
  for (size_t a(0); a < cph.size(); ++a) {
    if (entry[a] < 0) continue;
    const size_t b(entry[a]);
    // the ISR vertex shift taken off dv is helicity-blind
    const Complex shift(0.5*(nd[b].second - dv[a]), 0.);
    const Complex half(0.5*dv[a], 0.);
    if (b >= nh.size() ||
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
  // Everything above the double real lives here: NLO_Base::CalculateNLO() is not called.
  m_h.m_ev.m_nlo_rn = 0.;
  m_h.m_ev.m_nlo_vxr = 0.;
  // providers above NLO_MAX_PHOTONS exist only for ACRAIC: TREE_BETA_MAX
  const size_t nrmax(Min(m_h.p_nlo->MaxRealPhotons(),
                         YFS::NLO_Base::RequestedMaxRealPhotons()));
  for (size_t n(3); n <= nrmax; ++n) {
    InitNLO();
    m_h.m_ev.m_nlo_rn += m_h.p_nlo->CalculateRealMultiplicity(n);
  }
  /*
    YFS: VIRTUAL_COMBINE. Where v = V/B enters the squared-level weight.
    sum: 1 + v + sum_j delta_j. product (default): (1 + v)(1 + sum_j delta_j),
    the placement ACRAIC has (virtual inside sum_h |A_1 + (v/2) A_0|^2).
    Stands aside when a real-virtual is requested (legacy RV supplies v x real),
    except the RV_MODE 1 remainder, which is defined beyond (1 + v) x real.
    v x RR and v x R_n are O(alpha^3), left out with RR_MODE singles_product.
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
    YFS: BETA0_2 - squared-level v2 = |1 + delta + d2|^2 - 1 - v, delta = v1/2
    (NLO_Base::Beta02SquaredLevel). Placed where v is: on the Born, and on the
    real where v multiplies it.
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
  delta^(2) (ACRAIC: BETA0_2) of the external virtual v' = v - dIB; the resonant
  log under RESONANCE_EXPONENTIATION greco is in the propagators.
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

// ACRAIC: FIXED_ORDER_STREAM at ORDER 2: the column truncated to O(alpha^2)
double YFS_Weight_Builder::AcraicFixedOrder2(size_t nphot, double y, double v, double d2,
                                      double rvrem) const
{
  if (nphot == 0)
    return m_h.p_acraic->GetResult0()*(1. + y + v + 0.5*y*y + y*v + 0.25*v*v + 2.*d2);
  if (nphot == 1)
    return (1. + y)*m_h.p_acraic->GetResult() + v*m_h.p_acraic->GetResult01() + rvrem;
  if (nphot == 2) return m_h.p_acraic->GetResult();
  return 0.;
}

// the same, one order up: the column's O(alpha^3) coefficient
double YFS_Weight_Builder::AcraicFixedOrder3(size_t nphot, double y, double v, double d2,
                                      double rvrem) const
{
  const double t2(AcraicFixedOrder2(nphot, y, v, d2, rvrem));
  if (nphot == 0)
    return t2 + m_h.p_acraic->GetResult0()*(y*y*y/6. + 0.5*v*y*y + (0.25*v*v + 2.*d2)*y + v*d2);
  if (nphot == 1)
    return t2 + 0.5*y*y*m_h.p_acraic->GetResult() + y*(v*m_h.p_acraic->GetResult01() + rvrem)
              + 0.25*v*v*m_h.p_acraic->GetResultVV() + 2.*d2*m_h.p_acraic->GetResult0X();
  if (nphot == 2)
    return t2 + y*m_h.p_acraic->GetResult() + v*m_h.p_acraic->GetResult01() + rvrem;
  if (nphot == 3) return m_h.p_acraic->GetResult();
  return 0.;
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
  // Born-level YFS weight (crude x form factor, no NLO correction), what YFS.LO reproduces.
  double w_lo = m_h.m_ev.m_yfsweight * m_h.m_ev.m_formfactor * (1.-m_h.m_v);
  /*
    YFS.LO without IFI_Real: with IFI_Real on the IF exponent stops at IFIOmega()
    and wif (absent from LO) restores the rest, so LO takes the IF form factor at
    sqrt(s)/2, the LO of an IFI_Real 0 run. Same conditions as the ACRAIC IF cut.
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
    IFI_Real (wif) is EEX-only and must not reach the ACRAIC weight: ACRAIC's
    coherent partition sum already is the real IFI, so wif would double count.
  */
  const double corr_eex (m_h.m_ev.m_real + (wif - 1.));

  AcraicColumns c;
  c.corr_acraic = m_h.m_ev.m_acraicfactor;            // 0 if ACRAIC produced nothing
  // ORDER 2 also keeps the O(alpha) column of the same event (ACRAIC_O1).
  // GetResultVV() is |A_v|^2, A_v the amplitude (v/2) multiplies (A_0 unless
  // ORDER2_VIRTUAL_ON_BETA1).
  c.acraico2 = (m_h.p_acraic && m_h.p_acraic->Order() == 2);
  c.corr_acraic_o0 = AcraicO0Stream() && m_h.p_acraic && c.corr_acraic != 0.
                     && m_h.p_acraic->GetRhoCrude() > 0.
                     ? m_h.p_acraic->GetResult0()/m_h.p_acraic->GetRhoCrude() : 0.;
  c.corr_acraic_o0_ifcut = c.corr_acraic_o0;
  /*
    ACRAIC: BETA0_2 - delta^(2) on A_0 (Beta0_2.H), added after the virtual and
    before the IF cut; ACRAIC_NO_B02 keeps the column without it. ACRAIC_O1 never has it.
  */
  c.acraicb02 = (m_h.p_acraic && m_h.p_nlo && m_h.p_nlo->Beta02Acraic().On());
  /*
    YFS: DYSON_VP in ACRAIC: r(h) on A_0, one insertion e(h) = d2 + r(h) with
    delta^(2); vp_* are the parts beyond beta_0^(2) (ACRAIC_NO_VP omits them).
  */
  c.acraicvp = (m_h.p_acraic && m_h.p_nlo && m_h.p_nlo->DysonVPOn() &&
               m_h.p_nlo->DysonVPFactors().ok);
  if (c.acraico2 && c.corr_acraic != 0.) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    c.corr_acraic_o1 = rcr > 0. ? m_h.p_acraic->GetResult1()/rcr : 0.;
    c.corr_acraic_o1v1 = c.corr_acraic_o1;
    c.corr_acraic_v0 = c.corr_acraic_v1 = c.corr_acraic;
  }
  // ACRAIC_Virtual external: see ApplyExternalVirtual.
  if (c.corr_acraic != 0. && (m_h.m_acraicvirtsrc == acraicvirt::external ||
                          m_h.m_acraicvirtsrc == acraicvirt::automatic) && m_h.p_acraic &&
      m_h.p_nlo && m_h.p_nlo->HasVirtual() && m_h.m_ev.m_nlo_current && m_h.m_born != 0.)
    ApplyExternalVirtual(c);
  if (m_h.p_acraic) for (double &x : m_h.p_acraic->m_lhcol) x = 0.;
  // ACRAIC_Virtual helicity: see ApplyHelicityVirtual.
  if (c.corr_acraic != 0. && m_h.m_acraicvirtsrc == acraicvirt::helicity && m_h.p_acraic &&
      m_h.p_nlo && m_h.p_nlo->HasVirtual() && m_h.m_ev.m_nlo_current && m_h.m_born != 0.)
    ApplyHelicityVirtual(c);
  /*
    IF form factor cut where ACRAIC needs it: with IFI_Real off the exponent holds
    the whole soft IF integral, but ACRAIC's coherent sum carries it above the
    generation cutoff, so the ACRAIC column swaps exp(Y_IF(IFIOmega)) for
    exp(Y_IF(omega_gen)) (1 with IFI_Real on).
  */
  FinalizeAcraicBeta02VP(c);
  ApplyIFFormFactorCut(c);
  /*
    ACRAIC: FIXED_ORDER_STREAM - ACRAIC truncated to O(alpha), as ACRAIC_FixedOrder.
    See ComputeFixedOrderColumn.
  */
  c.corr_eex = corr_eex;
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
  m_h.m_ev.m_corr_acraic_fo2 = c.corr_acraic_fo2;
  m_h.m_ev.m_corr_acraic_fo3 = c.corr_acraic_fo3;
  m_h.m_ev.m_corr_acraic_o1s = c.corr_acraic_o1s;
  m_h.m_ev.m_corr_acraic_o1v1 = c.corr_acraic_o1v1;
  m_h.m_ev.m_corr_acraic_v0 = c.corr_acraic_v0;
  m_h.m_ev.m_corr_acraic_v1 = c.corr_acraic_v1;
  m_h.m_ev.m_corr_acraic_v2 = c.corr_acraic_v2;
  m_h.m_ev.m_corr_acraic_nob02 = c.corr_acraic_nob02;
  m_h.m_ev.m_corr_acraic_novp = c.corr_acraic_novp;
  m_h.m_ev.m_corr_acraic_rvnll = c.corr_acraic != 0. ? c.corr_acraic + c.rvnll_shift : 0.;
  m_h.m_ev.m_tnp1.swap(c.tnp1);
  m_h.m_ev.m_tnp2.swap(c.tnp2);
  m_h.m_ev.m_tnpifi.swap(c.tnpifi);
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
    ACRAIC: RESONANCE_EXPONENTIATION greco: the resonant IF log is in the
    propagators, so v' = v - dI_B (Acraic_Base::BuildResonanceBornDelta; 0 with
    none). Everything below takes v'.
  */
  const double v(m_h.m_ev.m_nlo_virtual/m_h.m_born - m_h.p_acraic->ResonanceBornDeltaV() + dvisr);
  c.tnp_vhalf = 0.5*v;
  c.tnp_vif = TnpVirtualIF();
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
    c.d2_fo = d2;
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
  if ((rvmode == acraicrv::averaged || rvmode == acraicrv::helicity
       || rvmode == acraicrv::hybrid) && rcr > 0.) {
    // the provider's per-photon loop has both vertices and the
    // initial-final boxes, which a one-sided MODE leaves out
    if (m_h.p_acraic->SingleSided())
      THROW(fatal_error, "ACRAIC: REAL_VIRTUAL averaged takes the loop "
            "provider's real-virtual, which contains both vertices and the "
            "initial-final loops; YFS: MODE ISR or FSR (one side only) needs ACRAIC: "
            "{REAL_VIRTUAL: factorisable} with ACRAIC's own virtual "
            "(YFS: ACRAIC_Virtual: acraic).");
    // RV_FRACTION f < 1: exact loops on a fraction f, isr_nll as control variate
    const double rvfrac(NLO_Base::RealVirtualFraction());
    double arv(0.);
    if (m_h.p_nlo && m_h.p_nlo->RealVirtualSkipped()) {
      arv = AcraicIsrNllRemainder()/rcr;
    }
    else {
      arv = (rvmode == acraicrv::helicity ? AcraicRealVirtualHelicityRemainder()
             : rvmode == acraicrv::hybrid ? AcraicHybridRemainder()
                                          : AcraicRealVirtualRemainder())/rcr;
      if (rvfrac < 1.) {
        const double nll(AcraicIsrNllRemainder()/rcr);
        arv = nll + (arv - nll)/rvfrac;
      }
    }
    if (!IsBad(arv)) { c.corr_acraic += arv; c.rvrem_fo = arv*rcr; }
    // ACRAIC_RVNLL column: the same event with the isr_nll remainder
    if (rvmode != acraicrv::helicity) {
      const double nll(AcraicIsrNllRemainder()/rcr);
      if (!IsBad(arv) && !IsBad(nll)) c.rvnll_shift = nll - arv;
    }
  }
  if (rvmode == acraicrv::isr_nll && rcr > 0.) {
    const double arv(AcraicIsrNllRemainder()/rcr);
    if (!IsBad(arv)) { c.corr_acraic += arv; c.rvrem_fo = arv*rcr; }
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
  const double rcr(m_h.p_acraic->GetRhoCrude());
  const double vfull(m_h.m_ev.m_nlo_virtual/m_h.m_born);
  std::vector<Complex> d;
  bool ok(m_h.p_nlo->m_lhelok &&
          m_h.p_acraic->LoopHelicityToAcraic(m_h.p_nlo->m_lhel, m_h.p_nlo->m_la0,
                                     m_h.p_nlo->m_ldelta, m_h.p_nlo->m_lhelmom,
                                     Complex(0.5*vfull, 0.), d));
  if (!ok) d.assign(Amplitude::NHel(), Complex(0.5*vfull, 0.));
  /*
    RESONANCE_EXPONENTIATION greco: the resonant log leaves delta^h per helicity
    (Acraic_Base::BuildResonanceBornDelta); v' = v - dI_B is its Born-weighted 2 Re.
  */
  const double v(vfull - m_h.p_acraic->ResonanceBornDeltaV());
  c.tnp_vhalf = 0.5*v;
  c.tnp_vif = TnpVirtualIF();
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
    ACRAIC: BETA0_2 per helicity: delta^h -> delta^h + delta^(2)(delta^h), complex,
    same phase convention as delta_f (Beta0_2.H). Events without helicity
    factors use the external, real form on v/2.
  */
  std::vector<Complex> db(d), dhb(dh);
  const bool hb02(c.acraicb02 && m_h.p_nlo->Beta02Ok() && rcr > 0.);
  if (hb02) {
    const double d2ext(AcraicBeta02External(*m_h.p_nlo, m_h.p_acraic->ResonanceBornDeltaV()));
    for (size_t f(0); f < db.size(); ++f)
      db[f] += ok ? Beta02Delta(m_h.p_nlo->Beta02Factors(), d[f], m_h.p_nlo->Beta02Acraic())
                  : Complex(d2ext, 0.);
    // BETA0_2 two_loop: photon-exchange two-loop boxes per helicity
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
    ACRAIC: EW_HIGHER_ORDER_HELICITY - EW_HIGHER_ORDER is otherwise missing from
    delta^h: averaged puts v_HO/2 on every helicity, helicity uses the provider's
    per-helicity factor. Applied after delta^(2).
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
    Columns ACRAIC_EXTV, ACRAIC_L/R and ACRAIC_EXTV_L/R (beam 1 fully left/right
    polarised, beam 2 unpolarised); absolute factors before the IF cut.
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
    ORDER 2: helicity analogues of the external compositions, delta^h for v/2:
    ACRAIC_O1 |A_1 + d A_0|^2, O1V1 |A_1 + d A_1|^2, O2V0/V1 |A_2 + d A_0|^2 and
    |A_2 + d A_1|^2, O2V2 |A_2 + d A_2|^2.
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
  // ACRAIC: TNP_IFI: [0] O(alpha), [1] ORDER 2, [2] the odd stream
  if (m_h.p_acraic && m_h.p_acraic->TnpIfiOn() && c.corr_acraic != 0. && m_h.p_acraic->GetRhoCrude() > 0.) {
    const double rcr(m_h.p_acraic->GetRhoCrude());
    double yif(0.);
    if (m_h.m_ifisub && m_h.m_fullform >= 1 && m_h.p_dipoles) {
      static const bool ifcutyif(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                 ["IF_FORMFACTOR_CUT"].SetDefault(true).Get<bool>());
      yif = ifcutyif ? m_h.p_dipoles->FormFactorSumIF(0.5*sqrt(m_h.m_s)*m_h.m_isrcut)
                     : m_h.p_dipoles->FormFactorSumIF();
    }
    const double eflip(IsBad(yif) ? 1. : exp(-2.*yif));
    double t1(0.), t2(0.), o1(0.);
    // [0], [1]: T_IFI at O(alpha), ORDER 2; [2]: Odd[w_1] (ACRAIC: TNP_ODD_STREAM)
    c.tnpifi.assign(3, 0.);
    if (m_h.p_acraic->TnpIfiTemplate(1, c.tnp_vhalf, c.tnp_vif, eflip, t1, &o1)) {
      c.tnpifi[0] = t1/rcr;
      c.tnpifi[2] = o1/rcr;
    }
    if (c.acraico2 && m_h.p_acraic->TnpIfiTemplate(2, c.tnp_vhalf, c.tnp_vif, eflip, t2))
      c.tnpifi[1] = t2/rcr;
  }
}

void YFS_Weight_Builder::ApplyIFFormFactorCut(AcraicColumns &c)
{
  if (c.corr_acraic != 0. && m_h.m_ifisub && m_h.m_fullform >= 1 &&
      m_h.p_dipoles) {
    const double wgen(0.5*sqrt(m_h.m_s)*m_h.m_isrcut);
    const double dy(m_h.p_dipoles->FormFactorSumIF(wgen)
                    - m_h.p_dipoles->FormFactorSumIF());
    if (!IsBad(dy)) { c.dyif_fo = dy; c.corr_acraic_o1s *= exp(dy);
                      c.corr_acraic *= exp(dy); c.corr_acraic_o1 *= exp(dy); c.corr_acraic_o1v1 *= exp(dy);
                      c.corr_acraic_o0_ifcut *= exp(dy);
                      c.corr_acraic_v0 *= exp(dy); c.corr_acraic_v1 *= exp(dy); c.corr_acraic_v2 *= exp(dy);
                      c.corr_acraic_nob02 *= exp(dy); c.corr_acraic_novp *= exp(dy);
                      c.rvnll_shift *= exp(dy);
                      for (double &x : c.tnp1) x *= exp(dy);
                      for (double &x : c.tnp2) x *= exp(dy);
                      for (double &x : c.tnpifi) x *= exp(dy); }
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
    // ORDER 2 needs the external virtual
    static const bool kkflux(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                             ["KKMC_FLUX_EMULATION"].Get<int>() != 0);
    if (c.acraico2 && (m_h.m_acraicvirtsrc == acraicvirt::acraic || c.acraicvp || kkflux))
      THROW(not_implemented, "ACRAIC: FIXED_ORDER_STREAM at ORDER 2 needs an external "
            "virtual (YFS: ACRAIC_Virtual external), no YFS: DYSON_VP and no "
            "ACRAIC: KKMC_FLUX_EMULATION.");
    const size_t nphot(m_h.p_acraic->NPhot());
    double yfo(m_h.m_ev.m_formfactor_sum + c.dyif_fo
               + m_h.p_dipoles->FormFactorShiftToRegion(0.5*sqrt(m_h.m_s)*m_h.m_isrcut));
    // ACRAIC: FIXED_ORDER_EVENT_EXPONENT: one-photon events divided by their own exponent,
    static const bool foevtexp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                               ["FIXED_ORDER_EVENT_EXPONENT"].SetDefault(false).Get<bool>());
    const double yfo_born(yfo);
    // only for a single FSR photon: its exponent replaces the final-state part of the Born-point shift
    if (foevtexp && nphot == 1 && m_h.p_fsr && m_h.p_fsr->YFSIR() > 0. && m_h.m_ev.m_FSRPhotons.size() == 1)
      yfo = yfo_born - m_h.p_dipoles->FormFactorShiftToRegionFinal(0.5*sqrt(m_h.m_s)*m_h.m_isrcut)
            + log(m_h.p_fsr->YFSIR());
    double rhofo(0.);
    if (nphot == 0)
      rhofo = m_h.p_acraic->GetResult1() + (c.vext_fo + yfo)*m_h.p_acraic->GetResult0();
    else if (nphot == 1)
      rhofo = m_h.p_acraic->GetResult1();
    c.corr_acraic_fo = exp(c.dyif_fo - yfo)*rhofo/m_h.p_acraic->GetRhoCrude();
    if (IsBad(c.corr_acraic_fo)) c.corr_acraic_fo = 0.;
    ProbeN1(c, nphot, yfo);
    if (c.acraico2) {
      const double norm(exp(c.dyif_fo - yfo_born)/m_h.p_acraic->GetRhoCrude());
      c.corr_acraic_fo2 = AcraicFixedOrder2(nphot, yfo_born, c.vext_fo, c.d2_fo, c.rvrem_fo)*norm;
      c.corr_acraic_fo3 = AcraicFixedOrder3(nphot, yfo_born, c.vext_fo, c.d2_fo, c.rvrem_fo)*norm;
    }
    if (IsBad(c.corr_acraic_fo2)) c.corr_acraic_fo2 = 0.;
    if (IsBad(c.corr_acraic_fo3)) c.corr_acraic_fo3 = 0.;
  }
}

// YFS: N1_PROBE: per-event dump of the one-photon weight pieces
void YFS_Weight_Builder::ProbeN1(const AcraicColumns &c, size_t nphot, double yfo)
{
  static const int n1probe(ATOOLS::Settings::GetMainSettings()["YFS"]["N1_PROBE"].SetDefault(0).Get<int>());
  static const double n1emin(ATOOLS::Settings::GetMainSettings()["YFS"]["N1_PROBE_EMIN"].SetDefault(0.02).Get<double>());
  static int nprobe(0);
  if (!(n1probe > 0 && nphot == 1 && nprobe < n1probe && m_h.p_fsr && m_h.m_mode == yfsmode::fsr)) return;
  const bool acraic_nom_probe(m_h.m_acraic_weight && c.corr_acraic != 0.);
  const std::string pl(m_h.p_acraic->ProbeN1Line());
  const size_t kp(pl.find(" k0="));
  double e0(0.);
  if (kp != std::string::npos) e0 = atof(pl.c_str() + kp + 4);
  if (e0 <= n1emin) return;
  ++nprobe;
  std::cerr << std::setprecision(17) << "@@@ N1PROBE" << pl << m_h.p_fsr->ProbeLine()
            << " yfo=" << yfo << " dyif=" << c.dyif_fo << " expfo=" << exp(c.dyif_fo - yfo)
            << " vext=" << c.vext_fo << " corrfo=" << c.corr_acraic_fo
            << " mformfactor=" << m_h.m_ev.m_formfactor << " yfsweight_in=" << m_h.m_ev.m_yfsweight
            << " mv=" << m_h.m_v << " ffsum=" << m_h.m_ev.m_formfactor_sum
            << " shift=" << m_h.p_dipoles->FormFactorShiftToRegion(0.5*sqrt(m_h.m_s)*m_h.m_isrcut)
            << " corrnom=" << c.corr_acraic << " mcorrnom=" << (acraic_nom_probe ? c.corr_acraic : c.corr_eex)
            << " fsrW=" << m_h.m_fsrWeight << std::endl;
}
// The YFS.* named weights: ratios turning the nominal weight into truncated or
// reordered matchings.
// w_lo   : Born-level YFS weight, before any NLO correction
// w_full : weight the event carries, before the clamps below
void YFS_Weight_Builder::BuildNamedWeights(double w_lo, double w_full) {
  m_h.m_ev.m_nlo_weightsmap = Weights_Map{1.0};
  // YFS: NLO_Weight_Breakdown false: the nominal weight only, no YFS.* columns
  if (!m_h.m_nlo_weight_breakdown) {
    RefuseStreamsWithoutBreakdown();
    return;
  }
  Weights wyfs{1.0};
  bool any(false);

  std::vector<std::string> names, tnpnames;
  TnpIfiColumns tnpifi;
  RegisterColumnNames(names, tnpnames, tnpifi, wyfs);
  if (!names.empty()) any = true;

  FillAcraicColumns(wyfs);
  FillTnpColumns(wyfs, tnpnames);
  FillTnpIfiColumns(wyfs, tnpifi);
  FillDiffBudgetColumns(wyfs);
  FillWithoutColumns(wyfs);
  if (FillNLOBreakdownColumns(wyfs, w_lo, w_full)) any = true;

  if (BuildLadderWeights(wyfs)) any = true;

  if (any) m_h.m_ev.m_nlo_weightsmap["YFS"] = wyfs;
}

//! A YFS.* stream asked for by name is not dropped silently when
//! YFS: NLO_Weight_Breakdown is false: the run stops instead.
void YFS_Weight_Builder::RefuseStreamsWithoutBreakdown() const
{
  static bool checked(false);
  if (checked) return;
  checked = true;
  std::vector<std::string> asked;
  if (m_h.m_nlo_columns) asked.push_back("YFS: NLO_Columns");
  if (m_h.m_ladder_weights) asked.push_back("YFS: Ladder_Weights");
  if (DiffBudgetStreamsOn()) asked.push_back("YFS: DIFF_BUDGET_STREAMS");
  if (m_h.m_useacraic) {
    if (AcraicO0Stream()) asked.push_back("ACRAIC: WEIGHT_O0");
    if (AcraicFixedOrderStream()) asked.push_back("ACRAIC: FIXED_ORDER_STREAM");
    if (AcraicStrictO1Stream()) asked.push_back("ACRAIC: STRICT_O1_STREAM");
    if (m_h.p_acraic && m_h.p_acraic->TnpOn()) asked.push_back("ACRAIC: TNP_STREAMS");
    if (AcraicTnpBandStreams()) asked.push_back("ACRAIC: TNP_BAND_STREAMS");
  }
  if (asked.empty()) return;
  std::string list;
  for (const std::string &a : asked) list += (list.empty() ? "" : ", ") + a;
  THROW(fatal_error, list + " write YFS.* named weights, which need "
        "YFS: NLO_Weight_Breakdown true.");
}

//! The YFS.* column names this configuration emits. Fixed once per event
//! from configuration alone, and the only names emit() will write -- a
//! per-event code path must not be able to conjure a column that other
//! events lack, since HepMC binds the set from event one.
void YFS_Weight_Builder::RegisterColumnNames(std::vector<std::string> &names,
                                             std::vector<std::string> &tnpnames,
                                             TnpIfiColumns &tnpifi,
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
    if (hasnlo)
      for (const char *n : {"Real","Virtual","BR","BV","NLO_1g",
                            "NLO_2g"})
        names.push_back(n);
    // EEX only when the beta expansion ran (m_eex is set only for m_betaorder > 0)
    if (hasnlo && m_h.m_betaorder > 0) names.push_back("EEX");
    if (hasnnlo)
      for (const char *n : {"RealVirtual","RealReal","NLO+RR","NLO+RV",
                            "NNLO_1g","NNLO_2g"})
        names.push_back(n);
    // the EEX estimate of the double virtual, retired by YFS: BETA0_2
    // (its LL part, without the NLL and the 100% envelope); never both
    if (hasnnlo && !m_h.p_nlo->Beta02Nlo().On())
      for (const char *n : {"VV_EEX","NNLO_VV","NNLO_VV_up","NNLO_VV_down"})
        names.push_back(n);
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
  if (m_h.m_useacraic && m_h.m_acraicvirtsrc != acraicvirt::acraic &&
      (YFS::Acraic_Base::RealVirtualMode() == acraicrv::hybrid ||
       YFS::Acraic_Base::RealVirtualMode() == acraicrv::averaged))
    names.push_back("ACRAIC_RVNLL");
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
  if (m_h.m_useacraic && AcraicFixedOrderStream()) {
    names.push_back("ACRAIC_FixedOrder");
    if (ATOOLS::Settings::GetMainSettings()["ACRAIC"]["ORDER"].SetDefault(1).Get<int>() == 2) {
      names.push_back("ACRAIC_FixedOrder2");
      names.push_back("ACRAIC_FixedOrder3");
    }
  }
  if (m_h.m_useacraic && AcraicStrictO1Stream()) names.push_back("ACRAIC_O1S");
  // ACRAIC: ORDER 2 also carries the O(alpha^1) ACRAIC column of each event
  static const bool acraicorder2(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["ORDER"]
                                 .SetDefault(1).Get<int>() == 2);
  if (m_h.m_useacraic && acraicorder2) {
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
  // ACRAIC: TNP_IFI - the initial-final interference template(s), each name
  // with its index in m_ev.m_tnpifi
  if (m_h.m_useacraic && m_h.p_acraic && m_h.p_acraic->TnpIfiOn()) {
    tnpifi.push_back({"ACRAIC_TNP1_IFI", 0});
    if (m_h.p_acraic->Order() == 2) tnpifi.push_back({"ACRAIC_TNP2_IFI", 1});
    // ACRAIC: TNP_ODD_STREAM: Odd[w_1], a low-variance F - B for A_FB
    static const bool oddstream(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["TNP_ODD_STREAM"].SetDefault(false).Get<bool>());
    if (oddstream) tnpifi.push_back({"ACRAIC_TNP1_ODD", 2});
    for (const auto &t : tnpifi) names.push_back(t.first);
  }
  if (m_h.m_ladder_weights) {
    if (m_h.m_coulomb && m_h.p_coulomb) names.push_back("NoCoulomb");
    if (m_h.m_ifisub && m_h.m_fullform >= 1 && m_h.p_dipoles)
      names.push_back("NoIFI");
  }
  for (const std::string &n : names) wyfs[n] = 1.;
  for (const std::string &n : tnpnames) wyfs[n] = 0.;
  for (const auto &t : tnpifi) wyfs[t.first] = 0.;
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
    if (m_h.m_wnames.count("ACRAIC_FixedOrder2")) {
      const double rfo2(m_h.m_ev.m_corr_acraic_fo2/m_h.m_ev.m_corr_nominal);
      if (!IsBad(rfo2)) wyfs["ACRAIC_FixedOrder2"] = rfo2;
    }
    if (m_h.m_wnames.count("ACRAIC_FixedOrder3")) {
      const double rfo3(m_h.m_ev.m_corr_acraic_fo3/m_h.m_ev.m_corr_nominal);
      if (!IsBad(rfo3)) wyfs["ACRAIC_FixedOrder3"] = rfo3;
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
  YFS: DIFF_BUDGET_STREAMS (diagnostic): the pieces of the two O(alpha)
  factors, ACRAIC C = a1 (1 + v'/2)^2 and YFS.NLO_EW N = (1 + v)(1 + r),
  v' = v - x (x = dI_B under RESONANCE_EXPONENTIATION greco, else 0), as
  columns DB_<q> = w_crude q. VSQ = a1 v^2/4, REST = a1 (1 + v) - N and
  NPV = N + VSQ are true weights. X1, X2, NPH, SXI, SXF are photon payloads
  (x = 2E/sqrt(s), + ISR / - FSR), not weights. 0 without an ACRAIC factor.
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
      if (m_h.m_ifisub && m_h.m_fullform >= 1 && m_h.p_dipoles) {
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
  if (m_h.m_wnames.count("ACRAIC_RVNLL") && m_h.m_ev.m_corr_acraic != 0. &&
      !IsZero(m_h.m_ev.m_corr_nominal)) {
    const double r(m_h.m_ev.m_corr_acraic_rvnll/m_h.m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC_RVNLL"] = r;
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
    The breakdown columns are x/m_real; with ACRAIC as nominal the denominator is
    rescaled by corr_nominal/corr_eex so each column keeps its absolute weight.
  */
  const double nomden(IsZero(m_h.m_ev.m_corr_eex) ? m_h.m_ev.m_real
                      : m_h.m_ev.m_real*(m_h.m_ev.m_corr_nominal/m_h.m_ev.m_corr_eex));

  // operator[] would create a column: refuse unregistered names, and leave
  // non-finite values at the registered 1.0.
  auto emit = [this, &wyfs](const std::string &name, double value) {
    if (!m_h.m_wnames.count(name)) return;
    if (!IsBad(value)) wyfs[name] = value;
  };

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
      // LO = w_lo/w_full, so nominal * YFS.LO == w_lo.
      if (!IsZero(w_full)) emit("LO", w_lo/w_full);
      emit("EEX", ratio(m_h.m_ev.m_eex + ifisoft/m_h.m_born, nomden));
      // Matching truncated to the 1 or 2 hardest photons; Virtual and the IF
      // restoration (part of the exponent) stay full. NLO_1g is NLO without m_nlo_vxr.
      const double nlo_1g = (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_hardest
                             + m_h.m_ev.m_nlo_virtual)/m_h.m_born;
      const double nlo_2g = (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_2hardest
                             + m_h.m_ev.m_nlo_virtual)/m_h.m_born;
      emit("NLO_1g", ratio(nlo_1g, nomden));
      emit("NLO_2g", ratio(nlo_2g, nomden));
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
      NLO+RV carries m_nlo_vxr like YFS.NLO and NLO+RR (0 with the legacy RV).
    */
    const double RV_total = (m_h.m_born + m_h.m_ev.m_nlo_real + m_h.m_ev.m_nlo_virtual + m_h.m_ev.m_nlo_rv
                             + m_h.m_ev.m_nlo_vxr)/m_h.m_born;
    // Every column here is x/m_ev.m_real.
    if (!IsZero(m_h.m_ev.m_real)) {
      emit("RealVirtual", ratio((m_h.m_ev.m_nlo_rv)/m_h.m_born, nomden));
      emit("RealReal", ratio((m_h.m_ev.m_nlo_rr)/m_h.m_born, nomden));
      emit("NLO+RR", ratio(RR_total, nomden));
      emit("NLO+RV", ratio(RV_total, nomden));
      emit("NNLO", ratio(nnlo_total, nomden));
      // NNLO analogue of NLO_1g/NLO_2g: Real + RealVirtual on the 1 or 2 hardest
      // photons, RealReal only for the pair; Virtual full.
      const double nnlo_1g =
          (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_hardest + m_h.m_ev.m_nlo_virtual
           + m_h.m_ev.m_nlo_rv_hardest)/m_h.m_born;
      const double nnlo_2g =
          (m_h.m_born + ifisoft + m_h.m_ev.m_nlo_real_2hardest + m_h.m_ev.m_nlo_virtual
           + m_h.m_ev.m_nlo_rv_2hardest + m_h.m_ev.m_nlo_rr_2hardest)/m_h.m_born;
      emit("NNLO_1g", ratio(nnlo_1g, nomden));
      emit("NNLO_2g", ratio(nnlo_2g, nomden));

      // Approximate double virtual: the O(alpha^2) term of the EEX virtual series,
      //     VV_EEX = prod(1 + 0.5g + 0.125g^2) - prod(1 + 0.5g)  (II and FF dipoles),
      // a leading-log estimate, not a calculation. Shipped as an envelope
      // 1 +/- VV_Approx_Uncertainty; the nominal NNLO stays VV-free, NNLO_VV includes it.
      // p_nlo's flag: the handler's own m_vvtool is a YFS_Base member that
      // nothing on the handler side ever sets (see YFS_Base.H)
      if (!m_h.p_nlo->HasVV()) {
        const double vv = m_h.p_dipoles->CalculateEEXVirtual(2)
                        - m_h.p_dipoles->CalculateEEXVirtual(1);
        if (!IsBad(vv)) {
          // Flat envelope (1 +/- d) x estimate: d > 1 would flip the sign of VV.
          // Not the dominant NNLO uncertainty (VV_EEX is far below RR).
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
void YFS_Weight_Builder::FillTnpIfiColumns(ATOOLS::Weights &wyfs, const TnpIfiColumns &tnpifi)
{
  if (IsZero(m_h.m_ev.m_corr_nominal)) return;
  for (const auto &t : tnpifi)
    if (t.second < m_h.m_ev.m_tnpifi.size() && !IsBad(m_h.m_ev.m_tnpifi[t.second]))
      wyfs[t.first] = m_h.m_ev.m_tnpifi[t.second]/m_h.m_ev.m_corr_nominal;
}
bool YFS_Weight_Builder::BuildLadderWeights(ATOOLS::Weights &w) {
  if (!m_h.m_ladder_weights) return false;
  bool any(false);

  if (m_h.m_coulomb && m_h.p_coulomb) {
    const double wc(m_h.p_coulomb->GetWeight());
    if (!IsZero(wc) && !IsBad(wc)) { w["NoCoulomb"] = 1./wc; any = true; }
  }

  if (m_h.m_ifisub && m_h.m_fullform >= 1 && m_h.p_dipoles) {
    const double fif(m_h.p_dipoles->FormFactorSumIF());
    if (!IsBad(fif)) { w["NoIFI"] = exp(-fif); any = true; }
  }

  return any;
}
