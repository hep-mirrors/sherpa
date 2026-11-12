#include "ATOOLS/Org/Message.H"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include "YFS/Main/YFS_Handler.H"
#include "YFS/Main/YFS_Weight_Builder.H"
#include "BEAM/Main/Beam_Base.H"
#include "YFS/Main/ISR.H"

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;


YFS_Handler::YFS_Handler()
{
  if(Mode()!=YFS::yfsmode::off){
    p_dipoles = std::make_unique<Define_Dipoles>();
    p_weights = std::make_unique<YFS_Weight_Builder>(*this);
    p_coulomb = std::make_unique<Coulomb>(m_coulomb);
    p_fsr = std::make_unique<FSR>();
    p_debug = std::make_unique<Debug>();
    p_yfsFormFact = std::make_unique<YFS::YFS_Form_Factor>();
    m_setparticles = false;
    p_isr = std::make_unique<YFS::ISR>();
    m_isrinital = true;
    if (m_photon_split != 0)
      THROW(not_implemented, "YFS: PHOTON_SPLITTER_MODE is not supported with"
            " YFS: its photons are in the YFS_Initial and QED_Radiation blobs,"
            " which the photon splitter does not handle.");
    RefuseUnweightedWithoutIntegratedNLO();
    m_rmode = 0;
    m_negskip = 0;
    rpa->gen.AddCitation(1,"The automation of YFS ISR is published in  \\cite{Krauss:2022ajk}. Which is based on \\cite{Jadach:1988gb}");
  }
}

void YFS_Handler::RefuseUnweightedWithoutIntegratedNLO() const
{
  const int weighted(0);   // EVENT_GENERATION_MODE: Weighted
  if (m_int_nlo) return;
  if (ToType<int>(rpa->gen.Variable("EVENT_GENERATION_MODE")) == weighted) return;
  THROW(fatal_error, "YFS: Integrate_NLO false needs EVENT_GENERATION_MODE: "
        "Weighted. Unweighted and PartiallyUnweighted events are unweighted "
        "against an integration without the beta, ACRAIC or YFS.NLO "
        "correction. Set YFS: Integrate_NLO: true, or generate weighted events.");
}

YFS_Handler::~YFS_Handler()
{
  if(Mode()!=YFS::yfsmode::off){
    if (m_acrv_n > 0)
      msg_Out()<<"ACRAIC: REAL_VIRTUAL 2 on "<<m_acrv_n<<" events (this rank): NLO RV photons "
               <<"without a ACRAIC match "<<m_acrv_nomatch<<", matched but without M_1 "
               <<m_acrv_nom1<<", matched to an initial-state photon (hybrid: not used) "
               <<m_acrv_isr<<std::endl;
    const Acraic_Stats &cs(m_acraicstats);
    if (cs.m_cmp_n > 0)
      msg_Out()<<"YFS ACRAIC_Compare over "<<cs.m_cmp_n<<" points, ACRAIC vs the "
               <<"external providers (relative difference, mean / worst):\n"
               <<"    virtual  "<<cs.m_vsum/cs.m_cmp_n<<" / "<<cs.m_vworst<<"\n"
               <<"    real     "<<cs.m_rsum/cs.m_cmp_n<<" / "<<cs.m_rworst<<"\n"
               <<"    total    "<<cs.m_tsum/cs.m_cmp_n<<" / "<<cs.m_tworst<<std::endl;
    if (!m_fsrfail.empty()) {
      msg_Out()<<"YFS FSR failure statistics (YFS: FSR_FAILURE_STATS), "
               <<m_fsrtrials<<" CalculateFSR trials:\n";
      for (const auto &f : m_fsrfail)
        msg_Out()<<"    "<<f.second<<"  "<<f.first<<"\n";
    }
    if (cs.m_oen > 0)
      msg_Debugging()<<"YFS: ACRAIC supplied the O(alpha) weight on "
               <<cs.m_oen<<" events; mean (ACRAIC factor)/(EEX factor) = "
               <<cs.m_oesum/cs.m_oen
               <<(cs.m_bad ? " ("+ATOOLS::ToString(cs.m_bad)
                             +" events fell back to EEX)" : "")
               <<std::endl;
    if(m_negskip!=0){
      msg_Debugging()<<"Total Events Skipped: "<<m_negskip<<std::endl;
    }
    if(m_regcheck_zero > 0)
      msg_Out()<<"YFS: SOFT_REGION_CHECK (this rank): "<<m_regcheck_zero
               <<" events skipped, ISR photon below K_min (weight zero)"<<std::endl;
    for (int w(0); w < 2; ++w)
      if(m_regcheck_n[w] > 0)
        msg_Out()<<"YFS: SOFT_REGION_CHECK (this rank, "<<(w ? "FSR" : "ISR")<<"): "
                 <<m_regcheck_bad[w]<<" of "<<m_regcheck_n[w]
                 <<" photons handed to ACRAIC lie inside the common soft region (E_P < IR_CUTOFF/2);"
                 <<" smallest E_P/Emin = "<<m_regcheck_min[w]<<std::endl;
    if(m_nonfinite!=0)
      msg_Debugging()<<"YFS (this rank): "<<m_nonfinite<<" events with a non-finite weight, "
               <<"set to zero"<<std::endl;
    if(m_ifireal && p_dipoles && p_dipoles->IFIClipped()!=0){
      msg_Debugging()<<"IFI real reweight: "<<p_dipoles->IFIClipped()
               <<" photon factors clamped to ["<<p_dipoles->IFIRClip()
               <<", "<<1./p_dipoles->IFIRClip()<<"]"<<std::endl;
    }
    if(m_ifireal && p_nlo && p_nlo->m_ifi_n>0){
      const double mean = p_nlo->m_ifi_sum/p_nlo->m_ifi_n;
      const double var  = p_nlo->m_ifi_sum2/p_nlo->m_ifi_n - mean*mean;
      msg_Debugging()<<"IFI real restoration: n="<<p_nlo->m_ifi_n
               <<"  subloc/subb mean="<<mean
               <<" rms="<<(var>0.?sqrt(var):0.)
               <<" min="<<p_nlo->m_ifi_min
               <<" max="<<p_nlo->m_ifi_max<<std::endl;
      msg_Debugging()<<"  profiled in x = E_gamma/sqrt(s):"<<std::endl;
      for(int i=0;i<5;++i){
        if(p_nlo->m_ifi_x_n[i]==0) continue;
        const double rr = p_nlo->m_ifi_x_r[i]/p_nlo->m_ifi_x_n[i];
        const double ee = p_nlo->m_ifi_x_e[i]/p_nlo->m_ifi_x_n[i];
        msg_Debugging()<<"    x="<<0.1*i<<"-"<<0.1*(i+1)
                 <<"  n="<<p_nlo->m_ifi_x_n[i]
                 <<"  applied="<<rr<<"  exact="<<ee
                 <<"  residue="<<(rr-ee)<<std::endl;
      }
    }
    if(m_ifireal && p_nlo && p_nlo->m_ifi_eligible>0){
      msg_Debugging()<<"IFI real restoration skipped on "
               <<(p_nlo->m_ifi_skip_emin+p_nlo->m_ifi_skip_momcons
                  +p_nlo->m_ifi_skip_unusable+p_nlo->m_ifi_skip_other)
               <<" of "<<p_nlo->m_ifi_eligible<<" photons above omega:"
               <<" below the NLO photon energy cut "<<p_nlo->m_ifi_skip_emin
               <<", momentum conservation "<<p_nlo->m_ifi_skip_momcons
               <<", zero or bad real ME "<<p_nlo->m_ifi_skip_unusable
               <<", ratio refused "<<p_nlo->m_ifi_skip_other<<std::endl;
    }
  }
}

NLO_Base *YFS_Handler::EnsureNLO()
{
  if (!p_nlo) p_nlo = std::make_unique<YFS::NLO_Base>();
  if (m_acraicborn) p_nlo->SetBornProc(m_acraicborn);
  return p_nlo.get();
}


void YFS_Handler::SetBeam(BEAM::Beam_Spectra_Handler *beam)
{
  p_beams = beam;
  m_beam1 = p_beams->GetBeam(0)->OutMomentum();
  m_beam2 = p_beams->GetBeam(1)->OutMomentum();
  if(m_beam1 != -m_beam2) m_asymbeams = true;
  else m_asymbeams = false;
}

void YFS_Handler::SetLimits(const double &smin, const double &s) {
  // s must be this event's actual s' (from the real, possibly beam-spread-
  // sampled, incoming momenta), not sqr(rpa->gen.Ecms())
  p_yfsFormFact->SetCharge(1);
  p_coulomb->SetAlphaQED(m_alpha);
  double maxV = 1. - smin / s;
  if (m_vmax > maxV && !m_asymbeams) {
    msg_Error() << "Warning: vmax to large in YFS integration reseting to " << maxV << std::endl;
    m_vmax = maxV;
  }
}

void YFS_Handler::SetChannelFlavours(const ATOOLS::Flavour_Vector &flavs) {
  if (m_flavs == flavs) return;
  if (flavs.size() != m_flavs.size())
    THROW(fatal_error, "YFS: a channel must have the same legs as its group.");
  if (m_mode != yfsmode::isr)
    THROW(not_implemented, "YFS: a container final state (e.g. 93 93) only "
          "works with MODE: ISR, the final-state masses enter the photon "
          "kinematics otherwise. List the processes one by one.");
  for (size_t i(0); i < 2; ++i)
    if (flavs[i] != m_flavs[i])
      THROW(fatal_error, "YFS: a channel cannot change the initial state.");
  for (size_t i(2); i < flavs.size(); ++i) {
    m_flavs[i] = flavs[i];
    m_mass[i] = flavs[i].Mass();
  }
  // the channel's Born and real processes, which the YFS.NLO providers need
  // with or without ACRAIC (InitializeACRAIC sets them as well)
  if (m_useacraic) InitializeACRAIC(m_flavs);
  else SetAcraicChannelProcs(m_flavs);
}

void YFS_Handler::SetFlavours(const ATOOLS::Flavour_Vector &flavs) {
  if(m_setparticles && m_flavs == flavs) return;
  if (m_mode_from_acraic) {
    bool charged(false);
    for (size_t i(2); i < flavs.size(); ++i) if (flavs[i].Charge() != 0.) charged = true;
    m_mode = charged ? yfsmode::isrfsr : yfsmode::isr;
  }
  m_particle_store.clear();
  m_particles.clear();
  m_flavs.clear();
  m_mass.clear();
  bool qed(false);
  for(size_t i = 0; i < flavs.size(); ++i) {
    m_flavs.push_back(flavs[i]);
    if (i < 2) {
      if (m_flavs[i].Mass() == 0 && m_mode!=yfsmode::fsr) {
        THROW(fatal_error, "Inital states must be massive for YFS");
      }
    }
    m_mass.push_back(m_flavs[i].Mass());
      m_particle_store.push_back(std::make_unique<ATOOLS::Particle>(
          i, m_flavs[i], ATOOLS::Vec4D{0, 0, 0, 0}, i < 2 ? 'i' : 'f'));
      m_particles.push_back(m_particle_store.back().get());
      m_particles[i]->ResetCounter();
    if (i >= 2) {
      if (flavs[i].IsQED()) qed = true;
    }
  }
  m_setparticles = true;
  if (m_useacraic && !m_hidephotons)
    THROW(fatal_error, "YFS: HIDE_PHOTONS false is incompatible with ACRAIC: the "
          "coherent sum needs one common soft region for all photons.");
  // a container (e.g. 93 93) is no process: ACRAIC is built per channel by
  // SetChannelFlavours
  bool container(false);
  for (size_t i(2); i < m_flavs.size(); ++i)
    if (m_flavs[i].Kfcode() >= 90 && m_flavs[i].Kfcode() <= 99) container = true;
  if (m_useacraic && !container) InitializeACRAIC(m_flavs);
}

void YFS_Handler::SetBornMomenta(const ATOOLS::Vec4D_Vector &p) {
  StartEvent();
  m_ev.m_bornMomenta.clear();
  for(size_t i = 0; i < p.size(); ++i) {
    m_ev.m_bornMomenta.push_back(p[i]);
  }
  // detect asymmetric beams from the original lab momenta, before any boost
  if(m_ev.m_bornMomenta[0] != -m_ev.m_bornMomenta[1]) m_asymbeams = true;
  else m_asymbeams = false;
  Vec4D Q(m_ev.m_bornMomenta[0] + m_ev.m_bornMomenta[1]);
  if (m_mode != yfsmode::fsr && !IsZero(Q.PSpat() / Q[0], 1e-10)) {
    m_ev.m_cmsboost = Poincare(Q);
    for (size_t i = 0; i < m_ev.m_bornMomenta.size(); ++i)
      m_ev.m_cmsboost.Boost(m_ev.m_bornMomenta[i]);
  } else {
    m_ev.m_cmsboost = Poincare();
  }
}

void YFS_Handler::SetMomenta(const ATOOLS::Vec4D_Vector &p) {
  m_ev.m_plab.clear();
  for(size_t i = 0; i < p.size(); ++i) {
    Vec4D pi(p[i]);
    m_ev.m_cmsboost.Boost(pi);
    m_ev.m_plab.push_back(pi);
  }
}

void YFS_Handler::CreatMomentumMap() {
  m_ev.m_outparticles.clear();
  for(size_t i = 0; i < 2; ++i)
  {
    m_particles[i]->SetMomentum(m_ev.m_bornMomenta[i]);
  }
  if(m_mode!=yfsmode::isr){
    for(size_t i = 2; i < m_flavs.size(); ++i)
    {
      m_ev.m_outparticles[m_particles[i]] = m_ev.m_bornMomenta[i];
      m_particles[i]->SetMomentum(m_ev.m_bornMomenta[i]);
    }
  }
}

void YFS_Handler::InitializeACRAIC(const ATOOLS::Flavour_Vector &fl) {
  std::string key;
  for (const Flavour &f : fl) key += ATOOLS::ToString((long)f) + ",";
  if (p_acraic && key == m_acraickey) { SetAcraicChannelProcs(fl); return; }
  if (p_acraic) m_acraicstore[m_acraickey] = std::move(p_acraic);
  auto it(m_acraicstore.find(key));
  if (it != m_acraicstore.end()) { p_acraic = std::move(it->second); m_acraicstore.erase(it); }
  else p_acraic = std::make_unique<Acraic_Base>(fl);
  m_acraickey = key;
  p_acraic->SetBornMomenta(m_ev.m_bornMomenta);
  SetAcraicChannelProcs(fl);
}

namespace {
  /*
    The member of a process group whose flavours are those of the channel
    fl plus nextra photons: the Comix amplitudes of a group exist only for its
    members. A single process, or a group of one, is returned as it is.
  */
  PHASIC::Process_Base *ChannelOf(PHASIC::Process_Base *group,
                                  const ATOOLS::Flavour_Vector &fl,
                                  const size_t nextra)
  {
    if (group == nullptr || !group->IsGroup() || group->Size() < 2) return group;
    for (size_t i(0); i < group->Size(); ++i) {
      PHASIC::Process_Base *c(ChannelOf((*group)[i], fl, nextra));
      if (c->IsGroup()) continue;
      const ATOOLS::Flavour_Vector &cf(c->Flavours());
      if (cf.size() != fl.size() + nextra) continue;
      bool same(true);
      std::multiset<long> a, b;
      for (size_t j(0); j < fl.size(); ++j) {
        if (j < 2 && cf[j] != fl[j]) same = false;
        if (j >= 2) { a.insert((long)fl[j]); b.insert((long)cf[j]); }
      }
      for (size_t j(0); j < nextra; ++j) a.insert((long)ATOOLS::Flavour(kf_photon));
      for (size_t j(0); j < nextra; ++j) b.insert((long)cf[fl.size() + j]);
      if (same && a == b) return c;
    }
    std::string name;
    for (const ATOOLS::Flavour &f : fl) name += ATOOLS::ToString((long)f) + " ";
    THROW(fatal_error, "YFS: no channel of " + group->Name()
                       + " has the flavours " + name);
    return nullptr;
  }
}

void YFS_Handler::SetAcraicChannelProcs(const ATOOLS::Flavour_Vector &fl) {
  m_acraicborn = ChannelOf(m_borngroup, fl, 0);
  m_acraicreal = ChannelOf(m_realgroup, fl, 1);
  if (p_acraic) {
    p_acraic->SetBornProc(m_acraicborn);
    p_acraic->SetRealProc(m_acraicreal);
  }
  if (p_nlo && m_acraicborn) p_nlo->SetBornProc(m_acraicborn);
  if (!p_nlo) return;
  // every n-photon real provider (the real, the double real, ACRAIC's tree
  // residuals) on its channel; the group it was built with is remembered the
  // first time it is seen, before it is re-pointed
  for (size_t n(1); n <= p_nlo->MaxRealPhotons(); ++n) {
    YFS::Real_Correction *prov(p_nlo->RealProvider(n));
    if (prov == nullptr) continue;
    auto it(m_realprovgroups.find(prov));
    if (it == m_realprovgroups.end())
      it = m_realprovgroups.emplace(prov, prov->p_proc).first;
    PHASIC::Process_Base *group(it->second);
    if (group == nullptr || !group->IsGroup() || group->Size() < 2) continue;
    // the order YFS hands the real its momenta in: the channel, then the photons
    ATOOLS::Flavour_Vector rfl(fl);
    for (size_t j(0); j < n; ++j) rfl.push_back(ATOOLS::Flavour(kf_photon));
    prov->SetChannelProc(ChannelOf(group, fl, n), rfl);
  }
}


void YFS_Handler::SetAcraicProcs(PHASIC::Process_Base *born,
                               PHASIC::Process_Base *real) {
  m_borngroup = born;
  m_realgroup = real;
  m_acraicborn = born;
  m_acraicreal = real;
  if (p_acraic) { p_acraic->SetBornProc(born); p_acraic->SetRealProc(real); }
  if (p_nlo)  p_nlo->SetBornProc(born);
}

bool YFS_Handler::MakeYFS(){
  return MakeYFS(m_ev.m_bornMomenta);
}

void YFS_Handler::SetEventS(const double s) {
  SetS(s);
  p_isr->SetS(s);
  p_fsr->SetS(s);
  p_dipoles->SetS(s);
  p_yfsFormFact->SetS(s);
  p_coulomb->m_s = s;
  if (p_nlo) p_nlo->SetS(s);
  if (p_acraic) p_acraic->SetS(s);
}

bool YFS_Handler::IIDipoleIsCurrent() {
  if (!p_dipoles->HasDipoleII()) return false;
  Dipole &D(p_dipoles->GetDipoleII());
  for (size_t i(0); i < 2; ++i)
    for (int mu(0); mu < 4; ++mu)
      if (D.GetBornMomenta(i)[mu] != m_ev.m_bornMomenta[i][mu]) return false;
  return true;
}

bool YFS_Handler::MakeYFS(ATOOLS::Vec4D_Vector &p)
{
  Reset();
  if (!ExchangeLinesPassTMin(m_ev.m_bornMomenta)) {
    CountFSRFailure("MakeYFS: reduced Born |t| < EXCHANGE_T_MIN");
    m_ev.m_yfsweight = 0.0;
    return false;
  }
  SetEventS((p[0] + p[1]).Abs2());
  if (m_isrinital || !IIDipoleIsCurrent()) {
    p_dipoles->MakeDipolesII(m_flavs, m_ev.m_plab, m_ev.m_bornMomenta);
    m_isrinital = true;
  }
  m_ev.m_ww_formfact = 1;
  m_fsrWeight = m_isrWeight = 1.0;
  CreatMomentumMap();
  if (m_mode == yfsmode::fsr) m_sp = m_s;
  m_v = 1. - m_sp / m_s;
  if ( m_v > m_vmax ) {
    CountFSRFailure("MakeYFS: v > vmax");
    m_ev.m_yfsweight = 0.0;
    return false;
  }
  p_isr->SetV(m_v);
  if (m_v <= m_deltacut && m_mode!=yfsmode::fsr) { // correction weight included in Generate photon
    CountFSRFailure("MakeYFS: v <= deltacut");
    Reset();
    return false;
  }
  if (!CalculateISR()) { CountFSRFailure("CalculateISR false"); return 0; }
  m_ev.m_FSRPhotons.clear();
  CalculateWWForm();
  CalculateCoulomb();
  p = m_ev.m_plab;
  return true;
}



void YFS_Handler::MakeACRAIC() {
  if (m_useacraic) {
    Vec4D_Vector vv;
    if (p_nlo) {
      for (size_t n(1); n <= p_nlo->MaxRealPhotons(); ++n)
        p_acraic->SetRealProvider(n, p_nlo->RealProvider(n));
      p_acraic->SetHasVirtual(p_nlo->HasVirtual() ||
                            (p_nlo->VirtualRequested() &&
                             p_nlo->AcraicSuppliesVirtual()));
      const bool autoown(m_acraicvirtsrc == acraicvirt::automatic &&
                         p_nlo->AcraicSuppliesVirtual());
      if (autoown && !p_acraic->IsFFbar())
        THROW(fatal_error, "YFS: ACRAIC_Virtual auto with V in NLO_Part and no"
              " Loop_Generator falls back to ACRAIC's own virtual, which is"
              " e+e- -> f fbar only. Name a Loop_Generator (OpenLoops, Recola)"
              " for this process, or drop V from NLO_Part.");
      // ACRAIC's own virtual has no beta_0^(2) or O2V streams
      if (p_nlo->AcraicSuppliesVirtual() &&
          (p_acraic->Order() >= 2 || p_nlo->Beta02Acraic().On()))
        THROW(fatal_error, "ACRAIC: ORDER 2 (or ACRAIC: BETA0_2) with ACRAIC's own"
              " virtual (YFS: ACRAIC_Virtual acraic, or auto with Loop_Generator"
              " Internal or none) has no beta_0^(2) and no external virtual for the"
              " O2V streams. Use YFS: ACRAIC_Virtual external with an external"
              " Loop_Generator (OpenLoops, Recola).");
      p_acraic->SetAutoOwnVirtual(autoown);
      p_acraic->SetNLO(p_nlo.get());
      if (m_acraicborn) p_nlo->SetBornProc(m_acraicborn);
    }
    p_acraic->SetBorn(m_born);
    for(size_t i = 0; i < m_ev.m_plab.size(); ++i) vv.push_back(m_ev.m_bornMomenta[i]);
    const Vec4D_Vector &finlegs(m_ev.m_acraiclegs.size() == m_ev.m_plab.size()
                                ? m_ev.m_acraiclegs : m_ev.m_plab);
    for(size_t i = 2; i < m_ev.m_plab.size(); ++i)
      vv.push_back(finlegs[i]);
    p_acraic->Init(vv);

    if (m_ev.m_sprimeBeams.size() == 2) {
      Vec4D_Vector sp(m_ev.m_sprimeBeams);
      for (size_t i = 2; i < m_ev.m_plab.size(); ++i) sp.push_back(finlegs[i]);
      p_acraic->SetLabMomenta(sp);
    }
    p_acraic->SetPreFSRMomenta(m_ev.m_reallab);
    p_acraic->SetISRPhotons(m_ev.m_ISRPhotons);
    if (HasFSR()) p_acraic->SetFSRPhotons(m_ev.m_FSRPhotons);

    static const bool regcheck(Settings::GetMainSettings()["YFS"]["SOFT_REGION_CHECK"]
                               .SetDefault(false).Get<bool>());
    
    if (regcheck && p_isr && p_isr->GetCut() == 0) ++m_regcheck_zero;
    else if (regcheck && m_ev.m_bornMomenta.size() >= 2) {
      const SoftRegion common(SoftRegion::RestOf(m_ev.m_bornMomenta[0] + m_ev.m_bornMomenta[1],
                                                 0.5*sqrt(m_s)*m_isrcut, "P"));
      auto check = [&](const Vec4D_Vector &ks, int w) {
        for (const Vec4D &k : ks) {
          const double r(common.PhotonEnergy(k)/common.Emin());
          ++m_regcheck_n[w];
          m_regcheck_min[w] = std::min(m_regcheck_min[w], r);
          if (r < 1. - 1e-9) ++m_regcheck_bad[w];
        }
      };
      check(m_ev.m_ISRPhotons, 0);
      if (HasFSR()) check(m_ev.m_FSRPhotons, 1);
    }
    p_acraic->SetBornMomenta(m_ev.m_bornMomenta);
    p_acraic->SetISRFormFactor(m_ev.m_formfactor);

    std::vector<std::vector<int> > groups;
    for (auto &D : p_dipoles->GetDipoleFF())
      if (D.IsResonance()) groups.push_back({D.Left(), D.Right()});
    // the generator's own stages, for the crude (ACRAIC: CRUDE_FROM_GENERATOR)
    p_acraic->SetCrudeStageGroups(groups);
    if (!SetACRAICWStageGroups()) p_acraic->SetStageGroups(groups);
    p_acraic->Calculate();
  }

}

bool YFS_Handler::SetACRAICWStageGroups() {
  if (!p_acraic || !p_acraic->WStagesRequested()) return false;
  const DipoleSet::WWLegs w(DipoleSet::FindWW(m_flavs, m_ev.m_plab));
  if (!w.ok) return false;
  if (p_acraic->WStagesMode() == tristate::automatic) {
    // auto: the pole scheme's own window, both W's near the pole
    const double MW(Flavour(kf_Wplus).Mass()), GW(Flavour(kf_Wplus).Width());
    if (GW > 0. && (fabs(w.wm.Mass() - MW)/GW > m_resonace_max
                    || fabs(w.wp.Mass() - MW)/GW > m_resonace_max)) return false;
  }
  std::vector<Acraic_Base::StageGroup> groups(2);
  groups[0].externalLegs       = {(int)w.lp};
  groups[0].resonanceDaughters = {(int)w.lp, (int)w.np};
  groups[0].resonance          = Flavour(kf_Wplus);
  groups[1].externalLegs       = {(int)w.lm};
  groups[1].resonanceDaughters = {(int)w.lm, (int)w.nm};
  groups[1].resonance          = Flavour(kf_Wplus).Bar();
  p_acraic->SetStageGroups(groups);
  return true;
}

void YFS_Handler::CalculateWWForm() {
  if (m_formWW) {
    if (!WPairMomenta(m_ev.m_plab, m_ev.m_Wp, m_ev.m_Wm)) {
      m_ev.m_Wp *= 0.; m_ev.m_Wm *= 0.;
      THROW(fatal_error, "YFS: WW_FORM needs a W+ W- final state or its four daughters");
    }
    m_ev.m_ww_formfact = p_yfsFormFact->BVV_WW(m_ev.m_plab, m_ev.m_ISRPhotons, m_ev.m_Wp, m_ev.m_Wm,
                                          m_photonMass, sqrt(m_sp) / 2.);
    if (IsBad(m_ev.m_ww_formfact) || m_ev.m_ww_formfact < 0) {
      msg_Error() << METHOD << ": BVV_WW returned " << m_ev.m_ww_formfact
                  << "; setting it to 1. Use WW_Scheme: pole instead.\n";
      m_ev.m_ww_formfact = 1.;
    }
  }
}

bool YFS_Handler::CalculateISR() {
  if (m_mode==yfsmode::fsr) return true;
  if (p_dipoles->GetDipoleII().size() != 2) {
    THROW(fatal_error, "Wrong dipole size for ISR");
  }
  if (m_isrinital) p_isr->SetIncoming(&p_dipoles->GetDipoleII());
  m_isrinital = false;
  m_g=p_dipoles->GetDipoleII().m_gamma;
  m_gp=p_dipoles->GetDipoleII().m_gammap;
  Vec4D_Vector me_acc;   // ISR photons are not hidden, so nothing accumulates here
  const YFS::EmissionResult res(
      p_dipoles->GetDipoleII().GenerateEmissions(p_isr.get(), p_fsr.get(), m_born, m_v, me_acc));
  m_ev.m_photonSumISR = res.photon_sum;
  m_ev.m_ISRPhotons.clear();
  for (const YFS::Photon &k : res.photons) m_ev.m_ISRPhotons.push_back(k.K());
  m_isrWeight = res.weight;
  m_ev.m_photons = res.photons;
  m_ev.m_sprimeBeams.assign(2, Vec4D());
  for(size_t i = 0; i < 2; ++i) {
    m_ev.m_plab[i] = p_dipoles->GetDipoleII().GetNewMomenta(i); 
    m_ev.m_sprimeBeams[i] = m_ev.m_plab[i];   // before ToLab: the pair frame
    ToLab(m_ev.m_plab[i]);
  }
  double sp = (m_ev.m_plab[0] + m_ev.m_plab[1]).Abs2();
  if (!IsEqual(sp, m_sp, 1e-4)) {
    msg_Error() << "Boost failed, sprime"
                << " is " << sp << " and should be "
                << m_sp << std::endl << "Diff = " <<
                m_sp - sp << std::endl << " Event with "
                << " N=" << p_dipoles->GetDipoleII().GetPhotons().size() << " photons" << std::endl
                << " V = " << m_v << std::endl
                << " Vmin = " << m_isrcut << std::endl
                << "ISR NPHotons = " << m_ev.m_ISRPhotons.size() << std::endl;
  }
  return true;
}



void YFS_Handler::AddFormFactor() {
  m_ev.m_formfactor_sum = 0.;
  if (m_fullform >= 1) {
    m_ev.m_formfactor_sum = p_dipoles->FormFactorSum();
    m_ev.m_formfactor = p_dipoles->FormFactor(m_ev.m_formfactor_sum);
  }
  else if (m_fullform == 2) {
    m_ev.m_formfactor = exp(m_g / 4.);//-m_alpha*M_PI);
  }
  else if (m_fullform == -1) {
    m_ev.m_formfactor = 1;
  }
  else m_ev.m_formfactor = exp(m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
}

bool YFS_Handler::CalculateFSR(){
  return CalculateFSR(m_ev.m_plab);
}

static void AbsorbHiddenPhotons(Vec4D &p1, Vec4D &p2, const double m1,
                                const double m2, const Vec4D &khidden)
{
  const Vec4D P(p1 + p2), Pnew(P + khidden);
  Vec4D q1(p1);
  Poincare rest(P);
  rest.Boost(q1);
  const Vec3D n(Vec3D(q1)/Vec3D(q1).Abs());
  const double snew(Pnew.Abs2());
  const double q(0.5*sqrt(YFS::Lambda(snew, m1*m1, m2*m2)/snew));
  p1 = Vec4D(sqrt(q*q + m1*m1), q*n);
  p2 = Vec4D(sqrt(q*q + m2*m2), -q*n);
  Poincare restnew(Pnew);
  restnew.BoostBack(p1);
  restnew.BoostBack(p2);
}

bool YFS_Handler::CalculateFSR(Vec4D_Vector & p) {
  m_ev.m_FSRPhotons.clear();
  m_ev.m_fsrphotonsforME.clear();
  m_ev.m_reallab = p;
  m_ev.m_plab=p;
  CreatMomentumMap();
  if (m_mode != yfsmode::fsr && !BornPhotonsPassReducedBeamAngle(p)) {
    CountFSRFailure("CalculateFSR: Born photon inside BORN_PHOTON_REDUCED_ANGLE");
    m_ev.m_yfsweight = 0.0;
    Reset();
    return false;
  }
  if(m_mode==yfsmode::isr) return true;
  m_fsrWeight=1;
  p_dipoles->MakeDipoles(m_flavs, m_ev.m_plab, m_ev.m_plab);
  CheckResonance();
  if(m_mode==yfsmode::isrfsr) {
    Vec4D_Vector ifmom(m_ev.m_plab);
    ifmom[0] = m_ev.m_bornMomenta[0];
    ifmom[1] = m_ev.m_bornMomenta[1];
    if (!ExchangeLinesPassTMin(ifmom)) {
      CountFSRFailure("CalculateFSR: initial-final |t| < EXCHANGE_T_MIN");
      m_ev.m_yfsweight = 0.0;
      Reset();
      return false;
    }
    p_dipoles->MakeDipolesIF(m_flavs, ifmom, ifmom);
  }
  Vec4D_Vector polemom(m_ev.m_plab);
  if (polemom.size() > 1) {
    polemom[0] = m_ev.m_bornMomenta[0];
    polemom[1] = m_ev.m_bornMomenta[1];
  }
  p_dipoles->MakeDipolesPole(m_flavs, polemom, polemom);
  YFS::DipoleView ffdip(p_dipoles->GetDipoleFF());
  ++m_fsrtrials;
  std::map<size_t, Vec4D> recordonly;
  for (auto Dip = ffdip.begin(); Dip != ffdip.end(); ++Dip) {
    if(!Dip->IsResonance()) continue;
    const YFS::EmissionResult res(
        Dip->GenerateEmissions(p_isr.get(), p_fsr.get(), m_born, m_v, m_ev.m_fsrphotonsforME));
    const std::string dname(Dip->GetFlav(0).IDName()+","+Dip->GetFlav(1).IDName());
    switch (res.fail) {
    case YFS::EmissionResult::Failure::initialize:
      CountFSRFailure("initialize ("+dname+")");
      Reset();
      return false;
    case YFS::EmissionResult::Failure::makefsr:
      CountFSRFailure("makefsr cut "+ToString((int)p_fsr->Cut())+" ("+dname+")");
      Reset();
      if (m_fsr_debug) p_debug->FillHist(m_ev.m_plab, p_isr.get(), p_fsr.get());
      return false;
    case YFS::EmissionResult::Failure::masswgt:
      CountFSRFailure("mass weight ("+dname+")");
      m_fsrWeight = 0;
      if (m_fsr_debug) p_debug->FillHist(m_ev.m_plab, p_isr.get(), p_fsr.get());
      return false;
    case YFS::EmissionResult::Failure::formfactor:
      CountFSRFailure("form factor ("+dname+")");
      return false;
    case YFS::EmissionResult::Failure::none:
      if (res.weight == 0.) CountFSRFailure("zero weight ("+dname+")");
      break;
    }
    m_ev.m_photonSumFSR = res.photon_sum;
    m_ev.m_FSRPhotons.clear();
    for (const YFS::Photon &k : res.photons) m_ev.m_FSRPhotons.push_back(k.K());
    m_fsrWeight *= res.weight;
    if (p_dipoles->PoleActive()) {
      // The radiating dipole is the W pair, and the W's are not entries in the
      // event record -- Left()/Right() point at the charged leptons they
      // decayed to. Writing the new W momenta there would put an 80 GeV
      // momentum in the muon's slot. Carry the recoil down to the four
      // fermions instead.
      if (!p_dipoles->ApplyPoleRecoil(m_ev.m_plab)) {
        CountFSRFailure("pole recoil");
        Reset();
        return false;
      }
      continue;
    }
    const Vec4D pairbefore(m_ev.m_plab[Dip->Left()] + m_ev.m_plab[Dip->Right()]);
    m_ev.m_plab[Dip->Left()]  =  Dip->GetNewMomenta(0);
    m_ev.m_plab[Dip->Right()] =  Dip->GetNewMomenta(1);
    if (res.photons.size() < res.me_photons.size()) {
      // photons were hidden: the pair conserved pairbefore against all of them
      Vec4D khidden(pairbefore - m_ev.m_plab[Dip->Left()] - m_ev.m_plab[Dip->Right()]);
      for (const YFS::Photon &k : res.photons) khidden -= k.K();
      Vec4D p1(m_ev.m_plab[Dip->Left()]), p2(m_ev.m_plab[Dip->Right()]);
      AbsorbHiddenPhotons(p1, p2, m_flavs[Dip->Left()].Mass(),
                          m_flavs[Dip->Right()].Mass(), khidden);
      recordonly[Dip->Left()]  = p1;
      recordonly[Dip->Right()] = p2;
    }
    if(!IsEqual(m_flavs[Dip->Left()].Mass(), m_ev.m_plab[Dip->Left()].Mass(),1e-5)){
      msg_Debugging()<<"Missmatch in Final state mass"<<std::endl
                 <<"Flavour = "<<m_flavs[Dip->Left()]<<std::endl
                 <<"Mass =   "<<m_flavs[Dip->Left()].Mass()<<std::endl
                 <<"Momentum =   "<<m_ev.m_plab[Dip->Left()]<<std::endl
                 <<"Mass =   "<<m_ev.m_plab[Dip->Left()].Mass()<<std::endl;
    }
    if(!IsEqual(m_flavs[Dip->Right()].Mass(), m_ev.m_plab[Dip->Right()].Mass(),1e-5)){
      msg_Debugging()<<"Missmatch in Final state mass"<<std::endl
                 <<"Flavour = "<<m_flavs[Dip->Right()]<<std::endl
                 <<"Mass =   "<<m_flavs[Dip->Right()].Mass()<<std::endl
                 <<"Momentum =   "<<m_ev.m_plab[Dip->Right()]<<std::endl
                 <<"Mass =   "<<m_ev.m_plab[Dip->Right()].Mass()<<std::endl;
    }
  }
  for(size_t i = 2; i < m_ev.m_plab.size(); ++i) {
    const auto rec(recordonly.find(i));
    m_ev.m_outparticles[m_particles[i]] =
      rec == recordonly.end() ? m_ev.m_plab[i] : rec->second;
  }
  m_ev.m_acraiclegs.clear();
  if (!recordonly.empty()) {
    m_ev.m_acraiclegs = m_ev.m_plab;
    for (const auto &r : recordonly) m_ev.m_acraiclegs[r.first] = r.second;
  }
  if (m_exchange_tmin_final > 0. && m_ev.m_bornMomenta.size() >= 2) {
    Vec4D_Vector phys(m_ev.m_acraiclegs.size() == m_ev.m_plab.size() ? m_ev.m_acraiclegs
                                                                    : m_ev.m_plab);
    phys[0] = m_ev.m_bornMomenta[0];
    phys[1] = m_ev.m_bornMomenta[1];
    if (!ExchangeLinesPassTMin(phys, m_exchange_tmin_final)) {
      CountFSRFailure("CalculateFSR: physical |t| < EXCHANGE_T_MIN_FINAL");
      m_ev.m_yfsweight = 0.0;
      Reset();
      return false;
    }
  }
  // get all photons
  m_ev.m_FSRPhotons.clear();
  m_ev.m_fsrphotonsforME.clear();
  // Rebuilt in full here rather than appended to, so a re-entered
  // CalculateFSR cannot leave last trial's photons behind.
  m_ev.m_photons.clear();
  m_ev.m_me_photons.clear();
  if (p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_ev.m_ISRPhotons)
      m_ev.m_photons.push_back(YFS::Photon(k, &p_dipoles->GetDipoleII()));
  YFS::DipoleView ffcollect(p_dipoles->GetDipoleFF());
  for (auto Dip = ffcollect.begin(); Dip != ffcollect.end(); ++Dip) {
    for(auto &k: Dip->GetPhotons()) {
      m_ev.m_FSRPhotons.push_back(k);
      m_ev.m_photons.push_back(YFS::Photon(k, &*Dip));
    }
    for(auto &k: Dip->GetMEPhotons()) {
      m_ev.m_fsrphotonsforME.push_back(k);
      m_ev.m_me_photons.push_back(YFS::Photon(k, &*Dip));
    }
  }
  if (!PhotonResolutionPasses()) {
    CountFSRFailure("CalculateFSR: photon resolution (RESOLVED_PHOTON_E/ANGLE)");
    m_ev.m_yfsweight = 0.0;
    Reset();
    return false;
  }
  CountFSRFailure("CalculateFSR completed");
  return true;
}


bool YFS_Handler::ExchangeLinesPassTMin(const Vec4D_Vector &p, double tmin) const {
  if (tmin <= 0.) return true;
  const size_t n(std::min(p.size(), m_flavs.size()));
  for (size_t i(0); i < 2 && i < n; ++i) {
    if (m_flavs[i].Charge() == 0.) continue;
    for (size_t j(2); j < n; ++j)
      if (m_flavs[j] == m_flavs[i] && -(p[i] - p[j]).Abs2() < tmin)
        return false;
  }
  return true;
}

// YFS: BORN_PHOTON_REDUCED_ANGLE: veto Born photons close to a reduced beam
bool YFS_Handler::BornPhotonsPassReducedBeamAngle(const Vec4D_Vector &p) const {
  if (m_bornphoton_reduced_angle <= 0.) return true;
  const double cosmax(cos(m_bornphoton_reduced_angle*M_PI/180.));
  Poincare rest(p[0] + p[1]);
  Vec4D beam(p[0]);
  rest.Boost(beam);
  const size_t n(std::min(p.size(), m_flavs.size()));
  for (size_t j(2); j < n; ++j) {
    if (!m_flavs[j].IsPhoton()) continue;
    Vec4D k(p[j]);
    rest.Boost(k);
    // the reduced beams are back to back here, so one cosine covers both
    if (fabs(k.CosTheta(beam)) > cosmax) return false;
  }
  return true;
}

// YFS: RESOLVED_PHOTON_E/ANGLE: all Born photons resolved, no YFS photon resolved
bool YFS_Handler::PhotonResolutionPasses() const {
  if (m_resolved_photon_e <= 0. && m_resolved_photon_angle <= 0.) return true;
  const Vec4D_Vector p(const_cast<YFS_Handler*>(this)->GetMomenta());
  const double cmin(cos(m_resolved_photon_angle*M_PI/180.));
  const Vec3D zaxis(0., 0., 1.);
  auto resolved = [&](const Vec4D &k) {
    if (k[0] <= m_resolved_photon_e) return false;
    const Vec3D kk(k);
    const double kabs(kk.Abs());
    if (!(kabs > 0.)) return false;
    if (fabs(kk*zaxis)/kabs > cmin) return false;          // near either beam
    for (size_t i(2); i < p.size() && i < m_flavs.size(); ++i) {
      if (m_flavs[i].Charge() == 0.) continue;
      const Vec3D q(p[i]);
      if (kk*q/(kabs*q.Abs()) > cmin) return false;       // near a charged leg
    }
    return true;
  };
  for (size_t i(2); i < p.size() && i < m_flavs.size(); ++i)
    if (m_flavs[i].IsPhoton() && !resolved(p[i])) return false;
  for (const Vec4D &k : GetISRPhotons()) if (resolved(k)) return false;
  for (const Vec4D &k : GetFSRPhotons()) if (resolved(k)) return false;
  return true;
}

void YFS_Handler::CountFSRFailure(const std::string &what) {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["FSR_FAILURE_STATS"].SetDefault(0).Get<int>() != 0);
  if (on) ++m_fsrfail[what];
}

bool YFS_Handler::WPairMomenta(const Vec4D_Vector &p, Vec4D &wp, Vec4D &wm) const {
  auto partner = [](const Flavour &a, const Flavour &b) {
    if (a.IsAnti() == b.IsAnti()) return false;
    const kf_code ka(a.Kfcode()), kb(b.Kfcode());
    const kf_code lo(std::min(ka, kb)), hi(std::max(ka, kb));
    const bool lep(lo >= 11 && lo <= 15 && lo % 2 == 1 && hi == lo + 1);
    const bool qrk((lo == 1 && hi == 2) || (lo == 3 && hi == 4) || (lo == 5 && hi == 6));
    return lep || qrk;
  };
  // an on-shell W+ W- final state
  int iwp(-1), iwm(-1);
  for (size_t i(2); i < p.size() && i < m_flavs.size(); ++i) {
    if (m_flavs[i].Kfcode() != kf_Wplus) continue;
    (m_flavs[i].Charge() > 0. ? iwp : iwm) = i;
  }
  if (iwp >= 0 && iwm >= 0) { wp = p[iwp]; wm = p[iwm]; return true; }
  std::vector<size_t> fin;
  for (size_t i(2); i < p.size() && i < m_flavs.size(); ++i)
    if (m_flavs[i].IsFermion()) fin.push_back(i);
  if (fin.size() != 4) return false;
  // the three ways to split four fermions into two pairs
  const int split[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
  for (const auto &sp : split) {
    const size_t a(fin[sp[0]]), b(fin[sp[1]]), c(fin[sp[2]]), d(fin[sp[3]]);
    if (!partner(m_flavs[a], m_flavs[b]) || !partner(m_flavs[c], m_flavs[d])) continue;
    const double q1(m_flavs[a].Charge() + m_flavs[b].Charge());
    const double q2(m_flavs[c].Charge() + m_flavs[d].Charge());
    if (std::abs(std::abs(q1) - 1.) > 1e-6 || std::abs(q1 + q2) > 1e-6) continue;
    const Vec4D P1(p[a] + p[b]), P2(p[c] + p[d]);
    wp = q1 > 0. ? P1 : P2;
    wm = q1 > 0. ? P2 : P1;
    return true;
  }
  return false;
}

void YFS_Handler::CalculateCoulomb() {
  if (m_coulombho != coulombho::none) {
    Vec4D wp, wm;
    if (WPairMomenta(m_ev.m_plab, wp, wm))
      m_ev.m_coulho = p_coulomb->HigherOrderWeight((wp + wm).Abs2(), wp.Abs2(), wm.Abs2(),
                                                   m_coulombho == coulombho::second ? 1 : 2);
    if (IsBad(m_ev.m_coulho)) m_ev.m_coulho = 1.;
  }
  if (!m_coulomb) return;
  if (!WPairMomenta(m_ev.m_plab, m_ev.m_Wp, m_ev.m_Wm))
    THROW(fatal_error, "YFS: COULOMB needs a W+ W- final state or its four daughters");
  p_coulomb->Calculate(m_ev.m_Wp, m_ev.m_Wm);
  if (m_formWW) {
    // need to Subtract the Coulomb loop from virtual form factor
    // double s  = (m_ev.m_Wp + m_ev.m_Wm).Abs2();
    double am1 = m_ev.m_Wp.Abs2();
    double am2 = m_ev.m_Wm.Abs2();
    double beta = sqrt(1. - 2.*(am1 + am2) / m_s + sqr((am1 - am2) / m_s));
    if (m_betatWW >= beta) {
      p_coulomb->Subtract();
    }
    else m_ev.m_coulSub = 0;
  }
}

void YFS_Handler::GenerateWeight() {
  p_weights->GenerateWeight();
}

void YFS_Handler::YFSDebug(double W){
  p_debug->FillHist(m_ev.m_plab, p_isr.get(), p_fsr.get(), W);
}

void YFS_Handler::Reset() {
  m_fsrWeight = 0;
  m_ev.m_yfsweight = 0;
  m_ev.m_ISRPhotons.clear();
  m_ev.m_FSRPhotons.clear();
  m_ev.m_photonSumISR *= 0;
  m_ev.m_photonSumFSR *= 0;
  m_ev.m_real = 1;
  m_ev.m_eex = 0.;
  // m_s = sqr(rpa->gen.Ecms());
}

Vec4D_Vector YFS_Handler::GetPhotons(){
  Vec4D_Vector k;
  for(auto p: m_ev.m_ISRPhotons) k.push_back(p);
  for(auto p: m_ev.m_FSRPhotons) k.push_back(p);
  return k;
}

void YFS_Handler::CheckResonance(){
  YFS::DipoleView ffres(p_dipoles->GetDipoleFF());
  for (auto D1 = ffres.begin(); D1 != ffres.end(); ++D1) {
    for (auto D2 = ffres.begin(); D2 != ffres.end(); ++D2) {
      if(D1==D2) continue;
      if(!D1->IsResonance() || !D2->IsResonance()) continue;
      if(D1->Right() == D2->Right() ||  D1->Right() == D2->Left()|| 
        D1->Left() == D2->Right()||  D1->Left() == D2->Left()){
        if(p_dipoles->ResonantDist(*D1) < p_dipoles->ResonantDist(*D2)) D2->SetResonance(false);
        else  D1->SetResonance(false);
        }
      }
    }
  }

