#include "ATOOLS/Org/Message.H"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include "YFS/Main/YFS_Handler.H"
#include "BEAM/Main/Beam_Base.H"
#include "YFS/Main/ISR.H"

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

YFS_Handler::YFS_Handler()
{
  if(Mode()!=YFS::yfsmode::off){
    p_dipoles = std::make_unique<Define_Dipoles>();
    p_coulomb = std::make_unique<Coulomb>(m_coulomb);
    p_fsr = std::make_unique<FSR>();
    p_debug = std::make_unique<Debug>();
    p_yfsFormFact = std::make_unique<YFS::YFS_Form_Factor>();
    m_setparticles = false;
    p_isr = std::make_unique<YFS::ISR>();
    if (m_fb_analysis)
      p_fb = std::make_unique<YFS::YFS_FB_Analysis>(std::vector<YFS::fbdef::code>{}, m_fb_kf);
    m_isrinital = true;
    // the YFS record puts its photons in YFS_Initial and QED_Radiation
    // blobs (Signal_Processes::FillYFSBlobs); the photon splitter works on
    // photons of the signal blob and has not been ported
    if (m_photon_split != 0)
      THROW(not_implemented, "YFS: PHOTON_SPLITTER_MODE is not supported with"
            " YFS: its photons are in the YFS_Initial and QED_Radiation blobs,"
            " which the photon splitter does not handle.");
    m_rmode = 0;
    m_negskip = 0;
    // m_ev needs nothing here: YFS_Event's default member initialisers are
    // the single definition of the per-event starting values, and every
    // event re-establishes them through StartEvent().
    rpa->gen.AddCitation(1,"The automation of YFS ISR is published in  \\cite{Krauss:2022ajk}.Which is based on \\cite{Jadach:1988gb}");
  }
}

YFS_Handler::~YFS_Handler()
{
  if(Mode()!=YFS::yfsmode::off){
    if (m_acrv_n > 0)
      msg_Out()<<"ACRAIC: REAL_VIRTUAL 2 on "<<m_acrv_n<<" events (this rank): NLO RV photons "
               <<"without a ACRAIC match "<<m_acrv_nomatch<<", matched but without M_1 "
               <<m_acrv_nom1<<std::endl;
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
      msg_Out()<<"YFS: ACRAIC supplied the O(alpha) weight on "
               <<cs.m_oen<<" events; mean (ACRAIC factor)/(EEX factor) = "
               <<cs.m_oesum/cs.m_oen
               <<(cs.m_bad ? " ("+ATOOLS::ToString(cs.m_bad)
                             +" events fell back to EEX)" : "")
               <<std::endl;
    // Everything this class owns is held by unique_ptr and released with it,
    // which also covers the two cases the hand written destructor got wrong:
    // p_acraic, new'd in InitializeACRAIC, was never released (so Acraic_Base's
    // end of run reporting never ran), and every member was deleted
    // unconditionally although the constructor only fills them when the mode
    // is not "off", so an "off" handler deleted uninitialised pointers.
    if(m_negskip!=0){
      msg_Out()<<"Total Events Skipped: "<<m_negskip<<std::endl;
    }
    // The emission-side IF reweight is only trustworthy while the clamp is
    // rarely hit; a large count here means photons are being handed an
    // interference factor the soft approximation cannot support, and the rate
    // is being shaped by IFI_RClip rather than by the physics.
    if(m_ifireal && p_dipoles && p_dipoles->IFIClipped()!=0){
      msg_Out()<<"IFI real reweight: "<<p_dipoles->IFIClipped()
               <<" photon factors clamped to ["<<p_dipoles->IFIRClip()
               <<", "<<1./p_dipoles->IFIRClip()<<"]"<<std::endl;
    }
    // The restoration is m_born*(subloc/subb - 1), so the spread of subloc/subb
    // is what sets both the shift and the MC error. A mean far from 1, or a
    // min/max spanning orders of magnitude, says the two eikonals are not the
    // matched pair the cancellation assumes - which is the thing to look at
    // before adjusting anything else.
    if(m_ifireal && p_nlo && p_nlo->m_ifi_n>0){
      const double mean = p_nlo->m_ifi_sum/p_nlo->m_ifi_n;
      const double var  = p_nlo->m_ifi_sum2/p_nlo->m_ifi_n - mean*mean;
      msg_Out()<<"IFI real restoration: n="<<p_nlo->m_ifi_n
               <<"  subloc/subb mean="<<mean
               <<" rms="<<(var>0.?sqrt(var):0.)
               <<" min="<<p_nlo->m_ifi_min
               <<" max="<<p_nlo->m_ifi_max<<std::endl;
      msg_Out()<<"  profiled in x = E_gamma/sqrt(s):"<<std::endl;
      for(int i=0;i<5;++i){
        if(p_nlo->m_ifi_x_n[i]==0) continue;
        const double rr = p_nlo->m_ifi_x_r[i]/p_nlo->m_ifi_x_n[i];
        const double ee = p_nlo->m_ifi_x_e[i]/p_nlo->m_ifi_x_n[i];
        msg_Out()<<"    x="<<0.1*i<<"-"<<0.1*(i+1)
                 <<"  n="<<p_nlo->m_ifi_x_n[i]
                 <<"  applied="<<rr<<"  exact="<<ee
                 <<"  residue="<<(rr-ee)<<std::endl;
      }
    }
    // Photons above omega whose IF ratio was never applied keep factor 1,
    // i.e. lose their share of the interference.
    if(m_ifireal && p_nlo && p_nlo->m_ifi_eligible>0){
      msg_Out()<<"IFI real restoration skipped on "
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
  // A raw observer on purpose: callers (YFS_Process) only use the NLO layer,
  // they never take it over.
  return p_nlo.get();
}


// bool YFS_Handler::On()
// {
//   return m_mode;
// }



void YFS_Handler::SetBeam(BEAM::Beam_Spectra_Handler *beam)
{
  p_beams = beam;
  // for(size_t i = 0; i < 2; ++i) m_beams.push_back(beam->GetBeam(i));
  m_beam1 = p_beams->GetBeam(0)->OutMomentum();
  m_beam2 = p_beams->GetBeam(1)->OutMomentum();
  if(m_beam1 != -m_beam2) m_asymbeams = true;
  else m_asymbeams = false;
}

void YFS_Handler::SetLimits(const double &smin, const double &s) {
  // s must be this event's actual s' (from the real, possibly beam-spread-
  // sampled, incoming momenta), not sqr(rpa->gen.Ecms()) -- that singleton is
  // the fixed nominal collider energy and does not track BEAM_SPECTRA
  // (Gaussian) event-by-event variation, which understated/overstated maxV
  // here let m_v exceed what the actual event could physically radiate.
  p_yfsFormFact->SetCharge(1);
  p_coulomb->SetAlphaQED(m_alpha);
  double maxV = 1. - smin / s;
  if (m_vmax > maxV && !m_asymbeams) {
    msg_Error() << "Warning: vmax to large in YFS integration reseting to " << maxV << std::endl;
    m_vmax = maxV;
  }
}

void YFS_Handler::SetFlavours(const ATOOLS::Flavour_Vector &flavs) {
  // One YFS_Handler is shared by every process in the run card, but m_flavs,
  // m_mass and m_particles describe a single process. Latching on
  // m_setparticles alone meant the first process to get here kept ownership of
  // those forever, so with more than one process the momenta of the active
  // process were combined with the flavours of a different one -- mismatched
  // multiplicities, which MakeDipoles used to index out of bounds. Rebuild
  // whenever the flavours actually change; the early return keeps the
  // per-event calls cheap when they do not.
  if(m_setparticles && m_flavs == flavs) return;
  // MODE implied by ACRAIC: 1 (YFS_Base::RegisterSettings): ISR alone when no
  // final-state particle is charged (nu nu, gamma gamma, H nu nu), else ISRFSR
  if (m_mode_from_acraic) {
    bool charged(false);
    for (size_t i(2); i < flavs.size(); ++i) if (flavs[i].Charge() != 0.) charged = true;
    m_mode = charged ? yfsmode::isrfsr : yfsmode::isr;
  }
  // Clearing the store frees the Particles; m_particles only observes them, so
  // it has to be emptied in step or it is left holding dangling pointers.
  // (The hand written version deleted through m_particles, which had to happen
  // before the clear -- an ordering the code had wrong once already.)
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
  if (m_useacraic) InitializeACRAIC(m_flavs);
}

void YFS_Handler::SetBornMomenta(const ATOOLS::Vec4D_Vector &p) {
  // The per-event boundary. Both the ISR/ISRFSR and the pure-FSR paths call
  // this first (Phase_Space_Point.C:229 and :241), so it is the one place
  // that sees every event; pure-FSR never reached MakeYFS()'s Reset().
  StartEvent();
  m_ev.m_bornMomenta.clear();
  for(size_t i = 0; i < p.size(); ++i) {
    m_ev.m_bornMomenta.push_back(p[i]);
  }
  // detect asymmetric beams from the original lab momenta, before any boost
  if(m_ev.m_bornMomenta[0] != -m_ev.m_bornMomenta[1]) m_asymbeams = true;
  else m_asymbeams = false;
  // NLO_Base::MapMomenta (and the ISR/FSR kinematics) assume the incoming pair
  // is at rest. For non-standard setups (fixed target, e.g. MUonE muon-e-;
  // asymmetric beams; beamstrahlung) it is not, so boost into the incoming-pair
  // rest frame here; the blob-facing getters (ToLab) undo it when handing the
  // event back. Pure-FSR mode feeds lab momenta straight into CalculateFSR(p)
  // and is left untouched.
  Vec4D Q(m_ev.m_bornMomenta[0] + m_ev.m_bornMomenta[1]);
  if (m_mode != yfsmode::fsr && !IsZero(Q.PSpat() / Q[0], 1e-10)) {
    m_ev.m_cmsboost = Poincare(Q);
    for (size_t i = 0; i < m_ev.m_bornMomenta.size(); ++i)
      m_ev.m_cmsboost.Boost(m_ev.m_bornMomenta[i]);
  } else {
    m_ev.m_cmsboost = Poincare();
  }
  if (m_formWW) MakeWWVecs(m_ev.m_bornMomenta);
  // AddFormFactor();
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
  m_ev.m_inparticles.clear();
  m_ev.m_outparticles.clear();
  for(size_t i = 0; i < 2; ++i)
  {
    m_ev.m_inparticles[m_particles[i]] = m_ev.m_bornMomenta[i];
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
  /*
    One ACRAIC object per flavour set. This returned as soon as ANY p_acraic
    existed, so with several processes in one card every process after the
    first ran ACRAIC with the first one's flavours, stages and helicity maps
    (YFS_Process::MakeActive re-points the flavours and providers of the
    shared handler, but ACRAIC kept its own). Now the current object is parked
    under its key and the one for fl is taken out, or created.
  */
  std::string key;
  for (const Flavour &f : fl) key += ATOOLS::ToString((long)f) + ",";
  if (p_acraic && key == m_acraickey) return;
  if (p_acraic) m_acraicstore[m_acraickey] = std::move(p_acraic);
  auto it(m_acraicstore.find(key));
  if (it != m_acraicstore.end()) { p_acraic = std::move(it->second); m_acraicstore.erase(it); }
  else p_acraic = std::make_unique<Acraic_Base>(fl);
  m_acraickey = key;
  p_acraic->SetBornMomenta(m_ev.m_bornMomenta);
  p_acraic->SetBornProc(m_acraicborn);
  p_acraic->SetRealProc(m_acraicreal);
}


void YFS_Handler::SetAcraicProcs(PHASIC::Process_Base *born,
                               PHASIC::Process_Base *real) {
  m_acraicborn = born;
  m_acraicreal = real;
  if (p_acraic) { p_acraic->SetBornProc(born); p_acraic->SetRealProc(real); }
  if (p_nlo)  p_nlo->SetBornProc(born);
}

bool YFS_Handler::MakeYFS(){
  return MakeYFS(m_ev.m_bornMomenta);
}

/*
  Every YFS object holds its own m_s, set to sqr(Ecms) when it was built. With
  a beam spectrum that is neither this event's s nor its mean, and the ISR
  photon energies, the form factors' soft limits, the FSR cutoff and ACRAIC's
  flux s/X^2 all read it. Without a spectrum s is the same every event and
  this rewrites the same value.
*/
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
  // p_dipoles->CreateAllDipoles(m_flavs, m_ev.m_plab, m_ev.m_bornMomenta);
  /*
    The initial-initial dipole is built once and kept, because ISR holds a
    pointer to it (CalculateISR re-points it while m_isrinital is set). With a
    fixed beam energy its beams are every event's. With a beam spectrum they
    are not: the kept dipole put every event's ISR on the FIRST event's beams,
    so s' ignored the spread and the photons did not balance the event's
    beams. Rebuild it whenever the beams have moved.
  */
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
    // The real providers, one per photon multiplicity, so ACRAIC: COMIX_REAL
    // can ask for Comix's helicity amplitudes at whatever multiplicity the
    // event came out with. Re-pointed every event because the handler is
    // shared by every process in the run card and the providers belong to
    // whichever one is currently selected.
    if (p_nlo) {
      for (size_t n(1); n <= p_nlo->MaxRealPhotons(); ++n)
        p_acraic->SetRealProvider(n, p_nlo->RealProvider(n));
      // Whether a virtual was requested at all, so ACRAIC's weight carries the
      // same perturbative content the card asked for.
      p_acraic->SetHasVirtual(p_nlo->HasVirtual() ||
                            (p_nlo->VirtualRequested() &&
                             p_nlo->AcraicSuppliesVirtual()));
      p_acraic->SetNLO(p_nlo.get());
      if (m_acraicborn) p_nlo->SetBornProc(m_acraicborn);
    }
    p_acraic->SetBorn(m_born);
    for(size_t i = 0; i < m_ev.m_plab.size(); ++i) vv.push_back(m_ev.m_bornMomenta[i]);
    /*
      Every lab final leg, not just two. BuildAcraicMomenta wants
      {beams, Born final state, lab final state} and decides it has the lab
      set by counting - so appending only legs 2 and 3 of a longer final state
      left it one short, and it fell back to the BORN legs, which do not
      balance against the photons. That surfaced as Comix refusing the real
      with "momenta do not balance", several layers away.
    */
    /*
      With hidden FSR photons the final legs of m_plab balance against ALL the
      FSR photons, while ACRAIC is handed only the visible ones: the momentum
      repair then moved the softest visible photons, and the coherent IF soft
      factor of photons at x ~ 1e-7 came out wrong by up to 77% on ~20% of them
      (Z-pole mu mu, 2026-10-02). m_acraiclegs has the hidden photons absorbed
      into the pairs, so legs and visible photons balance.
    */
    const Vec4D_Vector &finlegs(m_ev.m_acraiclegs.size() == m_ev.m_plab.size()
                                ? m_ev.m_acraiclegs : m_ev.m_plab);
    for(size_t i = 2; i < m_ev.m_plab.size(); ++i)
      vv.push_back(finlegs[i]);
    p_acraic->Init(vv);
    /*
      The phase-space point that actually realises s': the ISR-reduced beams
      in the SAME frame as the recoiled final legs. m_plab[0..1] have been
      boosted to the lab while m_plab[2..] have not, so m_plab itself does not
      balance; m_sprimeBeams is the pre-boost pair.
    */
    if (m_ev.m_sprimeBeams.size() == 2) {
      Vec4D_Vector sp(m_ev.m_sprimeBeams);
      for (size_t i = 2; i < m_ev.m_plab.size(); ++i) sp.push_back(finlegs[i]);
      p_acraic->SetLabMomenta(sp);
    }
    p_acraic->SetPreFSRMomenta(m_ev.m_reallab);
    p_acraic->SetISRPhotons(m_ev.m_ISRPhotons);
    if (HasFSR()) p_acraic->SetFSRPhotons(m_ev.m_FSRPhotons);
    p_acraic->SetBornMomenta(m_ev.m_bornMomenta);
    p_acraic->SetISRFormFactor(m_ev.m_formfactor);
    /*
      The radiating final-state dipoles are ACRAIC's final stages: one per
      resonance pair, so that each resonance propagator carries its own
      partition-shifted invariant. Leg labels are positions in the process
      flavour list.
    */
    std::vector<std::vector<int> > groups;
    for (auto &D : p_dipoles->GetDipoleFF())
      if (D.IsResonance()) groups.push_back({D.Left(), D.Right()});
    // the generator's own stages, for the crude (ACRAIC: CRUDE_FROM_GENERATOR)
    p_acraic->SetCrudeStageGroups(groups);
    if (!SetACRAICWStageGroups()) p_acraic->SetStageGroups(groups);
    p_acraic->Calculate();
  }

}

/*
  ACRAIC's W stages for a W+W- -> 4f final state (ACRAIC: W_STAGES, see
  NOTES-w-stages-2026-09-27.md): one decay group per W, holding the W as the
  incoming leg and its charged daughter, from the same pairing the pole
  scheme uses (DipoleSet::FindWW: each charged lepton with the neutrino of
  its own generation). ACRAIC puts the W's, outgoing, into the production
  stage itself. Independent of YFS: WW_Scheme - the generator's dipoles are
  whatever they are, the crude has to be theirs (section 4.8 of the notes).
  Returns false when W stages are off or the state is not W+W- -> 4f, and
  the flat radiating-dipole groups are handed over as before.
*/
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
    MakeWWVecs(m_ev.m_bornMomenta);
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
  // Address-of is deliberate: ISR keeps this pointer until MakeYFS rebuilds
  // the dipole (see m_isrinital), which is why DipoleSet holds its dipoles by
  // unique_ptr so they never move.
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
  m_ev.m_isrphotonsforME = m_ev.m_ISRPhotons;
  m_isrWeight = res.weight;
  m_ev.m_photons = res.photons;
  m_ev.m_sprimeBeams.assign(2, Vec4D());
  for(size_t i = 0; i < 2; ++i) {
    m_ev.m_plab[i] = p_dipoles->GetDipoleII().GetNewMomenta(i); 
    m_ev.m_sprimeBeams[i] = m_ev.m_plab[i];   // before ToLab: the pair frame
    ToLab(m_ev.m_plab[i]);
  }
  double sp = (m_ev.m_plab[0] + m_ev.m_plab[1]).Abs2();
  // An invariant, so it holds for boosted (asymmetric, beam-spread) pairs as
  // well; exempting them hid the stale initial-initial dipole (MakeYFS).
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
  if (m_CalForm) return;
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
  else {
    if(FixedOrder()==fixed_order::nlo){
      m_ev.m_formfactor = 1 + m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5);
    }
    else m_ev.m_formfactor = exp(m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
  }
}

bool YFS_Handler::CalculateFSR(){
  return CalculateFSR(m_ev.m_plab);
}

// Event-record copy of an FSR pair that has to absorb the photons
// HIDE_PHOTONS dropped. The pair was built recoiling against every generated
// photon, and the weights are computed at that point (m_plab together with the
// pre-hide m_me_photons), so only the copy written to the record is changed:
// the pair keeps its rest-frame direction and takes the total momentum
// p1 + p2 + k_hidden. Without it the record misses k_hidden ~ Emin =
// IR_CUTOFF/2 and fails the four-momentum check (relative accuracy
// sqrt(NUM_ACCURACY) = 1e-5) once Emin/sqrt(s) is above it: Beam_Remnants:None
// then rejects the event and the trials lower the generated cross section
// (Bhabha/Moller at sqrt(s) = 0.106 GeV, IR_CUTOFF 1e-5: 6.8% of events,
// sigma -7.9%; 2026-10-02).
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
  // update NLO momenta from PHASIC
  // m_ev.m_reallab should be used for 
  // fixed order corrections.
  // Final state eikonals should be constructed
  // for the final state momenta before emissions
  // of photons. 
  m_ev.m_FSRPhotons.clear();
  m_ev.m_fsrphotonsforME.clear();
  m_ev.m_reallab = p;
  m_ev.m_plab=p;
  // Pure-FSR mode never goes through MakeYFS, so CreatMomentumMap() (the only
  // place m_ev.m_inparticles/m_ev.m_outparticles get cleared) would otherwise never run
  // for this path, letting stale entries from earlier trials survive under
  // reused Particle* keys and leak into Signal_Processes::FillBlob via
  // GetOutParticles(). Reset it here on every call so it always starts from
  // the current born momenta.
  CreatMomentumMap();
  if(FixedOrder()==fixed_order::nlo && m_ev.m_ISRPhotons.size()!=0) {
    for(size_t i = 2; i < m_ev.m_plab.size(); ++i) m_ev.m_outparticles[m_particles[i]] = m_ev.m_plab[i];
    return true;
  }
  if(m_mode==yfsmode::isr) {
    // if(m_ev.m_ISRPhotons.size() < m_mingammaN){
    //   m_isrWeight=0;
    //   return false;
    // }
    return true;
  }
  m_fsrWeight=1;
  p_dipoles->MakeDipoles(m_flavs, m_ev.m_plab, m_ev.m_plab);
  // p_dipoles->CreateAllDipoles(m_flavs, m_ev.m_plab, m_ev.m_plab);
  CheckResonance();
  // p_dipoles->CreateAllDipoles(m_flavs, m_ev.m_plab, m_ev.m_plab);
  if(m_mode==yfsmode::isrfsr) {
    // Initial legs are the BORN beams, not the ISR-reduced ones. The
    // interference is between radiation off the incoming particles and off the
    // outgoing ones, so the initial leg of an initial-final pair is the
    // physical beam - which is also what KKMC's Yint uses (m_p1, m_p2 in
    // KKceex.cxx:315, the same momenta its Yisr = SForFac(alfpini, m_p1, m_p2)
    // uses). Final legs stay at m_ev.m_plab, i.e. after the ISR recoil and before
    // FSR emission. This does NOT match KKMC: its m_p3, m_p4 are the fermions
    // AFTER FSR (KKceex.cxx:197-198 builds them from m_Event->m_Qf1/m_Qf2,
    // which KKarFin.cxx:97-98 overwrites with the post-YFSfin m_q1, m_q2), and
    // Yint uses those (KKceex.cxx:334-337).
    //
    // Passing m_ev.m_plab for both also mixed frames once the beams were asymmetric:
    // CalculateISR() writes m_ev.m_plab[0..1] back through ToLab() while
    // m_ev.m_plab[2..] stay in the incoming-pair rest frame. m_ev.m_bornMomenta is in
    // that rest frame throughout, so the pair is now built in one frame - which
    // matters because Btilda depends on the leg energies, not just invariants.
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
  // event-record momenta that differ from m_plab: the legs of pairs that
  // absorbed hidden photons (AbsorbHiddenPhotons)
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
    ProbeFSRDensity(*Dip);
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
  // if(!CheckMomentumConservation()) return false;
  if(FixedOrder()==fixed_order::nlo){
    int totk = m_ev.m_ISRPhotons.size();
    if(m_nlo_fsr_photons) totk += m_ev.m_FSRPhotons.size();
    if(totk != 1) {
      CountFSRFailure("NLO photon count");
      if(totk > 1)
        msg_Error()<<"Wrong photon multiplicity at Fixed Order: "<<totk<<std::endl;
      return false;
    }
  }
  // if((m_ev.m_ISRPhotons.size() +  m_ev.m_FSRPhotons.size()) < m_mingammaN) {
  //   m_fsrWeight=0;
  //   return false;
  // }
  // CheckMasses();
  CountFSRFailure("CalculateFSR completed");
  return true;
}


/*
  YFS: EXCHANGE_T_MIN. A lower bound on |t| = |(p_in - p_out)^2| of every
  same-flavour in-out pair (an exchange line) of p, legs 0 and 1 incoming.
  The selectors only see the radiative event, so a lab-frame energy or angle
  cut does not keep YFS's own pairs away from t = 0. At fixed target (MUonE)
  an outgoing muon can come out equal to the beam after the ISR recoil while
  the electron passes a lab E_e cut: the initial-final form factor's
  ln(-t) times the photon-mass log then blows up exp(Y) (LO) or is NaN (the
  virtual subtraction). Checked on the reduced Born (the Born weight) and on
  the initial-final dipoles (beams against the outgoing legs after ISR).
  A generation cut: set it below the fiducial |t|.
*/
namespace {
  //! ACRAIC: WEIGHT_O0 - the O(alpha^0) ACRAIC column (rho_0/rho_crude) as a
  //! named weight, the like-for-like partner of KKMC's CEEX0 (diagnostic)
  bool AcraicO0Stream()
  {
    static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["WEIGHT_O0"].SetDefault(false).Get<bool>());
    return on;
  }
}

bool YFS_Handler::ExchangeLinesPassTMin(const Vec4D_Vector &p) const {
  if (m_exchange_tmin <= 0.) return true;
  const size_t n(std::min(p.size(), m_flavs.size()));
  for (size_t i(0); i < 2 && i < n; ++i) {
    if (m_flavs[i].Charge() == 0.) continue;
    for (size_t j(2); j < n; ++j)
      if (m_flavs[j] == m_flavs[i] && -(p[i] - p[j]).Abs2() < m_exchange_tmin)
        return false;
  }
  return true;
}

void YFS_Handler::CountFSRFailure(const std::string &what) {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["FSR_FAILURE_STATS"].SetDefault(0).Get<int>() != 0);
  if (on) ++m_fsrfail[what];
}

void YFS_Handler::MakeWWVecs(ATOOLS::Vec4D_Vector p) {
  m_ev.m_Wm *= 0;
  m_ev.m_Wp *= 0;
  Flavour_Vector wp, wm;
  for(size_t i = 2; i < p.size(); ++i)
  {
    if (m_flavs[i].IsAnti() && m_flavs[i].IntCharge()) {
      m_ev.m_Wp += m_ev.m_plab[i];
      wp.push_back(m_flavs[i]);
    }
    if (!m_flavs[i].IsAnti() && m_flavs[i].IntCharge()) {
      m_ev.m_Wm += m_ev.m_plab[i];
      wm.push_back(m_flavs[i]);
    }
    if (!m_flavs[i].IntCharge()) {
      if (m_flavs[i].IsAnti()) {
        m_ev.m_Wm += m_ev.m_plab[i];
        wm.push_back(m_flavs[i]);
      }
      else {
        m_ev.m_Wp += m_ev.m_plab[i];
        wp.push_back(m_flavs[i]);
      }
    }
  }
}


void YFS_Handler::CalculateCoulomb() {
  if (!m_coulomb) return;
  MakeWWVecs(m_ev.m_bornMomenta);
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

void YFS_Handler::CalculateBeta() {
  // Invalidate last event's NLO pieces first, so an early return cannot leave
  // them looking current. Zeroed as well as flagged: a stale value that is
  // never read is still a trap for the next person to add a weight here.
  m_ev.m_nlo_current = false;
  m_ev.m_nlo_real = m_ev.m_nlo_virtual = m_ev.m_nlo_rv = m_ev.m_nlo_rr = 0.;
  m_ev.m_nlo_vxr = 0.;
  m_ev.m_nlo_real_ifisoft = 0.;
  if(!m_rmode && !m_int_nlo) return;
  double realISR(0), realFSR(0);
  if (m_betaorder > 0) {
    if(m_real_only) {
      if(!m_no_born) m_ev.m_real = p_dipoles->CalculateEEX()+1;
      else m_ev.m_real = p_dipoles->CalculateEEX();
    }
    else if(m_virtual_only) {
      if(!m_no_born) m_ev.m_real = p_dipoles->CalculateEEXVirtual();
      else m_ev.m_real = p_dipoles->CalculateEEXVirtual()-1;
    }
    else {
      if(!m_no_born) m_ev.m_real = p_dipoles->CalculateEEX()+p_dipoles->CalculateEEXVirtual();
      else m_ev.m_real = p_dipoles->CalculateEEX()/m_born+p_dipoles->CalculateEEXVirtual()/m_born;
    }
    m_ev.m_eex = m_ev.m_real;
    if(IsNan(m_ev.m_eex)) m_ev.m_eex=0;
    // if(m_ev.m_real < 0) m_ev.m_real = 0;
    // m_ev.m_real /= m_born;
  }
  /*
    ACRAIC FIRST. Two things downstream need its result:
      - CalculateNLO() asks NLO_Base for the virtual, and when no
        Loop_Generator was named ACRAIC is what supplies it,
      - the nominal weight itself, if ACRAIC_WEIGHT is on.
    Both read a number ACRAIC has to have produced already.
  */
  double acraicfac(1.);
  bool   haveacraic(false);
  if (m_useacraic) {
    /*
      ACRAIC: TRACE_FACTOR_ABOVE. The beta_1 trace writes to std::cerr as the
      partition loop runs, before the factor is known; so cerr is captured
      into a buffer for the duration of MakeACRAIC and the buffer is printed
      only for the events whose factor exceeds the threshold.
    */
    static const double heavy(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                              ["TRACE_FACTOR_ABOVE"].Get<double>());
    std::ostringstream tracebuf;
    std::streambuf *cerrbuf(NULL);
    if (heavy > 0.) cerrbuf = std::cerr.rdbuf(tracebuf.rdbuf());
    MakeACRAIC();
    if (cerrbuf) std::cerr.rdbuf(cerrbuf);
    /*
      YFS: REAL_BORN_PHOTON_MULTICHANNEL, the ACRAIC half: rho_crude becomes the
      sum of the crudes of every assignment of photons to Born pair and ISR
      that the generator could have produced this final state with
      (NLO_Base::BornPhotonChannelSum, all Born-pair choices). on applies
      it, report only computes it for WEIGHT_PROBE. G = 1 exactly without
      Born photons or when no ISR photon passes the Born cuts.
    */
    m_bpmc_evtG = 1.; m_bpmc_evtn = 0;
    static const bornphotonmc::code bpmc(ATOOLS::Settings::GetMainSettings()["YFS"]
      ["REAL_BORN_PHOTON_MULTICHANNEL"].SetDefault(bornphotonmc::on)
      .Get<bornphotonmc::code>());
    if (bpmc != bornphotonmc::off && p_acraic && p_nlo && p_dipoles->HasDipoleII()
        && !m_ev.m_ISRPhotons.empty()) {
      InitNLO();
      const std::vector<Vec4D> isr(m_ev.m_ISRPhotons.begin(), m_ev.m_ISRPhotons.end());
      m_bpmc_evtG = p_nlo->BornPhotonChannelSum(m_ev.m_reallab, isr,
                                                p_dipoles->GetDipoleII(), -1,
                                                &m_bpmc_evtn);
      if (bpmc == bornphotonmc::on && m_bpmc_evtG != 1.)
        p_acraic->ScaleRhoCrude(m_bpmc_evtG);
    }
    if (heavy > 0. && p_acraic) {
      const double r0(p_acraic->GetRhoCrude()), r1(p_acraic->GetResult());
      if (r0 > 0. && r1/r0 > heavy) {
        std::cerr<<"@@@ HEAVY factor="<<r1/r0<<" nphot="<<p_acraic->NPhot()
                 <<" BRfactor="<<m_ev.m_real<<std::endl;
        p_acraic->PrintKinematics(std::cerr);
        std::cerr<<tracebuf.str()<<"@@@ HEAVY end"<<std::endl;
      }
    }
    if (p_acraic) {
      /*
        Denominator = the INCOHERENT partition sum (KKMC's RhoCrud), not rho0.

        rho0 is the coherent sum, so the partitions can cancel and it can come
        arbitrarily close to zero, while the crude weight this factor
        multiplies is a factorised eikonal times Born and does not. Measured at
        250 GeV with rho0 as the denominator: one event produced 1107 pb on a
        4.86 pb cross section, unweighting efficiency 2e-06. The incoherent sum
        is a sum of positive terms and cannot vanish.

        With a single partition - no FSR, so nothing to interfere - the two
        denominators are identical, which is why the Z-peak ISR results and the
        KKMC cross-check are unaffected.
      */
      const double r0(p_acraic->GetRhoCrude()), r1(p_acraic->GetResult());
      if (r0 > 0. && !IsBad(r1/r0)) { acraicfac = r1/r0; haveacraic = true; }
      else ++m_acraicstats.m_bad;
      // Published immediately: AcraicCompare() below reads it, and so does
      // GenerateWeight(). Assigning it only at the end of this function left
      // the comparison reading the PREVIOUS event's value.
      m_ev.m_acraicfactor = haveacraic ? acraicfac : 0.;
      ProbeGolden(r0, r1);
      ProbeAcraicWeight();
      if (p_nlo) p_nlo->SetAcraicVirtual(p_acraic->VirtualFactor());
    }
  }

  if(m_nlotype!=nlo_type::born && p_nlo) {
    if(m_no_born) m_ev.m_real=CalculateNLO()/m_born;
    else m_ev.m_real=(m_born+CalculateNLO())/m_born;
    m_ev.m_nlo_current = true;
    if (m_acraic_compare && haveacraic) AcraicCompare();
    ProbeAcraicCompare(haveacraic);
    ProbeReference(haveacraic);
  }

  /*
    The ACRAIC O(alpha) factor, kept whether or not it drives the nominal.

    m_ev.m_real is the same object from the EEX/fixed-order side - 1 + sum(beta)/Born
    there, rho1/rho0 here - so m_ev.m_acraicfactor/m_ev.m_real is the ratio that turns the
    nominal weight into the ACRAIC one, which is exactly what a named weight is.
    With ACRAIC_WEIGHT on, ACRAIC becomes the nominal and the column is 1.
  */
  if (haveacraic) {
    m_acraicstats.AddOverEex(m_ev.m_real != 0. ? acraicfac/m_ev.m_real : 0.);
    // m_ev.m_real is NOT overwritten here. Which correction drives the weight is
    // decided in GenerateWeight(), because only there is the IFI_Real term
    // known - and that term belongs to the EEX correction alone.
  }
}

/*!
  ACRAIC against the EXTERNAL providers, event by event.

  Both sides are O(alpha) corrections relative to the SAME Born, so they are
  directly comparable without any normalisation being matched by hand:

      external virtual / Born   <->   rho(Born+virtual)/rho(Born) - 1
      external real    / Born   <->   rho(Born+real)   /rho(Born) - 1

  The TOTAL row is the meaningful one: both sides are the O(alpha) correction
  factor to the same resummed Born, and it is what drives the weight.

  The split rows are INDICATIVE ONLY, and the real row (marked real*) is not a
  like-for-like comparison at all. ACRAIC decomposes at AMPLITUDE level -
  rho(B+R)/rho(B) - 1 is 2Re(B*R)/|B|^2 + |R|^2/|B|^2 - whereas NLO_Base's
  CalculateReal() returns the YFS-subtracted real ME, a squared object. They
  coincide only in the soft limit. Read the real row as "are these even the
  same size", not as a discrepancy.

  The ACRAIC side excludes IFI_Real by construction - its partition sum already
  contains the real initial-final interference - so the external side is taken
  without it too, or the comparison is against a different quantity.
*/
void YFS_Handler::AcraicCompare() {
  if (!p_acraic || !p_nlo || m_born == 0.) return;
  const double rc(p_acraic->RealFactor() - 1.);
  const double tc(m_ev.m_acraicfactor - 1.);
  /*
    The virtual BY DIFFERENCE, total minus real - KKMC's convention
    (Rho1(full) - Rho1(Born+real))/Rho0 - and the only one that is comparable
    with the external provider's additive decomposition.

    rho(B+V)/rho(B) - 1 is NOT the same object once the real emission is hard:
    it drops the V-R cross term. Measured on a dumped point with real/Born =
    -0.40, that definition gave 0.0920 where KKMC gave 0.0629; total - real
    gives 0.06287084 against KKMC's 0.0628708363144.
  */
  const double vc(tc - rc);
  const double ve(m_ev.m_nlo_virtual/m_born);
  const double re(m_ev.m_nlo_real   /m_born);
  const double te(ve + re);
  auto rel = [](double a, double b) {
    const double s(std::abs(a)+std::abs(b));
    return s > 0. ? std::abs(a-b)/s : 0.;
  };
  const double dv(rel(vc,ve)), dr(rel(rc,re)), dt(rel(tc,te));
  m_acraicstats.AddCompare(dv, dr, dt);
  const bool dumped(p_acraic->JustDumped());
  if (m_acraicstats.m_cmp_n <= (long)m_acraic_compare || dumped) {
    size_t ng(m_ev.m_ISRPhotons.size() + m_ev.m_FSRPhotons.size());
    msg_Out()<<std::setprecision(8)
             <<"@@@ ACRAICCMP"<<(dumped?"-DUMPED":"")<<" n="<<m_acraicstats.m_cmp_n<<" ngam="<<ng
             <<" born="<<m_born<<"\n"
             <<"    virtual  acraic="<<vc<<"  ext="<<ve<<"  reldiff="<<dv<<"\n"
             <<"    real*    acraic="<<rc<<"  ext="<<re<<"  reldiff="<<dr<<"\n"
             <<"    total    acraic="<<tc<<"  ext="<<te<<"  reldiff="<<dt<<"\n"
             <<"      virt pieces: raw/born="<<(m_born!=0.?p_nlo->m_virt_raw/m_born:0.)
             <<"  sub/born="<<(m_born!=0.?p_nlo->m_virt_subval/m_born:0.)
             <<"  raw-sub="<<(m_born!=0.?(p_nlo->m_virt_raw-p_nlo->m_virt_subval)/m_born:0.)
             <<"\n      formfactor="<<m_ev.m_formfactor<<"  log(FF)="<<(m_ev.m_formfactor>0.?log(m_ev.m_formfactor):0.)
             <<"\n";
    // The Born point of the virtual: cos theta of the first final-state leg
    // in the beam CM frame, and the virtual subtraction split by dipole type
    // (taken in NLO_Base::CalculateVirtual; the default, non-eps, s-channel
    // branch only).
    if (p_nlo->m_plab.size() > 2) {
      const Vec4D P(p_nlo->m_plab[0] + p_nlo->m_plab[1]);
      Poincare cms(P);
      Vec4D pa(p_nlo->m_plab[0]), pf(p_nlo->m_plab[2]);
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
        for (size_t i(2); i < m_ev.m_plab.size(); ++i) {
          if (m_flavs.size() > i && m_flavs[i].Charge() == 0.) continue;
          Vec4D pc(m_ev.m_plab[i]); cms.Boost(pc);
          thf = Min(thf, 1. - kc.CosTheta(pc));
        }
      };
      for (const Vec4D &k : m_ev.m_ISRPhotons) hardest(k);
      for (const Vec4D &k : m_ev.m_FSRPhotons) hardest(k);
      msg_Out()<<"      @@@ VPT cos="<<pf.CosTheta(pa)
               <<" sqrts="<<P.Mass()<<" raw/born="<<p_nlo->m_virt_raw/m_born
               <<" sub/born="<<p_nlo->m_virt_subval/m_born
               <<" subII="<<p_nlo->m_virt_sub_ii/m_born<<" subFF="<<p_nlo->m_virt_sub_ff/m_born
               <<" subIF="<<p_nlo->m_virt_sub_if/m_born
               <<" v_used="<<ve<<" dIB="<<p_acraic->ResonanceBornDeltaV()<<" ngam="<<ng<<" xmax="<<xmax<<" th1c="<<th1c<<" thf="<<thf
               <<" nISR="<<m_ev.m_ISRPhotons.size()<<" nFSR="<<m_ev.m_FSRPhotons.size()
               <<" xisr="<<[&]{ double x(0.); for (const Vec4D &k : m_ev.m_ISRPhotons) x = Max(x, 2.*k.E()/P.Mass()); return x; }()
               <<" xfsr="<<[&]{ double x(0.); for (const Vec4D &k : m_ev.m_FSRPhotons) x = Max(x, 2.*k.E()/P.Mass()); return x; }()
               <<" mff="<<(m_ev.m_plab.size() > 3 ? (m_ev.m_plab[2] + m_ev.m_plab[3]).Mass() : -1.)
               // the pair plus one photon (any label), the mass closest to M_Z: an
               // ISR-labelled photon that would be resonant FSR off the Z return
               <<" mffk="<<[&]{ double best(-1.), mz(Flavour(kf_Z).Mass());
                   if (m_ev.m_plab.size() < 4) return best;
                   const Vec4D q(m_ev.m_plab[2] + m_ev.m_plab[3]);
                   for (const auto *ks : {&m_ev.m_ISRPhotons, &m_ev.m_FSRPhotons})
                     for (const Vec4D &k : *ks) {
                       const double mk((q + k).Mass());
                       if (best < 0. || std::abs(mk - mz) < std::abs(best - mz)) best = mk;
                     }
                   return best; }()
               <<" acraicfac="<<m_ev.m_acraicfactor<<" eexfac="<<m_ev.m_real
               <<" re="<<re<<" rho0crude="<<(p_acraic->GetRhoCrude() > 0. ?
                                            p_acraic->GetResult0()/p_acraic->GetRhoCrude() : 0.)
               <<" rho01crude="<<(p_acraic->GetRhoCrude() > 0. ?
                                  p_acraic->GetResult01()/p_acraic->GetRhoCrude() : 0.)
               <<" rho1="<<p_acraic->GetResult()<<" rho0="<<p_acraic->GetResult0()
               <<" rhocrude="<<p_acraic->GetRhoCrude()<<" born="<<m_born
               <<" y_r="<<p_nlo->m_dbg_r<<" y_flux="<<p_nlo->m_dbg_flux
               <<" y_subloc="<<p_nlo->m_dbg_subloc<<" y_subb="<<p_nlo->m_dbg_subb
               <<" y_tot="<<p_nlo->m_dbg_tot<<" y_dk="<<p_nlo->m_dbg_dk
               <<" y_dp="<<p_nlo->m_dbg_dp
               <<" x2pt="<<p_nlo->m_dbg_x2pt<<" x2gen="<<p_nlo->m_dbg_x2gen
               <<" x2isr="<<p_nlo->m_dbg_x2isr<<" facisr="<<p_nlo->m_dbg_facisr
               <<" facfsr="<<p_nlo->m_dbg_facfsr<<"\n";
      p_nlo->m_dbg_x2pt = p_nlo->m_dbg_x2gen = p_nlo->m_dbg_x2isr = 0.;
      p_nlo->m_dbg_facisr = p_nlo->m_dbg_facfsr = 0.;
      // one-photon events: the momenta, for an offline soft-factor check
      if (ng == 1) {
        const Vec4D k(m_ev.m_ISRPhotons.size() ? m_ev.m_ISRPhotons[0] : m_ev.m_FSRPhotons[0]);
        msg_Out()<<std::setprecision(16)<<"      @@@ VPTMOM k="<<k[0]<<","<<k[1]<<","<<k[2]<<","<<k[3];
        for (size_t i(0); i < m_ev.m_plab.size(); ++i)
          msg_Out()<<" lab"<<i<<"="<<m_ev.m_plab[i][0]<<","<<m_ev.m_plab[i][1]<<","<<m_ev.m_plab[i][2]<<","<<m_ev.m_plab[i][3];
        for (size_t i(0); i < p_nlo->m_plab.size(); ++i)
          msg_Out()<<" nlo"<<i<<"="<<p_nlo->m_plab[i][0]<<","<<p_nlo->m_plab[i][1]<<","<<p_nlo->m_plab[i][2]<<","<<p_nlo->m_plab[i][3];
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
    const Vec4D_Vector &cph(p_acraic->AllPhotonsLab());
    const std::vector<std::pair<Vec4D,double> > &eph(m_ev.m_extrealphot);
    /*
      The two sides are NOT handed the same photons. ACRAIC gets m_ev.m_ISRPhotons +
      m_ev.m_FSRPhotons; p_nlo->m_photons is m_ev.m_ISRPhotons + m_ev.m_me_photons, and the
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
    msg_Out()<<"\n        nISR="<<m_ev.m_ISRPhotons.size()
             <<" nFSR="<<m_ev.m_FSRPhotons.size()
             <<" nFSRforME="<<m_ev.m_fsrphotonsforME.size()<<"\n";
    /*
      Which list belongs to THIS event? Momentum conservation decides it with
      no reference to either matrix element:  p_a + p_b - q_c - q_d - sum k
      must vanish. A list captured in the wrong frame, or at the wrong stage
      of the dipole's boost, cannot balance.
    */
    auto bal = [&](const Vec4D &sum) {
      Vec4D b(m_ev.m_plab[0] + m_ev.m_plab[1] - sum);
      for (size_t i(2); i < m_ev.m_plab.size(); ++i) b -= m_ev.m_plab[i];
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
      const double rce(p_acraic->RealFactorPhoton(a));
      const double rex(ok && m_born != 0. ? eph[match].second/m_born : 0.);
      msg_Out()<<"      gam["<<a<<"] E="<<cph[a].E()
               <<(a < m_ev.m_ISRPhotons.size() ? " ISR" : " FSR")
               <<"  acraic="<<rce
               <<(ok ? "  ext=" : "  ext=(no match) ")<<rex
               <<"  reldiff="<<(ok && std::abs(rce)+std::abs(rex) > 0. ?
                                std::abs(rce-rex)/(std::abs(rce)+std::abs(rex)) : -1.)
               <<"\n";
    }
  }
}

void YFS_Handler::InitNLO(){
  p_nlo->Init(m_flavs,m_ev.m_reallab,m_ev.m_bornMomenta);
  p_nlo->m_bpmc_vmax = m_vmax; p_nlo->m_bpmc_s = m_s;
  p_nlo->SetPostEmissionMomenta(m_ev.m_plab);
  ProbeRealStab();
  p_nlo->p_dipoles = p_dipoles.get();
  p_nlo->SetBorn(m_born);
  p_nlo->SetFSR(p_fsr.get());
  p_nlo->m_ISRPhotons = m_ev.m_ISRPhotons;
  if (m_nlo_fsr_photons)
    p_nlo->m_FSRPhotons = m_nlo_fsr_from_event ? m_ev.m_FSRPhotons : m_ev.m_fsrphotonsforME;
  else
    p_nlo->m_FSRPhotons.clear();
  // Mirror the two lines above: same photons, now carrying their dipole.
  p_nlo->m_photons.clear();
  if (p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_ev.m_ISRPhotons)
      p_nlo->m_photons.push_back(YFS::Photon(k, &p_dipoles->GetDipoleII()));
  if (m_nlo_fsr_photons) {
    // m_ev.m_photons holds the event-record photons WITH their dipoles: ISR first,
    // then the FF ones, so skip the ISR head to take only the final-state tail.
    if (m_nlo_fsr_from_event) {
      const size_t skip(p_dipoles->HasDipoleII() ? m_ev.m_ISRPhotons.size() : 0);
      for (size_t i(skip); i < m_ev.m_photons.size(); ++i)
        p_nlo->m_photons.push_back(m_ev.m_photons[i]);
    } else {
      for (const YFS::Photon &k : m_ev.m_me_photons) p_nlo->m_photons.push_back(k);
    }
  }
}

/*
  ACRAIC: REAL_VIRTUAL 2. The ACRAIC photon list (m_ev.m_ISRPhotons +
  m_ev.m_FSRPhotons, AllPhotonsLab) and YFS.NLO's (m_ISRPhotons +
  m_me_photons) are built independently, so each ACRAIC photon takes the dv of
  the NLO photon with the same lab momentum (to 1e-6 of its energy); photons
  without a match or below RV_SOFT_CUT keep dv = 0 (their remainder -> 0 in
  the soft limit anyway), and are counted. YFS: RV_PROBE 1 prints per event
  "@@@ ACRV": ACRAIC's addition over rho_crude next to YFS.NLO's rv/B.
*/
double YFS_Handler::AcraicRealVirtualRemainder()
{
  if (!p_acraic || !p_nlo || !p_nlo->RVRemainder()) return 0.;
  const Vec4D_Vector &cph(p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(p_nlo->m_rvdv);
  std::vector<double> dv(cph.size(), 0.);
  int nmatch(0);
  for (size_t b(0); b < nd.size(); ++b) {
    long a(-1); double best(1e30);
    for (size_t c(0); c < cph.size(); ++c) {
      const Vec4D d(cph[c] - nd[b].first);
      const double m(Max(Max(dabs(d[0]),dabs(d[1])),Max(dabs(d[2]),dabs(d[3]))));
      if (m < best) { best = m; a = (long)c; }
    }
    if (a < 0 || !(best < 1e-6*Max(nd[b].first.E(), 1e-30))) { ++m_acrv_nomatch; continue; }
    if (!p_acraic->PhotonHasM1(a)) { ++m_acrv_nom1; continue; }
    dv[a] = nd[b].second;
    ++nmatch;
  }
  ++m_acrv_n;
  const double r(p_acraic->RealVirtualRemainderRho(dv));
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  if (probe) {
    const double rcr(p_acraic->GetRhoCrude());
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ ACRV nph="<<cph.size()<<" nnlo="<<nd.size()
     <<" matched="<<nmatch<<" acrv="<<(rcr > 0. ? r/rcr : 0.)
     <<" nlorv="<<(m_born != 0. ? m_ev.m_nlo_rv/m_born : 0.)
     <<" rho1/rcr="<<(rcr > 0. ? p_acraic->GetResult()/rcr : 0.)
     <<" nloreal="<<(m_born != 0. ? 1. + m_ev.m_nlo_real/m_born : 0.)
     <<" dv="<<(nd.size() == 1 ? nd[0].second : 0.)
     <<" x="<<(nd.size() == 1 ? 2.*nd[0].first.E()/sqrt(m_s) : 0.)
     <<" nisr="<<m_ev.m_ISRPhotons.size()<<"\n";
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
double YFS_Handler::AcraicRealVirtualHelicityRemainder()
{
  if (!p_acraic || !p_nlo || !p_nlo->RVRemainder()) return 0.;
  static const bool force(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["RV_HEL_FORCE_AVERAGE"].SetDefault(false).Get<bool>());
  const Vec4D_Vector &cph(p_acraic->AllPhotonsLab());
  const std::vector<std::pair<Vec4D,double> > &nd(p_nlo->m_rvdv);
  const std::vector<RVHelicityFactors> &nh(p_nlo->m_rvhel);
  const int nhel(Amplitude::NHel());
  std::vector<double> dv(cph.size(), 0.);
  std::vector<std::vector<Complex> > r(cph.size(),
                                       std::vector<Complex>(nhel, Complex(0., 0.)));
  bool anyhel(false);
  for (size_t b(0); b < nd.size(); ++b) {
    long a(-1); double best(1e30);
    for (size_t c(0); c < cph.size(); ++c) {
      const Vec4D d(cph[c] - nd[b].first);
      const double m(Max(Max(dabs(d[0]),dabs(d[1])),Max(dabs(d[2]),dabs(d[3]))));
      if (m < best) { best = m; a = (long)c; }
    }
    if (a < 0 || !(best < 1e-6*Max(nd[b].first.E(), 1e-30))) { ++m_acrv_nomatch; continue; }
    if (!p_acraic->PhotonHasM1(a)) { ++m_acrv_nom1; continue; }
    dv[a] = nd[b].second;
    const Complex half(0.5*dv[a], 0.);
    if (force || b >= nh.size() ||
        !p_acraic->RealVirtualHelicityToAcraic(nh[b], p_acraic->PhotonHelicity(a),
                                               half, r[a]))
      r[a].assign(nhel, half);
    else anyhel = true;
  }
  ++m_acrv_n;
  const double rhel(p_acraic->HelicityRealVirtualRemainderRho(r));
  const double ravg(p_acraic->RealVirtualRemainderRho(dv));
  const double rcr(p_acraic->GetRhoCrude());
  if (rcr > 0.) p_acraic->AccumulateRealVirtualHelicity(anyhel, rhel/rcr, ravg/rcr);
  static const bool probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_PROBE"]
                          .SetDefault(0).Get<int>() != 0);
  if (probe) {
    std::ostringstream o;
    o<<std::setprecision(12)<<"@@@ ACRVH nph="<<cph.size()<<" nnlo="<<nd.size()
     <<" anyhel="<<anyhel<<" hel="<<(rcr > 0. ? rhel/rcr : 0.)
     <<" avg="<<(rcr > 0. ? ravg/rcr : 0.)
     <<" nlorv="<<(m_born != 0. ? m_ev.m_nlo_rv/m_born : 0.)
     <<" rho1/rcr="<<(rcr > 0. ? p_acraic->GetResult()/rcr : 0.)<<"\n";
    std::cerr<<o.str();
  }
  return rhel;
}

double YFS_Handler::CalculateNLO(){
// CheckMomentumConservation();
  InitNLO();
  // one-shot fixed-point dump for the KKMC CEEX comparison (YFS: ACRAIC_Compare)
  p_nlo->ACRAICComparePoint();
  InitNLO();
  m_ev.m_nlo_real = p_nlo->CalculateReal();
  if (m_acraic_compare) {
    m_ev.m_extrealphot.clear();
    for (const YFS::Photon &g : p_nlo->m_photons)
      m_ev.m_extrealphot.push_back(std::make_pair(g.K(), g.beta10()));
  }
  // Hardest-photon-only contributions are captured as a side effect of the
  // nominal sums above (see NLO_Base::CalculateReal/CalculateRealVirtual/
  // CalculateRealReal) - no extra ME evaluation needed here.
  m_ev.m_nlo_real_hardest = p_nlo->m_real_hard1;
  m_ev.m_nlo_real_2hardest = p_nlo->m_real_hard2;
  // read now: the double real calls CalculateReal(k) again, moving m_ifi_prod
  m_ev.m_nlo_real_ifisoft = p_nlo->m_real_ifi_soft;
  InitNLO();
  m_ev.m_nlo_virtual = p_nlo->CalculateVirtual();
  // YFS: VV_POLE_CHECK with a double-virtual provider (diagnostic)
  if (p_nlo->HasVV()) p_nlo->CheckVVPoles();
  InitNLO();
  m_ev.m_nlo_rv = p_nlo->CalculateRealVirtual();
  m_ev.m_nlo_rv_hardest = p_nlo->m_rv_hard1;
  m_ev.m_nlo_rv_2hardest = p_nlo->m_rv_hard2;
  InitNLO();
  m_ev.m_nlo_rr = p_nlo->CalculateRealReal();
  m_ev.m_nlo_rr_2hardest = p_nlo->m_rr_hard2;
  // Everything above the double real. This has to live HERE, not in
  // NLO_Base::CalculateNLO(): that function is never called - this is the
  // driver - so anything added to it silently does nothing.
  // CalculateRealMultiplicity returns 0 for any n without a provider, so a run
  // that leaves YFS: NLO_MAX_PHOTONS at its default of 2 pays nothing.
  m_ev.m_nlo_rn = 0.;
  m_ev.m_nlo_vxr = 0.;
  // NLO_MAX_PHOTONS bounds it: providers above it exist only for ACRAIC's
  // tree residuals (ACRAIC: TREE_BETA_MAX)
  const size_t nrmax(Min(p_nlo->MaxRealPhotons(),
                         YFS::NLO_Base::RequestedMaxRealPhotons()));
  for (size_t n(3); n <= nrmax; ++n) {
    InitNLO();
    m_ev.m_nlo_rn += p_nlo->CalculateRealMultiplicity(n);
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
  const bool vxr_real_only(p_nlo->RRSinglesProduct());
  const bool addvxr((vcomb == virtualcombine::product && !p_nlo->HasRealVirtual())
                    || p_nlo->RVRemainder());
  if (addvxr && m_born != 0.) {
    const double beyond_real(vxr_real_only ? 0. : m_ev.m_nlo_rr + m_ev.m_nlo_rn);
    const double vxr((m_ev.m_nlo_virtual/m_born)*(m_ev.m_nlo_real + beyond_real));
    if (!IsBad(vxr)) m_ev.m_nlo_vxr = vxr;
  }
  /*
    YFS: BETA0_2 - the squared-level beta_0^(2), v2 = |1 + delta + d2|^2
    - 1 - v with delta = v1/2 (NLO_Base::Beta02SquaredLevel): it carries
    the |m_1|^2 = (v1/2)^2 YFS.NLO lacks and 2 Re delta^(2). It goes where
    v goes: on the Born, and on the real where v multiplies the real, so
    that 1 + v + v2 is the Born-level virtual factor of the product.
  */
  m_ev.m_nlo_beta02 = 0.;
  if (p_nlo->Beta02Nlo().On() && p_nlo->Beta02Ok() && m_born != 0.) {
    const double v2(p_nlo->Beta02SquaredLevel());
    const double beyond_real(vxr_real_only ? 0. : m_ev.m_nlo_rr + m_ev.m_nlo_rn);
    const double b02(v2*(m_born + (addvxr ? m_ev.m_nlo_real + beyond_real : 0.)));
    if (!IsBad(b02)) m_ev.m_nlo_beta02 = b02;
  }
  // YFS: DYSON_VP - sum_h w_h (|1 + r_h|^2 - 1), placed like v2
  m_ev.m_nlo_vp = 0.;
  if (p_nlo->DysonVPOn() && p_nlo->DysonVPFactors().ok && m_born != 0.) {
    const double beyond_real(vxr_real_only ? 0. : m_ev.m_nlo_rr + m_ev.m_nlo_rn);
    const double vp(p_nlo->DysonVPFactors().v2
                    *(m_born + (addvxr ? m_ev.m_nlo_real + beyond_real : 0.)));
    if (!IsBad(vp)) m_ev.m_nlo_vp = vp;
  }
  ProbeNNLOEvent();
  const double total(m_ev.m_nlo_real + m_ev.m_nlo_virtual + m_ev.m_nlo_rv
                     + m_ev.m_nlo_rr + m_ev.m_nlo_rn + m_ev.m_nlo_vxr);
  const double with02(p_nlo->Beta02Nlo().On() ? total + m_ev.m_nlo_beta02 : total);
  return p_nlo->DysonVPOn() ? with02 + m_ev.m_nlo_vp : with02;
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

double YFS_Handler::AcraicBeta02Rho(double d2, double v, double r0x, double r0v) const
{
  return 2.*d2*(r0x + 0.5*v*r0v) + d2*d2*p_acraic->GetResult0();
}

/*
  YFS: BETA0_2_PROBE with ACRAIC: CHECK_XS - on the event written to
  acraic_point.dat, the photonic factors in KKMC's convention for the KKMC point
  harness (acraic-kkmc-point-harness): ISR at the beams' s and FSR at the lab
  final state's s', the pairs KKMC dresses (same_side), next to the event's
  own factors at the loop's Born point. KKMC's CEEX2 - CEEX1 Born factor is
  (F1ini2 - F1ini1) + (F1fin2 - F1fin1) = sum_d [f_d^2/2 + NLL_d] (no cross).
*/
void YFS_Handler::Beta02DumpProbe(double v) const
{
  static const int probe(ATOOLS::Settings::GetMainSettings()["YFS"]["BETA0_2_PROBE"]
                         .Get<int>());
  if (!probe || !p_acraic || !p_acraic->JustDumped() || !p_nlo) return;
  Beta02_Config kk(p_nlo->Beta02Acraic());
  kk.pairs = beta02pairs::same_side;
  kk.cross = false;
  kk.weak  = false;
  kk.mode  = beta02::llnll;
  const Photonic_Factors f(PhotonicFactors(m_flavs, p_acraic->m_pacraic, kk, 2));
  if (!f.ok) return;
  std::ostringstream o;
  o<<std::setprecision(16)<<"@@@ B02DUMP alpha="<<kk.alpha<<" c="<<kk.imc
   <<" cnll="<<Beta02NLLCoefficient()
   <<" sini="<<(p_acraic->m_pacraic[0] + p_acraic->m_pacraic[1]).Abs2()
   <<" sfin="<<(p_acraic->m_pacraic[2] + p_acraic->m_pacraic[3]).Abs2()<<"\n";
  for (size_t d(0); d < f.f.size(); ++d)
    o<<"@@@ B02DUMP pair "<<d<<" f="<<f.f[d].real()<<","<<f.f[d].imag()
     <<" nll="<<f.n[d].real()<<","<<f.n[d].imag()
     <<" step="<<(0.5*f.f[d]*f.f[d] + f.n[d]).real()<<","
     <<(0.5*f.f[d]*f.f[d] + f.n[d]).imag()<<"\n";
  const Complex step(Beta02Delta(f, f.deltaf, kk));
  kk.cross = true;
  const Complex full(Beta02Delta(f, f.deltaf, kk));
  o<<"@@@ B02DUMP kkstep="<<step.real()<<","<<step.imag()
   <<" withcross="<<full.real()<<","<<full.imag()
   <<" event_deltaf="<<p_nlo->Beta02Factors().deltaf.real()<<","
   <<p_nlo->Beta02Factors().deltaf.imag()
   <<" event_v1="<<p_nlo->VOneLoop()<<" v="<<v
   <<" event_d2ext="<<p_nlo->Beta02External()<<"\n";
  std::cerr<<o.str();
}


void YFS_Handler::GenerateWeight() {
  if (m_dump_dipoles) p_dipoles->DumpDipoles();
  AddFormFactor();
  if (m_mode == yfsmode::isrfsr) m_ev.m_yfsweight = m_isrWeight * m_fsrWeight;
  else if (m_mode == yfsmode::fsr) m_ev.m_yfsweight = m_fsrWeight;
  else m_ev.m_yfsweight = m_isrWeight;
  if (m_coulomb) m_ev.m_yfsweight *= p_coulomb->GetWeight();
  if (m_formWW) m_ev.m_yfsweight *= m_ev.m_ww_formfact; //*exp(m_ev.m_coulSub);
  // YFS: LO_PROBE - the LO factors before CalculateBeta, to see whether the
  // NLO pieces modify them (w_lo is formed after it)
  const double pre_yfsw(m_ev.m_yfsweight), pre_ff(m_ev.m_formfactor), pre_v(m_v),
               pre_born(m_born);
  CalculateBeta();
  ProbeLOPre(pre_yfsw, pre_ff, pre_v, pre_born);

  double wif = 1.;
  if (m_ifireal && m_mode == yfsmode::isrfsr && m_nlotype == nlo_type::born &&
      p_nlo && p_nlo->HasReal()) {
    Vec4D_Vector allphotons(m_ev.m_ISRPhotons);
    allphotons.insert(allphotons.end(), m_ev.m_FSRPhotons.begin(), m_ev.m_FSRPhotons.end());
    // Same eikonal ratios the real correction subtracted with (see
    // NLO_Base::m_wifterms); the event-dipole ratio only when the fixed-order
    // side did not see every photon.
    const std::vector<double> &wt(p_nlo->WIFTerms());
    static const bool wifnlo(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["IFI_REAL_FROM_NLO"].SetDefault(true).Get<bool>());
    if (wifnlo && wt.size() == allphotons.size() && !wt.empty()) {
      wif = 1.; for (double t : wt) wif *= t;
    } else
      wif = p_dipoles->RealIFWeight(allphotons);
  }
  // The Born-level YFS weight: ISR x FSR crude (plus Coulomb/WW if on) times
  // the form factor, with NO NLO correction applied. This is what YFS.LO has
  // to reproduce -- built here directly rather than recovered downstream as
  // 1/m_ev.m_real, so the LO column cannot inherit anything m_ev.m_real does.
  double w_lo = m_ev.m_yfsweight * m_ev.m_formfactor * (1.-m_v);
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
  ProbeLO();
  static const bool lonoifi(ATOOLS::Settings::GetMainSettings()["YFS"]
                            ["LO_WITHOUT_IFI"].SetDefault(true).Get<bool>());
  if (lonoifi && m_ifireal && m_ifisub && m_fullform >= 1 &&
      p_dipoles) {
    const double dy(p_dipoles->FormFactorSumIF(0.5*sqrt(m_s))
                    - p_dipoles->FormFactorSumIF());
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
  const double corr_eex (m_ev.m_real + (wif - 1.));
  double corr_acraic(m_ev.m_acraicfactor);            // 0 if ACRAIC produced nothing
  /*
    ACRAIC_Virtual: external - the provider's virtual in the ACRAIC column.

    ACRAIC's own virtual and KKMC's O(alpha^1) put the IR-subtracted virtual
    factor delta on the Born-level amplitude INSIDE the coherent sum:
        rho = sum_h |A_0 (1 + delta) + (A_1 - A_0)|^2
            = rho_1 + 2 Re(delta) Re<A_0, A_1> + |delta|^2 rho_0,
    with A_0 the Born-level partition sum and A_1 the full O(alpha^1) one.
    A squared-level provider gives v = V_sub/Born = 2 Re(delta) at O(alpha),
    helicity-summed, so with delta = v/2 (common to helicities and partitions,
    Im(delta) dropped) the column becomes
        (rho_1 + v Re<A_0,A_1> + (v/2)^2 rho_0) / rho_crude.

    Two simpler forms were tried and are wrong beyond O(alpha) by v * real:
    multiplying rho_1/rho_crude by (1 + v) put e+e- -> mu mu at 250 GeV
    (66-116 GeV, radiative return, real part of the factor about -0.3) 2.7%
    below ACRAIC's own virtual and KKMC; adding v rho_0/rho_crude put it 6.2%
    above. Re<A_0,A_1> lies between rho_0 and rho_1, where ACRAIC's own virtual
    landed (1.009 of KKMC). On photon-free events the provider's v and ACRAIC's
    own virtual agree to 0.34%.

    RV and RR are NOT added: ACRAIC here is O(alpha^1). Applied before the IF
    form-factor cut below, which multiplies the whole column. `acraic` (2 -> 2
    only) keeps ACRAIC's own amplitude-level virtual; `none` leaves the column
    at Born + real. Until 2026-09-25 `external` was read and never consumed.
  */
  /*
    ACRAIC: ORDER 2 keeps the O(alpha^1) column of the same event alongside
    (named weight ACRAIC_O1), composed with the external virtual exactly as an
    ORDER 1 run composes it, so ORDER2/ORDER1 is a same-event ratio.
    GetResultVV() is |A_v|^2 with A_v the amplitude (v/2) multiplies: A_0
    (= GetResult0(), bit for bit) unless ORDER2_VIRTUAL_ON_BETA1.
  */
  const bool acraico2(p_acraic && p_acraic->Order() == 2);
  double corr_acraic_o1(0.), corr_acraic_v0(0.), corr_acraic_v1(0.);
  // ACRAIC_O1V1: the O(alpha^1) column with (1+v/2) on the whole of A_1, the
  // O(alpha) composition of ACRAIC: REAL_VIRTUAL factorisable, on the same event
  double corr_acraic_o1v1(0.);
  const double corr_acraic_o0(AcraicO0Stream() && p_acraic && corr_acraic != 0.
                              && p_acraic->GetRhoCrude() > 0.
                              ? p_acraic->GetResult0()/p_acraic->GetRhoCrude() : 0.);
  double corr_acraic_o0_ifcut(corr_acraic_o0);
  /*
    ACRAIC_O2V2 (ORDER 2, external virtual): the nominal column with the virtual
    on the whole O(alpha^2) amplitude, sum_h |(1 + v/2) A_2|^2, instead of
    sum_h |A_2 + (v/2) A_v|^2. dv2_acraic is the difference, so
    ACRAIC_O2V2/Nominal - 1 is exactly the (v/2) beta_2 terms the nominal drops:
    O(alpha^3), not the full one-loop two-photon real, so a truncation
    estimate, not a correction. Everything else (beta_0^(2), VP, IF cut) is
    shared with the nominal.
  */
  double dv2_acraic(0.);
  /*
    ACRAIC: BETA0_2 - delta^(2) on A_0 (YFS/Main/Beta0_2.H), the increments of
    the nominal and the two ORDER 2 compositions, added after the virtual and
    before the IF cut; corr_acraic_nob02 keeps the column without it
    (ACRAIC_NO_B02). The O(alpha^1) column ACRAIC_O1 never carries it.
  */
  const bool acraicb02(p_acraic && p_nlo && p_nlo->Beta02Acraic().On());
  double b02_acraic(0.), b02_v0(0.), b02_v1(0.);
  /*
    YFS: DYSON_VP in ACRAIC: r(h) on A_0 alone, with delta^(2) when that is
    on too (one amplitude-level insertion e(h) = d2 + r(h)); vp_* are the
    parts beyond beta_0^(2) (ACRAIC_NO_VP keeps the column without them).
  */
  const bool acraicvp(p_acraic && p_nlo && p_nlo->DysonVPOn() &&
                    p_nlo->DysonVPFactors().ok);
  double vp_acraic(0.), vp_v0(0.), vp_v1(0.);
  if (acraico2 && corr_acraic != 0.) {
    const double rcr(p_acraic->GetRhoCrude());
    corr_acraic_o1 = rcr > 0. ? p_acraic->GetResult1()/rcr : 0.;
    corr_acraic_o1v1 = corr_acraic_o1;
    corr_acraic_v0 = corr_acraic_v1 = corr_acraic;
  }
  if (corr_acraic != 0. && (m_acraicvirtsrc == acraicvirt::external ||
                          m_acraicvirtsrc == acraicvirt::automatic) && p_acraic &&
      p_nlo && p_nlo->HasVirtual() && m_ev.m_nlo_current && m_born != 0.) {
    const double rcr(p_acraic->GetRhoCrude());
    /*
      ACRAIC: RESONANCE_EXPONENTIATION greco: the resonant initial-final log is
      in the resonance propagators of the ACRAIC amplitudes now, so it leaves
      the provider's virtual: v' = v - dI_B, dI_B its exact Born-weighted
      share at the Born point (Acraic_Base::BuildResonanceBornDelta; 0, and v
      bit for bit, with none). Everything below - beta_0^(2), Dyson VP, the
      ORDER 2 streams - takes v'.
    */
    const double v(m_ev.m_nlo_virtual/m_born - p_acraic->ResonanceBornDeltaV());
    const double add(rcr > 0. ? (v*p_acraic->GetResult01()
                                 + 0.25*v*v*p_acraic->GetResultVV())/rcr : 0.);
    if (!IsBad(add)) corr_acraic += add;
    else ++m_acraicstats.m_bad;
    if (acraicb02 && p_nlo->Beta02Ok() && rcr > 0.) {
      // delta^(2) of v' (its delta_f w must not see the resonant log either)
      const double d2(AcraicBeta02External(*p_nlo, p_acraic->ResonanceBornDeltaV()));
      b02_acraic = AcraicBeta02Rho(d2, v, p_acraic->GetResult0X(), p_acraic->GetResult0V())/rcr;
      if (acraico2) {
        b02_v0 = AcraicBeta02Rho(d2, v, p_acraic->GetResult02(), p_acraic->GetResult0())/rcr;
        b02_v1 = AcraicBeta02Rho(d2, v, p_acraic->GetResult02(), p_acraic->GetResult01Order1())/rcr;
      }
      Beta02DumpProbe(v);
    }
    if (acraicvp && rcr > 0.) {
      const double d2(acraicb02 && p_nlo->Beta02Ok()
                      ? AcraicBeta02External(*p_nlo, p_acraic->ResonanceBornDeltaV()) : 0.);
      std::vector<Complex> e(p_acraic->HelicityShares(p_nlo->DysonVPFactors()));
      for (Complex &x : e) x += d2;
      const std::vector<Complex> e2(e.size(), Complex(d2, 0.));
      const bool vb1(p_acraic->VirtualOnBeta1());
      // the label check: ACRAIC's A_0-weighted 2 Re r against the Born's
      static const int vpprobe(ATOOLS::Settings::GetMainSettings()["YFS"]
                               ["DYSON_VP_PROBE"].Get<int>());
      static long nvp(0);
      if (vpprobe > 0 && nvp < vpprobe) {
        ++nvp;
        const std::vector<Complex> r(p_acraic->HelicityShares(p_nlo->DysonVPFactors()));
        double lin(0.);
        for (const VP_Helicity &h : p_nlo->DysonVPFactors().hel) lin += h.w*2.*h.r.real();
        std::ostringstream o;
        o<<std::setprecision(10)<<"@@@ DYSONVPACRAIC acraic2Re="<<p_acraic->A0WeightedTwoRe(r)
         <<" born2Re="<<lin<<" rbar2Re="<<2.*p_nlo->DysonVPFactors().rbar.real()<<"\n";
        std::cerr<<o.str();
      }
      vp_acraic = (p_acraic->A0InsertionRho(e, 0.5*v, vb1)
                 - p_acraic->A0InsertionRho(e2, 0.5*v, vb1))/rcr;
      if (acraico2) {
        vp_v0 = (p_acraic->A0InsertionRho(e, 0.5*v, false)
                 - p_acraic->A0InsertionRho(e2, 0.5*v, false))/rcr;
        vp_v1 = (p_acraic->A0InsertionRho(e, 0.5*v, true)
                 - p_acraic->A0InsertionRho(e2, 0.5*v, true))/rcr;
      }
    }
    // ACRAIC: REAL_VIRTUAL 2 - the non-factorisable real-virtual, YFS.NLO's
    // per-photon v_{n+1} - v_B on each photon's M_1 (Acraic_Base.C)
    const acraicrv::code rvmode(YFS::Acraic_Base::RealVirtualMode());
    if ((rvmode == acraicrv::averaged || rvmode == acraicrv::helicity) && rcr > 0.) {
      // the provider's per-photon loop has the final-state vertex and the
      // initial-final boxes, which ISR only leaves out
      if (p_acraic->IsrOnly())
        THROW(fatal_error, "ACRAIC: REAL_VIRTUAL averaged takes the loop "
              "provider's real-virtual, which contains final-state and "
              "initial-final loops; YFS: MODE ISR (ISR only) needs ACRAIC: "
              "{REAL_VIRTUAL: factorisable} with ACRAIC's own virtual "
              "(YFS: ACRAIC_Virtual: acraic).");
      const double arv((rvmode == acraicrv::helicity ? AcraicRealVirtualHelicityRemainder()
                                                    : AcraicRealVirtualRemainder())/rcr);
      if (!IsBad(arv)) corr_acraic += arv;
    }
    if (acraico2) {
      const double add1(rcr > 0. ? (v*p_acraic->GetResult01Order1()
                                    + 0.25*v*v*p_acraic->GetResult0())/rcr : 0.);
      if (!IsBad(add1)) corr_acraic_o1 += add1;
      const double add1v1(rcr > 0. ? (v + 0.25*v*v)*p_acraic->GetResult1()/rcr : 0.);
      if (!IsBad(add1v1)) corr_acraic_o1v1 += add1v1;
      const double a0(rcr > 0. ? (v*p_acraic->GetResult02()
                                  + 0.25*v*v*p_acraic->GetResult0())/rcr : 0.);
      const double a1(rcr > 0. ? (v*p_acraic->GetResult12()
                                  + 0.25*v*v*p_acraic->GetResult1())/rcr : 0.);
      if (!IsBad(a0)) corr_acraic_v0 += a0;
      if (!IsBad(a1)) corr_acraic_v1 += a1;
      const double rho2(p_acraic->GetResult());   // |A_2|^2 at ORDER 2
      const double d2(rcr > 0. ? ((v + 0.25*v*v)*rho2
                                  - v*p_acraic->GetResult01()
                                  - 0.25*v*v*p_acraic->GetResultVV())/rcr : 0.);
      if (!IsBad(d2)) dv2_acraic = d2;
    }
  }
  /*
    ACRAIC_Virtual: helicity (ACRAIC) - the same insertion with a complex,
    helicity-resolved factor from the loop provider's helicity amplitudes
    (Recola), delta^h in place of v/2:
        rho_V = 1/4 sum_h |A_1^h + delta^h A_0^h|^2 ,
    added as (rho_V - rho_1)/rho_crude. delta^h comes from the loop call of
    this event (NLO_Base::BuildLoopHelicityFactors, Born point, one per
    event, applied to every partition); its Born-weighted 2 Re is v by
    construction. Helicities the provider cannot resolve, and whole events
    without helicity amplitudes, keep v/2, i.e. external.
    ACRAIC: VIRT_HEL_FORCE_HALF_V: 1 sets delta^h = v/2 everywhere: the
    regression against external, through this code path.
  */
  // With ACRAIC: ORDER 2 the stored amplitudes are the O(alpha^2) ones; the
  // helicity virtual is built and validated for ORDER 1 only.
  if (p_acraic) for (double &c : p_acraic->m_lhcol) c = 0.;
  if (corr_acraic != 0. && m_acraicvirtsrc == acraicvirt::helicity && p_acraic &&
      p_nlo && p_nlo->HasVirtual() && m_ev.m_nlo_current && m_born != 0.) {
    static const bool forcehalf(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["VIRT_HEL_FORCE_HALF_V"].SetDefault(0).Get<int>() != 0);
    const double rcr(p_acraic->GetRhoCrude());
    const double vfull(m_ev.m_nlo_virtual/m_born);
    std::vector<Complex> d;
    bool ok(!forcehalf && p_nlo->m_lhelok &&
            p_acraic->LoopHelicityToAcraic(p_nlo->m_lhel, p_nlo->m_la0,
                                       p_nlo->m_ldelta, p_nlo->m_lhelmom,
                                       Complex(0.5*vfull, 0.), d));
    if (!ok) d.assign(Amplitude::NHel(), Complex(0.5*vfull, 0.));
    /*
      ACRAIC: RESONANCE_EXPONENTIATION greco: the resonant log leaves delta^h
      helicity by helicity, exactly (Acraic_Base::BuildResonanceBornDelta),
      and v' = v - dI_B is its Born-weighted 2 Re; with none both are zero.
    */
    const double v(vfull - p_acraic->ResonanceBornDeltaV());
    if (p_acraic->ResonanceExponentiation()) {
      const std::vector<Complex> &dI(p_acraic->ResonanceBornDeltaHelicity());
      for (size_t f(0); f < d.size() && f < dI.size(); ++f) d[f] -= dI[f];
    }
    const double rho1(p_acraic->GetResult());
    const double add(rcr > 0. ? (p_acraic->HelicityVirtualRho(d) - rho1)/rcr : 0.);
    // the external form on the same event, and rho_V(v/2) against it
    const double ext(v*p_acraic->GetResult01() + 0.25*v*v*p_acraic->GetResult0());
    const std::vector<Complex> dh(Amplitude::NHel(), Complex(0.5*v, 0.));
    /*
      ACRAIC: BETA0_2 per helicity: delta^h -> delta^h + delta^(2)(delta^h),
      complex, with delta_f in the same phase convention (Beta0_2.H), so that
      the O(alpha^2) sum does not depend on VIRT_HEL_IM_C. Events without
      helicity factors use the external, real form on v/2.
    */
    std::vector<Complex> db(d), dhb(dh);
    const bool hb02(acraicb02 && p_nlo->Beta02Ok() && rcr > 0.);
    if (hb02) {
      const double d2ext(AcraicBeta02External(*p_nlo, p_acraic->ResonanceBornDeltaV()));
      for (size_t f(0); f < db.size(); ++f)
        db[f] += ok ? Beta02Delta(p_nlo->Beta02Factors(), d[f], p_nlo->Beta02Acraic())
                    : Complex(d2ext, 0.);
      // BETA0_2 two_loop: the photon-exchange two-loop boxes per helicity
      // (the external fallback above has them Born-weighted already)
      if (ok && p_nlo->Beta02Acraic().mode == beta02::two_loop &&
          !p_nlo->Beta02Factors().boxhel.empty()) {
        const std::vector<Complex> bx(p_acraic->PhysicalLabelsToAcraic
          (p_nlo->Beta02Factors().boxhel, p_nlo->Beta02Factors().boxd2, Complex(0., 0.)));
        for (size_t f(0); f < db.size() && f < bx.size(); ++f) db[f] += bx[f];
      }
      for (Complex &x : dhb) x += Complex(d2ext, 0.);
      b02_acraic = (p_acraic->HelicityVirtualRho(db) - p_acraic->HelicityVirtualRho(d))/rcr;
      Beta02DumpProbe(v);
    }
    /*
      ACRAIC: EW_HIGHER_ORDER_HELICITY - EW_HIGHER_ORDER (in v, so in the v/2
      fallback) is otherwise missing from delta^h: averaged puts v_HO/2 on
      every helicity, helicity GRIFFIN's own per-helicity factor. After
      delta^(2), whose delta_f w must not see this O(alpha^2) term.
    */
    double gho_acraic(0.);
    if (ok && rcr > 0. && p_nlo->EWHOHelicityMode() != ewhohel::none &&
        p_nlo->VEWHigherOrder() != 0.) {
      const Complex half(0.5*p_nlo->VEWHigherOrder(), 0.);
      std::vector<Complex> g(Amplitude::NHel(), half);
      if (p_nlo->EWHOHelicityMode() == ewhohel::helicity && p_nlo->m_ghelok)
        g = p_acraic->PhysicalLabelsToAcraic(p_nlo->m_ghel, p_nlo->m_gdelta, half);
      std::vector<Complex> dg(db);
      for (size_t f(0); f < dg.size() && f < g.size(); ++f) dg[f] += g[f];
      gho_acraic = (p_acraic->HelicityVirtualRho(dg) - p_acraic->HelicityVirtualRho(db))/rcr;
      db.swap(dg);
    }
    // YFS: DYSON_VP per helicity, on top of delta^h (+ delta^(2))
    const bool hvp(acraicvp && rcr > 0.);
    if (hvp) {
      const std::vector<Complex> r(p_acraic->HelicityShares(p_nlo->DysonVPFactors()));
      std::vector<Complex> dv(db);
      for (size_t f(0); f < dv.size() && f < r.size(); ++f) {
        dv[f] += r[f];
        dhb[f] += r[f];
      }
      vp_acraic = (p_acraic->HelicityVirtualRho(dv) - p_acraic->HelicityVirtualRho(db))/rcr;
      db.swap(dv);
    }
    const double reg(std::abs(p_acraic->HelicityVirtualRho(dh) - rho1 - ext)
                     /Max(std::abs(rho1), 1e-300));
    double split[3] = {0., 0., 0.};
    if (ok && rcr > 0.) {
      p_acraic->HelicityVirtualSplit(d, 0.5*v, split);
      for (double &x : split) x /= rcr;
    }
    p_acraic->AccumulateLoopHelicity(ok, add, rcr > 0. ? ext/rcr : 0., reg, split);
    // the factors in ACRAIC's index on the point dumped for the KKMC harness
    // (ACRAIC: CHECK_XS), whose KKPART/KKAMP lines use the same j1j2j3j4 labels
    if (p_acraic->JustDumped()) {
      std::ostringstream o;
      o<<std::setprecision(10)<<"@@@ VHELDUMP ok="<<ok<<" v="<<v<<"\n";
      for (int f(0); f < Amplitude::NHel() && f < (int)d.size(); ++f) {
        o<<"@@@ VHELDUMP ";
        for (int l(0); l < Amplitude::s_nlegs; ++l) o<<LegBitOf(f, l);
        o<<" delta="<<d[f].real()<<","<<d[f].imag()<<"\n";
      }
      std::cerr<<o.str();
    }
    /*
      The same event for the extra named weights: ACRAIC_EXTV (external
      virtual), ACRAIC_L/R and ACRAIC_EXTV_L/R (beam 1 fully left/right
      polarised, beam 2 unpolarised; the unpolarised column is their mean).
      Stored as absolute factors before the IF cut, which multiplies all.
    */
    const bool hextra(hb02 || hvp || gho_acraic != 0.);
    if (rcr > 0. && !IsBad(add) && !hextra) {
      p_acraic->m_lhcol[0] = corr_acraic + add;
      p_acraic->m_lhcol[1] = corr_acraic + ext/rcr;
      p_acraic->m_lhcol[2] = p_acraic->HelicityVirtualRho(d,  -1., 0.)/rcr;
      p_acraic->m_lhcol[3] = p_acraic->HelicityVirtualRho(d,  +1., 0.)/rcr;
      p_acraic->m_lhcol[4] = p_acraic->HelicityVirtualRho(dh, -1., 0.)/rcr;
      p_acraic->m_lhcol[5] = p_acraic->HelicityVirtualRho(dh, +1., 0.)/rcr;
    }
    // the same columns with beta_0^(2) in each (the external one on v/2)
    if (rcr > 0. && !IsBad(add) && hextra) {
      p_acraic->m_lhcol[0] = corr_acraic + add + b02_acraic + vp_acraic + gho_acraic;
      p_acraic->m_lhcol[1] = corr_acraic
        + (p_acraic->HelicityVirtualRho(dhb) - rho1)/rcr;
      p_acraic->m_lhcol[2] = p_acraic->HelicityVirtualRho(db,  -1., 0.)/rcr;
      p_acraic->m_lhcol[3] = p_acraic->HelicityVirtualRho(db,  +1., 0.)/rcr;
      p_acraic->m_lhcol[4] = p_acraic->HelicityVirtualRho(dhb, -1., 0.)/rcr;
      p_acraic->m_lhcol[5] = p_acraic->HelicityVirtualRho(dhb, +1., 0.)/rcr;
    }
    if (!IsBad(add)) corr_acraic += add;
    else ++m_acraicstats.m_bad;
    if (!IsBad(gho_acraic)) corr_acraic += gho_acraic;
  }
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
  double corr_acraic_nob02(corr_acraic);
  if (acraicb02 && corr_acraic != 0.) {
    if (!IsBad(b02_acraic)) corr_acraic += b02_acraic;
    if (acraico2 && !IsBad(b02_v0)) corr_acraic_v0 += b02_v0;
    if (acraico2 && !IsBad(b02_v1)) corr_acraic_v1 += b02_v1;
  }
  double corr_acraic_novp(corr_acraic);
  if (acraicvp && corr_acraic != 0.) {
    if (!IsBad(vp_acraic)) { corr_acraic += vp_acraic; corr_acraic_nob02 += vp_acraic; }
    if (acraico2 && !IsBad(vp_v0)) corr_acraic_v0 += vp_v0;
    if (acraico2 && !IsBad(vp_v1)) corr_acraic_v1 += vp_v1;
  }
  double corr_acraic_v2(acraico2 && corr_acraic != 0. ? corr_acraic + dv2_acraic : 0.);
  static const bool ifcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["IF_FORMFACTOR_CUT"].SetDefault(true).Get<bool>());
  if (ifcut && corr_acraic != 0. && m_ifisub && m_fullform >= 1 &&
      p_dipoles) {
    const double wgen(0.5*sqrt(m_s)*m_isrcut);
    const double dy(p_dipoles->FormFactorSumIF(wgen)
                    - p_dipoles->FormFactorSumIF());
    if (!IsBad(dy)) { corr_acraic *= exp(dy); corr_acraic_o1 *= exp(dy); corr_acraic_o1v1 *= exp(dy);
                      corr_acraic_o0_ifcut *= exp(dy);
                      corr_acraic_v0 *= exp(dy); corr_acraic_v1 *= exp(dy); corr_acraic_v2 *= exp(dy);
                      corr_acraic_nob02 *= exp(dy); corr_acraic_novp *= exp(dy); }
  }
  const bool   acraic_nom (m_acraic_weight && corr_acraic != 0.);
  ProbeFinalWeights(corr_eex, corr_acraic, wif, acraico2, corr_acraic_o1, corr_acraic_v0, corr_acraic_v1);
  m_ev.m_yfsweight *= acraic_nom ? corr_acraic : corr_eex;
  m_ev.m_yfsweight *= m_ev.m_formfactor*(1.-m_v);
  // What the named ACRAIC column has to divide by to become a ratio.
  m_ev.m_corr_nominal = acraic_nom ? corr_acraic : corr_eex;
  m_ev.m_corr_acraic    = corr_acraic;
  m_ev.m_corr_eex     = corr_eex;
  m_ev.m_corr_acraic_o1 = corr_acraic_o1;
  m_ev.m_corr_acraic_o0 = corr_acraic_o0_ifcut;
  m_ev.m_corr_acraic_o1v1 = corr_acraic_o1v1;
  m_ev.m_corr_acraic_v0 = corr_acraic_v0;
  m_ev.m_corr_acraic_v1 = corr_acraic_v1;
  m_ev.m_corr_acraic_v2 = corr_acraic_v2;
  m_ev.m_corr_acraic_nob02 = corr_acraic_nob02;
  m_ev.m_corr_acraic_novp = corr_acraic_novp;
  CheckInvariants();
  // Captured before the IsBad/negative-weight clamps below, since the named
  // weights are ratios against the weight the event actually carries.
  const double w_full = m_ev.m_yfsweight;
  ProbeLOFinal(w_lo, w_full, acraic_nom, corr_acraic, corr_eex);
  if(m_isr_debug) {
    Vec4D ele;
    for (int i = 2; i < m_flavs.size(); ++i)
    {
      if(IsEqual(m_flavs[i],kf_e)) {
        ele = m_ev.m_plab[p_dipoles->m_flav_label[m_flavs[i]]];
        p_beams->BoostBackLab(ele);
        p_debug->FillHist("Form_Factor_FS_Angle", ele.Theta()*1000,m_ev.m_formfactor,1);
      }
    }
  }
  DEBUG_FUNC("\nISR Weight = " << m_isrWeight << "\n" <<
             "  FSR Weight = " << m_fsrWeight << "\n" <<
             "  WW form Weight = " << m_ev.m_ww_formfact << "\n" <<
             "  Total form Weight = " << m_ev.m_formfactor << "\n" <<
             "  Coulomb Weight = " << p_coulomb->GetWeight() << "\n" <<
             " Coulomb Subtraction Weight = " << exp(m_ev.m_coulSub) << "\n" <<
             "Total Weight = " << m_ev.m_yfsweight << "\n");
  if(IsBad(m_ev.m_yfsweight)){
    msg_Error()<<"\nISR Weight = " << m_isrWeight << "\n" <<
             "  FSR Weight = " << m_fsrWeight << "\n" <<
             "  Form Factor = " << m_ev.m_formfactor << "\n" <<
             "  NLO  Correction = " << m_ev.m_real << "\n" <<
             "Total Weight = " << m_ev.m_yfsweight << "\n";
    m_ev.m_yfsweight = 0;
  }
  if(m_ev.m_yfsweight < 0 && m_skipNegWeights){
    msg_Debugging()<<"Skipping negative Weight in YFS"<<std::endl;
    m_ev.m_yfsweight=0;
    m_negskip++;
  }

  ProbeSign(w_full);
  BuildNamedWeights(w_lo, w_full);
}

// The YFS.* named weights: each is a ratio that turns the nominal weight into
// one of the truncated or reordered matchings. Split out of GenerateWeight,
// which was 216 lines with more than half of them here.
//
// w_lo   : the Born-level YFS weight, before any NLO correction
// w_full : the weight the event actually carries, before the clamps below
void YFS_Handler::BuildNamedWeights(double w_lo, double w_full) {
  // Build named NLO sub-weights. base_weight=1 so the nominal is unchanged.
  // YFSNLO  — Real + Virtual only (NLO denominator).
  // YFSNNLO — Real + Virtual + RealVirtual + RealReal (full NNLO denominator).
  m_ev.m_nlo_weightsmap = Weights_Map{1.0};
  Weights wyfs{1.0};
  bool any(false);

  std::vector<std::string> names;
  const bool hasnlo (p_nlo && p_nlo->HasNLO());
  const bool hasnnlo(p_nlo && p_nlo->HasNNLO());
  if (m_nlotype != nlo_type::born && (hasnlo || hasnnlo)) {
    if (hasnlo)  names.push_back("LO"), names.push_back("NLO");
    if (hasnnlo) names.push_back("NNLO");
    if (m_nlo_weight_breakdown) {
      if (hasnlo)
        for (const char *n : {"Real","Virtual","BR","BV","NLO_1g",
                              "NLO_2g","NLO_FixedOrder"})
          names.push_back(n);
      // EEX only when the beta expansion ran: m_ev.m_eex is assigned solely inside
      // CalculateBeta's `if (m_betaorder > 0)`, so at BETA:0 the column would
      // be a no-op 1.0 masquerading as a measurement. m_betaorder is
      // configuration, so gating the name on it keeps the set constant.
      if (hasnlo && m_betaorder > 0) names.push_back("EEX");
      if (hasnnlo)
        for (const char *n : {"RealVirtual","RealReal","NLO+RR","NLO+RV",
                              "NNLO_1g","NNLO_2g","NNLO_FixedOrder"})
          names.push_back(n);
      // the EEX estimate of the double virtual, retired by YFS: BETA0_2
      // (its LL part, without the NLL and the 100% envelope); never both
      if (hasnnlo && !p_nlo->Beta02Nlo().On())
        for (const char *n : {"VV_EEX","NNLO_VV","NNLO_VV_up","NNLO_VV_down"})
          names.push_back(n);
    }
  }
  /*
    ACRAIC gets its own column whenever ACRAIC is on, independently of the NLO
    names above: it is defined at Born level too, where none of those exist.
  */
  if (m_useacraic) names.push_back("ACRAIC");
  // ACRAIC nominal (YFS: ACRAIC_WEIGHT, default with ACRAIC on): the squared-level
  // YFS.NLO weight the event would otherwise carry, as a named column
  if (m_useacraic && m_acraic_weight) names.push_back("NLO_EW");
  // YFS: DIFF_BUDGET_STREAMS: the per-event ingredients (filled further down)
  if (DiffBudgetStreamsOn())
    for (const char *n : s_diffbudgetnames) names.push_back(n);
  // BETA0_2: the column without beta_0^(2) on the same event, per column
  if (m_useacraic && p_nlo && p_nlo->Beta02Acraic().On() &&
      m_acraicvirtsrc != acraicvirt::acraic)
    names.push_back("ACRAIC_NO_B02");
  if (p_nlo && p_nlo->Beta02Nlo().On()) names.push_back("NLO_EW_NO_B02");
  if (p_nlo && p_nlo->DysonVPOn()) {
    if (m_useacraic && m_acraicvirtsrc != acraicvirt::acraic) names.push_back("ACRAIC_NO_VP");
    names.push_back("NLO_EW_NO_VP");
  }
  // ACRAIC_Virtual: helicity - the external-virtual ACRAIC column on the same
  // events, and the beam-1-polarised columns, see GenerateWeight
  if (m_useacraic && m_acraicvirtsrc == acraicvirt::helicity)
    for (const char *n : {"ACRAIC_EXTV","ACRAIC_L","ACRAIC_R","ACRAIC_EXTV_L","ACRAIC_EXTV_R"})
      names.push_back(n);
  if (m_useacraic && AcraicO0Stream()) names.push_back("ACRAIC_O0");
  // ACRAIC: ORDER 2 also carries the O(alpha^1) ACRAIC column of each event
  if (m_useacraic && ATOOLS::Settings::GetMainSettings()["ACRAIC"]["ORDER"]
      .SetDefault(1).Get<int>() == 2) {
    names.push_back("ACRAIC_O1");
    names.push_back("ACRAIC_O1V1");
    names.push_back("ACRAIC_O2V0");
    names.push_back("ACRAIC_O2V1");
    names.push_back("ACRAIC_O2V2");
  }
  if (m_ladder_weights) {
    if (m_coulomb && p_coulomb) names.push_back("NoCoulomb");
    if (m_ifisub && m_fullform >= 1 &&
        FixedOrder() != fixed_order::nlo && p_dipoles)
      names.push_back("NoIFI");
  }
  for (const std::string &n : names) wyfs[n] = 1.;
  if (!names.empty()) any = true;
  m_wnames.clear();
  m_wnames.insert(names.begin(), names.end());
  if (m_useacraic && m_wnames.count("ACRAIC") && m_ev.m_corr_acraic != 0. &&
      !IsZero(m_ev.m_corr_nominal)) {
    // The ratio that turns the nominal weight into the ACRAIC one. Both are
    // O(alpha) correction FACTORS against the same crude, so everything else
    // in the weight cancels - including IFI_Real, which is in the denominator
    // and deliberately not in the numerator.
    const double r(m_ev.m_corr_acraic/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC"] = r;
    if (m_wnames.count("ACRAIC_EXTV") && p_acraic) {
      const char *nm[5] = {"ACRAIC_EXTV","ACRAIC_L","ACRAIC_R","ACRAIC_EXTV_L","ACRAIC_EXTV_R"};
      const double c0(p_acraic->m_lhcol[0]);
      for (int i(0); i < 5; ++i) {
        const double ri(c0 != 0. ? r*p_acraic->m_lhcol[i+1]/c0 : r);
        if (!IsBad(ri)) wyfs[nm[i]] = ri;
      }
    }
    if (m_wnames.count("ACRAIC_O0")) {
      const double r0(m_ev.m_corr_acraic_o0/m_ev.m_corr_nominal);
      if (!IsBad(r0)) wyfs["ACRAIC_O0"] = r0;
    }
    if (m_wnames.count("ACRAIC_O1")) {
      const double r1(m_ev.m_corr_acraic_o1/m_ev.m_corr_nominal);
      if (!IsBad(r1)) wyfs["ACRAIC_O1"] = r1;
      const double r1v1(m_ev.m_corr_acraic_o1v1/m_ev.m_corr_nominal);
      if (!IsBad(r1v1)) wyfs["ACRAIC_O1V1"] = r1v1;
      const double rv0(m_ev.m_corr_acraic_v0/m_ev.m_corr_nominal);
      const double rv1(m_ev.m_corr_acraic_v1/m_ev.m_corr_nominal);
      if (!IsBad(rv0)) wyfs["ACRAIC_O2V0"] = rv0;
      if (!IsBad(rv1)) wyfs["ACRAIC_O2V1"] = rv1;
      const double rv2(m_ev.m_corr_acraic_v2/m_ev.m_corr_nominal);
      if (!IsBad(rv2) && m_ev.m_corr_acraic_v2 != 0.) wyfs["ACRAIC_O2V2"] = rv2;
    }
  }
  if (m_wnames.count("NLO_EW") && !IsZero(m_ev.m_corr_nominal)) {
    const double r(m_ev.m_corr_eex/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW"] = r;
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
  if (DiffBudgetStreamsOn()) {
    std::map<std::string, double> q;
    for (const char *n : s_diffbudgetnames) q[n] = 0.;
    const double cnom(m_ev.m_corr_nominal);
    if (p_acraic && m_ev.m_corr_acraic != 0. && !IsZero(cnom) && m_born != 0. &&
        p_acraic->GetRhoCrude() > 0.) {
      double ifcut(1.);
      static const bool ifcuton(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                ["IF_FORMFACTOR_CUT"].Get<bool>());
      if (ifcuton && m_ifisub && m_fullform >= 1 && p_dipoles) {
        const double dy(p_dipoles->FormFactorSumIF(0.5*sqrt(m_s)*m_isrcut)
                        - p_dipoles->FormFactorSumIF());
        if (!IsBad(dy)) ifcut = exp(dy);
      }
      const double rcr(p_acraic->GetRhoCrude());
      std::vector<double> xs;   // signed photon x, + ISR / - FSR
      double sxi(0.), sxf(0.);
      for (const Vec4D &k : m_ev.m_ISRPhotons) { xs.push_back( 2.*k[0]/sqrt(m_s)); sxi += xs.back(); }
      for (const Vec4D &k : m_ev.m_FSRPhotons) { xs.push_back(-2.*k[0]/sqrt(m_s)); sxf -= xs.back(); }
      std::sort(xs.begin(), xs.end(),
                [](double a, double b) { return std::abs(a) > std::abs(b); });
      q["DB_CRUDE"] = 1.;
      q["DB_A1"]  = p_acraic->GetResult()/rcr*ifcut;
      q["DB_A0"]  = p_acraic->GetResult0()/rcr*ifcut;
      q["DB_V"]   = m_ev.m_nlo_virtual/m_born;
      q["DB_X"]   = p_acraic->ResonanceBornDeltaV();
      q["DB_R"]   = m_ev.m_nlo_real/m_born;
      q["DB_RIF"] = 1. + m_ev.m_nlo_real_ifisoft/m_born;
      q["DB_X1"]  = xs.size() > 0 ? xs[0] : 0.;
      q["DB_X2"]  = xs.size() > 1 ? xs[1] : 0.;
      q["DB_NPH"] = xs.size();
      q["DB_SXI"] = sxi;
      q["DB_SXF"] = sxf;
      q["DB_VSQ"] = q["DB_A1"]*q["DB_V"]*q["DB_V"]/4.;
      q["DB_REST"] = q["DB_A1"]*(1. + q["DB_V"]) - m_ev.m_corr_eex;
      q["DB_NPV"] = m_ev.m_corr_eex + q["DB_VSQ"];
      for (auto &kv : q) kv.second /= cnom;
    }
    for (const auto &kv : q) wyfs[kv.first] = IsBad(kv.second) ? 0. : kv.second;
  }
  if (m_wnames.count("ACRAIC_NO_B02") && m_ev.m_corr_acraic != 0. &&
      !IsZero(m_ev.m_corr_nominal)) {
    const double r(m_ev.m_corr_acraic_nob02/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC_NO_B02"] = r;
  }
  if (m_wnames.count("ACRAIC_NO_VP") && m_ev.m_corr_acraic != 0. &&
      !IsZero(m_ev.m_corr_nominal)) {
    const double r(m_ev.m_corr_acraic_novp/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["ACRAIC_NO_VP"] = r;
  }
  if (m_wnames.count("NLO_EW_NO_VP") && !IsZero(m_ev.m_corr_nominal) &&
      m_born != 0.) {
    const double r((m_ev.m_corr_eex - m_ev.m_nlo_vp/m_born)/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW_NO_VP"] = r;
  }
  if (m_wnames.count("NLO_EW_NO_B02") && !IsZero(m_ev.m_corr_nominal) &&
      m_born != 0.) {
    const double r((m_ev.m_corr_eex - m_ev.m_nlo_beta02/m_born)/m_ev.m_corr_nominal);
    if (!IsBad(r)) wyfs["NLO_EW_NO_B02"] = r;
  }
  if (m_ev.m_nlo_current && m_nlotype != nlo_type::born && !IsZero(m_ev.m_real) &&
      (p_nlo->HasNLO() || p_nlo->HasNNLO())) {
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
    const double nomden(IsZero(m_ev.m_corr_eex) ? m_ev.m_real
                        : m_ev.m_real*(m_ev.m_corr_nominal/m_ev.m_corr_eex));

    // Values only. Writing through operator[] CREATES a column, which is how
    // YFS.EEX kept appearing at BETA:0 after it was dropped from the
    // registration list -- and a column that is zero on every event comes back
    // as NaN once Rivet divides by its own sumOfWeights. So refuse any name the
    // configuration did not register, and leave a non-finite value at the
    // registered 1.0 rather than poisoning the column.
    auto emit = [this, &wyfs](const std::string &name, double value) {
      if (!m_wnames.count(name)) return;
      if (!IsBad(value)) wyfs[name] = value;
    };

    const bool have_fixed_order_ff =
        m_fullform >= 1 && !IsZero(m_ev.m_formfactor);
    // The cached exponent, not a fresh FormFactorSum(): the two must be the
    // same Y or this is not a truncation of anything.
    const double ff_fixedorder_ratio =
        have_fixed_order_ff ? (1. + m_ev.m_formfactor_sum) / m_ev.m_formfactor : 1.;

    /*
      YFS: IFI_Real true cuts the IF part of the exponent at IFIOmega()
      instead of sqrt(s)/2 and gives every photon above it the IF eikonal
      ratio r_j; the soft (beta_0) level paired with that exponent is
      B prod_j r_j. m_nlo_real carries B (prod_j r_j - 1); every column NOT
      built from m_nlo_real adds it back, or its exponent has lost the IF
      interference above omega (A_FB +0.15 at the Z pole). 0 with IFI_Real
      false, where all columns are unchanged.
    */
    const double ifisoft = m_ev.m_nlo_real_ifisoft;

    // NLO: Real + Virtual
    if (p_nlo->HasNLO()) {
      // m_nlo_vxr: YFS: VIRTUAL_COMBINE 1 (zero otherwise), so the named
      // column equals the nominal weight it is compared against
      const double nlo_sum   = (m_born + m_ev.m_nlo_real + m_ev.m_nlo_virtual
                                + m_ev.m_nlo_vxr)/m_born;
      const double real_sum  = (m_born + m_ev.m_nlo_real)/m_born;
      const double virt_sum  = (m_born + ifisoft + m_ev.m_nlo_virtual)/m_born;
      if (!IsZero(m_ev.m_real)) {
        emit("Real", ratio((m_ev.m_nlo_real)/m_born, nomden));
        emit("Virtual", ratio((m_ev.m_nlo_virtual + m_ev.m_nlo_vxr)/m_born, nomden));
        emit("NLO", ratio(nlo_sum, nomden));
        emit("BR", ratio(real_sum, nomden));
        emit("BV", ratio(virt_sum, nomden));
        // LO = (Born-level YFS weight) / (full weight), so that
        // nominal * YFS.LO == w_lo identically. Algebraically 1/m_ev.m_real, but
        // built from the two weights themselves.
        if (!IsZero(w_full)) emit("LO", w_lo/w_full);
        emit("EEX", ratio(m_ev.m_eex + ifisoft/m_born, nomden));
        // Matching truncated to a fixed real-photon multiplicity, to see the
        // result "as if" only the 1 or 2 hardest photons were used in the
        // matching (full "NLO" above keeps all generated photons). Real is
        // summed over the 1 / 2 hardest photons; Virtual is always full, and so
        // is the IF restoration, which belongs to the exponent, not to beta_1.
        // With one photon NLO_1g is NLO without m_nlo_vxr (either REAL_COMBINE).
        const double nlo_1g = (m_born + ifisoft + m_ev.m_nlo_real_hardest
                               + m_ev.m_nlo_virtual)/m_born;
        const double nlo_2g = (m_born + ifisoft + m_ev.m_nlo_real_2hardest
                               + m_ev.m_nlo_virtual)/m_born;
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
    if (p_nlo->HasNNLO()) {
      // YFS: BETA0_2's v2 belongs to the O(alpha^2) column (0 when off)
      const double nnlo_total = (m_born + m_ev.m_nlo_real + m_ev.m_nlo_virtual + m_ev.m_nlo_rv + m_ev.m_nlo_rr
                                 + m_ev.m_nlo_vxr + m_ev.m_nlo_beta02 + m_ev.m_nlo_vp)/m_born;
      const double RR_total = (m_born + m_ev.m_nlo_real + m_ev.m_nlo_virtual + m_ev.m_nlo_rr
                               + m_ev.m_nlo_vxr)/m_born;
      /*
        NLO+RV is YFS.NLO plus the real-virtual, so it carries m_nlo_vxr like
        YFS.NLO and NLO+RR. Without it the column was YFS.NLO + rv - v x real
        whenever vxr is set (YFS: RV_MODE remainder, or VIRTUAL_COMBINE product
        without an RV): Z-pole mu mu NLO+RV/NLO 1.003 with IFI_Real false and
        0.991 with a 2% cos theta tilt with true, where the real carries the IF
        restoration B (prod r_IF - 1) (2026-10-04). With the legacy RV vxr is 0,
        so nothing changes there.
      */
      const double RV_total = (m_born + m_ev.m_nlo_real + m_ev.m_nlo_virtual + m_ev.m_nlo_rv
                               + m_ev.m_nlo_vxr)/m_born;
      // No separate denominator: every column here is x/m_ev.m_real. NNLO used to be
      // make_ratio(nnlo_total, m_born + nnlo_total), which is
      // (m_born + nnlo_total)/(m_born + nnlo_total) -- identically 1, whatever
      // the physics did. It also added a dimensionful m_born to a ratio.
      if (!IsZero(m_ev.m_real)) {
        emit("RealVirtual", ratio((m_ev.m_nlo_rv)/m_born, nomden));
        emit("RealReal", ratio((m_ev.m_nlo_rr)/m_born, nomden));
        emit("NLO+RR", ratio(RR_total, nomden));
        emit("NLO+RV", ratio(RV_total, nomden));
        emit("NNLO", ratio(nnlo_total, nomden));
        // Matching truncated to a fixed real-photon multiplicity, the NNLO
        // analogue of NLO_1g/NLO_2g. 1 photon: Real + RealVirtual on the
        // single hardest, RealReal = 0 (a pair needs two photons). 2 photons:
        // Real + RealVirtual summed over the two hardest, plus the RealReal
        // pair formed by them. Virtual is always full.
        const double nnlo_1g =
            (m_born + ifisoft + m_ev.m_nlo_real_hardest + m_ev.m_nlo_virtual
             + m_ev.m_nlo_rv_hardest)/m_born;
        const double nnlo_2g =
            (m_born + ifisoft + m_ev.m_nlo_real_2hardest + m_ev.m_nlo_virtual
             + m_ev.m_nlo_rv_2hardest + m_ev.m_nlo_rr_2hardest)/m_born;
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
        if (!p_nlo->HasVV()) {
          const double vv = p_dipoles->CalculateEEXVirtual(2)
                          - p_dipoles->CalculateEEXVirtual(1);
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
            const double d = m_vv_approx_unc;
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

    if (p_fb) p_fb->SplitWeights(wyfs, m_ev.m_plab, m_flavs);
    any = true;
  }

  if (BuildLadderWeights(wyfs)) any = true;

  if (any) m_ev.m_nlo_weightsmap["YFS"] = wyfs;
}


bool YFS_Handler::BuildLadderWeights(ATOOLS::Weights &w) {
  if (!m_ladder_weights) return false;
  // Names are registered by BuildNamedWeights; this fills values only. An
  // event that cannot supply one leaves it at 1.0.
  bool any(false);

  if (m_coulomb && p_coulomb) {
    const double wc(p_coulomb->GetWeight());
    if (!IsZero(wc) && !IsBad(wc)) { w["NoCoulomb"] = 1./wc; any = true; }
  }

  if (m_ifisub && m_fullform >= 1 &&
      FixedOrder() != fixed_order::nlo && p_dipoles) {
    const double fif(p_dipoles->FormFactorSumIF());
    if (!IsBad(fif)) { w["NoIFI"] = exp(-fif); any = true; }
  }

  return any;
}


void YFS_Handler::YFSDebug(double W){
  p_debug->FillHist(m_ev.m_plab, p_isr.get(), p_fsr.get(), W);
}


void YFS_Handler::Reset() {
  // NOT the event boundary - that is StartEvent(), called from
  // SetBornMomenta. This is the failure path: MakeYFS and CalculateFSR call
  // it when an event is abandoned part-way, to drop the photons and zero the
  // weight while leaving the kinematics that were already handed in.
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

bool YFS_Handler::CheckMomentumConservation(){
  Vec4D incoming = m_ev.m_bornMomenta[0]+m_ev.m_bornMomenta[1];
  Vec4D outgoing;
  for(auto k: m_ev.m_ISRPhotons)  outgoing+=k;
  for(auto kk: m_ev.m_FSRPhotons) outgoing+=kk;
  for(size_t i = 2; i < m_ev.m_plab.size(); ++i)
  {
    outgoing+=m_ev.m_plab[i];
  }
  Vec4D diff = incoming - outgoing;
  if(!IsEqual(incoming,outgoing, 1e-5)){
    msg_Error()<<"Momentum not conserverd in YFS"<<std::endl
               <<"Incoming momentum = "<<incoming<<std::endl
               <<"Outgoing momentum = "<<outgoing<<std::endl
               <<"Difference = "<<diff<<std::endl
               <<"ISR Photons = "<<m_ev.m_ISRPhotons<<std::endl
               <<"FSR Photons = "<<m_ev.m_FSRPhotons<<std::endl;
  return false;
  }
  return true;
}


void YFS_Handler::CheckMasses(){
  bool allonshell=true;
  std::vector<double> mass;
  Vec4D_Vector p = m_ev.m_plab;
  for(auto k: m_ev.m_ISRPhotons) p.push_back(k);
  for(auto kk: m_ev.m_FSRPhotons) p.push_back(kk);

  for(size_t i = 0; i < p.size(); ++i)
  {
    if(i<m_ev.m_plab.size()){
      mass.push_back(m_flavs[i].Mass());
      if(!IsEqual(p[i].Mass(),m_flavs[i].Mass(),1e-5)){
        msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
                       <<"Flavour = "<<m_flavs[i]<<", with mass = "<<m_flavs[i].Mass()<<std::endl
                       <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
        allonshell = false;

      }
    }
    else{
      mass.push_back(0);
      if(!IsEqual(p[i].Mass(),0,1e-5)){
        msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
                       <<"Flavour = "<<Flavour(22)<<", with mass = "<<Flavour(22).Mass()<<std::endl
                       <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
        allonshell = false;
      }
    }
  }
  if(!allonshell) {
    m_stretcher.StretchMomenta(p, mass);
    for(size_t i = 0; i < m_ev.m_plab.size(); ++i)
    {
      msg_Debugging()<<"Mass after Mometum strechting"<<std::endl;
      if(i<m_ev.m_plab.size()){
         msg_Debugging()<<"Flavour = "<<m_flavs[i]<<", with mass = "<<m_flavs[i].Mass()<<std::endl
                       <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
      }
      else{
         msg_Debugging()<<"Flavour = "<<Flavour(22)<<", with mass = "<<Flavour(22).Mass()<<std::endl
                        <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
      }
      m_ev.m_plab[i] = p[i];
    }
  }
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

void YFS_Handler::CheckInvariants() const {
  if (!m_check_invariants) return;

  // Every photon knows the dipole it came from. A null one means a Photon was
  // built without it, and IsISR()/IsFSR() would dereference null.
  for (const YFS::Photon &k : m_ev.m_photons)
    if (!k.Dip())
      msg_Error() << METHOD << ": photon with no dipole; its origin cannot be "
                  << "determined." << std::endl;

  // E^2 >= m^2 for every dipole leg. Violated when an energy is paired with a
  // mass from a different source, which is what produced the WW NaN in
  // YFS_Form_Factor::A4 via sqrt(E^2 - m^2).
  YFS::DipoleView ffchk(p_dipoles->GetDipoleFF());
  for (auto D = ffchk.begin(); D != ffchk.end(); ++D)
    for (int i(0); i < 2; ++i) {
      const ATOOLS::Vec4D p(D->GetBornMomenta(i));
      if (p[0]*p[0] + 1e-9 < p.Abs2())
        msg_Error() << METHOD << ": dipole leg with E^2 < m^2, E=" << p[0]
                    << " m^2=" << p.Abs2() << std::endl;
    }

  // The weight the event carries must be a number. Catching it here names the
  // event; downstream it only shows up as a NaN cross section.
  if (ATOOLS::IsBad(m_ev.m_yfsweight))
    msg_Error() << METHOD << ": YFS weight is " << m_ev.m_yfsweight
                << " (isr=" << m_isrWeight << " fsr=" << m_fsrWeight
                << " form=" << m_ev.m_formfactor << " real=" << m_ev.m_real << ")"
                << std::endl;
  if (ATOOLS::IsBad(m_ev.m_formfactor))
    msg_Error() << METHOD << ": form factor is " << m_ev.m_formfactor << std::endl;
}
