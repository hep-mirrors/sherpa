#include "ATOOLS/Org/Message.H"
#include "YFS/Main/YFS_Handler.H"
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
    p_dipoles = new Define_Dipoles();
    p_coulomb = new Coulomb();
    p_fsr = new FSR();
    p_debug = new Debug();
    p_yfsFormFact = new YFS::YFS_Form_Factor();
    m_setparticles = false;
    p_isr = new YFS::ISR();
    p_nlo = nullptr;
    p_fb  = m_fb_analysis ? new YFS::YFS_FB_Analysis({}, m_fb_kf) : nullptr;
    m_formfactor = 1;
    m_isrinital = true;
    p_splitter = new PHOTONS::Photon_Splitter(m_photon_split);
    m_rmode = 0;
    m_real = 1;
    m_negskip = 0;
    m_nlo_real = 0.0;
    m_nlo_virtual = 0.0;
    m_nlo_rv = 0.0;
    m_nlo_rr = 0.0;
    m_nlo_real_hardest = 0.0;
    m_nlo_rv_hardest = 0.0;
    m_nlo_rr_2hardest = 0.0;
    m_nlo_real_2hardest = 0.0;
    m_nlo_rv_2hardest = 0.0;
    rpa->gen.AddCitation(1,"The automation of YFS ISR is published in  \\cite{Krauss:2022ajk}.Which is based on \\cite{Jadach:1988gb}");
  }
}

YFS_Handler::~YFS_Handler()
{
  if(Mode()!=YFS::yfsmode::off){
    if (p_isr) delete p_isr;
    if (p_fsr) delete p_fsr;
    if (p_coulomb) delete p_coulomb;
    if (p_debug)   delete p_debug;
    if (p_yfsFormFact) delete p_yfsFormFact;
    if (p_dipoles) delete p_dipoles;
    if (p_nlo) delete p_nlo;
    if (p_fb) delete p_fb;
    if (p_splitter) delete p_splitter;
    for (auto &p: m_particles){
      if(p) delete p;
    }
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
  }
}

NLO_Base *YFS_Handler::EnsureNLO()
{
  if (!p_nlo) p_nlo = new YFS::NLO_Base();
  return p_nlo;
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

void YFS_Handler::SetLimits(const double &smin) {
  double s = sqr(rpa->gen.Ecms());
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
  // delete before clearing: the old code cleared first, so the loop below ran
  // over an empty vector and leaked every Particle it was meant to free.
  for(auto particle : m_particles) {
    delete particle;
  }
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
      if (i < 2) m_particles.push_back(new ATOOLS::Particle(i, m_flavs[i], {0, 0, 0, 0}, 'i'));
      else    m_particles.push_back(new ATOOLS::Particle(i, m_flavs[i], {0, 0, 0, 0}, 'f'));
      m_particles[i]->ResetCounter();
    if (i >= 2) {
      if (flavs[i].IsQED()) qed = true;
    }
  }
  m_setparticles = true;
  if (m_useceex) InitializeCEEX(m_flavs);
}

void YFS_Handler::SetBornMomenta(const ATOOLS::Vec4D_Vector &p) {
  m_bornMomenta.clear();
  for(size_t i = 0; i < p.size(); ++i) {
    m_bornMomenta.push_back(p[i]);
  }
  // detect asymmetric beams from the original lab momenta, before any boost
  if(m_bornMomenta[0] != -m_bornMomenta[1]) m_asymbeams = true;
  else m_asymbeams = false;
  // NLO_Base::MapMomenta (and the ISR/FSR kinematics) assume the incoming pair
  // is at rest. For non-standard setups (fixed target, e.g. MUonE muon-e-;
  // asymmetric beams; beamstrahlung) it is not, so boost into the incoming-pair
  // rest frame here; the blob-facing getters (ToLab) undo it when handing the
  // event back. Pure-FSR mode feeds lab momenta straight into CalculateFSR(p)
  // and is left untouched.
  Vec4D Q(m_bornMomenta[0] + m_bornMomenta[1]);
  if (m_mode != yfsmode::fsr && !IsZero(Q.PSpat() / Q[0], 1e-10)) {
    m_cmsboost = Poincare(Q);
    for (size_t i = 0; i < m_bornMomenta.size(); ++i)
      m_cmsboost.Boost(m_bornMomenta[i]);
  } else {
    m_cmsboost = Poincare();
  }
  if (m_formWW) MakeWWVecs(m_bornMomenta);
  // AddFormFactor();
}

void YFS_Handler::SetMomenta(const ATOOLS::Vec4D_Vector &p) {
  m_plab.clear();
  for(size_t i = 0; i < p.size(); ++i) {
    Vec4D pi(p[i]);
    m_cmsboost.Boost(pi);
    m_plab.push_back(pi);
  }
}

void YFS_Handler::CreatMomentumMap() {
  m_inparticles.clear();
  m_outparticles.clear();
  for(size_t i = 0; i < 2; ++i)
  {
    m_inparticles[m_particles[i]] = m_bornMomenta[i];
    m_particles[i]->SetMomentum(m_bornMomenta[i]);
  }
  if(m_mode!=yfsmode::isr){
    for(size_t i = 2; i < m_flavs.size(); ++i)
    {
      m_outparticles[m_particles[i]] = m_bornMomenta[i];
      m_particles[i]->SetMomentum(m_bornMomenta[i]);
    }
  }
}

void YFS_Handler::InitializeCEEX(const ATOOLS::Flavour_Vector &fl) {
  if (p_ceex) return;
  p_ceex = new Ceex_Base(fl);
  p_ceex->SetBornMomenta(m_bornMomenta);
}


bool YFS_Handler::MakeYFS(){
  return MakeYFS(m_bornMomenta);
}

bool YFS_Handler::MakeYFS(ATOOLS::Vec4D_Vector &p)
{
  Reset();
  // p_dipoles->CreateAllDipoles(m_flavs, m_plab, m_bornMomenta);
  if (m_isrinital) {
    p_dipoles->MakeDipolesII(m_flavs, m_plab, m_bornMomenta);
  }
  m_ww_formfact = 1;
  m_fsrWeight = m_isrWeight = 1.0;
  CreatMomentumMap();
  if (m_mode == yfsmode::fsr) m_sp = m_s;
  m_v = 1. - m_sp / m_s;
  if ( m_v > m_vmax ) {
    m_yfsweight = 0.0;
    return false;
  }
  p_isr->SetV(m_v);
  if (m_v <= m_deltacut && m_mode!=yfsmode::fsr) { // correction weight included in Generate photon
    Reset();
    return false;
  }
  if (!CalculateISR()) return 0;
  m_FSRPhotons.clear();
  CalculateWWForm();
  CalculateCoulomb();
  p = m_plab;
  return true;
}



void YFS_Handler::MakeCEEX() {
  if (m_useceex) {
    Vec4D_Vector vv;
    p_ceex->SetBorn(m_born);
    for(size_t i = 0; i < m_plab.size(); ++i) vv.push_back(m_bornMomenta[i]);
    for(size_t i = 2; i < 4; ++i) vv.push_back(m_plab[i]);
    p_ceex->Init(vv);
    p_ceex->SetISRPhotons(m_ISRPhotons);
    p_ceex->SetBornMomenta(m_bornMomenta);
    p_ceex->SetISRFormFactor(m_formfactor);
    p_ceex->Calculate();
  }

}

void YFS_Handler::CalculateWWForm() {
  if (m_formWW) {
    MakeWWVecs(m_bornMomenta);
    m_ww_formfact = p_yfsFormFact->BVV_WW(m_plab, m_ISRPhotons, m_Wp, m_Wm, 1e-60, sqrt(m_sp) / 2.);
    if (m_ww_formfact < 0) PRINT_VAR(m_ww_formfact);
    if (IsBad(m_formfactor)) {
      THROW(fatal_error, "YFS Form Factor is NaN");
    }
  }
}

bool YFS_Handler::CalculateISR() {
  if (m_mode==yfsmode::fsr) return true;
  if (p_dipoles->GetDipoleII().size() != 2) {
    THROW(fatal_error, "Wrong dipole size for ISR");
  }
  // Address-of is deliberate: ISR keeps this pointer for the whole run (see
  // m_isrinital), which is why DipoleSet holds its dipoles by unique_ptr so
  // they never move.
  if (m_isrinital) p_isr->SetIncoming(&p_dipoles->GetDipoleII());
  m_isrinital = false;
  p_isr->NPhotons();
  p_isr->GeneratePhotonMomentum();
  p_isr->Weight();
  m_g=p_dipoles->GetDipoleII().m_gamma;
  m_gp=p_dipoles->GetDipoleII().m_gamma;
  p_dipoles->GetDipoleII().SetBorn(m_born);
  m_photonSumISR = p_isr->GetPhotonSum();
  m_ISRPhotons   = p_isr->GetPhotons();
  m_isrphotonsforME = m_ISRPhotons; 
  m_isrWeight = p_isr->GetWeight();
  m_photons.clear();
  for (const Vec4D &k : m_ISRPhotons)
    m_photons.push_back(YFS::Photon(k, &p_dipoles->GetDipoleII()));
  p_dipoles->GetDipoleII().AddPhotonsToDipole(m_ISRPhotons);
  p_dipoles->GetDipoleII().Boost();
  for(size_t i = 0; i < 2; ++i) {
    m_plab[i] = p_dipoles->GetDipoleII().GetNewMomenta(i); 
    ToLab(m_plab[i]);
  }
  double sp = (m_plab[0] + m_plab[1]).Abs2();
  if (!IsEqual(sp, m_sp, 1e-4) && !m_asymbeams) {
    msg_Error() << "Boost failed, sprime"
                << " is " << sp << " and should be "
                << m_sp << std::endl << "Diff = " <<
                m_sp - sp << std::endl << " Event with "
                << " N=" << p_dipoles->GetDipoleII().GetPhotons().size() << " photons" << std::endl
                << " V = " << m_v << std::endl
                << " Vmin = " << m_isrcut << std::endl
                << "ISR NPHotons = " << m_ISRPhotons.size() << std::endl;
  }
  return true;
}



void YFS_Handler::AddFormFactor() {
  if (m_CalForm) return;
  if (m_fullform >= 1) {
    if(m_tchannel!=0) m_formfactor = p_dipoles->TFormFactor();
    else {
      m_formfactor = p_dipoles->FormFactor();
    }
  }
  else if (m_fullform == 2) {
    m_formfactor = exp(m_g / 4.);//-m_alpha*M_PI);
  }
  else if (m_fullform == -1) {
    m_formfactor = 1;
  }
  else {
    if(FixedOrder()==fixed_order::nlo){
      m_formfactor = 1 + m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5);
    }
    else m_formfactor = exp(m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
  }
}

bool YFS_Handler::CalculateFSR(){
  return CalculateFSR(m_plab);
}

bool YFS_Handler::CalculateFSR(Vec4D_Vector & p) {
  // update NLO momenta from PHASIC
  // m_reallab should be used for 
  // fixed order corrections.
  // Final state eikonals should be constructed
  // for the final state momenta before emissions
  // of photons. 
  m_FSRPhotons.clear();
  m_fsrphotonsforME.clear();
  m_reallab = p;
  m_plab=p;
  // Pure-FSR mode never goes through MakeYFS, so CreatMomentumMap() (the only
  // place m_inparticles/m_outparticles get cleared) would otherwise never run
  // for this path, letting stale entries from earlier trials survive under
  // reused Particle* keys and leak into Signal_Processes::FillBlob via
  // GetOutParticles(). Reset it here on every call so it always starts from
  // the current born momenta.
  CreatMomentumMap();
  if(FixedOrder()==fixed_order::nlo && m_ISRPhotons.size()!=0) {
    for(size_t i = 2; i < m_plab.size(); ++i) m_outparticles[m_particles[i]] = m_plab[i];
    return true;
  }
  if(m_mode==yfsmode::isr) {
    // if(m_ISRPhotons.size() < m_mingammaN){
    //   m_isrWeight=0;
    //   return false;
    // }
    return true;
  }
  m_fsrWeight=1;
  p_dipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  // p_dipoles->CreateAllDipoles(m_flavs, m_plab, m_plab);
  CheckResonance();
  // p_dipoles->CreateAllDipoles(m_flavs, m_plab, m_plab);
  if(m_mode==yfsmode::isrfsr) {
    // Initial legs are the BORN beams, not the ISR-reduced ones. The
    // interference is between radiation off the incoming particles and off the
    // outgoing ones, so the initial leg of an initial-final pair is the
    // physical beam - which is also what KKMC's Yint uses (m_p1, m_p2 in
    // KKceex.cxx:315, the same momenta its Yisr = SForFac(alfpini, m_p1, m_p2)
    // uses). Final legs stay at m_plab, i.e. after the ISR recoil and before
    // FSR emission, matching KKMC's m_p3, m_p4.
    //
    // Passing m_plab for both also mixed frames once the beams were asymmetric:
    // CalculateISR() writes m_plab[0..1] back through ToLab() while
    // m_plab[2..] stay in the incoming-pair rest frame. m_bornMomenta is in
    // that rest frame throughout, so the pair is now built in one frame - which
    // matters because Btilda depends on the leg energies, not just invariants.
    Vec4D_Vector ifmom(m_plab);
    ifmom[0] = m_bornMomenta[0];
    ifmom[1] = m_bornMomenta[1];
    p_dipoles->MakeDipolesIF(m_flavs, ifmom, ifmom);
  }
  YFS::DipoleView ffdip(p_dipoles->GetDipoleFF());
  for (auto Dip = ffdip.begin(); Dip != ffdip.end(); ++Dip) {
    if(!Dip->IsResonance()) continue;
    p_fsr->Reset();
    Dip->BoostToQFM(0);
    Dip->SetBorn(m_born);
    p_fsr->SetV(m_v);
    if (!p_fsr->Initialize(*Dip)) {
      Reset();
      return false;
    }
    if (!p_fsr->MakeFSR()) {
      Reset();
      if (m_fsr_debug) p_debug->FillHist(m_plab, p_isr, p_fsr);
      return false;
    }
    m_photonSumFSR = p_fsr->GetPhotonSum();
    m_FSRPhotons   = p_fsr->GetPhotons();
    if (!p_fsr->F()) {
      m_fsrWeight = 0;
      if (m_fsr_debug) p_debug->FillHist(m_plab, p_isr, p_fsr);
      return false;
    } 

    // m_fsrphotonsforME = m_FSRPhotons;
    for(auto &k: m_FSRPhotons) m_fsrphotonsforME.push_back(k);
    Dip->AddPhotonsToDipole(m_FSRPhotons);
    Dip->SetMEPhotons(m_fsrphotonsforME);
    Dip->Boost();
    if(!p_fsr->YFS_FORM()) return false;
    p_fsr->HidePhotons();
    m_FSRPhotons = p_fsr->GetPhotons();
    Dip->AddPhotonsToDipole(m_FSRPhotons);
    p_fsr->Weight();
    m_fsrWeight *= p_fsr->GetWeight();
    m_plab[Dip->Left()]  =  Dip->GetNewMomenta(0);
    m_plab[Dip->Right()] =  Dip->GetNewMomenta(1);
    if(!IsEqual(m_flavs[Dip->Left()].Mass(), m_plab[Dip->Left()].Mass(),1e-5)){
      msg_Debugging()<<"Missmatch in Final state mass"<<std::endl
                 <<"Flavour = "<<m_flavs[Dip->Left()]<<std::endl
                 <<"Mass =   "<<m_flavs[Dip->Left()].Mass()<<std::endl
                 <<"Momentum =   "<<m_plab[Dip->Left()]<<std::endl
                 <<"Mass =   "<<m_plab[Dip->Left()].Mass()<<std::endl;
    }
    if(!IsEqual(m_flavs[Dip->Right()].Mass(), m_plab[Dip->Right()].Mass(),1e-5)){
      msg_Debugging()<<"Missmatch in Final state mass"<<std::endl
                 <<"Flavour = "<<m_flavs[Dip->Right()]<<std::endl
                 <<"Mass =   "<<m_flavs[Dip->Right()].Mass()<<std::endl
                 <<"Momentum =   "<<m_plab[Dip->Right()]<<std::endl
                 <<"Mass =   "<<m_plab[Dip->Right()].Mass()<<std::endl;
    }
  }
  for(size_t i = 2; i < m_plab.size(); ++i) {
    m_outparticles[m_particles[i]] = m_plab[i];
  }
  // get all photons
  m_FSRPhotons.clear();
  m_fsrphotonsforME.clear();
  // Rebuilt in full here rather than appended to, so a re-entered
  // CalculateFSR cannot leave last trial's photons behind.
  m_photons.clear();
  m_me_photons.clear();
  if (p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_ISRPhotons)
      m_photons.push_back(YFS::Photon(k, &p_dipoles->GetDipoleII()));
  YFS::DipoleView ffcollect(p_dipoles->GetDipoleFF());
  for (auto Dip = ffcollect.begin(); Dip != ffcollect.end(); ++Dip) {
    for(auto &k: Dip->GetPhotons()) {
      m_FSRPhotons.push_back(k);
      m_photons.push_back(YFS::Photon(k, &*Dip));
    }
    for(auto &k: Dip->GetMEPhotons()) {
      m_fsrphotonsforME.push_back(k);
      m_me_photons.push_back(YFS::Photon(k, &*Dip));
    }
  }
  // if(!CheckMomentumConservation()) return false;
  if(FixedOrder()==fixed_order::nlo){
    int totk = m_ISRPhotons.size();
    if(m_nlo_fsr_photons) totk += m_FSRPhotons.size();
    if(totk != 1) {
      if(totk > 1)
        msg_Error()<<"Wrong photon multiplicity at Fixed Order: "<<totk<<std::endl;
      return false;
    }
  }
  // if((m_ISRPhotons.size() +  m_FSRPhotons.size()) < m_mingammaN) {
  //   m_fsrWeight=0;
  //   return false;
  // }
  // CheckMasses();
  return true;
}


void YFS_Handler::MakeWWVecs(ATOOLS::Vec4D_Vector p) {
  m_Wm *= 0;
  m_Wp *= 0;
  Flavour_Vector wp, wm;
  for(size_t i = 2; i < p.size(); ++i)
  {
    if (m_flavs[i].IsAnti() && m_flavs[i].IntCharge()) {
      m_Wp += m_plab[i];
      wp.push_back(m_flavs[i]);
    }
    if (!m_flavs[i].IsAnti() && m_flavs[i].IntCharge()) {
      m_Wm += m_plab[i];
      wm.push_back(m_flavs[i]);
    }
    if (!m_flavs[i].IntCharge()) {
      if (m_flavs[i].IsAnti()) {
        m_Wm += m_plab[i];
        wm.push_back(m_flavs[i]);
      }
      else {
        m_Wp += m_plab[i];
        wp.push_back(m_flavs[i]);
      }
    }
  }
}


void YFS_Handler::CalculateCoulomb() {
  if (!m_coulomb) return;
  MakeWWVecs(m_bornMomenta);
  p_coulomb->Calculate(m_Wp, m_Wm);
  if (m_formWW) {
    // need to Subtract the Coulomb loop from virtual form factor
    // double s  = (m_Wp + m_Wm).Abs2();
    double am1 = m_Wp.Abs2();
    double am2 = m_Wm.Abs2();
    double beta = sqrt(1. - 2.*(am1 + am2) / m_s + sqr((am1 - am2) / m_s));
    if (m_betatWW >= beta) {
      p_coulomb->Subtract();
    }
    else m_coulSub = 0;
  }
}

void YFS_Handler::CalculateBeta() {
  // Invalidate last event's NLO pieces first, so an early return cannot leave
  // them looking current. Zeroed as well as flagged: a stale value that is
  // never read is still a trap for the next person to add a weight here.
  m_nlo_current = false;
  m_nlo_real = m_nlo_virtual = m_nlo_rv = m_nlo_rr = 0.;
  if(!m_rmode && !m_int_nlo) return;
  double realISR(0), realFSR(0);
  if (m_betaorder > 0) {
    if(m_real_only) {
      if(!m_no_born) m_real = p_dipoles->CalculateEEX()+1;
      else m_real = p_dipoles->CalculateEEX();
    }
    else if(m_virtual_only) {
      if(!m_no_born) m_real = p_dipoles->CalculateEEXVirtual();
      else m_real = p_dipoles->CalculateEEXVirtual()-1;
    }
    else {
      if(!m_no_born) m_real = p_dipoles->CalculateEEX()+p_dipoles->CalculateEEXVirtual();
      else m_real = p_dipoles->CalculateEEX()/m_born+p_dipoles->CalculateEEXVirtual()/m_born;
    }
    m_eex = m_real;
    if(IsNan(m_eex)) m_eex=0;
    // if(m_real < 0) m_real = 0;
    // m_real /= m_born;
  }
  if(m_nlotype!=nlo_type::born) {
    if(m_no_born) m_real=CalculateNLO()/m_born;
    else m_real=(m_born+CalculateNLO())/m_born;
    m_nlo_current = true;
  }
  if (m_useceex) MakeCEEX();
}

void YFS_Handler::InitNLO(){
  p_nlo->Init(m_flavs,m_reallab,m_bornMomenta);
  p_nlo->p_dipoles = p_dipoles;
  p_nlo->m_eikmom = m_plab;
  p_nlo->SetBorn(m_born);
  p_nlo->SetFSR(p_fsr);
  p_nlo->m_ISRPhotons = m_ISRPhotons;
  if (m_nlo_fsr_photons) p_nlo->m_FSRPhotons = m_fsrphotonsforME;
  else                   p_nlo->m_FSRPhotons.clear();
  // Mirror the two lines above: same photons, now carrying their dipole.
  p_nlo->m_photons.clear();
  if (p_dipoles->HasDipoleII())
    for (const Vec4D &k : m_ISRPhotons)
      p_nlo->m_photons.push_back(YFS::Photon(k, &p_dipoles->GetDipoleII()));
  if (m_nlo_fsr_photons)
    for (const YFS::Photon &k : m_me_photons) p_nlo->m_photons.push_back(k);
}

double YFS_Handler::CalculateNLO(){
// CheckMomentumConservation();
  InitNLO();
  // one-shot fixed-point dump for the KKMC CEEX comparison (YFS: CEEX_Compare)
  p_nlo->CEEXComparePoint();
  InitNLO();
  m_nlo_real = p_nlo->CalculateReal();
  // Hardest-photon-only contributions are captured as a side effect of the
  // nominal sums above (see NLO_Base::CalculateReal/CalculateRealVirtual/
  // CalculateRealReal) - no extra ME evaluation needed here.
  m_nlo_real_hardest = p_nlo->m_real_hard1;
  m_nlo_real_2hardest = p_nlo->m_real_hard2;
  InitNLO();
  m_nlo_virtual = p_nlo->CalculateVirtual();
  InitNLO();
  m_nlo_rv = p_nlo->CalculateRealVirtual();
  m_nlo_rv_hardest = p_nlo->m_rv_hard1;
  m_nlo_rv_2hardest = p_nlo->m_rv_hard2;
  InitNLO();
  m_nlo_rr = p_nlo->CalculateRealReal();
  m_nlo_rr_2hardest = p_nlo->m_rr_hard2;
  return m_nlo_real + m_nlo_virtual + m_nlo_rv + m_nlo_rr;
}


void YFS_Handler::GenerateWeight() {
  AddFormFactor();
  if (m_mode == yfsmode::isrfsr) m_yfsweight = m_isrWeight * m_fsrWeight;
  else if (m_mode == yfsmode::fsr) m_yfsweight = m_fsrWeight;
  else m_yfsweight = m_isrWeight;
  if (m_coulomb) m_yfsweight *= p_coulomb->GetWeight();
  if (m_formWW) m_yfsweight *= m_ww_formfact; //*exp(m_coulSub);
  CalculateBeta();

  double wif = 1.;
  if (m_ifireal && m_mode == yfsmode::isrfsr && m_nlotype == nlo_type::born &&
      p_nlo && p_nlo->HasReal()) {
    Vec4D_Vector allphotons(m_ISRPhotons);
    allphotons.insert(allphotons.end(), m_FSRPhotons.begin(), m_FSRPhotons.end());
    wif = p_dipoles->RealIFWeight(allphotons);
  }
  // The Born-level YFS weight: ISR x FSR crude (plus Coulomb/WW if on) times
  // the form factor, with NO NLO correction applied. This is what YFS.LO has
  // to reproduce -- built here directly rather than recovered downstream as
  // 1/m_real, so the LO column cannot inherit anything m_real does.
  const double w_lo = m_yfsweight * m_formfactor * (1.-m_v);
  m_yfsweight *= m_real + (wif - 1.);
  m_yfsweight *= m_formfactor*(1.-m_v);
  // Captured before the IsBad/negative-weight clamps below, since the named
  // weights are ratios against the weight the event actually carries.
  const double w_full = m_yfsweight;
  if(m_isr_debug) {
    Vec4D ele;
    for (int i = 2; i < m_flavs.size(); ++i)
    {
      if(IsEqual(m_flavs[i],kf_e)) {
        ele = m_plab[p_dipoles->m_flav_label[m_flavs[i]]];
        p_beams->BoostBackLab(ele);
        p_debug->FillHist("Form_Factor_FS_Angle", ele.Theta()*1000,m_formfactor,1);
      }
    }
  }
  DEBUG_FUNC("\nISR Weight = " << m_isrWeight << "\n" <<
             "  FSR Weight = " << m_fsrWeight << "\n" <<
             "  WW form Weight = " << m_ww_formfact << "\n" <<
             "  Total form Weight = " << m_formfactor << "\n" <<
             "  Coulomb Weight = " << p_coulomb->GetWeight() << "\n" <<
             " Coulomb Subtraction Weight = " << exp(m_coulSub) << "\n" <<
             "Total Weight = " << m_yfsweight << "\n");
  if(IsBad(m_yfsweight)){
    msg_Error()<<"\nISR Weight = " << m_isrWeight << "\n" <<
             "  FSR Weight = " << m_fsrWeight << "\n" <<
             "  Form Factor = " << m_formfactor << "\n" <<
             "  NLO  Correction = " << m_real << "\n" <<
             "Total Weight = " << m_yfsweight << "\n";
    m_yfsweight = 0;
  }
  if(m_yfsweight < 0 && m_skipNegWeights){
    msg_Debugging()<<"Skipping negative Weight in YFS"<<std::endl;
    m_yfsweight=0;
    m_negskip++;
  }

  // Build named NLO sub-weights. base_weight=1 so the nominal is unchanged.
  // YFSNLO  — Real + Virtual only (NLO denominator).
  // YFSNNLO — Real + Virtual + RealVirtual + RealReal (full NNLO denominator).
  m_nlo_weightsmap = Weights_Map{1.0};
  if (m_nlo_current && m_nlotype != nlo_type::born && !IsZero(m_real) &&
      (p_nlo->HasNLO() || p_nlo->HasNNLO())) {
    auto make_ratio = [this](double term, double denom) -> double {
      return m_no_born ? term / denom : (m_born + term) / denom;
    };
    auto ratio = [this](double term, double denom) -> double {
      return term / denom;
    };

    Weights wyfsnlo{1.0};

    auto emit = [this, &wyfsnlo](const std::string &name, double value) {
      if (!m_nlo_weight_breakdown &&
          name != "LO" && name != "NLO" && name != "NNLO")
        return;
      wyfsnlo[name] = value;
    };

    const bool have_fixed_order_ff =
        m_fullform >= 1 && m_tchannel == 0 && !IsZero(m_formfactor);
    const double ff_fixedorder_ratio =
        have_fixed_order_ff ? (1. + p_dipoles->FormFactorSum()) / m_formfactor
                             : 1.;

    // NLO: Real + Virtual
    if (p_nlo->HasNLO()) {
      const double nlo_sum   = (m_born + m_nlo_real + m_nlo_virtual)/m_born;
      const double eex_sum   = (m_born + m_eex)/m_born;
      const double real_sum  = (m_born + m_nlo_real)/m_born;
      const double virt_sum  = (m_born + m_nlo_virtual)/m_born;
      const double nlo_denom = m_no_born ? nlo_sum : (m_born + nlo_sum);
      if (!IsZero(nlo_denom)) {
        emit("Real", ratio((m_nlo_real)/m_born, m_real));
        emit("Virtual", ratio((m_nlo_virtual)/m_born, m_real));
        emit("NLO", ratio(nlo_sum, m_real));
        emit("BR", ratio(real_sum, m_real));
        emit("BV", ratio(virt_sum, m_real));
        // LO = (Born-level YFS weight) / (full weight), so that
        // nominal * YFS.LO == w_lo identically. Algebraically 1/m_real, but
        // built from the two weights themselves.
        // if (!IsZero(w_full)) emit("LO", w_lo/w_full);
        emit("EEX", ratio(m_eex, m_real));
        // Matching truncated to a fixed real-photon multiplicity, to see the
        // result "as if" only the 1 or 2 hardest photons were used in the
        // matching (full "NLO" above keeps all generated photons). Real is
        // summed over the 1 / 2 hardest photons; Virtual is always full.
        const double nlo_1g = (m_born + m_nlo_real_hardest  + m_nlo_virtual)/m_born;
        const double nlo_2g = (m_born + m_nlo_real_2hardest + m_nlo_virtual)/m_born;
        emit("NLO_1g", ratio(nlo_1g, m_real));
        emit("NLO_2g", ratio(nlo_2g, m_real));
        // Fixed-order comparison point: 1-photon NLO correction with the
        // resummed exp(form) form factor undone in favour of its 1+form
        // fixed-order truncation - matches a plain (non-YFS-resummed) NLO EW
        // calculation, which only ever has at most one real photon.
        if (have_fixed_order_ff)
          emit("NLO_FixedOrder", ratio(nlo_1g, m_real) * ff_fixedorder_ratio);
      }
    }

    // NNLO: RealVirtual + RealReal
    if (p_nlo->HasNNLO()) {
      const double nnlo_total = (m_born + m_nlo_real + m_nlo_virtual + m_nlo_rv + m_nlo_rr)/m_born;
      const double RR_total = (m_born + m_nlo_real + m_nlo_virtual + m_nlo_rr)/m_born;
      const double RV_total = (m_born + m_nlo_real + m_nlo_virtual + m_nlo_rv)/m_born;
      const double nnlo_denom = m_no_born ? nnlo_total : (m_born + nnlo_total);
      if (!IsZero(nnlo_denom)) {
        emit("RealVirtual", ratio((m_nlo_rv)/m_born, m_real));
        emit("RealReal", ratio((m_nlo_rr)/m_born, m_real));
        emit("NLO+RR", ratio(RR_total, m_real));
        emit("NLO+RV", ratio(RV_total, m_real));
        emit("NNLO", make_ratio(nnlo_total, nnlo_denom));
        // Matching truncated to a fixed real-photon multiplicity, the NNLO
        // analogue of NLO_1g/NLO_2g. 1 photon: Real + RealVirtual on the
        // single hardest, RealReal = 0 (a pair needs two photons). 2 photons:
        // Real + RealVirtual summed over the two hardest, plus the RealReal
        // pair formed by them. Virtual is always full.
        const double nnlo_1g =
            (m_born + m_nlo_real_hardest  + m_nlo_virtual + m_nlo_rv_hardest)/m_born;
        const double nnlo_2g =
            (m_born + m_nlo_real_2hardest + m_nlo_virtual + m_nlo_rv_2hardest +
             m_nlo_rr_2hardest)/m_born;
        emit("NNLO_1g", ratio(nnlo_1g, m_real));
        emit("NNLO_2g", ratio(nnlo_2g, m_real));
        // Fixed-order NNLO comparison: the 2-photon truncation (fixed-order
        // NNLO EW allows up to two real photons) with the resummed form factor
        // undone to its 1+form truncation.
        if (have_fixed_order_ff)
          emit("NNLO_FixedOrder", ratio(nnlo_2g, m_real) * ff_fixedorder_ratio);

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
        if (!m_vvtool) {
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
            emit("VV_EEX", ratio(vv, m_real));
            emit("NNLO_VV", ratio(nnlo_total + vv, m_real));
            emit("NNLO_VV_up", ratio(nnlo_total + (1.+d)*vv, m_real));
            emit("NNLO_VV_down", ratio(nnlo_total + (1.-d)*vv, m_real));
          } else {
            msg_Error() << METHOD << ": EEX double-virtual estimate is "
                        << vv << ", skipping the VV weights\n";
          }
        }
      }
    }

    // Emit forward/backward variants of each named contribution so they flow
    // through the normal multiweight machinery (Single_Process multiplies in
    // the universal ME x PDF x flux factors). Must run before the assignment.
    if (p_fb) p_fb->SplitWeights(wyfsnlo, m_plab, m_flavs);

    m_nlo_weightsmap["YFS"] = wyfsnlo;
  }
}

void YFS_Handler::YFSDebug(double W){
  p_debug->FillHist(m_plab, p_isr, p_fsr, W);
}


void YFS_Handler::Reset() {
  m_fsrWeight = 0;
  m_yfsweight = 0;
  m_ISRPhotons.clear();
  m_FSRPhotons.clear();
  m_photonSumISR *= 0;
  m_photonSumFSR *= 0;
  m_real = 1;
}

bool YFS_Handler::CheckMomentumConservation(){
  Vec4D incoming = m_bornMomenta[0]+m_bornMomenta[1];
  Vec4D outgoing;
  for(auto k: m_ISRPhotons)  outgoing+=k;
  for(auto kk: m_FSRPhotons) outgoing+=kk;
  for(size_t i = 2; i < m_plab.size(); ++i)
  {
    outgoing+=m_plab[i];
  }
  Vec4D diff = incoming - outgoing;
  if(!IsEqual(incoming,outgoing, 1e-5)){
    msg_Error()<<"Momentum not conserverd in YFS"<<std::endl
               <<"Incoming momentum = "<<incoming<<std::endl
               <<"Outgoing momentum = "<<outgoing<<std::endl
               <<"Difference = "<<diff<<std::endl
               <<"ISR Photons = "<<m_ISRPhotons<<std::endl
               <<"FSR Photons = "<<m_FSRPhotons<<std::endl;
  return false;
  }
  return true;
}


void YFS_Handler::CheckMasses(){
  bool allonshell=true;
  std::vector<double> mass;
  Vec4D_Vector p = m_plab;
  for(auto k: m_ISRPhotons) p.push_back(k);
  for(auto kk: m_FSRPhotons) p.push_back(kk);

  for(size_t i = 0; i < p.size(); ++i)
  {
    if(i<m_plab.size()){
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
    for(size_t i = 0; i < m_plab.size(); ++i)
    {
      msg_Debugging()<<"Mass after Mometum strechting"<<std::endl;
      if(i<m_plab.size()){
         msg_Debugging()<<"Flavour = "<<m_flavs[i]<<", with mass = "<<m_flavs[i].Mass()<<std::endl
                       <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
      }
      else{
         msg_Debugging()<<"Flavour = "<<Flavour(22)<<", with mass = "<<Flavour(22).Mass()<<std::endl
                        <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
      }
      m_plab[i] = p[i];
    }
  }
}

void YFS_Handler::SplitPhotons(ATOOLS::Blob * blob){
  if(IsEqual(m_photon_split,0)) return;
  p_splitter->SplitPhotons(blob);
}

Vec4D_Vector YFS_Handler::GetPhotons(){
  Vec4D_Vector k;
  for(auto p: m_ISRPhotons) k.push_back(p);
  for(auto p: m_FSRPhotons) k.push_back(p);
  return k;
}

void YFS_Handler::CheckResonance(){
  YFS::DipoleView ffres(p_dipoles->GetDipoleFF());
  for (auto D1 = ffres.begin(); D1 != ffres.end(); ++D1) {
    for (auto D2 = ffres.begin(); D2 != ffres.end(); ++D2) {
      if(D1==D2) continue;
      // Only two dipoles that both still radiate can double count a leg.
      // Without this, an overlapping non-resonant pair could switch off a
      // dipole that Define_Dipoles::SelectResonantDipoles picked.
      if(!D1->IsResonance() || !D2->IsResonance()) continue;
      if(D1->Right() == D2->Right() ||  D1->Right() == D2->Left()|| 
        D1->Left() == D2->Right()||  D1->Left() == D2->Left()){
        if(p_dipoles->ResonantDist(*D1) < p_dipoles->ResonantDist(*D2)) D2->SetResonance(false);
        else  D1->SetResonance(false);
        }
      }
    }
  }