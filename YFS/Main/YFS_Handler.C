#include "YFS/Main/YFS_Handler.H"

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/Main/ISR.H"
#include "YFS/NLO/EEX.H"
#include "YFS/NLO/Real.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "METOOLS/Loops/Master_Integrals.H"

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;

double evtcount = 0;
double failed = 0;
double maxv = 0;

YFS_Handler::YFS_Handler():
  p_realff(NULL), p_virtff(NULL)
{
  p_dipoles = new Define_Dipoles();
  p_coulomb = new Coulomb();
  p_realff = new Real_ff(m_betaorder);
  p_fsr = new FSR();
  p_debug = new Debug();
  m_born_set = false;
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  m_setparticles = false;
  p_isr = new YFS::ISR();
  p_nlo = new YFS::NLO_Base();
  m_formfactor = 1;
  m_isrinital = true;
  p_splitter = new PHOTONS::Photon_Splitter(m_photon_split);
  m_rmode = 0;
}

YFS_Handler::~YFS_Handler()
{
  if (p_isr) delete p_isr;
  if (p_fsr) delete p_fsr;
  if (p_realff) delete p_realff;
  if (p_coulomb) delete p_coulomb;
  if (p_debug)   delete p_debug;
  if (p_yfsFormFact) delete p_yfsFormFact;
  if (p_dipoles) delete p_dipoles;
  if (p_nlo) delete p_nlo;
  if (p_splitter) delete p_splitter;
  for (auto &p: m_particles){
    if(p) delete p;
  }
}


bool YFS_Handler::On()
{
  return m_mode;
}



void YFS_Handler::SetBeam(BEAM::Beam_Spectra_Handler *beam)
{
  p_beams = beam;
  // for(size_t i = 0; i < 2; ++i) m_beams.push_back(beam->GetBeam(i));
  m_beam1 = p_beams->GetBeam(0)->OutMomentum();
  m_beam2 = p_beams->GetBeam(1)->OutMomentum();
  if(m_beam1 != -m_beam2) m_asymbeams = true;
    else m_asymbeams = false;
}

void YFS_Handler::SetSprimeLimits(std::vector<double> &splimits) {
  if (m_flavs.size() == 0) return;
  double s = sqr(rpa->gen.Ecms());
  p_yfsFormFact->SetCharge(1);
  p_coulomb->SetAlphaQED(m_alpha);
  double maxV = 1. - m_smin / s;
  if (m_vmax > maxV && !m_asymbeams) {
    msg_Error() << "Warning: vmax to large in YFS integration reseting to " << maxV << std::endl;
    m_vmax = maxV;
  }
  splimits[0] = m_smin;
  splimits[1] = s;
  splimits[2] = s;
  if(sqrt(m_smin) > 91.2+10 ) m_resonance_mode = 0;
  else m_resonance_mode = 1;
  // m_resonance_mode = 1;
  p_nlo->SetResonanceMode(m_resonance_mode);
  // p_dipoles->m_resonance_mode=m_resonance_mode;
  for (int i = 0; i < splimits.size(); ++i) m_splimits[i] = splimits[i];
}

void YFS_Handler::SetFlavours(const ATOOLS::Flavour_Vector &flavs) {
  if(m_setparticles) return;
  m_flavs.clear();
  m_mass.clear();
  bool qed(false);
  for (int i = 0; i < flavs.size(); ++i) {
    m_flavs.push_back(flavs[i]);
    if (i < 2) {
      if (m_flavs[i].Mass() == 0 && m_fsrmode!=2) {
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
  if (!qed) m_fsrmode = 0;
  if (m_useceex) InitializeCEEX(m_flavs);
}

void YFS_Handler::SetBornMomenta(const ATOOLS::Vec4D_Vector &p) {
  m_bornMomenta.clear();
  m_momeikonal.clear();
  for (int i = 0; i < p.size(); ++i) {
    m_bornMomenta.push_back(p[i]);
    m_momeikonal.push_back(p[i]);
  }
  if (m_formWW) MakeWWVecs(m_bornMomenta);
  if(m_bornMomenta[0] != -m_bornMomenta[1]) m_asymbeams = true;
    else m_asymbeams = false;
  // AddFormFactor();
}

void YFS_Handler::SetMomenta(const ATOOLS::Vec4D_Vector &p) {
  m_plab.clear();
  for (int i = 0; i < p.size(); ++i) {
    m_plab.push_back(p[i]);
  }
}

void YFS_Handler::CreatMomentumMap() {
  m_inparticles.clear();
  m_outparticles.clear();
  for (int i = 0; i < 2; ++i)
  {
    m_inparticles[m_particles[i]] = m_bornMomenta[i];
    m_particles[i]->SetMomentum(m_bornMomenta[i]);
    if (m_particles[i]->Flav().IsAnti()) m_beamP = m_bornMomenta[i];
    else m_beamM = m_bornMomenta[i];

  }
  if(m_fsrmode!=0){
    for (int i = 2; i < m_flavs.size(); ++i)
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
  if (m_beamspread || m_isrinital ) {
    p_dipoles->MakeDipolesII(m_flavs, m_plab, m_bornMomenta);
  }
  m_ww_formfact = 1;
  m_fsrWeight = m_isrWeight = 1.0;
  CreatMomentumMap();
  if (m_fsrmode == 2) m_sp = m_s;
  m_v = 1. - m_sp / m_s;
  if ( m_v > m_vmax ) {
    m_yfsweight = 0.0;
    return false;
  }
  p_isr->SetV(m_v);
  if (m_v <= m_deltacut && m_fsrmode != 2  ) { // correction weight included in Generate photon
    m_yfsweight = 0.0;
    return false;
  }
  if (!CalculateISR()) return 0;
  m_FSRPhotons.clear();
  CalculateWWForm();
  CalculateCoulomb();
  p = m_plab;
  return true;
}



void YFS_Handler::MakeCEEX(const Vec4D_Vector & p) {
  if (m_useceex) {
    Vec4D_Vector vv;
    p_ceex->SetBorn(m_born);
    for (int i = 0; i < m_plab.size(); ++i) vv.push_back(m_bornMomenta[i]);
    for (int i = 2; i < 4; ++i) vv.push_back(m_plab[i]);
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
    Vec4D_Vector born;
    born.push_back(m_beamM);
    born.push_back(m_beamP);
    m_ww_formfact = p_yfsFormFact->BVV_WW(m_plab, m_ISRPhotons, m_Wp, m_Wm, 1e-60, sqrt(m_sp) / 2.);
    if (m_ww_formfact < 0) PRINT_VAR(m_ww_formfact);
    if (IsBad(m_formfactor)) {
      THROW(fatal_error, "YFS Form Factor is NaN");
    }
  }
}

bool YFS_Handler::CalculateISR() {
  if (m_fsrmode == 2) return true;
  if (p_dipoles->GetDipoleII()->size() != 2) {
    THROW(fatal_error, "Wrong dipole size for ISR");
  }
  if (m_isrinital) p_isr->SetIncoming(p_dipoles->GetDipoleII());
  m_isrinital = false;
  p_isr->NPhotons();
  m_N = p_isr->GetN();
  p_isr->GeneratePhotonMomentum();
  p_isr->Weight();
  m_g=p_dipoles->GetDipoleII()->m_gamma;
  m_gp=p_dipoles->GetDipoleII()->m_gamma;
  m_photonSumISR = p_isr->GetPhotonSum();
  m_ISRPhotons   = p_isr->GetPhotons();
  m_isrWeight = p_isr->GetWeight();
  // if(m_isrWeight==0) return 0;
  p_dipoles->GetDipoleII()->AddPhotonsToDipole(m_ISRPhotons);
  p_dipoles->GetDipoleII()->Boost();
  for (int i = 0; i < 2; ++i) m_plab[i] = p_dipoles->GetDipoleII()->GetNewMomenta(i);
  double sp = (m_plab[0] + m_plab[1]).Abs2();
  if (!IsEqual(sp, m_sp, 1e-4) && !m_asymbeams) {
    msg_Error() << "Boost failed, sprime"
                << " is " << sp << " and should be "
                << m_sp << std::endl << "Diff = " <<
                m_sp - sp << std::endl << " Event with "
                << " N=" << p_dipoles->GetDipoleII()->GetPhotons().size() << " photons" << std::endl
                << " V = " << m_v << std::endl
                << " Vmin = " << m_isrcut << std::endl
                << "ISR NPHotons = " << p_isr->m_N << std::endl;
  }
  m_isrinital = false;
  YFSDebug(1.0);
  if(m_photonSumISR.E() > sqrt(m_s)){
    return false;
  //   // PRINT_VAR(1-m_sp/m_s);
  //   // PRINT_VAR(p_isr->m_v);
  //   // PRINT_VAR(p_isr->m_v*sqrt(m_s)/2.);
  //   // msg_Error()<<"ISR photons are too energetic!" << std::endl;
  //   //            // <<"Photons "<< m_ISRPhotons<< std::endl
  //   //            // <<"Photon sum is "<< m_photonSumISR<< std::endl;
  //   // PRINT_VAR(m_ISRPhotons);
  //   // PRINT_VAR(m_photonSumISR);
  //   // PRINT_VAR(m_plab);
  }
  return true;

}



void YFS_Handler::AddFormFactor() {
  if (m_CalForm) return;
  // if(m_fsrmode==2) return; // Calculated in FSR.C for fsr dipoles
  if (m_fullform == 1) {
    // m_formfactor = p_yfsFormFact->BVV_full(m_bornMomenta[0], m_bornMomenta[1], m_photonMass, sqrt(m_s) / 2., 0);
    if(m_tchannel) m_formfactor = p_dipoles->TFormFactor();
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
    // high energy limit
    m_formfactor = exp(m_g / 4. + m_alpha / M_PI * (pow(M_PI, 2.) / 3. - 0.5));
    // m_CalForm = true;
  }
}

bool YFS_Handler::CalculateFSR(){
  // CheckMasses();
  // CheckMomentumConservation();
  return CalculateFSR(m_plab);
}

bool YFS_Handler::CalculateFSR(Vec4D_Vector & p) {
  // update NLO momenta from PHASIC
  // m_reallab should be used for 
  // fixed order corrections.
  // Final state eikonals should be constructed
  // for the final state momenta before emissions
  // of photons. 
  m_reallab = p;
  // m_reallab[0] = m_bornMomenta[0];
  // m_reallab[1] = m_bornMomenta[1];
  if (m_fsrmode == 0) return true;
  m_plab=p;
  
  p_dipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  if(m_mode==yfsmode::isrfsr)  p_dipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  for (int i = 2; i < m_momeikonal.size(); ++i)
  {
    m_momeikonal[i] = m_plab[i];
  }
  m_fsrWeight = 1;
  m_finalFSRPhotons.clear();
  m_newdipoles.clear();
  m_real = 1;
  Vec4D_Vector oldplab = m_plab;
  m_FSRPhotons.clear();
  if (m_fsrmode >= 1) {
    if (p_dipoles->GetDipoleFF()->size() == 0) {
      msg_Error() << "No dipoles found for YFS FSR " << std::endl;
      return true;
    }
    for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
         Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
      if(!Dip->IsResonance()) continue;
      p_fsr->Reset();
      Dip->BoostToQFM(0);
      Dip->SetBeams(0, m_plab[0]);
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
      if (!p_fsr->F(m_FSRPhotons)) {
        m_fsrWeight = 0;
        if (m_fsr_debug) p_debug->FillHist(m_plab, p_isr, p_fsr);
        return false;
      } 

      m_fsrphotonsforME = m_FSRPhotons;
      Dip->AddPhotonsToDipole(m_FSRPhotons);
      Dip->Boost();
      if(!p_fsr->YFS_FORM()) return false;
      p_fsr->HidePhotons();
      m_FSRPhotons   = p_fsr->GetPhotons();
      Dip->AddPhotonsToDipole(m_FSRPhotons);
      // Dip->Boost();
      // p_fsr->HidePhotons();
      // p_fsr->HidePhotons(m_FSRPhotons);
      // Dip->AddPhotonsToDipole(m_FSRPhotons);
      // m_fsrphotonsforME = Dip->GetPhotons();
      // if(m_photonSumFSR.E() > sqrt(m_s)/2) {
      //   // PRINT_VAR(m_photonSumFSR);
      //   m_fsrWeight = 0;
      //   return false;
      // }
      p_fsr->Weight();
      m_fsrWeight = p_fsr->GetWeight();
      int i(0);
      for (auto f : Dip->m_flavs) {
        m_plab[p_dipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
        i++;
      }
    }
    Vec4D_Vector oldplab = m_plab;
    for (int i = 2; i < m_plab.size(); ++i) {
      m_outparticles[m_particles[i]] = m_plab[i];
    }
  }
  // if (m_fsr_debug) p_debug->FillHist(m_plab, p_isr, p_fsr);
  // if (m_looptool) {
  //   CalculateVirtual(m_bornMomenta, m_born);
  // }
  // get all photons
  m_FSRPhotons.clear();
  for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
         Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
    for(auto &k: Dip->GetPhotons()) m_FSRPhotons.push_back(k);
      Dip->Clean();
  }
  // CheckMasses();
  CheckMomentumConservation();
  return true;
}


void YFS_Handler::MakeWWVecs(ATOOLS::Vec4D_Vector p) {
  m_Wm *= 0;
  m_Wp *= 0;
  Flavour_Vector wp, wm;
  for (int i = 2; i < p.size(); ++i)
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
  m_real=1;
  if(!m_rmode && !m_int_nlo) return;
  double realISR(0), realFSR(0);
  if (m_betaorder > 0) {
    p_realff->SetBorn(m_born);
    p_realff->SetMode(m_fsrmode);
    if (m_fsrmode != 2) {
      p_realff->SetMode(0);
      p_realff->SetIncoming(p_dipoles->GetDipoleII());
      p_realff->CalculateVirt();
      // p_realff->Calculate();
      // realISR = p_realff->GetReal();
      for(auto const &k: p_dipoles->GetDipoleII()->GetPhotons()){
          realISR += p_dipoles->GetDipoleII()->EEX(m_betaorder)*m_born/p_dipoles->GetDipoleII()->Eikonal(k)-m_born;
        }
    }
    if (m_fsrmode >= 1) {
      p_realff->SetMode(m_fsrmode);
      for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
           Dip != p_dipoles->GetDipoleFF()->end(); ++Dip)  {
        p_realff->SetIncoming(Dip, m_reallab, m_FSRPhotons, p_fsr->m_yini, p_fsr->m_zini);
        p_realff->CalculateVirt();
        // p_realff->Calculate();
        for(auto const &k: m_fsrphotonsforME){
          realFSR += Dip->EEX(m_betaorder)*m_born/Dip->Eikonal(k)-m_born;
        }
        // PRINT_VAR(realFSR);
        // PRINT_VAR(p_realff->GetReal());
      }
      if (m_use_fsr_beta == 0) realFSR = 0;
    }
    m_real = (realISR + realFSR + p_realff->AddVirtual())/m_born;
    if (m_virtual_only) {
      if(m_no_born) m_real = p_realff->AddVirtual()/m_born-1;
      else m_real = p_realff->AddVirtual()/m_born;
    }
    if (m_real_only) {
      if (m_looptool) m_real = realISR + realFSR + ((m_betaorder > 1 ? p_realff->AddVirtual(m_betaorder) - p_realff->AddVirtual(1) : 0));
      else {
        if(m_no_born) m_real = (realISR + realFSR)/m_born;
        else m_real = (m_born+realISR + realFSR)/m_born;
      }
    }
    if (m_useint) m_real += p_realff->IntIF(m_ISRPhotons, m_fsrphotonsforME, p_isr->m_yini, p_isr->m_zini, p_fsr->m_yini, p_fsr->m_zini);
  }
  // PRINT_VAR(m_nlotype);
  if(m_nlotype==nlo_type::loop || m_nlotype==nlo_type::real) {
    if(m_no_born) m_real=CalculateNLO()/m_born;
    else m_real=(m_born+CalculateNLO())/m_born;
  }
}



double YFS_Handler::CalculateNLO(){
  // CheckMomentumConservation();
  p_nlo->Init(m_flavs,m_reallab,m_bornMomenta);
  p_nlo->p_dipoles = p_dipoles;
  p_nlo->m_eikmom = m_plab;
  p_nlo->SetBorn(m_born);
  p_nlo->m_ISRPhotons = m_ISRPhotons;
  p_nlo->m_FSRPhotons = m_fsrphotonsforME;
  return p_nlo->CalculateNLO();
}


void YFS_Handler::GenerateWeight() {
  AddFormFactor();
  if (m_fixed_weight != 0) {
    m_yfsweight = m_fixed_weight;
    return;
  }
  if (m_semiyfs == 0) {
    if (m_constfsrW) m_fsrWeight = 1.;
    if (m_fsrmode == 1) m_yfsweight = m_isrWeight * m_fsrWeight;
    else if (m_fsrmode == 2) m_yfsweight = m_fsrWeight;
    else m_yfsweight = m_isrWeight;
  }
  else {
    m_yfsweight = SY->Weight(m_v, m_alpha, m_semiyfs);
  }
  if (m_coulomb) m_yfsweight *= p_coulomb->GetWeight();
  if (m_formWW) m_yfsweight *= m_ww_formfact; //*exp(m_coulSub);
  CalculateBeta();
  m_yfsweight*=m_real;
  m_yfsweight *= m_formfactor*(1-m_v);
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
}

void YFS_Handler::YFSDebug(double W){
  p_debug->FillHist(m_plab, p_isr, p_fsr, W);
}

int YFS_Handler::NHardPhotons(const Vec4D_Vector &k) {
  int N = 0;
  for (auto p : k) {
    if (p.E() > m_hardmin) N += 1;
  }
  return N;
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
  for (int i = 2; i < m_plab.size(); ++i)
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

  for (int i = 0; i < p.size(); ++i)
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
    for (int i = 0; i < m_plab.size(); ++i)
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

double YFS_Handler::Flux(const Vec4D& p1, const Vec4D& p2)
{
  return 0.25 / sqrt(sqr(p1 * p2) - p1.Abs2() * p2.Abs2());
}

double YFS_Handler::Flux()
{
  return 0.25 / sqrt(sqr(m_beam1 * m_beam2) - m_beam1.Abs2() * m_beam2.Abs2());
}



double YFS_Handler::Eikonal(Vec4D k) {
  double frac = (m_beam1 / (k * m_beam1) - m_beam2 / (k * m_beam2)).Abs2(); //*(m_beam1/(k*m_beam1) - m_beam2/(k*m_beam2));
  double s = -m_alpha / (4 * M_PI * M_PI) * frac ;
  return s;
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