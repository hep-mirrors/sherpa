#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <algorithm>
#include <limits>
#include <set>
#ifdef USING__LOOPTOOLS
  #include "clooptools.h"
#endif

using namespace YFS;
using namespace ATOOLS;
using namespace std;


Define_Dipoles::Define_Dipoles() {
  m_in = 2; // This is fine in YFS. It will not work for any other inital state multiplicity
  m_softphotonSum = {0., 0., 0., 0.};
  p_yfsFormFact = new YFS::YFS_Form_Factor();
}

Define_Dipoles::~Define_Dipoles() {
  if(p_yfsFormFact) delete p_yfsFormFact;
}


void Define_Dipoles::MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  m_N_born_Gamma=1;
  if(!HasISR()) return;
  if ((mom.size() < 2 || fl.size() < 2) ) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::initial<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for initial-initial dipole");
  }
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesII.clear();
  m_bornmomenta = born;
  Dipole_II(fl, mom);
  for(auto f: fl) if(f.IsPhoton()) m_N_born_Gamma+=1;  
}


void Define_Dipoles::MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born) {
  if(m_mode==yfsmode::fsr) return;
  if ((mom.size() != fl.size())) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::ifi<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for initial-final dipole");
  }
  if (!HasFSR() ) return;
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_out = fl.size() - m_in;
  m_dipolesIF.clear();
  Dipole_IF(fl, mom, born);
}

void Define_Dipoles::MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  if ((mom.size() != fl.size())) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::ifi<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for final-final dipole");
  }
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesFF.clear();
  m_bornmomenta = born;
  Dipole_FF(fl, mom);
}

void Define_Dipoles::MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born ) {
  if ((mom.size() != fl.size()) || (born.size() != fl.size())) {
    msg_Out()<<"Dipole type is  =  "<<dipoletype::final<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for final-final dipole");
  }
  ATOOLS::Flavour_Vector dipoleFlav;
  ATOOLS::Vec4D_Vector dipoleMom;
  Dipole_Vector dipoles;
  m_test_dip.clear();
  m_flav_label.clear();
  m_softphotonSum *= 0;
  m_bornmomenta = born;
  m_out = fl.size() - m_in;
  m_olddipoles.clear();
  m_dipolesFF.clear();
  int fsrc=0;
  for (int i = 2; i < fl.size(); ++i) if(fl[i].Charge()!=0) fsrc++;
  // m_dipolesIF.clear();
  for(size_t i = 0; i < fl.size(); ++i)
  {
    m_flav_label[fl[i]] = i;
  }
  if (!HasFSR() ) return;
  if (fl.size() == 4) {
    Flavour_Vector ff;
    Vec4D_Vector mm, bm;
    m_flav_label[fl[2]] = 2;
    m_flav_label[fl[3]] = 3;
    for(size_t i = 2; i < fl.size(); i++) {
      if(fl[i].IntCharge()!=0){
        ff.push_back(fl[i]);
        mm.push_back(mom[i]);
        bm.push_back(m_bornmomenta[i]);
      }
    }
    if(ff.size()==0) return;
    Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
    D.SetResonance(true);
    D.SetFlavLab(2,3);
    Dipole_FF(ff, mm);
    m_dipolesFF.push_back(D);
    return;
  }
  map<ATOOLS::Flavour, ATOOLS::Vec4D>::iterator itr;
  if (m_dip.size() != 0) {
    for (auto a : m_dip) {
      Get4Mom(fl, mom); // makes map for flavour momentum
      Flavour_Vector ff;
      Flavour f;
      Vec4D_Vector mm, bm;
      for(size_t i = 2; i < fl.size(); ++i)
      {
        f = fl[i];
        if (f.Charge()!=0) {
          ff.push_back(f);
          mm.push_back(m_test_dip[f]);
          bm.push_back(m_born_dip[f]);
          m_flav_label[f] = i;
          if (!IsEqual(f.Mass(), m_test_dip[f].Mass(), 1e-5)) {
            msg_Error() << "Incorrect mass mapping in dipole" << std::endl
                        << "Flavour mass is " << f.Mass() << std::endl
                        << "Four-Momentum mass is " << m_test_dip[f].Mass() << std::endl;
          }
          if (ff.size() == 2) break;
        }
      }
      Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
      Dipole_FF(ff, mm);
      if(fsrc==2)  {
        D.SetResonance(true);
        D.SetFlavLab(2,3);
      }
      else IsResonant(D);
      m_dipolesFF.push_back(D);
      msg_Debugging() << "Added " << ff << " to dipole " << a << std::endl;
    }
  }
  else {
    Get4Mom(fl, mom);
    int N = 0; // number of leptons minus the inital state
    for(int i=2; i < fl.size(); i++){
      if (fl[i].Charge()!=0) N += 1;
    }
    if (N == 2) {
      //only two leptons in final state
      // one dipole
      Flavour_Vector ff;
      Vec4D_Vector mm, bm;
      std::vector<int> id;
      for(size_t i = 2; i < fl.size(); i++) {
        if (fl[i].Charge()!=0) {
          ff.push_back(fl[i]);
          mm.push_back(mom[i]);
          bm.push_back(m_bornmomenta[i]);
          id.push_back(i);
        }
      }
      Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
      if(fsrc==2)  D.SetResonance(true);
      else IsResonant(D);
      Dipole_FF(ff, mm);
      D.SetFlavLab(id[0],id[1]);
      m_dipolesFF.push_back(D);
      return;
    }
    // N > 2 charged final-state particles: one dipole per unique pair, each
    // built exactly once. The labels are positions in fl/mom, so neutral
    // final-state particles cannot shift the mapping, and an odd N is treated
    // like any other (the previous perfect-matching enumeration returned no
    // pairings at all for odd N, and duplicated every pair three times for
    // N = 6).
    std::vector<int> charged_id;
    for(int i = 2; i < fl.size(); i++) {
      if (fl[i].Charge()!=0) charged_id.push_back(i);
    }
    for(size_t a = 0; a < charged_id.size(); a++) {
      for(size_t b = a+1; b < charged_id.size(); b++) {
        const int d1(charged_id[a]), d2(charged_id[b]);
        Flavour_Vector ff{fl[d1], fl[d2]};
        Vec4D_Vector mm{mom[d1], mom[d2]};
        Vec4D_Vector bm{m_bornmomenta[d1], m_bornmomenta[d2]};
        Dipole D(ff, mm, bm, dipoletype::final,m_alpha);
        D.SetFlavLab(d1,d2);
        m_dipolesFF.push_back(D);
        msg_Debugging() << "Added " << ff << " to dipole (" << d1 << "," << d2 << ")" << std::endl;
      }
    }
    // Every pair is needed for the virtual/form-factor sums, but a charged
    // particle may enter photon generation and the real eikonal current only
    // once. The radiating subset is flagged here.
    SelectResonantDipoles();
    // The charged/neutral out-particle lists describe the final state, not an
    // individual pair, so fill them once. Dipole_FF clears them on entry, so
    // the old per-dipole call left only the last pair behind.
    Flavour_Vector fsflav;
    Vec4D_Vector fsmom;
    for(int i = 2; i < fl.size(); i++) {
      fsflav.push_back(fl[i]);
      fsmom.push_back(mom[i]);
    }
    Dipole_FF(fsflav, fsmom);
  }
}

void Define_Dipoles::CreateAllDipoles(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta) {
    
    if (momenta.size() != flavors.size() || momenta.size() != born_momenta.size()) {
        THROW(fatal_error, "Inconsistent vector sizes in CreateAllDipoles");
    }
    
    ResetAllDipoleState();
    m_bornmomenta = born_momenta;
    m_out = flavors.size() - INITIAL_STATE_PARTICLES;
    
    for (size_t i = 0; i < flavors.size(); ++i) {
        m_flav_label[flavors[i]] = i;
    }
    
    auto initial_particles = ExtractInitialStateCharged(flavors, momenta, born_momenta);
    auto final_particles = ExtractFinalStateCharged(flavors, momenta, born_momenta);
    
    if (HasISR() && initial_particles.size() >= 2) {
        CreateInitialDipoles(initial_particles);
    }
    
    if (HasFSR() && final_particles.size() >= 2) {
        CreateFinalDipoles(final_particles);
    }
    
    if (HasISR() && HasFSR()) {
        CreateInitialFinalDipoles(initial_particles, final_particles);
    }
}

void Define_Dipoles::CreateInitialDipoles(
    const std::vector<ParticleInfo>& initial_particles) {
    
    CleanInParticles();
    m_dipolesII.clear();
    
    if (initial_particles.size() != 2) {
        msg_Error() << "Expected exactly 2 initial particles for II dipole, got "
                   << initial_particles.size() << std::endl;
        return;
    }
    
    Dipole dipole = CreateDipole(
        initial_particles[0],
        initial_particles[1],
        dipoletype::initial);
    
    m_g = dipole.m_gamma;
    m_gp = dipole.m_gammap;
    
    m_dipolesII.push_back(dipole);
    m_olddipoles.push_back(dipole);
}

void Define_Dipoles::CreateFinalDipoles(
    const std::vector<ParticleInfo>& final_particles) {
    
    CleanOutParticles();
    m_dipolesFF.clear();
    
    // Special handling for exactly 2 final-state particles
    bool is_two_body = (final_particles.size() == 2);
    
    // Create all unique pairs
    for (size_t i = 0; i < final_particles.size(); ++i) {
        for (size_t j = i + 1; j < final_particles.size(); ++j) {
            
            Dipole dipole = CreateDipole(
                final_particles[i],
                final_particles[j],
                dipoletype::final);
            
            // Set resonance
            if (is_two_body) {
                dipole.SetResonance(true);
            } else {
                IsResonant(dipole);
            }
            
            // Update bookkeeping
            ATOOLS::Flavour_Vector flavors = {
                final_particles[i].flavor,
                final_particles[j].flavor
            };
            ATOOLS::Vec4D_Vector momenta = {
                final_particles[i].momentum,
                final_particles[j].momentum
            };
            // Dipole_FF(flavors, momenta);
            
            m_dipolesFF.push_back(dipole);
        }
    }
}

void Define_Dipoles::CreateInitialFinalDipoles(
    const std::vector<ParticleInfo>& initial_particles,
    const std::vector<ParticleInfo>& final_particles) {
    
    CleanInParticles();
    m_dipolesIF.clear();
    
    // Create all initial-final combinations
    for (const auto& initial : initial_particles) {
        for (const auto& final : final_particles) {
            
            Dipole dipole = CreateDipole(
                initial,
                final,
                dipoletype::ifi);
            
            dipole.SetResonance(false);  // IF dipoles are never resonant
            m_dipolesIF.push_back(dipole);
        }
    }
}

void Define_Dipoles::ResetAllDipoleState() {
    m_test_dip.clear();
    m_flav_label.clear();
    m_softphotonSum *= 0;
    m_olddipoles.clear();
    m_dipolesII.clear();
    m_dipolesFF.clear();
    m_dipolesIF.clear();
}

void Define_Dipoles::Get4Mom(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector mom) {
  Vec4D_Vector P;
  for(size_t i = 2; i < fl.size(); ++i)
  {
    if (fl[i].IntCharge()!=0) {
      m_test_dip[fl[i]] = mom[i];
      P.push_back(mom[i]);
      if (P.size() == 2) break;
    }
  }
  if (P.size() != 2) {
    PRINT_VAR(P.size());
    THROW(fatal_error, "Wrong size dipole");
  }
}


void Define_Dipoles::Dipole_II(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  CleanInParticles();
  Flavour_Vector ff;
  Vec4D_Vector mm, bm;
  for(size_t i = 0; i < 2; i++) {
    ff.push_back(fl[i]);
    mm.push_back(mom[i]);
    bm.push_back(m_bornmomenta[i]);
  }
  Dipole D(ff, mm, bm, dipoletype::initial,m_alpha);
  m_g=D.m_gamma;
  m_gp=D.m_gammap;
  D.SetFlavLab(0, 1);
  m_olddipoles.push_back(D);
  m_dipolesII.push_back(D);
}


void Define_Dipoles::Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  CleanOutParticles();
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_FF Momenta");
  }
  for(size_t i = 0; i < fl.size(); ++i)
  {
    if (fl[i].IsQED()) {
      m_chargedoutparticles.push_back(mom[i]);
      m_massOutC.push_back(mom[i].Mass());

    }
    else {
      m_neutraloutparticles.push_back(mom[i]);
      m_massOutN.push_back(mom[i].Mass());
    }
  }
}



void Define_Dipoles::Dipole_IF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  CleanInParticles();
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_IF Momenta");
  }
  Flavour_Vector ff;
  Vec4D_Vector mm, bm;
  //create IF dipoles
    for(size_t i = 0; i < 2; ++i)
    {
      for(size_t j = 2; j < fl.size(); ++j)
      {
        if(fl[i].IntCharge()==0) continue;
        if(fl[j].IntCharge()==0) continue;
        ff.clear();
        mm.clear();
        bm.clear();
        ff.push_back(fl[i]);
        ff.push_back(fl[j]);

        mm.push_back(mom[i]);
        mm.push_back(mom[j]);


        bm.push_back(born[i]);
        bm.push_back(born[j]);
        Dipole D(ff, mm, bm, dipoletype::ifi,m_alpha);
        D.SetResonance(false);
        m_dipolesIF.push_back(D);
      }
    }
}



double Define_Dipoles::CalculateRealSub(const Vec4D &k) {
  double sub(0);
  // if(FixedOrder()!=fixed_order::full) return sub;
  Vec4D eik{0.,0.,0.,0.};
  for (auto &D : m_dipolesII) {
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
       Vec4D p = D.GetMomenta(i);
      eik += D.m_Q[i]*p/(p*k);
    }
  }
  for (auto &D : m_dipolesFF) {
     if(!D.IsResonance()) continue;
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
      Vec4D p = D.GetMomenta(i);
      eik += -D.m_Q[i]*p/(p*k);
    }
  }
  sub = -m_alpha / (4 * M_PI * M_PI)*eik*eik;
  return sub/(m_N_born_Gamma!=0?m_N_born_Gamma:1.0);
}

double Define_Dipoles::CalculateRealSubIF(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_dipolesIF){
    if(m_massless_sub) sub += D.EikonalMassless(k, D.GetMomenta(0), D.GetMomenta(1));
    else sub +=  D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  return sub;
}


double Define_Dipoles::CalculateVirtualSub() {
  double sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  if(m_dim_reg==1) return CalculateVirtualSubEps();
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  for (auto &D : m_dipolesFF) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  return sub;
}

double Define_Dipoles::CalculateVirtualSubEps() {
  DivArrD sub(0);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  for (auto &D : m_dipolesFF) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }

  for (auto &D : m_dipolesIF){
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateVVSubEps() {
  DivArrD sub(0);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  for (auto &D : m_dipolesFF) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }

  for (auto &D : m_dipolesIF){
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  m_vvSub=0.5*sub*sub;
  return (0.5*sub*sub).Finite();
}

double Define_Dipoles::CalculateRealVirtualSubEps(const Vec4D &k) {
  // Virtual YFS B-hat in dim-reg at the reduced (photon-removed) kinematics.
  // Must use the same BVV function as CalculateVirtualSubEps: the RV finite
  // remainder is defined by subtracting 2*alpha*B_fin off the one-loop real
  // emission ME, with the identical finite convention as the Born virtual —
  // otherwise beta_1^1 does not vanish in the soft limit.
  DivArrD sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_dipolesFF) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }
  m_virtSub=sub;
  return sub.Finite();
}



double Define_Dipoles::FormFactorSum(){
  double form = 0;

  for(auto &D: m_dipolesII){
    form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D, sqrt(m_s)/2);
  }
  // if(!m_hidephotons){
      for(auto &D: m_dipolesFF){
        form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D, sqrt(D.Sprime())/2); 
      }
    // }
  if(m_ifisub==1){
    // IFForFac = Btilda + t-channel virtual, i.e. KKMC's TForFac, and
    // ChargeNorm() = -QiQj*thetaij reproduces KKMC's +/- pattern across its
    // four TForFac calls (KKceex.cxx:315). +ChargeNorm is also what
    // TFormFactor() and every CalculateVirtualSub*() use on these same
    // dipoles, so at NLO the exponentiated IF form factor and the IF virtual
    // subtracted from the one-loop ME now carry the same sign.
    //
    // omega is sqrt(s)/2 -- one cutoff shared by all four pairs, as in KKMC,
    // and the same maximal choice the II term above makes. It was
    // sqrt(D.Sprime())/2, which for an initial-final pair is
    // (s/2)(1 -+ beta cos theta) and so is ANGLE-DEPENDENT: a soft cutoff odd
    // in cos(theta) manufactures A_FB by construction. That, not the
    // interference, was supplying essentially all of Sherpa's asymmetry.
    // Measured standalone at 0.7 GeV, |cos| < 0.55 (YFS/Tools/IFI_Budget.C):
    //
    //   assembly                              A_FB(mu+)
    //   -ChargeNorm, omega = sqrt(Sprime)/2    -0.01614   <- was live here
    //   +ChargeNorm, omega = sqrt(Sprime)/2    +0.01614
    //   +ChargeNorm, omega = sqrt(s)/2         -0.00056   <- this line
    //   KKMC CEEX2 / BabaYaga measured         -0.0103 / -0.0110
    //   Sherpa measured, before this change    -0.01605
    //
    // The -0.01614 reproduces the measured -0.01605, i.e. the whole of the old
    // asymmetry came from this one line, and it was a cutoff artifact whose
    // sign had been flipped to make it land near KKMC.
    //
    // A_FB is exactly linear in log(omega) and crosses KKMC's -0.0103 at
    // omega ~ 56 MeV, which is not a scale in the problem: no common cutoff is
    // a prediction. The remaining gap to KKMC is the real interference, which
    // KKMC gets from summing over photon partitions in the CEEX amplitude and
    // which Sherpa can only get on the emission side -- see RealIFWeight() and
    // the IFI_Real switch.
    for(auto &D: m_dipolesIF){
      form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega());
    }
  }
  return form;
}

double Define_Dipoles::IFIOmega() const {
  // Soft cutoff for the IF form factor. It must be ONE scale for all four
  // pairs, and it must be the boundary above which the interference is carried
  // by explicit photons instead.
  //
  // With IFI_Real off nothing carries it, so the only consistent choice is the
  // kinematic maximum: the exponent then holds the whole soft integral and the
  // result is cutoff-independent by construction, matching what the II term
  // does with sqrt(s)/2.
  //
  // With IFI_Real on, RealIFWeight() supplies the region above the generation
  // cutoff, so the exponent must stop there or the two double-count.
  if (m_ifireal && m_ifiomega > 0.) return m_ifiomega;
  return sqrt(m_s)/2.;
}


double Define_Dipoles::RealIFWeight(const ATOOLS::Vec4D_Vector &photons) {
  // The interference the IF form factor no longer carries once IFIOmega() has
  // been lowered to the generation cutoff. Photons are generated from the
  // dipole-diagonal density S_II + S_FF; the YFS radiation function is
  // S_II + S_FF + S_IF, so each generated photon is reweighted by
  //
  //     1 + S_IF(k) / (S_II(k) + S_FF(k))
  //
  // whose average over the crude exponentiates to the soft integral of S_IF
  // above the cutoff -- exactly what came out of the exponent. Verified
  // numerically in YFS/Tools/IFI_Budget.C: integrating S_IF over photon phase
  // space between two cutoffs reproduces the Btilda difference to 4+ digits.
  //
  // CalculateRealSubIF() is S_IF and CalculateRealSubEEX() is S_II + S_FF;
  // both go through Dipole::Eikonal(k,p1,p2), so they share one convention,
  // and NLO_Base already uses CalculateRealSubEEX() as the crude eikonal
  // (NLO_Base.C:524), so this is the same "crude" the generation defines.
  //
  // Numerator and denominator now share the initial legs: both the II dipole
  // and the IF dipoles are built from the Born beams (YFS_Handler.C, in
  // CalculateFSR), so the ratio is not spoiled by the ISR recoil.
  //
  // The ISR and FSR generation cutoffs differ, but FSR::YFS_FORM()'s Piatek
  // term m_DelYFS puts the FSR bookkeeping back onto m_Emin, which is what
  // IFIOmega() defaults to - the same single-Emin arrangement KKMC uses for
  // Yisr, Yfsr and Yint. Still scan IFI_Omega and check the answer is flat
  // before believing a number: that is the observable statement that the
  // exponent and these photons are cancelling each other.
  if (!m_ifireal || m_dipolesIF.empty()) return 1.;
  const double omega = IFIOmega();
  double w(1.);
  for (const auto &k : photons) {
    if (IsZero(k.E())) continue;
    // Only photons ABOVE the cutoff. Everything below omega is already held by
    // the IF form factor, so reweighting it here counts that region twice.
    //
    // ISR photons cannot trip this - their generation cutoff is
    // (sqrt(s)/2)*m_isrcut = IR_CUTOFF/2, which is exactly IFIOmega()'s default
    // - but FSR photons can: FSR_CUT defaults to 1e-2*IR_CUTOFF, so the FSR
    // generation reaches about two decades lower. Integrating S_IF over that
    // band is worth roughly dY_IF/dlog(omega) * log(100), a percent-level shift
    // on the rate, and it is one-signed, so it does not average away.
    if (k.E() <= omega) continue;
    const double crude = CalculateRealSubEEX(k);
    if (IsZero(crude) || IsBad(crude)) continue;
    double r = 1. + CalculateRealSubIF(k)/crude;
    // S_IF is not sign-definite, so a photon in a region where the
    // interference dominates the diagonal can drive the ratio negative or very
    // large. Both are outside the soft approximation this weight is built on -
    // |S_IF| << S_II + S_FF is what makes the eikonal reweighting valid, and a
    // photon violating it is one whose interference belongs to the matrix
    // element, not here.
    //
    // CLAMP to the boundary, do not drop the factor. Dropping means using 1,
    // which for r below the floor throws away a suppression and for r above the
    // ceiling throws away an enhancement - and since the low side is the more
    // common one, that is a net upward bias on the rate, not a neutral guard.
    // Anything clamped is counted so the rate is reportable rather than silent.
    if (r < m_ifi_rclip)      { r = m_ifi_rclip;    ++m_ifi_clipped; }
    else if (r > 1./m_ifi_rclip) { r = 1./m_ifi_rclip; ++m_ifi_clipped; }
    w *= r;
  }
  if (IsBad(w)) return 1.;
  return w;
}

double Define_Dipoles::FormFactor(){
  double form = FormFactorSum();
  if(FixedOrder()==fixed_order::nlo){
    return 1.+form;
  }
  return exp(form);
}


double Define_Dipoles::TFormFactor(){
  double form = 0;
  for(auto &D: m_dipolesII){
    form+= D.ChargeNorm()*p_yfsFormFact->R1(D);
  }
    for(auto &D: m_dipolesFF){
      form += D.ChargeNorm()*p_yfsFormFact->R1(D);
  }
  if(m_ifisub==1){
    // IF dipoles use IFForFac here too, NOT R1. An initial-final pair is
    // t-channel-like by construction, so its form factor does not depend on the
    // TChannel setting -- that flag is about how the II/FF (s-channel) dipoles
    // are treated. Routing IF through R1 made the interference term change
    // whenever TChannel changed, which is why the IF treatment appeared entangled
    // with a switch that has nothing to do with it. Both FormFactorSum() and
    // TFormFactor() now give the IF dipoles the same object.
    for(auto &D: m_dipolesIF){
      form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega());
    }
  }
  if(FixedOrder()==fixed_order::nlo){
    return 1.+form;
  }
  return exp(form);
}

double Define_Dipoles::CalculateVirtualSubTchannel(){
   // YFSij = 2.d0*B0ij - B0ii - B0jj
   //   .         + 4.d0 * mi2 * C0singular(mi2,phmass)
   //   .         + 4.d0 * mj2 * C0singular(mj2,phmass)
   //   .         + 8.d0*pi_pj * C0ij
  if(m_dim_reg) return CalculateVirtualSubTchannelEps();
  double sub(0);
  // Vec4D_Vector pvirt;
  // std::vector<double> z,th;
  // pvirt.push_back(m_dipolesII[0].GetNewMomenta(0));
  // pvirt.push_back(m_dipolesII[0].GetNewMomenta(1));
  // pvirt.push_back(m_dipolesFF[0].GetBornMomenta(0));
  // pvirt.push_back(m_dipolesFF[0].GetBornMomenta(1));
  // th.push_back(1);
  // th.push_back(1);
  // th.push_back(-1);
  // th.push_back(-1);
  // z.push_back(m_dipolesII[0].m_Qi);
  // z.push_back(m_dipolesII[0].m_Qj);
  // z.push_back(m_dipolesFF[0].m_Qi);
  // z.push_back(m_dipolesFF[0].m_Qj);
  // // double m1 = pvirt[0].Mass();
  // // double m2 = pvirt[1].Mass();
  // // double m3 = pvirt[2].Mass();
  // // double m4 = pvirt[3].Mass();
  // for(size_t i = 0; i < pvirt.size(); ++i)
  // {
  //   for(size_t j=i; j<pvirt.size(); ++j ){
  //     double etaij = z[i]*z[j]*th[i]*th[j];
  //     Complex YFSij = 0.;
  //     double mi = pvirt[i].Mass();
  //     double mj = pvirt[j].Mass();
  //     double s = (pvirt[i]-pvirt[j]).Abs2();

  //     Complex bii = B0(0,mi*mi,mi*mi);
  //     Complex bjj = B0(0,mj*mj,mj*mj);
  //     Complex bij = B0(0,mi*mi,mj*mj);
  //     Complex cij = C0(mi*mi,(th[i]*pvirt[i]+th[j]*pvirt[j]).Abs2(),mj*mj,
  //                      0.,mi*mi,mj*mj);
  //     Complex cii = C0(mi*mi, 0., mi*mi, 0.0, mi*mi, mi*mi);
  //     // YFSij = 8*pvirt[i]*pvirt[j]*cij;
  //     if(i==j){
  //       YFSij = bii-4.*mi*mi*0.5*log(m_photonMass*m_photonMass/mi/mi);
  //     }
  //     else{
  //       YFSij = 2.*pvirt[i]*pvirt[j]*cij+0.5*bij;
  //       }
  //     if(IsBad(YFSij)){
  //       msg_Error()<<"YFS Virtual Sub is NaN"<<endl
  //                  <<"bii = "<<bii<<endl
  //                  <<"bij = "<<bij<<endl
  //                  <<"cii = "<<cii<<endl
  //                  <<"cij = "<<cij<<endl;
  //     }
  //     sub+=m_alpi*etaij*YFSij.real();
  // //     // PRINT_VAR(etaij*YFSij);
  //   }
  // }
  // PRINT_VAR(count);
  for (auto &D : m_dipolesII){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_dipolesFF){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_dipolesIF){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  // clearcache();
  return sub;
}

double Define_Dipoles::CalculateVirtualSubTchannelEps() {
  DivArrD sub(0);
  DivArrD massph(0,-1.,0,0,0,0);
  for (auto &D : m_dipolesII) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
    #endif
  }
  for (auto &D : m_dipolesFF) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
    #endif
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateRealVirtualSub(const Vec4D & k) {
  double sub(0);
  for (auto &D : m_dipolesII) {
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetNewMomenta(0), D.GetNewMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_dipolesFF) {
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetOldMomenta(0), D.GetOldMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);

  }

  for (auto &D : m_dipolesIF){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.m_QiQj*p_yfsFormFact->BVV_full(D.GetBornMomenta(0), D.GetBornMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  return sub;
}


double Define_Dipoles::CalculateEEX(){
  double eex=0;
  for (auto &D: m_dipolesII){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_dipolesFF){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_dipolesIF){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  return eex;
}

double Define_Dipoles::CalculateEEXVirtual(){
  return CalculateEEXVirtual(m_betaorder);
}

// Same, at an EXPLICIT order rather than the runcard's BETA. Differencing two
// orders isolates one term of the EEX virtual series
// (Dipole::VirtualEEX: 0.5*gamma at order 1, +0.125*gamma^2 at order 2), which
// is how the approximate double-virtual is built in
// YFS_Handler::GenerateWeight when no exact VV provider exists.
double Define_Dipoles::CalculateEEXVirtual(int betaorder){
  double vint{1.}, vfin{1};
  for (auto &D: m_dipolesII){
    vint*=1+D.VirtualEEX(betaorder);
  }
  for (auto &D: m_dipolesFF){
    vfin*=1+D.VirtualEEX(betaorder);
  }
  return vint*vfin;
}

double Define_Dipoles::EEXRealVirtual(const Vec4D &k){
  double eex = 0;
  for(auto &D: m_dipolesII){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  for(auto &D: m_dipolesFF){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  // for(auto &D: m_dipolesIF){
  //   D.m_betaorder = 2;
  //   eex += D.Beta1(k)/D.Eikonal(k);
  // }
  if(IsNan(eex)) return 0;
  return eex;
}

double Define_Dipoles::CalculateRealSubEEX(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_dipolesII) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  for (auto &D : m_dipolesFF) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  // for (auto &D : m_dipolesIF) {
  //   sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  // }

  return sub;
}


void Define_Dipoles::CleanInParticles() {
  m_chargedinparticles.clear();
  m_neutralinparticles.clear();
  m_massInC.clear();
  m_massInN.clear();
}

void Define_Dipoles::CleanOutParticles() {
  m_chargedoutparticles.clear();
  m_neutraloutparticles.clear();
  m_massOutC.clear();
  m_massOutN.clear();
}

void Define_Dipoles::CleanUp() {
  m_dipoles.clear();
}

double Define_Dipoles::CalculateFlux(const Vec4D &k){
  if(!HasISR()) return 1;
  double sq, sx;
  double flux = 1;
  dipoletype::code fluxtype;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  if(HasISR()&&HasFSR()){
    fluxtype = dipoletype::initial;
  }
  else if(!HasFSR()){
    fluxtype = dipoletype::initial;
  }
  else if(!HasISR()){
    fluxtype = dipoletype::final;
  }
  else{
    msg_Error()<<"Unknown dipole type in "<<METHOD<<std::endl;
  }
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetBornMomenta(0)+D.GetBornMomenta(1);
      sq = (QX).Abs2(); 
      sx = (QX-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_dipolesFF) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_dipolesFF.size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2(); 
      sx = (Q-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_dipolesFF) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_dipolesFF.size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, const Vec4D &kk){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  dipoletype::code fluxtype1, fluxtype2;
  if(m_noflux==1) return 1;
  fluxtype1 = WhichResonant(k);
  fluxtype2 = WhichResonant(kk);
  if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k-kk).Abs2();
      flux = sx/sq;
      return flux;
    }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::final){
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k+kk).Abs2();
      flux = sq/sx;
    }
    return flux;
  }
  else if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::final){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+kk).Abs2();
      flux*= sq/sx;
   }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::initial){
    for (auto &D : m_dipolesII) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-kk).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_dipolesFF) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux*= sq/sx;
    }
  }
  else{
    msg_Error()<<"Unknown flux type in "<<METHOD<<" with \nfluxtype1 = "<<fluxtype1<<"\n and fluxtype2 = "<< fluxtype2 <<std::endl;
  }
  return flux;
}


double Define_Dipoles::Propagator(const double &s, int width){
  Flavour fl;
  Complex Prop = Complex(0.,0.);///Complex(s,0.0);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      fl = v->in[0];
      if(IsZero(fl.Mass())) continue;
      Prop += Complex(1.,0.)/Complex(s-sqr(fl.Mass()),fl.Width()*fl.Mass());
    }
  }
  return ((Prop*conj(Prop)).real());
}

double Define_Dipoles::ResonanceWidthDistance(YFS::Dipole &D) {
  // |m_ij - M| / Gamma, minimised over the resonances of the process. Returns a
  // large sentinel if the process has no usable resonance, so that a pair stays
  // selectable and the ordering simply falls back to the order of the pairs.
  const double mass_d((D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass());
  double mdist(std::numeric_limits<double>::max());
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if (IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width())) continue;
      mdist = std::min(mdist, std::abs(mass_d - v->in[0].Mass()) / v->in[0].Width());
    }
  }
  return mdist;
}

void Define_Dipoles::SelectResonantDipoles() {
  // m_dipolesFF holds every unique pair of charged final-state particles, which
  // is what the virtual and form-factor sums want. Photon generation
  // (YFS_Handler::CalculateFSR) and the real eikonal current
  // (CalculateRealSub) must instead see each charged particle exactly once,
  // otherwise a leg radiates -- and is boosted -- more than once. So the
  // radiating subset has to be a matching of the charged final state, and that
  // is what IsResonance() flags.
  //
  // Selection is greedy over two passes. Same-flavour opposite-charge pairs are
  // the physical dipoles and are matched first, most resonant first, so that
  // with several candidates of one flavour the pairing closest to a resonance
  // of the process wins. A second pass pairs whatever is left by opposite
  // charge alone, so no charged leg is dropped when the flavours do not pair up.
  for (auto &D : m_dipolesFF) D.SetResonance(false);
  std::set<int> used;
  for (int pass(0); pass < 2; ++pass) {
    std::vector<std::pair<double, size_t> > cand;
    for (size_t i(0); i < m_dipolesFF.size(); ++i) {
      Dipole &D(m_dipolesFF[i]);
      if (D.m_QiQj >= 0) continue;                      // opposite charges only
      if (pass == 0 && !D.IsDecayAllowed()) continue;    // same flavour first
      if (used.count(D.Left()) || used.count(D.Right())) continue;
      cand.push_back(std::make_pair(ResonanceWidthDistance(D), i));
    }
    std::stable_sort(cand.begin(), cand.end());
    for (size_t c(0); c < cand.size(); ++c) {
      Dipole &D(m_dipolesFF[cand[c].second]);
      if (used.count(D.Left()) || used.count(D.Right())) continue;
      D.SetResonance(true);
      used.insert(D.Left());
      used.insert(D.Right());
      msg_Debugging() << METHOD << "(): radiating dipole (" << D.Left() << ","
                      << D.Right() << ") " << D.m_flavs[0] << D.m_flavs[1]
                      << " pass=" << pass << " |m-M|/Gamma=" << cand[c].first
                      << std::endl;
    }
  }
  std::set<int> legs;
  for (auto &D : m_dipolesFF) {
    legs.insert(D.Left());
    legs.insert(D.Right());
  }
  for (auto l : legs) {
    if (used.count(l)) continue;
    static bool warned(false);
    if (!warned) {
      msg_Error() << METHOD << "(): charged final-state particle at position "
                  << l << " enters no radiating dipole, so the final-state "
                  << "eikonal current does not conserve charge. This is "
                  << "expected for an odd number of charged final-state "
                  << "particles. Further warnings suppressed." << std::endl;
      warned = true;
    }
    msg_Debugging() << METHOD << "(): unpaired charged leg at " << l << std::endl;
  }
}

void Define_Dipoles::IsResonant(YFS::Dipole &D) {
double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      // if(D.IsDecayAllowed()){
      //   D.SetResonance(true);
      //   continue;
      // }
      if(D.m_QiQj==1 || !D.IsDecayAllowed()){
        D.SetResonance(false);
        continue;
        }   
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<m_resonace_max) {
        D.SetResonance(true);
        return;
      }
      else D.SetResonance(false);
    }
    D.SetResonance(false);
  }
  if(D.IsDecayAllowed())  D.SetResonance(true);
  else  D.SetResonance(false);
}

bool Define_Dipoles::CheckResonant(YFS::Dipole &D) {
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<5) {
        return true;
      }
    }
  }
  return false;
}

bool Define_Dipoles::IsResonant() {
  for(auto &D: m_dipolesFF){
    double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
    double mdist;
    for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
      for (auto *v : it->second) {
        mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
        if(mdist<5) {
          return true;
        }
      }
    }
  }
  return false;
}

bool Define_Dipoles::CheckResonant(){
  bool isres = false;
  for(auto &D: m_dipolesII){
    if(CheckResonant(D)) isres=true;
  }
  for(auto &D: m_dipolesFF){
    if(CheckResonant(D)) isres=true;
  }
  return isres;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D, const Vec4D &k){
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)-k).Mass();
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_d- mass_i);
      if(mcheck < mdist) {
        mdist = mcheck;
        D.SetResonaceFlavour(v->in[0]);
      }
    }
  }
  return mdist;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D){
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_i - v->in[0].Mass());
      if(mcheck < mdist) mdist = mcheck;
    }
  }
  return mdist;
}

dipoletype::code Define_Dipoles::WhichResonant(const Vec4D &k){
  if(!HasFSR()) return dipoletype::initial;
  if(!HasISR()) return dipoletype::final;
  double mdistisr(10000), mdistfsr(100000),mdistifi(100000);
  double mindis(10000);
  dipoletype::code min(dipoletype::initial);
  for(auto &D: m_dipolesII){
    mdistisr = ResonantDist(D,k);
    mindis = mdistisr;
    min = dipoletype::initial;
  }
  for(auto &D: m_dipolesFF){
    mdistfsr = ResonantDist(D,k);  
    if(mdistfsr < mdistisr){
      min = dipoletype::final;
    }
  }
  return min;
}

void Define_Dipoles::generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums) {
  if (available_nums.empty()) {
    pairings.push_back(curr_pairing);
    return;
  }
  int curr_num = available_nums[0];
  available_nums.erase(available_nums.begin());
  for(size_t i = 0; i < available_nums.size(); i++) {
    int next_num = available_nums[i];
    available_nums.erase(available_nums.begin() + i);
    curr_pairing.push_back(curr_num);
    curr_pairing.push_back(next_num);
    generate_pairings(pairings, curr_pairing, available_nums);
    curr_pairing.pop_back();
    curr_pairing.pop_back();
    available_nums.insert(available_nums.begin() + i, next_num);
  }
  available_nums.insert(available_nums.begin(), curr_num);
}

std::ostream& Define_Dipoles::operator<<(std::ostream &out) {
  out << "N_in = " << m_in << "\n m_out = " << m_out <<
      "Number of Charged incoming particles = " << m_chargedinparticles.size() << std::endl <<
      "Number of Charged outgoing particles = " << m_chargedoutparticles.size() << std::endl <<
      "Number of Neutral incoming particles = " << m_neutralinparticles.size() << std::endl <<
      "Number of Neutral outgoing particles = " << m_neutraloutparticles.size() << std::endl;
  return out;
}

std::vector<ParticleInfo> Define_Dipoles::ExtractChargedParticles(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta,
    size_t start_index,
    size_t end_index,
    bool is_initial_state) const {
    
    std::vector<ParticleInfo> charged_particles;
    
    for (size_t i = start_index; i < end_index; ++i) {
        if (flavors[i].IntCharge() != 0) {
            charged_particles.emplace_back(
                flavors[i], momenta[i], born_momenta[i], i, is_initial_state);
        }
    }
    
    return charged_particles;
}

std::vector<ParticleInfo> Define_Dipoles::ExtractInitialStateCharged(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta) const {
    
    constexpr size_t INITIAL_STATE_PARTICLES = 2;
    return ExtractChargedParticles(flavors, momenta, born_momenta, 
                                   0, INITIAL_STATE_PARTICLES, true);
}

std::vector<ParticleInfo> Define_Dipoles::ExtractFinalStateCharged(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta) const {
    
    constexpr size_t INITIAL_STATE_PARTICLES = 2;
    return ExtractChargedParticles(flavors, momenta, born_momenta, 
                                   INITIAL_STATE_PARTICLES, flavors.size(), false);
}

Dipole Define_Dipoles::CreateDipole(
    const ParticleInfo& particle1,
    const ParticleInfo& particle2,
    dipoletype::code type) const {
    
    ATOOLS::Flavour_Vector flavors = {particle1.flavor, particle2.flavor};
    ATOOLS::Vec4D_Vector momenta = {particle1.momentum, particle2.momentum};
    ATOOLS::Vec4D_Vector born_momenta = {particle1.born_momentum, particle2.born_momentum};
    
    Dipole dipole(flavors, momenta, born_momenta, type, m_alpha);
    dipole.SetFlavLab(particle1.index, particle2.index);

    return dipole;
}
