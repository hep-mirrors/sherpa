#include "PHOTONS++/PhotonSplitter/Sudakov.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Blob.H"
#include "PHOTONS++/PhotonSplitter/Splitting_Functions.H"
#include <algorithm>
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Histogram_2D.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "AHADIC++/Main/Ahadic.H"

using namespace PHOTONS;
using namespace ATOOLS;

#ifdef PHOTONSPLITTER_DEBUG
std::string PHOTONS::Sudakov::s_histo_base_name;
Histogram PHOTONS::Sudakov::s_histo_dipole = Histogram(10,1e-18,1e3,200,"");
Histogram_2D PHOTONS::Sudakov::s_histo_tdR = Histogram_2D(1,-6.,2.,100,-5,log10(0.5),100);
#endif 

Sudakov::Sudakov() : m_mode(0), m_addedanything(false), m_NInP(0) {}

Sudakov::Sudakov(int mode) : m_mode(mode), m_addedanything(false), m_NInP(0), m_decay(0)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };

  m_masscutoff = s["PHOTON_SPLITTER_MAX_HADMASS"].Get<double>();

  // YFS_PHOTON_SPLITTER_ORDERING_SCHEME:
  // 0 = transverse momentum ordering 
  // 1 = virtuality ordering 
  // 2 = mixed ordering - kT for initial conditions, virt for photon splitting (default)
  m_virtualityOrdering = s["PHOTON_SPLITTER_ORDERING_SCHEME"].Get<int>();
  // YFS_PHOTON_SPLITTER_SPECTATOR_SCHEME:
  // 0 = all final-state charged particles that exist prior to this module being called (default)
  // 1 = only the final-state charged particle that the soft photon is calculated to be emitted off 
  m_spectatorScheme = s["PHOTON_SPLITTER_SPECTATOR_SCHEME"].Get<int>();

  m_debug_initProbabilistic = s["PHOTON_SPLITTER_STARTING_SCALE_SCHEME"].Get<int>();

  // replace by a proper flavour-dependent read-in, for now use the enhance
  // factor for all leptons only
  double enh(s["PHOTON_SPLITTER_ENHANCE_FACTOR"].Get<double>());
  m_enhancefac[kf_e]=enh;
  if (m_enhancefac[kf_e]!=1.) msg_Info()<<"METHOD(): Enhancing P->ee splittings by factor "<<m_enhancefac[kf_e]<<std::endl;
  m_enhancefac[kf_mu]=enh;
  if (m_enhancefac[kf_mu]!=1.) msg_Info()<<"METHOD(): Enhancing P->mumu splittings by factor "<<m_enhancefac[kf_mu]<<std::endl;
  m_enhancefac[kf_tau]=enh;
  if (m_enhancefac[kf_tau]!=1.) msg_Info()<<"METHOD(): Enhancing P->tautau splittings by factor "<<m_enhancefac[kf_tau]<<std::endl;

  #ifdef PHOTONSPLITTER_DEBUG
  s_histo_base_name = s["PHOTON_SPLITTER_HISTO_BASE_NAME"].Get<std::string>();
  #endif
}

void Sudakov::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };

  s["PHOTON_SPLITTER_MAX_HADMASS"].SetDefault(0.5);
  s["PHOTON_SPLITTER_ORDERING_SCHEME"].SetDefault(2);
  s["PHOTON_SPLITTER_SPECTATOR_SCHEME"].SetDefault(0);
  s["PHOTON_SPLITTER_STARTING_SCALE_SCHEME"].SetDefault(1);
  s["PHOTON_SPLITTER_HISTO_BASE_NAME"].SetDefault("histos/");
  s["PHOTON_SPLITTER_ENHANCE_FACTOR"].SetDefault(1.);
}

Sudakov::~Sudakov()
{
  ClearAll();

  #ifdef PHOTONSPLITTER_DEBUG
  size_t pos(s_histo_base_name.find_last_of("/"));
  if (pos!=std::string::npos) MakeDir(s_histo_base_name.substr(0,pos));
  s_histo_dipole.Finalize();
  s_histo_dipole.Output(s_histo_base_name+"starting_scale.dat");
  s_histo_tdR.Finalize();
  s_histo_tdR.Output(s_histo_base_name+"tdR.dat");
  #endif
}

void Sudakov::SetNInParticles(int n)
{
  m_NInP = n;
  if (n == 1) { m_decay = 1; }
  else if (n == 2) { m_decay = 0; }
  else { msg_Out() << "More than 2 incoming particles...\n"; m_decay = 0; }
}

void Sudakov::AddSplitter(ATOOLS::Particle *softphoton, const size_t &id)
{
  // id is for blob->GetParticle()
  YFS_Particle * p = new YFS_Particle(softphoton->Momentum(),id,softphoton->Flav());
  m_remainingphotons.push_back(p);
  m_splitterIds.push_back(id);
  m_splittingphotons.push_back(p); // this list will stay invariant

  if (m_mode & 1) m_splitters.push_back(new SF_FF(p,kf_photon,-kf_e,kf_e,1,id,m_enhancefac[kf_e]));
  if (m_mode & 2) m_splitters.push_back(new SF_FF(p,kf_photon,-kf_mu,kf_mu,1,id,m_enhancefac[kf_mu]));
  if (m_mode & 4) m_splitters.push_back(new SF_FF(p,kf_photon,-kf_tau,kf_tau,1,id,m_enhancefac[kf_tau]));
  if (m_mode & 8) {
    for (KF_Table::const_iterator it(s_kftable.begin()); it!=s_kftable.end(); ++it) {
      if (it->second->m_hadron && it->second->m_icharge && it->second->m_mass<m_masscutoff) {
        kf_code kfc(it->second->m_kfc);
        double enh(m_enhancefac.find(kfc)!=m_enhancefac.end()?m_enhancefac[kfc]:1.);
        m_splitters.push_back(new SF_FF(p,kf_photon,kfc,-kfc,2*it->second->m_spin,id,enh));
      }
    }
  }
}

void Sudakov::AddChargedParticle(Particle* p, const size_t &id)
{
  m_spectators.push_back(new YFS_Particle(p->Momentum(),id,p->Flav()));
}

bool Sudakov::ClearAll()
{
  m_splitterIds.clear();
  m_NInP = 0;

  m_addedanything = false;

  for (size_t i=0; i<m_spectators.size(); i++) delete m_spectators[i];
  m_spectators.clear();

  for (size_t i=0; i<m_splitters.size(); i++) delete m_splitters[i];
  m_splitters.clear();

  for (size_t i=0; i<m_addedparticles.size(); i++) delete m_addedparticles[i];
  m_addedparticles.clear();

  m_remainingphotons.clear();

  for (size_t i=0; i<m_splittingphotons.size(); i++) delete m_splittingphotons[i];
  m_splittingphotons.clear();

  if (m_splitters.size() != 0 || m_spectators.size() != 0 
      || m_addedparticles.size() != 0 || m_splittingphotons.size() != 0) return false;
  return true;
}

void Sudakov::SetCutoff()
{
  // called from Photon_Splitter before Run()
  // cutoff is 4*mass^2 of lowest-mass fermion in SFs
  m_t0 = 1.;
  for (size_t i=0; i<m_splitters.size(); i++) {
    m_t0 = std::min(m_splitters[i]->Cutoff(),m_t0);
  }
}

YFS_Particle* Sudakov::DefineInitialConditions(double &t, Vec4D pphoton)
{
  // loop over possible emitters
  for (size_t i=0; i < m_spectators.size(); ++i) {
    // loop over possible spectators 
    for (size_t j=0; j < m_spectators.size(); ++j) {
      if (i == j) continue;
    }
  }
  return 0; // to do
}

bool Sudakov::Run()
{
  // set initial conditions 
  double tstart;
  m_t = 0.;
  for (size_t i = 0; i < m_splittingphotons.size(); ++i)
  {
    YFS_Particle* initial_emitter = DefineInitialConditions(tstart,m_splittingphotons[i]->Momentum());

    if (!initial_emitter) return true; // PhotonSplitter cannot act here
    
    for (size_t j=0; j<m_splitters.size(); ++j) {
      if (m_splitters[j]->Id() == i) {
        m_splitters[j]->SetStartScale(tstart);
        if (m_spectatorScheme == 1) { m_splitters[j]->AddSpec(initial_emitter); }
        else {
          for (YFS_Particle* s : m_spectators) {
            if (s->Id() >= m_NInP) { // for now, no W spectator for photon splittings 
              m_splitters[j]->AddSpec(s);
            }
          }
        }
      }
    }
    if (tstart > m_t) {
      m_t = tstart;
    }
    #ifdef PHOTONSPLITTER_DEBUG
    s_histo_dipole.Insert(tstart);
    #endif 
  }
  // run
  while (m_t > m_t0) 
  {
    if (!Generate()) return false;
  }
  return true;
}

bool Sudakov::Generate()
{
  double t, Q2, zmax, zmin, f, g, tgen, z, y, phi;
  int ind;
  ATOOLS::Vec4D pij, pi, pj, pk;
  while (m_t > m_t0)
  {
    t = m_t0; // comparing value
    for (size_t i=0; i<m_splitters.size(); i++)
    {
      if (!m_splitters[i]->On()) continue; // if photon no longer exists
      if (m_t > m_splitters[i]->StartScale()) continue; // if we're above the photon's starting scale
      if (m_t < m_splitters[i]->Cutoff()) continue; // if we're below the cutoff

      YFS_Particle *split = m_splitters[i]->GetSplitter();
      for (YFS_Particle* spectator : m_splitters[i]->GetSpecs())
      {
        // compute z boundaries
        if (m_virtualityOrdering) {
          // the z boundaries are complicated, we just reject if the generated z is not allowed later
          zmax = 1.;
          zmin = 1.-zmax;
          // double Q2tmp, mui2, muk2, viji, vijk, tQ;
          // Q2tmp = (split->Momentum() + spectator->Momentum()).Abs2();
          // mui2 = m_splitters[i]->Mass2B()/Q2tmp;
          // muk2 = spectator->Momentum().Abs2()/Q2tmp;
          // tQ = m_t/Q2tmp;

          // viji = sqrt(1. - 4*mui2/tQ);
          // vijk = sqrt(1. + 4*muk2*tQ/sqr(1-tQ-muk2));

          // zmax = 1./2. * (1.+viji*vijk);
          // zmin = 1./2. * (1.-viji*vijk);
        }
        else {
          //compute z boundaries
          // eq. 45/46 in arxiv:0709.1027. Note kT^2_MAX = Q2/4
          // note (split->Momentum() + spectator->Momentum()).Abs2() = Q2, but not assigned yet
          if ((split->Momentum() + spectator->Momentum()).Abs2() < 4*m_t0) continue;
          zmax = 0.5 * (1 + sqrt(1-4*m_t0/(split->Momentum() + spectator->Momentum()).Abs2()));
          zmin = 1 - zmax;
        }

        // trial emission
        g = m_splitters[i]->OverIntegrated(zmin, zmax) / (2*M_PI); // this also sets the limits
        tgen = m_t * pow(ran->Get(),1./g);

        if (tgen > t) {
          t = tgen;
          Q2 = (split->Momentum() + spectator->Momentum()).Abs2();
          pij = split->Momentum();
          msg_Debugging() << "Winner found with energy " << pij[0] << " and flb = " << m_splitters[i]->FlB() << "\n";
          pk = spectator->Momentum();
          m_splitters[i]->SetSpec(spectator);
          ind = i;
        }
      }
    }
    m_t = t;
    // if winner found 
    if (t > m_t0) {
      z = m_splitters[ind]->Z();
      if (m_virtualityOrdering) {
        // virtuality or mixed scheme 
        // what is passed in the second argument is q2, not t. t = q2 - m2 so we add the mass here. For photon splittings m2=0.
        y = m_splitters[ind]->KinFF()->GetYVirt(Q2,t+m_splitters[ind]->Mass2A(),m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
                            m_splitters[ind]->Mass2Spec());
      }
      else {
        y = m_splitters[ind]->KinFF()->GetY(Q2,t,z,m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
                            m_splitters[ind]->Mass2Spec()); // this t is kT2
      }

      if (y >= 1) break;
      if (y <= 0) break;

      f = (*m_splitters[ind])(t,z,y,Q2); // this t is only used for cutoff 
      g = m_splitters[ind]->OverEstimated(z);

      if (f/g > 1.) {
        msg_Debugging() << "f = " << f << ", g = " << g << "\n";
        msg_Error() << "Error: splitting function overestimate is not an overestimate!\n";
        return false;
      }

      // veto 
      if (f/g > ran->Get()) {
        msg_Debugging() << "A photon split!\n";
        msg_Debugging() << "Spectator flavour is " << m_splitters[ind]->GetSpectator()->GetFlavour() << "\n";
        if (m_splitters[ind]->FlB().Kfcode() == 15) { msg_Debugging() << "Split into taus!\n"; }
        else if (m_splitters[ind]->FlB().Kfcode() == 211) { msg_Debugging() << "Split into pions!\n"; }
        else if (m_splitters[ind]->FlB().Kfcode() == 321) { msg_Debugging() << "Split into kaons!\n"; }
        else if (m_splitters[ind]->FlB().Kfcode() == 13) { msg_Debugging() << "Split into muons\n"; }

        phi = 2*M_PI * ran->Get();
        bool madeKinematics = m_splitters[ind]->KinFF()->MakeKinematics(z,y,phi,pij,pk,pi,pj,
          m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
          m_splitters[ind]->Mass2Spec(), m_splitters[ind]->Mass2A());
        if (!madeKinematics) THROW(fatal_error, "Invalid kinematics");

        // set ID to -1 such that this will get initialised correctly by the blob
        YFS_Particle *newparticle = new YFS_Particle(pi,-1,m_splitters[ind]->FlB());
        YFS_Particle *newantiparticle = new YFS_Particle(pj,-1,m_splitters[ind]->FlC());
        m_addedparticles.push_back(newparticle);
        m_addedparticles.push_back(newantiparticle);

        m_splitters[ind]->GetSpectator()->SetMomentum(pk); 
        
        // remove photon as splitter
        YFS_Particle_List::iterator PLIt = std::find(m_remainingphotons.begin(),m_remainingphotons.end(),m_splitters[ind]->GetSplitter());
        m_remainingphotons.erase(PLIt);

        // turn off future splittings of this photon
        for (size_t i=0; i<m_splitters.size(); i++) {
          if (m_splitters[i]->Id() == m_splitters[ind]->Id()){
            m_splitters[i]->SetOn(false);
          }
        }

        #ifdef PHOTONSPLITTER_DEBUG
        double dtheta2 = sqr(pi.Theta()-pj.Theta());
        double dphi2 = sqr(pi.Phi()-pj.Phi());
        s_histo_tdR.Insert(log10(m_splitters[ind]->StartScale()),log10(sqrt(dtheta2+dphi2)));
        #endif

        m_addedanything = true;
        return true;
      }
    }
  }
 return true;
}
