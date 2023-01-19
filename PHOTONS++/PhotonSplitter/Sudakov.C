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
Histogram_2D PHOTONS::Sudakov::s_histo_zy = Histogram_2D(1,0.,1.,100,0.,1.,100);
Histogram_2D PHOTONS::Sudakov::s_histo_zy_llbar = Histogram_2D(1,0.,1.,100,0.,1.,100);
Histogram_2D PHOTONS::Sudakov::s_histo_tdR = Histogram_2D(1,-6.,2.,100,-5,log10(0.5),100);
#endif 

Sudakov::Sudakov() : m_mode(0), m_addedanything(false), p_kinematics(), p_FIkinematics(), m_NInP(0) {}

Sudakov::Sudakov(int mode) : m_mode(mode), m_addedanything(false), p_kinematics(), p_FIkinematics(), m_NInP(0)
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
}

Sudakov::~Sudakov()
{
  m_splitterIds.clear();
  for (size_t i=0; i<m_splitters.size(); i++) delete m_splitters[i];
  if (m_splitters.size()) m_splitters.clear();
  for (size_t i=0; i<m_spectators.size(); i++) delete m_spectators[i];
  if (m_spectators.size()) m_spectators.clear();

  #ifdef PHOTONSPLITTER_DEBUG
  size_t pos(s_histo_base_name.find_last_of("/"));
  if (pos!=std::string::npos) MakeDir(s_histo_base_name.substr(0,pos));
  s_histo_dipole.Finalize();
  s_histo_dipole.Output(s_histo_base_name+"starting_scale.dat");
  s_histo_zy.Finalize();
  s_histo_zy.Output(s_histo_base_name+"zy.dat");
  s_histo_zy_llbar.Finalize();
  s_histo_zy_llbar.Output(s_histo_base_name+"zy_llbar.dat");
  s_histo_tdR.Finalize();
  s_histo_tdR.Output(s_histo_base_name+"tdR.dat");
  #endif
}

void Sudakov::AddSplitter(ATOOLS::Particle *softphoton, const size_t &id)
{
  // id is for blob->GetParticle()
  Particle* p = new Particle();
  p->Copy(softphoton);
  m_remainingphotons.push_back(p);
  m_splitterIds.push_back(id);

  if (m_mode & 1) m_splitters.push_back(new Splitting_Function(p,22,-11,11,1,id));
  if (m_mode & 2) m_splitters.push_back(new Splitting_Function(p,22,-13,13,1,id));
  if (m_mode & 4) m_splitters.push_back(new Splitting_Function(p,22,-15,15,1,id));
  if (m_mode & 8) {
    for (KF_Table::const_iterator it(s_kftable.begin()); it!=s_kftable.end(); ++it) {
      if (it->second->m_hadron && it->second->m_icharge && it->second->m_mass<m_masscutoff) {
        m_splitters.push_back(new Splitting_Function(p,22,it->second->m_kfc,-it->second->m_kfc,2*it->second->m_spin,id));
      }
    }
  }
}

void Sudakov::AddChargedParticle(Particle* p, const size_t &id)
{
  m_spectators.push_back(new Spectator(p->Momentum(),id,p->Flav()));
}

bool Sudakov::ClearAll()
{
  m_splitterIds.clear();
  m_NInP = 0;

  if (!m_addedanything)
  {
    for (Part_Iterator pvit=m_remainingphotons.begin();pvit!=m_remainingphotons.end();++pvit)
    {
      delete *pvit;
    }
  }
  m_addedanything = false;
  for (size_t i=0; i<m_spectators.size(); i++) delete m_spectators[i];
  m_spectators.clear();
  for (size_t i=0; i<m_splitters.size(); i++) delete m_splitters[i];
  m_splitters.clear();
  m_addedparticles.clear();
  m_remainingphotons.clear();
  if (m_splitters.size() != 0 || m_spectators.size() != 0 
      || m_addedparticles.size() != 0 || m_remainingphotons.size() != 0) return false;
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

Spectator* Sudakov::DefineInitialConditions(double &t, Vec4D pphoton)
{
  // choose starting t based on photon and lepton momenta 

  std::vector<Vec4D> pis, pks;
  Spec_Vector initial_splitters = {};
  std::vector<bool> dipoleFI = {};
  std::vector<double> Kp, Kcumu, ts;
  Kcumu = {0.};
  Vec4D ptotal = Vec4D(0.,0.,0.,0.);
  
  // build momentum vectors
  for (int i = 0; i < m_spectators.size(); ++i) {
    if (m_spectators[i]->Id() < m_NInP) continue; // incoming particles cannot be emitters 
    
    ptotal += m_spectators[i]->Momentum();
    // choose spectator for splitter i 
    bool selected = false;
    for (int k = 0; k < m_spectators.size(); ++k) {
      if (!selected && i != k) {
        if ((m_spectators[k]->Id() < m_NInP && m_spectators[k]->GetFlavour().Charge() == m_spectators[i]->GetFlavour().Charge()) || 
            (m_spectators[k]->Id() >= m_NInP && m_spectators[k]->GetFlavour().Charge() == -m_spectators[i]->GetFlavour().Charge())) {
          initial_splitters.push_back(m_spectators[i]); // add as photon emitter
          pis.push_back(m_spectators[i]->Momentum());
          pks.push_back(m_spectators[k]->Momentum()); // set spectator in pks[i]
          dipoleFI.push_back(m_spectators[k]->Id() < m_NInP); // true if FI dipole, false if FF dipole 
          selected = true;
        }
      }
    }
  }

  // build splitting functions 
  for (int i = 0; i < pis.size(); ++i) 
  {
    double Q2, s12, s13, s23, z, y, mi2, mij2, mk2, vtijk, vijk, pipj;

    s12 = 2.*pis[i]*pks[i];
    s13 = 2.*pis[i]*pphoton;
    s23 = 2.*pks[i]*pphoton;

    // calculate masses
    mi2 = mij2 = pis[i].Abs2();
    mk2 = pks[i].Abs2();

    // calculate invariants 
    if (dipoleFI[i]) { 
      Q2 = (pis[i] - pks[i] + pphoton).Abs2();

      y = s13/(s12+s23-s13-2.*mi2);
      z = (s12-s13-mi2)/(s12+s23-2.*s13-mi2);

      pipj  = (mk2-mi2-Q2)*y/2.0;

      // splitting functions for photon emission off fermion
      // no need for coupling and other constants 
      double sf = 1./y * (2./(1.-z+z*y) * (1. + 2.*mi2/(mk2-mi2-Q2)) 
                    - (1.+z) - mi2/pipj 
                    -pipj*mk2/sqr(mk2-mi2-Q2) * 4./sqr(1.-z+z*y) );
      if (IsBad(sf)) { sf = 0; }
      Kp.push_back(sf);
    }
    else { 
      Q2 = (pis[i] + pks[i] + pphoton).Abs2();

      y = s13/(s12+s13+s23);
      z = s12/(s12+s23);

      vtijk = sqrt(Q2)*sqrt(Q2+sqr(mij2)+sqr(mk2)-2.*mij2-2.*mij2*mk2-2.*mk2)/(Q2-mij2-mk2);
      vijk = sqrt(Q2)*sqrt(sqr(2.*mk2+(Q2-mi2-mk2)*(1.-y))-4.*mk2*Q2)/((Q2-mij2-mk2)*(1.-y));

      pipj  = (Q2-mi2-mk2)*y/2.0;

      // splitting functions for photon emission off fermion
      // no need for coupling and other constants 
      double sf = 1./y * ( 2./(1.-z+z*y) - vtijk/vijk * (1.+z + mi2/pipj) );
      if (IsBad(sf)) { sf = 0; }
      Kp.push_back(sf);
    }

    if (i==0) Kcumu[0] = Kp[0]; // first time 
    else Kcumu.push_back(Kcumu[i-1] + Kp[i]); // every other time 

    if (m_virtualityOrdering%2==1) {
      // this t is q2-mij2 = virtuality 
      double q2;
      if (dipoleFI[i]) q2 = p_FIkinematics->GetVirt(Q2,y,z,mi2,0,mk2);
      else q2 = p_kinematics->GetVirt(Q2,y,z,mi2,0,mk2);
      ts.push_back(q2-mij2);
    }
    else {
      // mixed or kT scheme 
      // this t is kT2, Krauss & Schumann 07 eq.27, mj2 = 0
      double kT2;
      if (dipoleFI[i]) kT2 = p_FIkinematics->GetKT2(Q2,y,z,mi2,0,mk2);
      else kT2 = p_kinematics->GetKT2(Q2,y,z,mi2,0,mk2);
      ts.push_back(kT2);
    }
  }

  // select emitter particle and starting scale 
  int winnerIndex;
  if (!m_debug_initProbabilistic) {
    // select winner 
    winnerIndex = std::max_element(Kp.begin(),Kp.end()) - Kp.begin();
  }
  else {
    // select probabilistic
    double rand = ran->Get();
    if (rand < Kcumu[0]/Kcumu.back()) { winnerIndex = 0; }
    else {
      for (int i = 1; i < Kcumu.size(); ++i) {
        if (rand >= Kcumu[i-1]/Kcumu.back() && rand < Kcumu[i]/Kcumu.back()) {
          winnerIndex = i;
        }
      }
    }
  }
  
  if (winnerIndex < 0 || winnerIndex >= ts.size()) { 
    msg_Error() << "No max splitting function found, skipping blob\n"; 
    return 0;
  }

  t = ts[winnerIndex];

  return initial_splitters[winnerIndex];
}

bool Sudakov::Run(ATOOLS::Blob *blob)
{
  // set initial conditions 
  double tstart;
  m_t = 0.;
  for (size_t i : m_splitterIds)
  {
    Spectator* initial_emitter = DefineInitialConditions(tstart,blob->GetParticle(i)->Momentum());

    if (!initial_emitter) return true; // PhotonSplitter cannot act here
    
    for (size_t j=0; j<m_splitters.size(); ++j) {
      if (m_splitters[j]->Id() == i) {
        m_splitters[j]->SetStartScale(tstart);
        if (m_spectatorScheme == 1) { m_splitters[j]->AddSpec(initial_emitter); }
        else {
          for (Spectator* s : m_spectators) {
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
    if (!Generate(blob)) return false;
  }
  return true;
}

bool Sudakov::Generate(ATOOLS::Blob *blob)
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

      Particle *split = blob->GetParticle(m_splitters[i]->Id());
      for (Spectator* spectator : m_splitters[i]->GetSpecs())
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
        y = p_kinematics->GetYVirt(Q2,t+m_splitters[ind]->Mass2A(),m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
                            m_splitters[ind]->Mass2Spec());
      }
      else {
        y = p_kinematics->GetY(Q2,t,z,m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
                            m_splitters[ind]->Mass2Spec()); // this t is kT2
      }

      if (y >= 1) break;
      if (y <= 0) break;

      f = (*m_splitters[ind])(t,z,y,Q2); // this t is only used for cutoff 
      g = m_splitters[ind]->OverEstimated();

      if (f/g > 1.) {
        msg_Debugging() << "f = " << f << ", g = " << g << "\n";
        msg_Error() << "Error: splitting function overestimate is not an overestimate!\n";
        return false;
      }

      #ifdef PHOTONSPLITTER_DEBUG
      s_histo_zy_llbar.Insert(z,y);
      #endif 

      // veto 
      if (f/g > ran->Get()) {
        msg_Debugging() << "A photon split!\n";
        msg_Debugging() << "Spectator flavour is " << m_splitters[ind]->GetSpectator()->GetFlavour() << "\n";
        if (abs(m_splitters[ind]->FlB().Kfcode()) == 15) { msg_Debugging() << "Split into taus!\n"; }
        else if (abs(m_splitters[ind]->FlB().Kfcode()) == 211) { msg_Debugging() << "Split into pions!\n"; }
        else if (abs(m_splitters[ind]->FlB().Kfcode()) == 321) { msg_Debugging() << "Split into kaons!\n"; }
        else if (abs(m_splitters[ind]->FlB().Kfcode()) == 13) { msg_Debugging() << "Split into muons\n"; }

        phi = 2*M_PI * ran->Get();
        bool madeKinematics = p_kinematics->MakeKinematics(z,y,phi,pij,pk,pi,pj,
          m_splitters[ind]->Mass2B(),m_splitters[ind]->Mass2C(),
          m_splitters[ind]->Mass2Spec(), m_splitters[ind]->Mass2A());
        if (!madeKinematics) THROW(fatal_error, "Invalid kinematics");

        // create new particles 
        Particle *newparticle = new Particle(-1,m_splitters[ind]->FlB(),pi,'F');
        Particle *newantiparticle = new Particle(-1,m_splitters[ind]->FlC(),pj,'F');
        m_addedparticles.push_back(newparticle);
        m_addedparticles.push_back(newantiparticle);

        m_splitters[ind]->GetSpectator()->SetMomentum(pk); 
        
        // remove photon as splitter
        Part_Iterator PLIt = std::find(m_remainingphotons.begin(),m_remainingphotons.end(),m_splitters[ind]->GetSplitter());
        m_remainingphotons.erase(PLIt);
        delete *PLIt; // since it will not be added back

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
