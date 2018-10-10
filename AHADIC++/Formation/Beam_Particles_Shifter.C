#include "AHADIC++/Formation/Beam_Particles_Shifter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Beam_Particles_Shifter::Beam_Particles_Shifter(list<Singlet *> * singlets,
					       Soft_Cluster_Handler * softclusters) :
p_singlets(singlets), p_softclusters(softclusters)
{}

Beam_Particles_Shifter::~Beam_Particles_Shifter() {}

void Beam_Particles_Shifter::Init() {
  p_constituents = hadpars->GetConstituents();
}

void Beam_Particles_Shifter::Reset() {
  m_beamparts.clear();
}

bool Beam_Particles_Shifter::operator()() {
  RescueLightClusters();
  ExtractBeamParticles();
  return ShiftBeamParticles(); 
}

void Beam_Particles_Shifter::ExtractBeamParticles() {
  m_beamparts.clear();
  Singlet * singlet, * bsinglet;
  Vec4D  mom(0.,0.,0.,0.);
  double mass(0.);
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();sit++) {
    singlet = (*sit);
    for (list<Proto_Particle *>::iterator pit=singlet->begin();
	 pit!=singlet->end();pit++) {
      if ((*pit)->IsBeam()) {
	mom  += (*pit)->Momentum();
	mass += p_constituents->Mass((*pit)->Flavour());
	m_beamparts.push_back((*pit));
      }
    }
  }
  if (m_beamparts.size()==0 ||
      m_beamparts.size()==1 && dabs(mom.Abs2()-mass*mass)<1.e-6) return;
  if (m_beamparts.size()==1 || mom.Abs2()<sqr(mass+0.1)) {
    for (list<Singlet *>::iterator sit=p_singlets->begin();
	 sit!=p_singlets->end();sit++) {
      singlet = (*sit);
      for (list<Proto_Particle *>::iterator pit=singlet->begin();
	   pit!=singlet->end();pit++) {
	if ((*pit)->IsBeam()) continue;
	mom  += (*pit)->Momentum();
	m_beamparts.push_back(*pit);
	if (mom.Abs2()>sqr(mass)) return;
      }
    }
  } 
}

bool Beam_Particles_Shifter::ShiftBeamParticles() {
  size_t n = m_beamparts.size(), i(0);
  //msg_Out()<<METHOD<<"(n = "<<n<<").\n";
  if (n<=1) return true;
  Vec4D  * moms   = new Vec4D[n];
  double * masses = new double[n];
  
  for (list<Proto_Particle *>::iterator pit=m_beamparts.begin();
       pit!=m_beamparts.end();pit++,i++) {
    moms[i]   = (*pit)->Momentum();
    masses[i] = p_constituents->Mass((*pit)->Flavour());  
  }
  bool success = hadpars->AdjustMomenta(n,moms,masses);
  if (success) {
    i = 0;
    for (list<Proto_Particle *>::iterator pit=m_beamparts.begin();
	 pit!=m_beamparts.end();pit++,i++) {
      (*pit)->SetMomentum(moms[i]);
    }
  }
  delete[] moms;
  delete[] masses;
  return success;
}

bool Beam_Particles_Shifter::RescueLightClusters() {
  Singlet * sing;
  Flavour flav, trip, anti;
  bool    beam, decayed;
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();) {
    sing = (*sit);
    trip = (*sing->begin())->Flavour();
    anti = (*sing->rbegin())->Flavour();
    beam = decayed = false;
    for (list<Proto_Particle *>::iterator pit=sing->begin();
	 pit!=sing->end();pit++) {
      if ((*pit)->IsBeam()) { beam = true; break; }
    }
    if (beam) {
      double mass = sqrt(sing->Mass2());
      if (p_softclusters->MustPromptDecay(trip,anti,mass)) {
	//msg_Out()<<"Gotcha pair ("<<trip<<", "<<anti<<") -> mass = "<<mass<<"!\n";
	if (sing->size()>2) {
	  msg_Out()<<"   have to add gluons to trip/anti.\n";
	  exit(1);
	}
	Cluster cluster((*sing->begin()),(*sing->rbegin()));
	if (p_softclusters->Treat(&cluster,true)) {
	  list<Proto_Particle * > * hadrons = p_softclusters->GetHadrons();
	  //for (list<Proto_Particle * >::iterator hit=hadrons->begin();
	  //   hit!=hadrons->end();hit++) {
	  //msg_Out()<<(**hit);
	  //}
	  decayed = true;
	}
      }
    }
    if (decayed) sit = p_singlets->erase(sit);
    else sit++;
  }
}

