#include "AHADIC++/Formation/Beam_Particles_Shifter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Beam_Particles_Shifter::Beam_Particles_Shifter(list<Singlet *> * singlets) :
  p_singlets(singlets)
{}

Beam_Particles_Shifter::~Beam_Particles_Shifter() {}

void Beam_Particles_Shifter::Init() {
  p_constituents = hadpars->GetConstituents();
}

bool Beam_Particles_Shifter::operator()() {
  ExtractBeamParticles();
  return ShiftBeamParticles();
}

void Beam_Particles_Shifter::ExtractBeamParticles() {
  m_beamparts.clear();
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();sit++) {
    Singlet * singlet = (*sit);
    for (list<Proto_Particle *>::iterator pit=singlet->begin();
	 pit!=singlet->end();pit++) {
      if ((*pit)->IsBeam()) m_beamparts.push_back((*pit));
    }
  }
}

bool Beam_Particles_Shifter::ShiftBeamParticles() {
  size_t n = m_beamparts.size(), i(0);
  if (n==0) return true;
  Vec4D  * moms   = new Vec4D[n];
  double * masses = new double[n];
  
  for (list<Proto_Particle *>::iterator pit=m_beamparts.begin();
       pit!=m_beamparts.end();pit++,i++) {
    moms[i]   = (*pit)->Momentum();
    masses[i] = p_constituents->Mass((*pit)->Flavour());  
  }
  if (hadpars->AdjustMomenta(n,moms,masses)) {
    i = 0;
    for (list<Proto_Particle *>::iterator pit=m_beamparts.begin();
	 pit!=m_beamparts.end();pit++,i++) {
      (*pit)->SetMomentum(moms[i]);
    }
    return true;
  }
  return false;
}
