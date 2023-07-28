#include "AHADIC++/Formation/Beam_Particles_Shifter.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Poincare.H"

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
  //msg_Out()<<"      --- "<<METHOD<<" for "<<p_singlets->size()<<" singlets.\n";
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
      (m_beamparts.size()==1 && dabs(mom.Abs2()-mass*mass)<1.e-6)) return;
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
  //msg_Out()<<"     --- "<<METHOD<<" for "<<m_beamparts.size()<<" particles.\n";
  size_t n = m_beamparts.size(), i(0);
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

void Beam_Particles_Shifter::RescueLightClusters() {
  //msg_Out()<<"    --- "<<METHOD<<" for "<<p_singlets->size()<<" singlets.\n";
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
	if (sing->size()>2) {
	  sing->StripSingletOfGluons();
	}
	Cluster cluster((*sing->begin()),(*sing->rbegin()));
	if (p_softclusters->Treat(&cluster,true)==1) {
	  decayed = true;
	}
	else {
	  Flavour transition = p_softclusters->LowestTransition(trip,anti);
	  double  transmass  = transition.Mass();
	  Proto_Particle * recoiler = GetRecoilPartner(transmass,cluster.Momentum(),sing);
	  if (recoiler && ShuffleMomenta(recoiler,&cluster,transition,transmass)) {
	    decayed = true;
	  }
	}
      }
    }
    if (decayed) {
      delete (*sit);
      sit = p_singlets->erase(sit);
    }
    else sit++;
  }
}

bool Beam_Particles_Shifter::
ShuffleMomenta(Proto_Particle * recoiler,Cluster * cluster,const Flavour & target,
	       const double & targetmass) {
  Vec4D momR = recoiler->Momentum(), momC = cluster->Momentum(), lab = momR+momC;
  Poincare boost(lab);
  boost.Boost(momR); boost.Boost(momC);
  Poincare rotat(momR,Vec4D::ZVEC);
  double cmsE = momR[0] + momC[0];
  rotat.Rotate(momR); rotat.Rotate(momC);
  double mR2 = sqr(p_constituents->Mass(recoiler->Flavour()));
  double mC2 = sqr(targetmass);
  double ER  = (sqr(cmsE)+mR2-mC2)/(2.*cmsE);
  double EC  = (sqr(cmsE)+mC2-mR2)/(2.*cmsE);
  double p   = sqrt(Max(0.,sqr(ER)-mR2));
  momR = Vec4D(ER,0.,0.,p); momC = Vec4D(EC,0.,0.,-p);
  rotat.RotateBack(momR); rotat.RotateBack(momC);
  boost.BoostBack(momR); boost.BoostBack(momC);
  recoiler->SetMomentum(momR);
  Proto_Particle * part = new Proto_Particle(target,momC,false);
  p_softclusters->GetHadrons()->push_back(part);
  return true;
}

Proto_Particle *Beam_Particles_Shifter::
GetRecoilPartner(const double & targetmass,const ATOOLS::Vec4D & mom,
		 const Singlet * veto) {
  Proto_Particle * recoiler(NULL);
  double pt2max = 1.e6, pt2, mass;
  list<Proto_Particle *> * hadrons = p_softclusters->GetHadrons();
  if (hadrons->size()>0) {
    for (list<Proto_Particle *>::iterator hit=hadrons->begin();
	 hit!=hadrons->end();hit++) {
      mass = sqrt((mom+(*hit)->Momentum()).Abs2());
      pt2  = (*hit)->Momentum().PPerp2();
      if (mass>targetmass+(*hit)->Flavour().Mass() &&
	  pt2<pt2max) {
	recoiler = (*hit);
	pt2max = pt2; 
      }
    }
  }
  if (!recoiler) {
    Proto_Particle * lastresort(NULL);
    double pt2max_last = 1.e12;
    for (list<Singlet *>::iterator sit=p_singlets->begin();
	 sit!=p_singlets->end();sit++) {
      Singlet * sing = (*sit);
      if (sing==veto) continue;
      for (list<Proto_Particle *>::iterator pit=sing->begin();
	   pit!=sing->end();pit++) {
	mass = sqrt((mom+(*pit)->Momentum()).Abs2());
	pt2  = (*pit)->Momentum().PPerp2();
	if (mass>targetmass+(*pit)->Flavour().Mass()) {
	  if (pt2<pt2max && (*pit)->IsBeam()) {
	    recoiler = (*pit);
	    pt2max = pt2; 
	  }
	  else if (!recoiler && pt2<pt2max_last) {
	    lastresort  = (*pit);
	    pt2max_last = pt2; 
	  }
	}
      }
    }
    if (!recoiler) recoiler = lastresort;
  }
  return recoiler;
}

