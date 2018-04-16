#include "AHADIC++/Formation/Singlet_Checker.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

/*
Define: 

m_min = minimal constituent mass of quark (m_minQmass)
q     = quark or diquark
m_q   = constituent mass

Strategy:

if 2-particle singlet:
** gg: m_gg < 4.*m_min  -->  q qbar
       m_gg < 2.*m_min  -->  pi+ pi-/pi0 pi0/eta + gamma
                             in ration 60:30:10
       m_gg < 2.*m_pi   -->  pi0 + gamma
       m_gg < m_pi      -->  gamma + gamma
       
       must write code in Soft_Cluster_Handler

** q q(bar): possible problem for massless quarks ...
       m_qq < m_q+m_q --> hadron + hadron or hadron + gamma or gamma + gamma
                          (if q qbar is not charged)
                          by existing code in Soft_Cluster_Handler


if multi-particle singlet
       m_gg < 4.*m_min --> q qbar
              must split singlet:
	      -- all-gluon singlet: split & reorder
	      -- otherwise: split & add new singlet to end of list,
	                    kick respective particles out
       m_gg < 2.* m_min --> combine the gluons into one:
                 find third parton (spectator) to take recoil

       m_qg similar to m_gg 

 */



Singlet_Checker::Singlet_Checker(list<Singlet *> * singlets,
				 Soft_Cluster_Handler * softclusters) :
  Singlet_Tools(),
  p_singlets(singlets), p_softclusters(softclusters),
  p_hadrons(softclusters->GetHadrons()),
  m_direct_transitions(0)
{}

Singlet_Checker::~Singlet_Checker() {
  msg_Tracking()<<METHOD<<" with "<<m_direct_transitions
		<<" direct enforced transitions in total.\n";
}

void Singlet_Checker::Init() {
  Singlet_Tools::Init();
  m_splitter.Init();
}

void Singlet_Checker::Reset() {
  m_badones.clear();
  m_transitions.clear();
}

bool Singlet_Checker::operator()() {
  list<Singlet *>::iterator lsit(p_singlets->begin());
  while (lsit!=p_singlets->end()) {
    p_singlet = (*lsit);
    // check if singlet is too light
    // (mass smaller than summed constituent masses - may hint at problem) 
    if (!CheckSinglet()) {
      // there are only two partons in it - this will have to be fixed
      // we put all of those into a separate list to be deatl with in
      // rescue system
      if (p_singlet->size()==2) {
	m_badones.push_back(lsit);
	lsit++;
      }
      // more than two partons - fusing partons may help in retry.
      // if fusing does not work we will re-try the event
      else if (!FusePartonsInLowMassSinglet()) {
	msg_Error()<<METHOD<<" throws error - fusing didn't work out.\n"
		   <<(*p_singlet)<<"\n";
	return false;
      }
    }
    // everything is fine -- the singlet can go on as normal.
    else lsit++;
  }
  // invoking the rescue system, if neccessary.
  if (m_badones.size()>0) {
    if (!DealWithProblematicSinglets()) {
      msg_Error()<<METHOD<<" throws error - no rescue possible.\n";
      if (msg_LevelIsTracking()) {
	for (list<list<Singlet *>::iterator>::iterator bit=m_badones.begin();
	     bit!=m_badones.end();bit++) {
	  msg_Out()<<(***bit)<<"\n";
	}
      }
      return false;
    }
  }
  return true;
}

bool Singlet_Checker::CheckSinglet() {
  // Checking the mass for pairs of colour-connected particles
  list<Proto_Particle *>::iterator plit1(p_singlet->begin()), plit2(plit1);
  plit2++;
  if (p_singlet->size()==2) {
    // check if singlet must transit directly to hadron
    double mass    = sqrt(((*plit1)->Momentum()+(*plit2)->Momentum()).Abs2());
    double minmass = Min(p_softclusters->MinSingleMass((*plit1)->Flavour(),
						       (*plit2)->Flavour()),
			 p_softclusters->MinDoubleMass((*plit1)->Flavour(),
						       (*plit2)->Flavour()));
    if (mass < minmass) return false;
  }
  while (plit2!=p_singlet->end()) {
    // this is more or less a plausibility check driven by some external mass
    // parameter (minmass2) --- same below for "gluon ring"
    p_part1 = (*plit1);
    p_part2 = (*plit2);
    if (!CheckMass(p_part1,p_part2)) return false;
    plit2++;
    plit1++;
  }
  // is gluon "ring" must also check for pair made of first and last particle
  if (p_singlet->front()->Flavour().IsGluon() &&
      p_singlet->back()->Flavour().IsGluon()) {
    p_part1 = (*plit1);
    p_part2 = p_singlet->front();
    if (!CheckMass(p_part1,p_part2)) return false;
  }
  return true;
}

bool Singlet_Checker::FusePartonsInLowMassSinglet() {
  if (p_singlet->front()->Flavour().IsGluon() &&
      sqrt(m_mass) > 2.*m_minQmass && m_splitter(p_part1,p_part2)) {
    // gluons heavy enough to be replaced by two quarks, splits gluon ring
    // and necessitates reordering the singlet to have a quark as first
    // particle.
    p_singlet->Reorder();
    return true;
  }
  return p_singlet->Combine(p_part1,p_part2);
}

bool Singlet_Checker::DealWithProblematicSinglets() {
  m_transitions.clear();
  SortProblematicSinglets();
  if (m_transitions.size()>1) {
    if (!TransitProblematicSinglets()) {
      msg_Error()<<METHOD<<" throws error for more than one transition.\n";
      return false;
    }
  }
  else if (m_transitions.size()==1) {
    // have one transition, try to balance against other tricky singlets
    // to sort them out - two birds with one stone.
    if (FindOtherSingletToTransit()) {
      if (!TransitProblematicSinglets()) {
	msg_Error()<<METHOD<<" throws error for one transition (1).\n";
	return false;
      }
    }
    // if this does not work, we'll try to find a "regular" singlet ....
    else if (FindRecoilerForTransit()) {
      if (!TransitProblematicSingletWithRecoiler()) {
	msg_Error()<<METHOD<<" throws error for one transition (2).\n";
	return false;
      }
    }
    // or just leave the particle off-shell.
    else {
      Vec4D mom   = m_transitions.begin()->first->Momentum();
      bool isbeam = false;
      Proto_Particle * part = new Proto_Particle(m_transitions.begin()->second,
						 mom,false,isbeam);
      p_hadrons->push_back(part);
      m_direct_transitions++;
      msg_Tracking()<<METHOD<<" with a transition for "
		    <<"("<<p_singlets->size()<<" singlets).\n"
		    <<m_transitions.begin()->second<<" from "
		    <<(*m_transitions.begin()->first)<<"\n";
      return true;
    }
  }
  ForcedDecays();
  return (m_badones.size()==0);
}

void Singlet_Checker::SortProblematicSinglets() {
  list<list<Singlet *>::iterator>::iterator bit=m_badones.begin();
  while (bit!=m_badones.end()) {
    p_singlet = (**bit);
    Flavour flav1 = p_singlet->front()->Flavour();
    Flavour flav2 = p_singlet->back()->Flavour();
    if (!flav1.IsGluon()) {
      Flavour had = p_softclusters->LowestTransition(flav1,flav2);
      if (had.Mass()>sqrt(p_singlet->Mass2())) {
	m_transitions[p_singlet] = had;
	p_singlets->erase((*bit));
	bit = m_badones.erase(bit);
	continue;
      }
    }
    bit++;
  }
}

bool Singlet_Checker::FindOtherSingletToTransit() {
  if (m_badones.size()==0) return false;
  list<list<Singlet *>::iterator>::iterator bit=m_badones.begin();
  list<list<Singlet *>::iterator>::iterator hit=m_badones.end();
  Flavour hadron(kf_none);
  double  massdiff(1.e6), masshad;
  while (bit!=m_badones.end()) {
    p_singlet = (**bit);
    Flavour flav1 = p_singlet->front()->Flavour();
    Flavour flav2 = p_singlet->back()->Flavour();
    if (!flav1.IsGluon()) {
      Flavour hadtest = p_softclusters->LowestTransition(flav1,flav2);
      masshad = hadtest.Mass();
      if (dabs(masshad-sqrt(p_singlet->Mass2()))<massdiff) {
	hadron   = hadtest;
	hit      = bit;
	massdiff = dabs(masshad-sqrt(p_singlet->Mass2()));
      }
    }
    bit++;
  }
  if (hit!=m_badones.end() && hadron!=Flavour(kf_none)) {
    m_transitions[(**hit)] = hadron;
    p_singlets->erase(*hit);
    m_badones.erase(hit);
    return true;
  }
  msg_Error()<<METHOD<<" throws error: no partner found.\n";
  return false;
}

bool Singlet_Checker::FindRecoilerForTransit() {
  if (m_transitions.size()!=1 && m_badones.size()!=1) abort();
  m_singletmom = m_transitions.begin()->first->Momentum();
  m_targetmass = m_transitions.begin()->second.Mass();
  p_recoiler = NULL;
  double testmass2 = 0., mass2;
  bool   isbeam(false);
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();sit++) {
    p_singlet = (*sit);
    mass2 = p_singlet->Mass2();
    if (mass2>testmass2 && TestRecoiler()) {
      if (p_recoiler==NULL || !isbeam) {
	p_recoiler = p_singlet;
	testmass2  = mass2;
	isbeam     = (p_recoiler->front()->IsBeam() ||
		      p_recoiler->back()->IsBeam());
      }
      else {
	if (isbeam && (p_singlet->front()->IsBeam() ||
		       p_singlet->back()->IsBeam())) {
	  p_recoiler = p_singlet;
	  testmass2  = mass2;
	}
      }
    }
  }
  return (p_recoiler!=NULL);
}

bool Singlet_Checker::TestRecoiler() {
  return ((m_singletmom+p_singlet->Momentum()).Abs2() >
	  sqr(m_targetmass+sqrt(p_singlet->Mass2()) ) );
}

bool Singlet_Checker::TransitProblematicSinglets() {
  size_t   n      = m_transitions.size(), i=0;
  Vec4D *  moms   = new Vec4D[n],  totmom  = Vec4D(0.,0.,0.,0.);
  double * masses = new double[n], totmass = 0;
  for (map<Singlet *,Flavour>::iterator tit=m_transitions.begin();
       tit!=m_transitions.end();tit++,i++) {
    totmom  += moms[i]   = tit->first->Momentum();
    totmass += masses[i] = tit->second.Mass();
  }
  if (totmom.Abs2()<sqr(totmass)) {
    for (map<Singlet *,Flavour>::iterator tit=m_transitions.begin();
	 tit!=m_transitions.end();tit++,i++) {
      msg_Debugging()<<"Singlet with "<<tit->first->Momentum()<<" --> "
		     <<tit->second<<" ("<<tit->second.Mass()<<")\n";
    }
    delete[] moms;
    delete[] masses;
    return false;
  }
  bool success = hadpars->AdjustMomenta(n,moms,masses);
  if (success) {
    i = 0;
    for (map<Singlet *,Flavour>::iterator tit=m_transitions.begin();
	 tit!=m_transitions.end();tit++,i++) {
      bool isbeam = (tit->first->front()->IsBeam() ||
		     tit->first->back()->IsBeam());
      Proto_Particle * part = new Proto_Particle(tit->second,moms[i],
						 false,isbeam);
      p_hadrons->push_back(part);
      delete tit->first;
    }
    m_transitions.clear();
  }
  delete[] moms;
  delete[] masses;
  return success;
}

bool Singlet_Checker::TransitProblematicSingletWithRecoiler() {
  Vec4D *  moms   = new Vec4D[2];
  double * masses = new double[2];
  p_singlet       = m_transitions.begin()->first;
  Flavour hadron  = m_transitions.begin()->second;
  moms[0]   = p_singlet->Momentum();
  moms[1]   = p_recoiler->Momentum();
  masses[0] = hadron.Mass();
  masses[1] = sqrt(p_recoiler->Mass2());
  bool success = hadpars->AdjustMomenta(2,moms,masses);
  if (success) {
    bool isbeam = (p_singlet->front()->IsBeam() ||
		   p_singlet->back()->IsBeam());
    Proto_Particle * part = new Proto_Particle(hadron,moms[0],false,isbeam);
    p_hadrons->push_back(part);
    BoostRecoilerInNewSystem(moms[1]);
    delete p_singlet;
    m_transitions.clear();
  }
  delete[] moms;
  delete[] masses;
  return success;
}

bool Singlet_Checker::BoostRecoilerInNewSystem(const Vec4D & newmom) {
  Vec4D oldmom = p_recoiler->Momentum();
  Poincare intocms(p_recoiler->Momentum());
  Poincare intonew(newmom);
  for (list<Proto_Particle *>::const_iterator pliter=p_recoiler->begin();
       pliter!=p_recoiler->end();pliter++) {
    Vec4D mom = (*pliter)->Momentum();
    intocms.Boost(mom);
    intonew.BoostBack(mom);
    (*pliter)->SetMomentum(mom);
  }
  return true;
}


void Singlet_Checker::ForcedDecays() {
  list<list<Singlet *>::iterator>::iterator bit=m_badones.begin();
  while (bit!=m_badones.end()) {
    p_singlet = (**bit);
    if (ForcedDecayOfTwoPartonSinglet()) {
      p_singlets->erase((*bit));
      bit = m_badones.erase(bit);
    }
    else {
      Flavour flav1 = (**bit)->front()->Flavour(); 
      Flavour flav2 = (**bit)->back()->Flavour();
      Flavour had   = p_softclusters->LowestTransition(flav1,flav2);
      bit++;
    }
  }
}

bool Singlet_Checker::ForcedDecayOfTwoPartonSinglet() {
  if (!ExtractAndCheckFlavours()) abort();
  if ((p_part1->Flavour().IsGluon() && p_part2->Flavour().IsGluon() && 
       TwoGluonSingletToHadrons()) ||
      (!(p_part1->Flavour().IsGluon() && p_part2->Flavour().IsGluon()) && 
       TwoQuarkSingletToHadrons())) {
    delete p_singlet;
    return true;
  }
  return false;
}

bool Singlet_Checker::ExtractAndCheckFlavours() {
  p_part1 = p_singlet->front();
  p_part2 = p_singlet->back();
  m_mass  = sqrt((p_part1->Momentum()+p_part2->Momentum()).Abs2());
  // check that both are gluons or both are not gluons.
  return ((p_part1->Flavour().IsGluon() && p_part2->Flavour().IsGluon()) ||
	  (!p_part1->Flavour().IsGluon() && !p_part2->Flavour().IsGluon()));
}

bool Singlet_Checker::TwoGluonSingletToHadrons() {
  // The gluon pair is too light to allow standard treatment.
  // If it is a bit too light, we make two quarks out of it.
  if (m_mass > 2.*m_minQmass) {
    if (m_splitter(p_part1,p_part2)) {
      Cluster * cluster = new Cluster(p_part1,p_part2);
      if ((!p_softclusters->Treat(cluster,true))==1) {
	msg_Error()<<"Error in "<<METHOD<<": transformed two gluons into\n"
		   <<(*cluster)
		   <<"but did not decay further.  Insert into cluster list.\n";
      }
      else delete cluster;
      return true;
    }
  }
  // If it is way to light, we make two hadrons/photons.
  Cluster * cluster = new Cluster(p_part1,p_part2);
  if (p_softclusters->TreatTwoGluons(cluster)) {
    delete cluster;
    return true;
  }
  msg_Error()<<"Error in "<<METHOD<<": could not decay two-gluon cluster\n"
	     <<(*cluster);
  return false;
}

bool Singlet_Checker::TwoQuarkSingletToHadrons() {
  // Regular two-hadron decay, if singlet mass larger than lowest decay mass,
  // else force a radiative decay.  Return true if either successful.
  Cluster * cluster = new Cluster(p_part1,p_part2);
  if ((m_mass > p_softclusters->MinDoubleMass(p_part1->Flavour(),
					      p_part2->Flavour()) &&
       p_softclusters->Treat(cluster,true)) ||
      p_softclusters->RadiativeDecay(cluster)) {
    delete cluster;
    return true;
  }
  return false;
}

