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
  p_hadrons(softclusters->GetHadrons())
{}

Singlet_Checker::~Singlet_Checker() {}

void Singlet_Checker::Init() {
  Singlet_Tools::Init();
  m_splitter.Init();
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
      else if (!FusePartonsInLowMassSinglet()) return false;
    }
    // everything is fine -- the singlet can go on as normal.
    else lsit++;
  }
  // invoking the rescue system, if neccessary.  
  if (m_badones.size()>0) {
    //msg_Out()<<METHOD<<"  --> "<<m_badones.size()<<" problems.\n";
    //for (list<list<Singlet *>::iterator>::iterator sit=m_badones.begin();
    //	 sit!=m_badones.end();sit++)
    //  msg_Out()<<(***sit)<<"\n";
    //msg_Out()<<"---------------------------------------------------\n";
    if (!DealWithProblematicSinglets()) {
      msg_Error()<<METHOD<<" throw error - no rescue possible.\n";
      return false;
    }
  }
  return true;
}

bool Singlet_Checker::CheckSinglet() {
  // Checking the mass for pairs of colour-connected particles
  list<Proto_Particle *>::iterator plit1(p_singlet->begin()), plit2(plit1);
  plit2++;
  while (plit2!=p_singlet->end()) {
    p_part1 = (*plit1);
    p_part2 = (*plit2);
    if (!CheckMass(p_part1,p_part2)) return false;
    plit2++;
    plit1++;
  }
  // is gluon "ring" must also check for pair made of first and last particle
  if (m_isring) {
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
  if ((m_transitions.size()==1 && FindOtherSingletToTransit()) ||
      m_transitions.size()>1) {
    if (!TransitProblematicSinglets()) {
      msg_Error()<<METHOD<<" throws error for one transition.\n";
      return false;
    }
  }
  else if (m_transitions.size()==1 && FindRecoilerForTransit()) {
    if (!TransitProblematicSingletWithRecoiler()) {
      msg_Error()<<METHOD<<" throws error for one transition.\n";
      return false;
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
  double  massdiff(1.e-6);
  while (bit!=m_badones.end()) {
    p_singlet = (**bit);
    Flavour flav1 = p_singlet->front()->Flavour();
    Flavour flav2 = p_singlet->back()->Flavour();
    if (!flav1.IsGluon()) {
      Flavour hadtest = p_softclusters->LowestTransition(flav1,flav2);
      if (hadtest.Mass()-m_mass<massdiff) {
	hadron   = hadtest;
	hit      = bit;
	massdiff = hadtest.Mass()-m_mass;
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
  msg_Error()<<METHOD<<" throws error.\n";
  return false;
}

bool Singlet_Checker::FindRecoilerForTransit() {
  if (m_transitions.size()!=1 && m_badones.size()!=1) abort();
  m_singletmom = m_transitions.begin()->first->Momentum();
  m_targetmass = m_transitions.begin()->second.Mass();
  p_recoiler = NULL;
  for (list<Singlet *>::iterator sit=p_singlets->begin();
       sit!=p_singlets->end();sit++) {
    p_singlet = (*sit);
    if (TestSingletForRecoiler() && p_recoiler->IsBeam()) return true;
  }
  return (p_recoiler==NULL);
}

bool Singlet_Checker::TestSingletForRecoiler() {
  // Logic: prefer beam particles at ends of singlets and take them
  // as first choice recoilers for troublesome singlets.
  // If no beam particles in singlet, go for other particles.
  if (TestForRecoilerBeams()) return true;
  return TestForOtherRecoilers();
}

bool Singlet_Checker::TestForOtherRecoilers() {
  Proto_Particle * test = NULL;
  list<Proto_Particle *>::iterator pit1,pit2,hit1,hit2;
  pit1 = p_singlet->begin(); pit2=pit1; pit2++;
  hit1 = hit2 = p_singlet->end();
  double smin(0.),stest;
  int i(0);
  do {
    stest = (((*pit1)->Momentum()+(*pit2)->Momentum()).Abs2()-
	     sqr((*pit1)->Flavour().Mass()+(*pit2)->Flavour().Mass()));
    if (stest>smin && (TestRecoiler((*pit1)) || TestRecoiler((*pit2)))) {
      smin = stest;
      hit1 = pit1;
      hit2 = pit2;
    }
    pit1++;
    pit2++;
  } while (pit2!=p_singlet->end());
  if (hit1!=p_singlet->end()) {
    if ((*hit1)->Flavour().IsGluon()) {
      if ((*hit2)->Flavour().IsGluon())
	test = (((*hit1)->Momentum()+m_singletmom).Abs2()>
		((*hit2)->Momentum()+m_singletmom).Abs2()) ? (*hit1):(*hit2);
      else test = (*hit2);
    }
    else if ((*hit2)->Flavour().IsGluon())
      test = (*hit1);
    else if (p_recoiler==NULL || !p_recoiler->Flavour().IsGluon()) {
      test = (((*hit1)->Momentum()+m_singletmom).Abs2()>
	      ((*hit2)->Momentum()+m_singletmom).Abs2()) ? (*hit1):(*hit2);
    }
  }
  if (test!=NULL) {
    if (p_recoiler==NULL ||
	(test->Momentum()+m_singletmom).Abs2()>
	(p_recoiler->Momentum()+m_singletmom).Abs2()) {
      p_recoiler = test;
      return true;
    }
  }
  return false;
}

bool Singlet_Checker::TestForRecoilerBeams() {
  if (p_singlet->front()->IsBeam() &&
      TestRecoiler(p_singlet->front())) {
    p_recoiler = p_singlet->front();
    return true;
  }
  if (p_singlet->back()->IsBeam() &&
      TestRecoiler(p_singlet->back())) {
    p_recoiler = p_singlet->back();
    return true;
  }
  return false;
}

bool Singlet_Checker::TestRecoiler(Proto_Particle * part) {
  return ((m_singletmom+part->Momentum()).Abs2() >
	  sqr(m_targetmass+part->Flavour().Mass()));
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
  if (totmom.Abs2()<sqr(totmass)) return false;
  if (hadpars->AdjustMomenta(n,moms,masses)) {
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
    return true;
  }
  return false;
}

bool Singlet_Checker::TransitProblematicSingletWithRecoiler() {
  Vec4D *  moms   = new Vec4D[2];
  double * masses = new double[2];
  p_singlet       = m_transitions.begin()->first;
  Flavour hadron  = m_transitions.begin()->second;
  moms[0]   = p_singlet->Momentum();
  moms[1]   = p_recoiler->Momentum();
  masses[0] = hadron.Mass();
  masses[1] = p_recoiler->Flavour().Mass();
  if (hadpars->AdjustMomenta(2,moms,masses)) {
    bool isbeam = (p_singlet->front()->IsBeam() ||
		   p_singlet->back()->IsBeam());
    Proto_Particle * part = new Proto_Particle(hadron,moms[0],false,isbeam);
    p_hadrons->push_back(part);
    p_recoiler->SetMomentum(moms[1]);
    delete p_singlet;
    m_transitions.clear();
    return true;
  }
  return false;
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
      double  hmass = p_softclusters->MinSingleMass(flav1,flav2);
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
  if (m_mass > 2.*m_minQmass && m_splitter(p_part1,p_part2))
    return true;
  // If it is way to light, we make two hadrons/photons.
  Cluster * cluster = new Cluster(p_part1,p_part2);
  return p_softclusters->TreatTwoGluons(cluster);
}

bool Singlet_Checker::TwoQuarkSingletToHadrons() {
  // Regular two-hadron decay, if singlet mass larger than lowest decay mass,
  // else force a radiative decay.  Return true if either successful.
  Cluster cluster(p_part1,p_part2);
  return ((m_mass > p_softclusters->MinDoubleMass(p_part1->Flavour(),
						  p_part2->Flavour()) &&
	   p_softclusters->Treat(&cluster,true)) ||
	  p_softclusters->RadiativeDecay(&cluster));
}

