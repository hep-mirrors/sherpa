#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Gluon_Decayer::Gluon_Decayer(list<Cluster *> * cluster_list,
			     Soft_Cluster_Handler * softclusters) :
  Singlet_Tools(),
  p_cluster_list(cluster_list), p_softclusters(softclusters),
  m_splitter(Gluon_Splitter(cluster_list,softclusters))
{}

Gluon_Decayer::~Gluon_Decayer() {}

void Gluon_Decayer::Init() {
  Singlet_Tools::Init();
  m_splitter.Init();
  m_breaker.Init();
}

bool Gluon_Decayer::operator()(Singlet * singlet) {
  p_singlet = singlet;
  if (p_singlet->front()->Flavour().IsGluon() && !SplitGluonRing()) {
    msg_Error()<<"Couldn't split the gluon ring.\n"<<(*singlet)<<"\n";
    return false;
  }
  if (p_singlet->size()==2) {
    return Trivial(p_singlet->front(),p_singlet->back());
  }
  Proto_Particle * part1,* part2;
  size_t count(0);
  while (p_singlet->size()>3) {
    bool direction(ran->Get()>0.5);
    if ((*p_singlet->begin())->Flavour()==Flavour(kf_b) ||
	(*p_singlet->begin())->Flavour()==Flavour(kf_b).Bar() ||
	(*p_singlet->begin())->Flavour()==Flavour(kf_c) ||
	(*p_singlet->begin())->Flavour()==Flavour(kf_c).Bar())
      direction = true;
    else if ((*p_singlet->rbegin())->Flavour()==Flavour(kf_b) ||
	(*p_singlet->rbegin())->Flavour()==Flavour(kf_b).Bar() ||
	(*p_singlet->rbegin())->Flavour()==Flavour(kf_c) ||
	(*p_singlet->rbegin())->Flavour()==Flavour(kf_c).Bar())
      direction = false;
    if (direction) {
      list<Proto_Particle *>::iterator pliter=p_singlet->begin();
      part1 = (*pliter);
      pliter++;
      part2 = (*pliter);
    }
    else {
      list<Proto_Particle *>::reverse_iterator pliter=p_singlet->rbegin();
      part1 = (*pliter);
      pliter++;
      part2 = (*pliter);
    }
    switch (Step(part1,part2)) {
    case 1:
      if (direction) 
	p_singlet->pop_front();
      else 
	p_singlet->pop_back();
    case 0:
      break;
    case -1:
    default:
      msg_Error()<<METHOD<<" says: huh?\n";
      return false;
    }
  }
  if (LastStep()) {
    delete singlet;
    return true;
  }
  msg_Error()<<METHOD<<": LastStep failed.\n";
  return false;
}

bool Gluon_Decayer::SplitGluonRing() {
  // Reorder to make sure first & second gluon have highest combined mass
  p_singlet->Reorder(FirstGluon());
  return m_breaker(p_singlet);
}

Proto_Particle * Gluon_Decayer::FirstGluon() {
  double minm2(1.e12), m2thres(sqr(2.*m_breaker.MinMass())), m2;
  list<Proto_Particle *>::iterator ppiter1, ppiter2, winner;
  for (list<Proto_Particle *>::iterator ppiter1=p_singlet->begin();
       ppiter1!=p_singlet->end();ppiter1++) {
    Proto_Particle * part1(*ppiter1);
    ppiter2 = ppiter1;
    ppiter2++;
    if (ppiter2==p_singlet->end()) ppiter2=p_singlet->begin();
    Proto_Particle * part2(*ppiter2);
    m2 = (part1->Momentum()+part2->Momentum()).Abs2();
    if (m2<minm2 && m2>m2thres) {
      minm2  = m2;
      winner = ppiter1;
    }
  }
  return (*winner);
}

int Gluon_Decayer::Step(Proto_Particle * part1,Proto_Particle * part2,
			Proto_Particle * part3) {
  if (CheckMass(part1,part2) && m_splitter(part1,part2,part3)) return 1;
  return (p_singlet->Combine(part1,part2)?0:-1);
}

bool Gluon_Decayer::LastStep() {
  Proto_Particle * part[3];
  size_t i(0);
  for (list<Proto_Particle *>::iterator pliter=p_singlet->begin();
       pliter!=p_singlet->end();pliter++) part[i++] = (*pliter);
  size_t gluon(1),split(0),spect(2);
  if ((!part[0]->IsLeading() && part[2]->IsLeading()) ||
      (!part[0]->IsLeading() && !part[2]->IsLeading() &&
       (part[0]->Momentum()+part[1]->Momentum()).Abs2()>
       (part[2]->Momentum()+part[1]->Momentum()).Abs2())) {
    split = 2; spect = 0;
  }
  int stepres(Step(part[split],part[gluon],part[spect]));
  if (stepres==0) {
    return Trivial(p_singlet->front(),p_singlet->back());
  }
  if (split==0) p_singlet->pop_front();
           else p_singlet->pop_back();
  return Trivial(p_singlet->front(),p_singlet->back());
}

bool Gluon_Decayer::Trivial(Proto_Particle * part1,Proto_Particle * part2) {
  Cluster * cluster = new Cluster(part1,part2);
  p_singlet->pop_front();
  p_singlet->pop_back();
  if (p_softclusters->Treat(cluster)) delete cluster;
  else p_cluster_list->push_back(cluster);
  return true;
}
