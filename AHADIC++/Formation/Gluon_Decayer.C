#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <cassert>

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Gluon_Decayer::Gluon_Decayer(list<Cluster *> * cluster_list,
			     Soft_Cluster_Handler * softclusters) :
  Singlet_Tools(),
  p_cluster_list(cluster_list), p_softclusters(softclusters),
  m_splitter(Gluon_Splitter(cluster_list,softclusters)),
  m_analyse(true)
{
  if (m_analyse) {
    m_histos[string("N_primaries")] = new Histogram(0,0.,100.,100);
    m_histos[string("N_clusters")]  = new Histogram(0,0.,100.,100);
    m_histos[string("M_all")]       = new Histogram(0,0.,100.,1000);
    m_histos[string("M_light")]     = new Histogram(0,0.,100.,1000);
    m_histos[string("M_c")]         = new Histogram(0,0.,100.,1000);
    m_histos[string("M_b")]         = new Histogram(0,0.,100.,1000);
    m_histos[string("Y_asym_1")]    = new Histogram(0,0.,8.,32);
  }
}

Gluon_Decayer::~Gluon_Decayer() {
  if (m_analyse) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      histo = hit->second;
      name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histos.clear();
  }
}

void Gluon_Decayer::Init() {
  Singlet_Tools::Init();
  m_splitter.Init();
  m_breaker.Init();
}

void Gluon_Decayer::Reset() {}

bool Gluon_Decayer::operator()(Singlet * singlet) {
  p_singlet = singlet;
  if (p_singlet->front()->Flavour().IsGluon() && !SplitGluonRing()) {
    // protection for low-mass two-gluon systems.
    if (p_singlet->size()==2) {
      if (Trivial(p_singlet->front(),p_singlet->back())) return true;
    }
    msg_Error()<<"Couldn't split the gluon ring.\n"<<(*singlet)<<"\n";
    return false;
  }
  if (p_singlet->size()==2) {
    bool flag = Trivial(p_singlet->front(),p_singlet->back(),false);
    if (!flag) {
      msg_Error()<<(*singlet)<<"\n";
      THROW(fatal_error,"Couldn't deal with 2-parton singlet.");
    }
    delete p_singlet;
    return flag;
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
      msg_Out()<<METHOD<<" failed at intermediate step:\n"
	       <<(*part1)<<"\n"
	       <<(*part2)<<"\n";
      return false;
    }
  }
  if (LastStep()) {
    delete singlet;
    return true;
  }
  return false;
}

bool Gluon_Decayer::SplitGluonRing() {
  // Reorder to make sure first & second gluon have highest combined mass
  p_singlet->Reorder(FirstGluon());
  return m_breaker(p_singlet);
}

Proto_Particle * Gluon_Decayer::FirstGluon() {
  double minm2(1.e12), m2thres(sqr(2.*m_breaker.MinMass())), m2;
  list<Proto_Particle *>::iterator ppiter1, ppiter2, winner(p_singlet->end());
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
  return (winner==p_singlet->end()?(*p_singlet->begin()):(*winner));
}

int Gluon_Decayer::Step(Proto_Particle * part1,Proto_Particle * part2,
			Proto_Particle * part3) {
  assert(part1 != nullptr);
  assert(part2 != nullptr);
  Vec4D momsave1(part1->Momentum()), momsave2(part2->Momentum());
  if (CheckMass(part1,part2) && m_splitter(part1,part2,part3)) {
    if (m_analyse) {
      m_Nclusters++;
      double mass = 0;
      bool isB = false, isC = false;
      m_splitter.GetLast(mass,isB,isC);
      m_histos[string("M_all")]->Insert(mass);
      if (isB)      m_histos[string("M_b")]->Insert(mass);
      else if (isC) m_histos[string("M_c")]->Insert(mass);
      else          m_histos[string("M_light")]->Insert(mass);
    }
    return 1;
  }
  part1->SetMomentum(momsave1);
  part2->SetMomentum(momsave2);
  return (p_singlet->Combine(part1,part2)?0:-1);
}

bool Gluon_Decayer::LastStep() {
  // collect particles
  Proto_Particle* part[3] {nullptr, nullptr, nullptr};
  {
    size_t i {0};
    for (auto& p : *p_singlet) {
      part[i] = p;
      ++i;
    }
  }
  // assign roles
  size_t gluon(1),split(0),spect(2);
  if ( (!part[0]->IsLeading() && part[2]->IsLeading()) ||
       (!part[0]->IsLeading() && !part[2]->IsLeading() &&
	(((part[0]->Momentum()+part[1]->Momentum()).Abs2()-
	  sqr(p_constituents->Mass(part[0]->Flavour())))    <
	 ((part[2]->Momentum()+part[1]->Momentum()).Abs2()-
	  sqr(p_constituents->Mass(part[2]->Flavour())))))) {
    split = 2; spect = 0;
  }
  // perform step and return result
  int stepres = Step(part[split],part[gluon],part[spect]);
  if (stepres==0) {
    return Trivial(p_singlet->front(),p_singlet->back());
  }
  if (split==0) p_singlet->pop_front();
           else p_singlet->pop_back();
  return Trivial(p_singlet->front(),p_singlet->back(),false);
}

bool Gluon_Decayer::Trivial(Proto_Particle * part1,Proto_Particle * part2,
			    const bool & force) {
  Cluster * cluster = new Cluster(part1,part2);
  if (m_analyse) {
    m_Nclusters++;
    double mass = sqrt(dabs(cluster->Momentum().Abs2()));
    bool isB    = (part1->Flavour()==Flavour(kf_b) ||
		   part1->Flavour()==Flavour(kf_b).Bar() ||
		   part2->Flavour()==Flavour(kf_b) ||
		   part2->Flavour()==Flavour(kf_b).Bar());
    bool isC    = (!isB &&
		   (part1->Flavour()==Flavour(kf_c) ||
		    part1->Flavour()==Flavour(kf_c).Bar() ||
		    part2->Flavour()==Flavour(kf_c) ||
		    part2->Flavour()==Flavour(kf_c).Bar()));
    m_histos[string("M_all")]->Insert(mass);
    if (isB)      m_histos[string("M_b")]->Insert(mass);
    else if (isC) m_histos[string("M_c")]->Insert(mass);
    else          m_histos[string("M_light")]->Insert(mass);
    double y = cluster->Momentum().Y();
    m_histos[string("Y_asym_1")]->Insert(dabs(y),(y>0.?1.:-1.));
  }
  p_singlet->pop_front();
  p_singlet->pop_back();
  switch (p_softclusters->Treat(cluster,force)) {
  case 1:
    delete cluster;
    break;
  case -1:
    if (p_softclusters->Rescue(cluster)) {
      delete cluster;
      break;
    }
    delete cluster;
    return false;
  default:
    p_cluster_list->push_back(cluster);
    break;
  }
  return true;
}

void Gluon_Decayer::FillNs(const int & Nhad) {
  if (m_analyse) {
    m_histos[string("N_clusters")]->Insert(m_Nclusters+0.5);
    m_histos[string("N_primaries")]->Insert(Nhad+0.5);
  }
}

void Gluon_Decayer::AnalyseClusters() {}

