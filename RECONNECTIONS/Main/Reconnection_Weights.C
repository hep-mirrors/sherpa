#include "RECONNECTIONS/Main/Reconnection_Weights.H"
#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Weights::Reconnection_Weights(Reconnection_Handler * rhandler) :
  p_rhandler(rhandler)
{
  SetLists();
}

void Reconnection_Weights::SetLists() {
  for (size_t pos=0;pos<2;pos++) p_parts[pos] = p_rhandler->GetParts(pos);
  p_singlets = p_rhandler->GetSinglets();
}

void Reconnection_Weights::Initialize(Default_Reader *const defaultreader) {
  // Pmode is the mode for the distance measure in momentum space.
  // 0 - mode is "linear":    dist = log(1+sij/Q0^2)
  // 1 - mode is "power law": dist = exp[eta * log(1+sij/Q0^2) ] 
  m_Pmode     = defaultreader->GetValue<double>("RECONNECTIONS::PMODE",0);
  m_Q02       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::Q_0",0.25));
  m_eta       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::eta",0.16));
  m_R02       = sqr(defaultreader->GetValue<double>("RECONNECTIONS::R_0",1.));
  m_reshuffle = 1./defaultreader->GetValue<double>("RECONNECTIONS::RESHUFFLE",1./3.);
  m_restring  = 1./defaultreader->GetValue<double>("RECONNECTIONS::RESTRING",1./3.);
}

void Reconnection_Weights::Reset() {
  while (!m_distances.empty()){
    delete (m_distances.begin()->second);
    m_distances.erase(m_distances.begin());
  }
}

void Reconnection_Weights::FillTables() {
  Particle * trip, * anti;
  for (ParticleSet::iterator tit=p_parts[0]->begin();tit!=p_parts[0]->end();tit++) {
    trip = (*tit);
    m_distances[trip] = new map<Particle *, double>;
    for (ParticleSet::iterator ait=p_parts[1]->begin();ait!=p_parts[1]->end();ait++) {
      anti = (*ait);
      if (trip==anti) continue;
      double pdist = MomDistance(trip,anti);
      double xdist = PosDistance(trip,anti);
      double cdist = ColDistance(trip,anti);
      double dist  = pdist * xdist * cdist;
      (*m_distances[trip])[anti] = dist;
    }
  }
  //OutputWeightTable();
}

double Reconnection_Weights::MomDistance(Particle * part1,Particle * part2) {
  // Here we take a variant of the Lund lambda measure for the distance in momentum space
  double p1p2 = ((part1->Momentum()+part2->Momentum()).Abs2() -
		 (part1->Momentum().Abs2()+part2->Momentum().Abs2()));
  return exp(0.16*log(1.+p1p2/m_Q02));
}

double Reconnection_Weights::PosDistance(Particle * part1,Particle * part2) {
  double xdist2 = dabs((part1->XProd().Perp()-part2->XProd().Perp()).Abs2());
  return exp(sqrt(xdist2/m_R02));
}

double Reconnection_Weights::ColDistance(Particle * part1,Particle * part2) {
  // For colour connected partons there is no extra colour suppression, which would increase
  // the distance in our combined space.
  if ((part1->GetFlow(1)==part2->GetFlow(2) && part1->GetFlow(1)!=0) ||
      (part1->GetFlow(2)==part2->GetFlow(1) && part1->GetFlow(2)!=0)) return 1.;
  // We find out if the two partons are in the smae singlet or not
  Part_List * singlet1(NULL), * singlet2(NULL);
  Part_Iterator pit1, pit2;
  for (list<Part_List *>::iterator sit=p_singlets->begin();sit!=p_singlets->end();sit++) {
    pit1=find((*sit)->begin(),(*sit)->end(),part1);
    if (pit1!=(*sit)->end()) singlet1 = (*sit);
    pit2=find((*sit)->begin(),(*sit)->end(),part2);
    if (pit2!=(*sit)->end()) singlet2 = (*sit);
    if (singlet1!=NULL && singlet2!=0) break;
  }
  // If they are not in the same singlet, the distance in colour space is given by the
  // "flat" re-stringing probability m_restring;
  if (singlet1!=singlet2) return m_restring;
  // If they are in the same singlet, we first figure out how many steps in colour space
  // the two partons are apart from each other.  The reshuffling probability is given by
  // a weight m_reshuffle ^ (steps-1)
  // TODO: Deal with gluon rings here.
  int steps=0;
  pit2=pit1;
  while (pit2!=singlet1->end() && (*pit2)!=part2) {
    pit2++;
    steps++;
  }
  int minsteps=steps;
  steps = 0;
  pit2=pit1;
  while (pit2!=singlet1->end() && (*pit2)!=part2) {
    pit2--;
    steps++;
  }
  if (steps!=0 && steps<minsteps) minsteps = steps;
  return pow(m_reshuffle,Max(2,minsteps-1));
}

void Reconnection_Weights::OutputWeightTable() {
  for (map<Particle *,distances * >::iterator dit=m_distances.begin();
       dit!=m_distances.end();dit++) {
    msg_Out()<<"Distances for particle ["<<dit->first->Number()<<"]"
    	     <<"("<<dit->first->GetFlow(1)<<", "
    	     <<dit->first->GetFlow(2)<<"):\n";
    distances * dists = dit->second;
    for (distances::iterator distit=dists->begin();
	 distit!=dists->end();distit++) {
      msg_Out()<<"   ["<<distit->first->Number()<<"]"
	       <<"("<<distit->first->GetFlow(1)<<", "
	       <<distit->first->GetFlow(2)<<") = "<<distit->second<<"\n";
    }
  }
}

