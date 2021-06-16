#include "RECONNECTIONS/Main/Reconnection_Weights.H"
#include "RECONNECTIONS/Main/Reconnection_Base.H"
#include "ATOOLS/Org/Message.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Weights::Reconnection_Weights(Reconnection_Base * reconnector) :
  p_reconnector(reconnector)
{
  SetLists();
}

void Reconnection_Weights::SetLists() {
  for (size_t pos=0;pos<2;pos++) p_parts[pos] = p_reconnector->GetParts(pos);
}

void Reconnection_Weights::Reset() {
  while (!m_distances.empty()){
    delete (m_distances.begin()->second);
    m_distances.erase(m_distances.begin());
  }
}

void Reconnection_Weights::FillTables() {
  Particle * trip, * anti;
  for (ParticleSet::iterator tit=p_parts[0]->begin();
       tit!=p_parts[0]->end();tit++) {
    trip = (*tit);
    m_distances[trip] = new map<Particle *, double>;
    for (ParticleSet::iterator ait=p_parts[1]->begin();
	 ait!=p_parts[1]->end();ait++) {
      anti = (*ait);
      if (trip==anti) continue;
      (*m_distances[trip])[anti] = p_reconnector->Distance(trip,anti);
    }
  }
  //OutputWeightTable();
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

