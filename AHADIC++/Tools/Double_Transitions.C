#include "AHADIC++/Tools/Double_Transitions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Double_Transitions::Double_Transitions(Single_Transitions * singles) :
  m_wtthres(1.e-6)
{
  FillMap(singles);
  Normalise();
}

void Double_Transitions::FillMap(Single_Transitions * singletransitions)
{
  Constituents * constituents     = hadpars->GetConstituents();
  Single_Transition_Map * singles = singletransitions->GetMap();
  for (Single_Transition_Map::iterator stmit1=singles->begin();
       stmit1!=singles->end();stmit1++) {
    Flavour_Pair pair1   = stmit1->first;
    Flavour popped       = pair1.second;
    for (Single_Transition_Map::iterator stmit2=singles->begin();
	 stmit2!=singles->end();stmit2++) {
      Flavour_Pair pair2 = stmit2->first;
      if (pair2.first!=popped.Bar()) continue;
      Flavour_Pair pair;
      pair.first         = pair1.first;
      pair.second        = pair2.second; 
      double weight      = constituents->TotWeight(popped.Bar()); 
      if (weight<1.e-6) continue;
      if (m_transitions.find(pair)==m_transitions.end())
	m_transitions[pair] = new Double_Transition_List;
      for (Single_Transition_List::iterator hit1=stmit1->second->begin();
	   hit1!=stmit1->second->end();hit1++) {
	for (Single_Transition_List::iterator hit2=stmit2->second->begin();
	     hit2!=stmit2->second->end();hit2++) {
	  Flavour_Pair hads;
	  hads.first  = hit1->first;
	  hads.second = hit2->first;
	  double wt   = weight*hit1->second*hit2->second;
	  if (wt<m_wtthres) continue;
	  (*m_transitions[pair])[hads] = wt;
	}
      }
    }
  }
}

void Double_Transitions::Normalise() {
  for (Double_Transition_Map::iterator dtmit=m_transitions.begin();
       dtmit!=m_transitions.end();dtmit++) {
    double totweight = 0.;
    for (Double_Transition_List::iterator dtlit=dtmit->second->begin();
	 dtlit!=dtmit->second->end();dtlit++)
      totweight += dtlit->second;
    for (Double_Transition_List::iterator dtlit=dtmit->second->begin();
	 dtlit!=dtmit->second->end();dtlit++) 
      dtlit->second /= totweight;
  }
}

void Double_Transitions::Print(const bool & full) {
  map<Flavour,double> checkit;
  msg_Out()<<"---------------------------------------------------------\n"
	   <<METHOD<<":\n";
  for (Double_Transition_Map::iterator dtmit=m_transitions.begin();
       dtmit!=m_transitions.end();dtmit++) {
    if (full)
      msg_Out()<<"---------------------------------------------------------\n"
	       <<"*** ["<<dtmit->first.first<<" "<<dtmit->first.second<<"]:\n";
    for (Double_Transition_List::iterator dtlit=dtmit->second->begin();
	 dtlit!=dtmit->second->end();dtlit++) {
      if (checkit.find(dtlit->first.first)==checkit.end())
	checkit[dtlit->first.first] = 0.;
      if (checkit.find(dtlit->first.second)==checkit.end())
	checkit[dtlit->first.second] = 0.;
      checkit[dtlit->first.first]  += dtlit->second;
      checkit[dtlit->first.second] += dtlit->second;
      if (full)
	msg_Out()<<"  -> ["<<dtlit->first.first<<" + "
		 <<dtlit->first.second<<"], "
		 <<"weight = "<<dtlit->second<<"\n";
    }
  }
  msg_Out()<<"---------------------------------------------------------\n"
	   <<"---------------------------------------------------------\n";
  for (map<Flavour,double>::iterator cit=checkit.begin();
       cit!=checkit.end();cit++) {
    msg_Out()<<" --> "<<cit->first<<" with total = "<<cit->second<<".\n";
  }
}
      
Double_Transitions::~Double_Transitions() {
  while (!m_transitions.empty()) {
    delete (m_transitions.begin()->second); 
    m_transitions.erase(m_transitions.begin());
  }
  m_transitions.clear();
}

Double_Transition_List *
Double_Transitions::operator[](const Flavour_Pair & flavs) {
  if (m_transitions.find(flavs)==m_transitions.end()) {
    msg_Error()<<"Error in "<<METHOD<<"["<<m_transitions.size()<<"] for "
	       <<"["<<flavs.first<<", "<<flavs.second<<"]:\n"
	       <<"   Illegal flavour combination, will return 0.\n";
    return 0;
  }
  return m_transitions.find(flavs)->second;
}

Flavour_Pair Double_Transitions::
GetLightestTransition(const Flavour_Pair & fpair) {
  Flavour_Pair pair;
  pair.first = pair.second = Flavour(kf_none);
  Double_Transition_Map::iterator dtiter = m_transitions.find(fpair);
  if (dtiter==m_transitions.end())  return pair;
  Double_Transition_List * dtl  = dtiter->second;
  if (dtl->empty()) return pair;
  return (--dtl->end())->first;
}

Flavour_Pair Double_Transitions::
GetHeaviestTransition(const Flavour_Pair & fpair) {
  Flavour_Pair pair;
  pair.first = pair.second = Flavour(kf_none);
  Double_Transition_Map::iterator dtiter = m_transitions.find(fpair);
  if (dtiter!=m_transitions.end()) pair = dtiter->second->begin()->first;
  return pair;
}

double Double_Transitions::GetLightestMass(const Flavour_Pair & fpair) {
  Flavour_Pair pair = GetLightestTransition(fpair);
  if (pair.first==Flavour(kf_none) || pair.second==Flavour(kf_none)) return -1.;
  return pair.first.HadMass()+pair.second.HadMass();
}

double Double_Transitions::GetHeaviestMass(const Flavour_Pair & fpair) {
  Flavour_Pair pair = GetHeaviestTransition(fpair);
  if (pair.first==Flavour(kf_none) || pair.second==Flavour(kf_none)) return -1.;
  return pair.first.HadMass()+pair.second.HadMass();
}

