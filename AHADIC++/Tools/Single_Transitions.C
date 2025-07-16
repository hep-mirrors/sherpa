#include "AHADIC++/Tools/Single_Transitions.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Single_Transitions::Single_Transitions(Wave_Functions * wavefunctions)
{
  FillMap(wavefunctions);
  Normalise();
}

Single_Transitions::~Single_Transitions()
{
  for (Single_Transition_Map::iterator stiter=m_transitions.begin();
       stiter!=m_transitions.end();stiter++) {
    delete stiter->second;
  }
  m_transitions.clear();
}

void Single_Transitions::FillMap(Wave_Functions * wavefunctions) {
  // Go through the wavefunctions of all hadrons, extract
  // their components (Flavour_Pairs) and make a list of
  // all transitions for a single pair, consisting of the
  // hadron and the transition probability.
  for (Wave_Functions::iterator wfit=wavefunctions->begin();
       wfit!=wavefunctions->end();wfit++) {
    Flavour hadron = wfit->first;
    double  weight = (wfit->second->MultipletWeight() *
		      wfit->second->SpinWeight() *
		      wfit->second->ExtraWeight());
    if (weight<1.e-6) continue;
    WaveComponents * singlewaves = wfit->second->GetWaves();
    for (WaveComponents::iterator cit=singlewaves->begin();
	 cit!=singlewaves->end();cit++) {
      Flavour_Pair pair = (*cit->first);
      double wt = weight * sqr(cit->second);
      if (m_transitions.find(pair)==m_transitions.end()) {
	m_transitions[pair] = new Single_Transition_List;
      }
      (*m_transitions[pair])[hadron] = wt;
    }
  }
}

void Single_Transitions::Normalise() {
  double totwt;
  for (Single_Transition_Map::iterator stmit=m_transitions.begin();
       stmit!=m_transitions.end();stmit++) {
    totwt = 0.;
    for (Single_Transition_List::iterator stlit=stmit->second->begin();
	 stlit!=stmit->second->end();stlit++) totwt += stlit->second;
    for (Single_Transition_List::iterator stlit=stmit->second->begin();
	 stlit!=stmit->second->end();stlit++) stlit->second/=totwt;
  }
}

Single_Transition_List *
Single_Transitions::operator[](const Flavour_Pair & flavs) {
  if (m_transitions.find(flavs)==m_transitions.end()) {
    msg_Error()<<"Error in "<<METHOD<<" for "
	       <<"["<<flavs.first<<", "<<flavs.second<<"]:\n"
	       <<"   Illegal flavour combination, will return 0.\n";
    return 0;
  }
  return m_transitions.find(flavs)->second;
}

Flavour Single_Transitions::GetLightestTransition(const Flavour_Pair & fpair) {
  Single_Transition_Map::iterator stmit = m_transitions.find(fpair);
  if (stmit!=m_transitions.end()) return stmit->second->rbegin()->first;
  return Flavour(kf_none);
}

Flavour Single_Transitions::GetHeaviestTransition(const Flavour_Pair & fpair) {
  Single_Transition_Map::iterator stmit = m_transitions.find(fpair);
  if (stmit!=m_transitions.end()) return stmit->second->begin()->first;
  return Flavour(kf_none);
}

double Single_Transitions::GetLightestMass(const Flavour_Pair & fpair) {
  Flavour had = GetLightestTransition(fpair);
  if (had==Flavour(kf_none)) {
    return -(hadpars->GetConstituents()->Mass(fpair.first)+
	     hadpars->GetConstituents()->Mass(fpair.second));
  }
  return had.HadMass();
}

double Single_Transitions::GetHeaviestMass(const Flavour_Pair & fpair) {
  Flavour had = GetHeaviestTransition(fpair);
  if (had==Flavour(kf_none)) return -1.;
  return had.HadMass();
}

void Single_Transitions::Print() 
{
  double totwt;
  map<Flavour,double> checkit;
  for (Single_Transition_Map::iterator stmit=m_transitions.begin();
       stmit!=m_transitions.end();stmit++) {
    totwt = 0.;
    msg_Out()<<"----- "<<stmit->first.first<<" "<<stmit->first.second
	     <<" --------------------------\n";
    for (Single_Transition_List::iterator stlit=stmit->second->begin();
	 stlit!=stmit->second->end();stlit++) {
      msg_Out()<<"   "<<stlit->first<<" --> "<<stlit->second<<"\n";
      totwt += stlit->second;
      if (checkit.find(stlit->first)==checkit.end()) checkit[stlit->first] = 0.;
      checkit[stlit->first] += stlit->second;
    }
    msg_Out()<<"   Total weight = "<<totwt<<"\n\n";
  }
  msg_Out()<<"-------------------------------------------------------------\n";
}
