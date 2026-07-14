#include "AHADIC++/Tools/Single_Transitions.H"
#include "AHADIC++/Tools/Hadronisation_Reweighting.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Single_Transitions::Single_Transitions(Wave_Functions * wavefunctions,
					Hadronisation_Reweighting * reweighting)
{
  m_n_variations = reweighting->NumberOfVariations();
  FillMap(wavefunctions, reweighting);
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

void Single_Transitions::FillMap(Wave_Functions * wavefunctions,
				 Hadronisation_Reweighting * reweighting) {
  // Go through the wavefunctions of all hadrons, extract
  // their components (Flavour_Pairs) and make a list of
  // all transitions for a single pair, consisting of the
  // hadron and the transition probability.
  for (Wave_Functions::iterator wfit=wavefunctions->begin();
       wfit!=wavefunctions->end();wfit++) {
    Flavour hadron = wfit->first;
    const std::vector<double> & extrawts = wfit->second->ExtraWeights();
    const std::vector<double> & mpletwts = wfit->second->MultipletWeights();
    std::vector<double> weight(m_n_variations);
    for (size_t ivar=0;ivar<m_n_variations;ivar++) {
      weight[ivar] = (mpletwts[ivar] *
		   wfit->second->SpinWeight() *
		   extrawts[ivar]);
    }
    reweighting->CheckTransitionGuard(hadron, weight, 1.e-6);
    if (weight[0]<1.e-6) continue;
    WaveComponents * singlewaves = wfit->second->GetWaves();
    for (WaveComponents::iterator cit=singlewaves->begin();
	 cit!=singlewaves->end();cit++) {
      Flavour_Pair pair = (*cit->first);
      const std::vector<double> & amps = cit->second;
      std::vector<double> wt(m_n_variations);
      for (size_t ivar=0;ivar<m_n_variations;ivar++)
	wt[ivar] = weight[ivar] * sqr(amps[amps.size()>1 ? ivar : 0]);
      if (m_transitions.find(pair)==m_transitions.end()) {
	m_transitions[pair] = new Single_Transition_List;
      }
      (*m_transitions[pair])[hadron] = wt;
    }
  }
}

void Single_Transitions::Normalise() {
  for (Single_Transition_Map::iterator stmit=m_transitions.begin();
       stmit!=m_transitions.end();stmit++) {
    std::vector<double> totwt(m_n_variations,0.);
    for (Single_Transition_List::iterator stlit=stmit->second->begin();
	 stlit!=stmit->second->end();stlit++) {
      for (size_t ivar=0;ivar<m_n_variations;ivar++) totwt[ivar] += stlit->second[ivar];
    }
    for (Single_Transition_List::iterator stlit=stmit->second->begin();
	 stlit!=stmit->second->end();stlit++) {
      for (size_t ivar=0;ivar<m_n_variations;ivar++) stlit->second[ivar] /= totwt[ivar];
    }
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
      msg_Out()<<"   "<<stlit->first<<" --> "<<stlit->second[0]<<"\n";
      totwt += stlit->second[0];
      if (checkit.find(stlit->first)==checkit.end()) checkit[stlit->first] = 0.;
      checkit[stlit->first] += stlit->second[0];
    }
    msg_Out()<<"   Total weight = "<<totwt<<"\n\n";
  }
  msg_Out()<<"-------------------------------------------------------------\n";
}
