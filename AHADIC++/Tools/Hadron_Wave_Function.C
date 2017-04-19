#include "AHADIC++/Tools/Hadron_Wave_Function.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Hadron_Wave_Function::Hadron_Wave_Function() :
  m_hadron(Flavour(kf_none)), m_spin(-1), m_kfcode(0), m_mpletwt(1.), m_barrable(false)
{ }

Hadron_Wave_Function::Hadron_Wave_Function(const ATOOLS::Flavour & _hadron) :
  m_hadron(_hadron), m_spin(_hadron.IntSpin()), 
  m_kfcode(_hadron.Kfcode()), m_mpletwt(1.), m_barrable(false)
{ }

Hadron_Wave_Function::~Hadron_Wave_Function() 
{
  for (WFcompiter wf=m_waves.begin();wf!=m_waves.end();wf++) {
    if (wf->first) delete wf->first; 
  }
  m_waves.clear();
}

void Hadron_Wave_Function::AddToWaves(Flavour_Pair * pair,double weight)
{
  if (m_waves.find(pair)==m_waves.end()) m_waves[pair] = weight;
  else {
    msg_Error()<<"Potential error in Hadron_Wave_Function::AddToWaves"<<endl
	       <<"   Pair "<<pair->first<<"/"<<pair->second<<" already in map."<<endl
	       <<"   Will ignore this and continue."<<endl;
    return;
  }
  if (pair->first!=pair->second.Bar()) m_barrable = true;
}

Hadron_Wave_Function * Hadron_Wave_Function::Anti() {
  Hadron_Wave_Function * wf = NULL;
  if (m_barrable) {
    wf = new Hadron_Wave_Function(m_hadron.Bar());
    wf->SetSpin(m_spin);
    wf->SetKfCode(-m_kfcode);
    Flavour_Pair * pair;
    for (WFcompiter wfc=m_waves.begin();wfc!=m_waves.end();wfc++) {
      pair         = new Flavour_Pair;
      pair->first  = wfc->first->second.Bar();
      pair->second = wfc->first->first.Bar();
      wf->AddToWaves(pair,wfc->second);
    }
  }
  return wf;
}

double Hadron_Wave_Function::WaveWeight(ATOOLS::Flavour first,ATOOLS::Flavour second) 
{
  Flavour_Pair * fpair;
  for (WFcompiter wit=m_waves.begin();wit!=m_waves.end();wit++) {
    fpair = wit->first;
    if ((fpair->first==first && fpair->second==second) ||
	(fpair->first==second && fpair->second==first)) return wit->second;
  }
  return 0.;
} 

namespace AHADIC {
ostream & operator<<(ostream & s, Hadron_Wave_Function & wf) 
{
  WFcomponent * waves = wf.GetWaves();
  double wf2(0.);
  for (WFcompiter wfc=waves->begin();wfc!=waves->end();wfc++) 
    wf2 += wfc->second*wfc->second;
  s<<" "<<wf.m_hadron<<" ("<<wf.m_kfcode<<"), spin = "<<(wf.m_spin-1)/2.
   <<", weight = "<<wf2<<"."<<endl;
  for (WFcompiter wfc=waves->begin();wfc!=waves->end();wfc++) {
    s<<"     "<<wfc->first->first<<" "<<wfc->first->second
     <<" : "<<wfc->second<<" ---> "<<(1./(wfc->second*wfc->second))<<endl;
  }
  return s;
}
}

