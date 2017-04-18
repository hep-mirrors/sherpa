#include "AHADIC++/Tools/Wave_Function.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Wave_Function::Wave_Function() :
  m_hadron(Flavour(kf_none)),
  m_kfcode(0), m_spin2(0), m_mpletwt(1.), m_extrawt(1.),
  m_barrable(false)
{ }

Wave_Function::Wave_Function(const ATOOLS::Flavour & _hadron) :
  m_hadron(_hadron), m_kfcode(_hadron.Kfcode()),
  m_spin2(0), m_mpletwt(1.), m_extrawt(1.),
  m_barrable(false)
{ }

Wave_Function::~Wave_Function() 
{
  for (WaveComponents::iterator wf=m_waves.begin();wf!=m_waves.end();wf++) {
    if (wf->first) delete wf->first; 
  }
  m_waves.clear();
}

void Wave_Function::AddToWaves(Flavour_Pair * pair,double weight)
{
  if (m_waves.find(pair)==m_waves.end()) m_waves[pair] = weight;
  else {
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   "<<pair->first<<"/"<<pair->second<<" already in map.\n"
	       <<"   Will ignore this and continue."<<endl;
    exit(1);
  }
  if (pair->first!=pair->second.Bar()) m_barrable = true;
}

Wave_Function * Wave_Function::GetAnti() {
  Wave_Function * wf = new Wave_Function(m_hadron.Bar());
  wf->SetSpin(m_spin2);
  wf->SetKfCode(-m_kfcode);
  wf->SetMultipletWeight(m_mpletwt);
  wf->SetExtraWeight(m_extrawt);
  Flavour_Pair * pair;
  for (WaveComponents::iterator wfc=m_waves.begin();wfc!=m_waves.end();wfc++) {
    pair         = new Flavour_Pair;
    pair->first  = wfc->first->second.Bar();
    pair->second = wfc->first->first.Bar();
    wf->AddToWaves(pair,wfc->second);
  }
  return wf;
}

double Wave_Function::WaveWeight(ATOOLS::Flavour first,ATOOLS::Flavour second) 
{
  Flavour_Pair * fpair;
  for (WaveComponents::iterator wit=m_waves.begin();wit!=m_waves.end();wit++) {
    fpair = wit->first;
    if ((fpair->first==first && fpair->second==second) ||
	(fpair->first==second && fpair->second==first)) return wit->second;
  }
  return 0.;
} 

namespace AHADIC {
  ostream & operator<<(ostream & s, Wave_Function & wf) 
  {
    WaveComponents * waves = wf.GetWaves();
    double wf2(0.);
    for (WaveComponents::iterator wfc=waves->begin();wfc!=waves->end();wfc++) 
      wf2 += wfc->second*wfc->second;
    s<<" "<<wf.m_hadron<<" ("<<wf.m_kfcode<<"), spin = "<<((wf.m_spin2-1)/2.)
     <<", weight = "<<wf2<<"."<<endl;
    for (WaveComponents::iterator wfc=waves->begin();wfc!=waves->end();wfc++) {
      s<<"     "<<wfc->first->first<<" "<<wfc->first->second
       <<" : "<<wfc->second<<" ---> "<<(1./(wfc->second*wfc->second))<<endl;
    }
    return s;
  }
}

