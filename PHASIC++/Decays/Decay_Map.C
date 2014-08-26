#include "PHASIC++/Decays/Decay_Map.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Decay_Map::Decay_Map(const Mass_Selector* ms) :
  map<Flavour, std::vector<Decay_Table*>, FlavourComp>(FlavourComp(ms)), p_ms(ms)
{
}

Decay_Map::~Decay_Map()
{
  for (Decay_Map::iterator pos = this->begin(); pos != this->end(); ++pos) {
    for(size_t i=0; i<pos->second.size(); i++) {
      delete pos->second[i];
    }
  }
}

bool Decay_Map::Knows(const ATOOLS::Flavour & decayer)
{
  Decay_Map::iterator it = find(decayer);
  if(it==end()) {
    it = find(decayer.Bar());
  }
  if(it==end()) return false;
  else return true;
}

Decay_Table* Decay_Map::FindDecay(const ATOOLS::Flavour & decayer)
{
  Flavour tempdecayer=decayer;
  Decay_Map::iterator it = find(decayer);
  if(it==end()) {
    it = find(decayer.Bar());
    tempdecayer=decayer.Bar();
  }
  if(it==end()) return NULL;

  // there may be multiple decay tables for one flavour, so find the right one
  int count = it->second.size()-1; // default to last decay table available
  map<ATOOLS::Flavour,int>::iterator counterit = m_counters.find(tempdecayer);
  if(counterit!=m_counters.end() &&
     counterit->second < int(it->second.size()-1) )
  {
    count = counterit->second;
    counterit->second++;
  }
  return it->second[count];
}

pair<Decay_Table*, Decay_Channel*> Decay_Map::FindDecayChannel(string idcode,
                                                               bool create)
{
  stringstream ss(idcode);
  string item;
  Flavour_Vector flavs;
  while (getline(ss, item, ',')) {
    int kfc(ToType<int>(item));
    flavs.push_back(Flavour(abs(kfc), kfc<0));
  }
  Decay_Table* dt = FindDecay(flavs[0]);
  if (dt) {
    for (Decay_Table::iterator it=dt->begin(); it!=dt->end(); ++it) {
      if ((*it)->IDCode()==idcode) return make_pair(dt, (*it));
    }
    if (create) {
      Decay_Channel* dc = new Decay_Channel(flavs[0], p_ms);
      for (int j=1; j<flavs.size(); ++j) dc->AddDecayProduct(flavs[j]);
      dc->SetActive(0);
      dt->AddDecayChannel(dc);
      return make_pair(dt, dc);
    }
    else return make_pair(dt, (Decay_Channel*) NULL);
  }
  else return make_pair((Decay_Table*) NULL, (Decay_Channel*) NULL);
}

void Decay_Map::ResetCounters()
{
  map<ATOOLS::Flavour,int>::iterator it;
  for(it=m_counters.begin(); it!=m_counters.end(); it++) {
    it->second=0;
  }
}


bool FlavourComp::operator()(const ATOOLS::Flavour& fl1, const ATOOLS::Flavour& fl2 )
  const {
  if (p_ms->Mass(fl1)!=p_ms->Mass(fl2))
    return p_ms->Mass(fl1)<p_ms->Mass(fl2);
  else if (fl1.Kfcode()!=fl2.Kfcode()) return fl1.Kfcode()<fl2.Kfcode();
  else if (fl1.IsAnti()!=fl2.IsAnti()) return fl2.IsAnti();
  else return false;
}


namespace PHASIC {
  std::ostream &operator<<(std::ostream &os,const Decay_Map &dm)
  {
    for (Decay_Map::const_iterator it=dm.begin(); it!=dm.end(); ++it) {
      for (size_t i=0; i<it->second.size(); ++i) {
        os<<*(it->second.at(i))<<endl;
      }
    }
    return os;
  }
}
