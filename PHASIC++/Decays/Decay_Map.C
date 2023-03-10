#include "PHASIC++/Decays/Decay_Map.H"

#include <algorithm>

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Decay_Map::Decay_Map(const Mass_Selector* ms) :
  map<Flavour, Decay_Table*, FlavourComp>(FlavourComp(ms)), p_ms(ms)
{
}

Decay_Map::~Decay_Map()
{
  for (Decay_Map::iterator pos = this->begin(); pos != this->end(); ++pos) {
    delete pos->second;
  }
}

bool Decay_Map::Knows(const ATOOLS::Flavour & decayer)
{
  return (FindDecay(decayer)!=NULL);
}

Decay_Table* Decay_Map::FindDecay(const ATOOLS::Flavour & decayer)
{
  Decay_Map::iterator it = find(decayer);
  if(it==end()) {
    it = find(decayer.Bar());
  }
  if(it==end()) return NULL;
  else return it->second;
}

void Decay_Map::ResetCounters()
{
  for (Decay_Map::const_iterator it=begin(); it!=end(); ++it) {
    it->second->ResetCounter();
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
      os<<*(it->second)<<endl;
    }
    return os;
  }
}
