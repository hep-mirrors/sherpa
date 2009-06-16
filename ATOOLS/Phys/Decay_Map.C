#include "ATOOLS/Phys/Decay_Map.H"

using namespace ATOOLS;
using namespace std;

Decay_Map::Decay_Map() :
  map<Flavour, std::vector<Decay_Table*>, FlavourComp >()
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

std::ostream &ATOOLS::operator<<(std::ostream &os,const Decay_Map &dm)
{
  for (Decay_Map::const_iterator it=dm.begin(); it!=dm.end(); ++it) {
    for (size_t i=0; i<it->second.size(); ++i) {
      os<<*(it->second.at(i))<<endl;
    }
  }
  return os;
}
