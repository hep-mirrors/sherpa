#ifndef ATOOLS_Phys_Decay_Map_H
#define ATOOLS_Phys_Decay_Map_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Decay_Table.H"

namespace ATOOLS {

  struct FlavourComp {
    bool operator()(const Flavour& fl1, const Flavour& fl2 ) const {
      if (fl1.HadMass()!=fl2.HadMass()) return fl1.HadMass()<fl2.HadMass();
      else if (fl1.Kfcode()!=fl2.Kfcode()) return fl1.Kfcode()<fl2.Kfcode();
      else if (fl1.IsAnti()!=fl2.IsAnti()) return fl2.IsAnti();
      else return false;
    }
  };

  class Decay_Map : 
    public std::map<ATOOLS::Flavour, std::vector<Decay_Table *>, FlavourComp> {
  public :
    Decay_Map();
    virtual ~Decay_Map();

    friend std::ostream &operator<<(std::ostream &os, const Decay_Map &dm);
  };
}

#endif
