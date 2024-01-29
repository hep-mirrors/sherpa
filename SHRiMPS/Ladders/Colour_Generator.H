#ifndef SHRIMPS_Ladders_Colour_Generator_H
#define SHRIMPS_Ladders_Colour_Generator_H

#include "SHRiMPS/Ladders/Ladder.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>


namespace SHRIMPS {
  struct lcolours {
    enum code {
      random        = 0,
      fix_emissions = 1,
      fix_all       = 2,
      no_singlets   = 3
    };
  };
  
  class Colour_Generator {
  private:
    Ladder       * p_ladder;
    LadderMap    * p_emissions;
    TPropList    * p_props;
    int            m_propcolours[2];
    lcolours::code m_mode;
    size_t         m_keepindex;
    
    void PickStartColours();
    void IterateColours(LadderMap::iterator out,TPropList::iterator prop);
    void PickEndColours();
  public:
    Colour_Generator();
    ~Colour_Generator();

    bool operator()(Ladder * ladder);
  };
}

#endif
