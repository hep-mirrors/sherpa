#ifndef SHRIMPS_Ladders_Colour_Generator_H
#define SHRIMPS_Ladders_Colour_Generator_H

#include "SHRiMPS/Ladders/Ladder.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>


namespace SHRIMPS {
  class Colour_Generator {
  private:
    Ladder    * p_ladder;
    LadderMap * p_emissions;
    TPropList * p_props;
    int         m_propcolours[2];
    
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
