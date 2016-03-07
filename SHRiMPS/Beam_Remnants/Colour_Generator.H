#ifndef SHRIMPS_Beam_Remnants_Colour_Generator_H
#define SHRIMPS_Beam_Remnants_Colour_Generator_H

#include "SHRiMPS/Event_Generation/Ladder.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>


namespace SHRIMPS {
  class Colour_Generator {
  private:
    std::set<int> m_colours[2][2];
    Ladder      * p_ladder;
    LadderMap   * p_emissions;
    TPropList   * p_props;
    
    void PickStartColours();
    void IterateColours(LadderMap::iterator out,TPropList::iterator prop,
			int col1,int col2);
    void FinishColours(int col1,int col2);
    void UpdateColours();
    bool ReplaceColours(const size_t & pos);
  public:
    Colour_Generator();
    ~Colour_Generator();

    bool operator()(Ladder * ladder);

    void Reset(); 
    void OutputStack();
  };
}

#endif
