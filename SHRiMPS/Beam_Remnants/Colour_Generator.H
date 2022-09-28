#ifndef SHRIMPS_Beam_Remnants_Colour_Generator_H
#define SHRIMPS_Beam_Remnants_Colour_Generator_H

#include "SHRiMPS/Ladders/Ladder.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <set>


namespace SHRIMPS {
  class Colour_Generator {
  private:
    int           m_dir;
    Ladder      * p_ladder;
    LadderMap   * p_emissions;
    TPropList   * p_props;
    std::set<int> m_colours[2][2];
    int           m_propcolours[2];
    
    void PickStartColours();
    void IterateColours(LadderMap::iterator out,TPropList::iterator prop);
    void PickEndColours();
    bool ReplaceFSColour(const size_t & pos,const int & orig,const int & repl);
  public:
    Colour_Generator();
    ~Colour_Generator();

    bool operator()(Ladder * ladder);

    void Reset(); 
    void OutputStack();
    inline std::set<int> & GetColours(const size_t & beam, const size_t & index) {
      return m_colours[beam][index];
    }
  };
}

#endif
