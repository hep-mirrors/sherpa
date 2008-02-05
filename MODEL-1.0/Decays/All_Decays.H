#ifndef All_Decays_h
#define All_Decays_h

#include "Flavour.H"
#include "Decay_Table.H"
#include "Vertex.H"

namespace MODEL {
  class Model_Base;
  
  typedef std::map<ATOOLS::Flavour,ATOOLS::Decay_Table *> DecayMap;
  typedef DecayMap::iterator                              DMIterator;
  
  //! The container for all hard decays.
  class All_Decays {
  private :
    Vertex_Table * p_vertextable;
    DecayMap     * p_decays;

    std::map<ATOOLS::Flavour,double> m_maxestimates;
    
    bool CheckInVertex(const ATOOLS::Flavour &);
    void BinaryDecays();
    void ThreeBodyDecays();
    void ArrangeDecays();
    void PrintDecayTables();
  public:
    All_Decays(Model_Base *);
    bool AddToDecays(const ATOOLS::Flavour &);
    bool AddToDecays(ATOOLS::Decay_Channel *);
    bool InitializeDecayTables();
    double Width(const ATOOLS::Flavour & fl);
    
    DecayMap * const            GetDecayMap()    { return p_decays; }
    size_t                      NumberOfDecays() { return p_decays->size(); }
    ATOOLS::Decay_Table * const GetDecayTable(const ATOOLS::Flavour & flav) {
      if (p_decays->find(flav)==p_decays->end()) return NULL;
      return p_decays->find(flav)->second;
    }
  };
  
  inline double All_Decays::Width(const ATOOLS::Flavour & fl) {
    double width(-1.0);
    DMIterator dmit = p_decays->find(fl);
    if (dmit!=p_decays->end()) width = dmit->second->TotalWidth();
    return width;
  }
}
#endif
