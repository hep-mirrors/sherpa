#ifndef All_Decays_h
#define All_Decays_h

#include "Process_Group.H"
#include "Single_Process.H"

#include "Interaction_Model_Base.H"
#include "Topology.H"
#include "Flavour.H"
#include "Full_Decay_Table.H"
#include <set>
#include <deque>

namespace AMEGIC {
  typedef std::vector<Single_Vertex *>           VertexList;
  typedef std::map<ATOOLS::Flavour, Vertex_List> VertexTable;
  typedef VertexTable::iterator                  VTIterator;

  typedef std::map<ATOOLS::Flavour, Full_Decay_Table *> DecayMap;
  typedef DecayMap::iterator                            DMIterator;


  //! The container for a list of processes.
  class All_Decays {
  private :
    Interaction_Model_Base * p_model;
    Topology               * p_top;
    VertexTable              m_vertextable;
    ATOOLS::FlavourSet       m_particles;
    DecayMap                 m_decays;

    bool CheckInVertex(ATOOLS::Flavour);
    void BinaryDecays();
    void ThreeBodyDecays();
    void ArrangeDecays();
    bool InitializeDecays();
  public:
    All_Decays(Interaction_Model_Base *,Topology *);
    bool AddToDecays(const ATOOLS::Flavour &);
    void PrintDecayings();
    bool InitializeDecayTables();
    bool CalculateWidths(std::string=std::string(""));

    double             Width(ATOOLS::Flavour _fl);
    Full_Decay_Table * GetFullDecayTable(ATOOLS::Flavour _fl);
  };
}
#endif

