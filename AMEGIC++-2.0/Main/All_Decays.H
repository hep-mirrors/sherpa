#ifndef All_Decays_h
#define All_Decays_h

#include "Process_Group.H"
#include "Single_Process.H"

#include "Interaction_Model_Base.H"
#include "Topology.H"
#include "Flavour.H"
#include <set>
#include <deque>

namespace AMEGIC {
  typedef  std::vector<Single_Vertex *>              VertexList;
  typedef  std::map<APHYTOOLS::Flavour, Vertex_List> VertexTable;
  typedef  std::deque<APHYTOOLS::Flavour>            DecayingParticleList;
  typedef  DecayingParticleList::iterator            DecayingParticleIterator;
  typedef  std::vector<Process_Base *>               DecayList;
  typedef  DecayList::iterator                       DecayIterator;
  typedef  std::map<APHYTOOLS::Flavour, DecayList>   Decay_Table;
  typedef  Decay_Table::iterator                     DTIterator;

  //! The container for a list of processes.
  class All_Decays {
  private :
    Interaction_Model_Base     * p_model;
    Topology                   * p_top;
    VertexTable                  m_vertextable;
    DecayingParticleList         m_particles;
    Decay_Table                  m_decaytable;

    void FindUnstableParticles();
    bool CheckInVertex(APHYTOOLS::Flavour);
    void InitializeDecayTables();
    void BinaryDecays();
  public:
    All_Decays(Interaction_Model_Base *,Topology *);
    bool InitializeDecays();
    bool CalculateBranchingWidths(std::string=std::string(""));
  };
}


#endif
