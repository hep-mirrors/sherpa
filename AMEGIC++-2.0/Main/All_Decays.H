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
  typedef Vertex_Table::iterator                        VTIterator;

  typedef std::map<ATOOLS::Flavour, Full_Decay_Table *> DecayMap;
  typedef DecayMap::iterator                            DMIterator;


  //! The container for a list of processes.
  class All_Decays {
  private :
    Interaction_Model_Base * p_model;
    Topology               * p_top;
    Vertex_Table             m_vertextable;
    ATOOLS::FlavourSet       m_particles;
    DecayMap                 m_decays;
    Full_Decay_Channel     * p_decay;

    bool CheckInVertex(ATOOLS::Flavour);
    void BinaryDecays();
    void ThreeBodyDecays();
    void ArrangeDecays();
    bool InitializeDecays();
  public:
    All_Decays(Interaction_Model_Base *,Topology *);
    bool AddToDecays(const ATOOLS::Flavour &);
    bool AddToDecays(ATOOLS::Decay_Channel *);
    void PrintDecayings();
    bool InitializeDecayTables();
    bool CalculateWidths(std::string=std::string(""));
    double             Width(ATOOLS::Flavour _fl);
    Full_Decay_Table * GetFullDecayTable(ATOOLS::Flavour _fl);

    bool               UnweightedEvent(ATOOLS::Decay_Channel *,double);
    int                Nin();
    int                Nout();
    std::string        Name();
    ATOOLS::Flavour  * Flavs();            
    ATOOLS::Vec4D    * Momenta();          
    int                NumberOfDiagrams(); 
    AMEGIC::Point    * Diagram(int _diag);
    void               Reset();
  };
  inline int All_Decays::Nin() { 
    if (p_decay) return p_decay->Nin(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Nin(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline int All_Decays::Nout() { 
    if (p_decay) return p_decay->Nout(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Nout(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline std::string All_Decays::Name() { 
    if (p_decay) return p_decay->Name(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Name(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline ATOOLS::Flavour * All_Decays::Flavs() { 
    if (p_decay) return p_decay->Flavs(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Flavs(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline ATOOLS::Vec4D * All_Decays::Momenta() { 
    if (p_decay) return p_decay->Momenta(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Momenta(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline int All_Decays::NumberOfDiagrams() { 
    if (p_decay) return p_decay->NumberOfDiagrams(); 
    ATOOLS::msg.Error()<<"Error in All_Decays::NumberOfDiagrams(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline AMEGIC::Point * All_Decays::Diagram(int _diag) { 
    if (p_decay) return p_decay->Diagram(_diag); 
    ATOOLS::msg.Error()<<"Error in All_Decays::Diagram(): no decay selected. Abort run."<<std::endl;
    abort();
  }
  inline void All_Decays::Reset() { p_decay = NULL; }
}
#endif

