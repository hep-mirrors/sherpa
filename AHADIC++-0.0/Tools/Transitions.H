#ifndef Transitions_H
#define Transitions_H

#include "Hadron_Multiplet.H"
#include "Hadron_Wave_Function.H"
#include "Cluster.H"
#include "Flavour.H"
#include <map>


namespace AHADIC {

  class Flavour_Sorting_Mass {
  public :
    bool operator() (const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2) {
      if (fl1.Mass()<fl2.Mass()) return 0;
      return 1;
    }
  };

  typedef std::set<ATOOLS::Flavour,Flavour_Sorting_Mass> Single_Transition_List;
  typedef Single_Transition_List::iterator               Single_Transition_Siter;

  typedef std::map<FlavPair,Single_Transition_List *>    Single_Transition_Map;
  typedef Single_Transition_Map::iterator                Single_Transition_Miter;

  struct stm {
    enum code {
      simplemass = 1,
      massXwaves = 2
    };
  };

  class All_Single_Transitions {
  private:
    stm::code               m_stmode;
    All_Hadron_Multiplets * p_multiplets;
    Single_Transition_Map * p_transitions;
    void SimpleMassCriterion(Cluster *,Single_Transition_List *,ATOOLS::Flavour &);
    void MassTimesWavefunction(Cluster *,Single_Transition_List *,ATOOLS::Flavour &);
  public:
    All_Single_Transitions(All_Hadron_Multiplets *);
    ~All_Single_Transitions();

    Single_Transition_Map * GetSingleTRansitions() { return p_transitions; }
    bool                    MustTransit(Cluster *,ATOOLS::Flavour &,const double=0.);
    bool                    NextLightest(Cluster *,ATOOLS::Flavour &);
    void                    PrintSingleTransitions();
  };
}

#endif
