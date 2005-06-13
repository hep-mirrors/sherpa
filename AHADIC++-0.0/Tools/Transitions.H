#ifndef Transitions_H
#define Transitions_H

#include "Hadron_Multiplet.H"
#include "Hadron_Wave_Function.H"
#include "Constituents.H"
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


  typedef std::map<FlavPair,double>                      Double_Transition_List;
  typedef Double_Transition_List::iterator               Double_Transition_Siter;

  typedef std::map<FlavPair,Double_Transition_List *>    Double_Transition_Map;
  typedef Double_Transition_Map::iterator                Double_Transition_Miter;

  struct dtm {
    enum code {
      PS_only      = 0,
      waves_PS     = 1,
      PS_pop       = 10,
      waves_PS_pop = 11
    };
  };

  class All_Double_Transitions {
  private:
    dtm::code               m_dtmode;
    All_Hadron_Multiplets * p_multiplets;
    Double_Transition_Map * p_transitions;
    FlavCCMap             * p_alloweds;

    void FillDoubleTransitionMap();
  public:
    All_Double_Transitions(All_Hadron_Multiplets *);
    ~All_Double_Transitions();

    void                    PrintDoubleTransitions();
  };
}

#endif
