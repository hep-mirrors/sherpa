#ifndef Hadron_Part_H
#define Hadron_Part_H

#include "Transitions.H"
#include "Pair_Popper.H"
#include "Cluster.H"


namespace AHADIC {
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

  class Hadron_Part {
  private:
    dtm::code               m_dtmode;
    Double_Transition_Map * p_transitions;
    Pair_Popper           * p_popper;
    double                  m_fraction;
  public:
    Hadron_Part();
    ~Hadron_Part();

    bool FixHHDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    bool MustTransit(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,const double=0.);
      
    void PrintDoubleTransitions();
  };
}
#endif
