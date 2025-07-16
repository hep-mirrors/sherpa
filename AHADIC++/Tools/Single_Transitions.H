#ifndef AHADIC_Tools_Single_Transitions_H
#define AHADIC_Tools_Single_Transitions_H

#include "AHADIC++/Tools/Wave_Function.H"
#include "AHADIC++/Tools/Constituents.H"


namespace AHADIC {
  class Flavour_Sorting_Mass {
  public :
    bool operator() (const ATOOLS::Flavour & fl1,
		     const ATOOLS::Flavour & fl2) const {
      if (fl1.HadMass()>fl2.HadMass()) return true;
      return false;
    }
  };
  
  typedef std::map<ATOOLS::Flavour,double,Flavour_Sorting_Mass>
  Single_Transition_List;

  typedef std::map<Flavour_Pair,Single_Transition_List *>
  Single_Transition_Map;



  class Single_Transitions {
  private:
    Single_Transition_Map m_transitions;

    void FillMap(Wave_Functions * wavefunctions);
    void Normalise();
  public:
    Single_Transitions(Wave_Functions * wavefunctions);
    ~Single_Transitions();
    ATOOLS::Flavour          GetLightestTransition(const Flavour_Pair &);
    ATOOLS::Flavour          GetHeaviestTransition(const Flavour_Pair &);
    double                   GetLightestMass(const Flavour_Pair &);
    double                   GetHeaviestMass(const Flavour_Pair &);
    Single_Transition_Map *  GetMap() { return &m_transitions; }
    Single_Transition_List * operator[](const Flavour_Pair & flavs);
    void                     Print();
  };
}
#endif
