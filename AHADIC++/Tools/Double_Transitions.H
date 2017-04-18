#ifndef AHADIC_Tools_Double_Transitions_H
#define AHADIC_Tools_Double_Transitions_H

#include "AHADIC++/Tools/Single_Transitions.H"


namespace AHADIC {
  class Flavour_Pair_Sorting_Mass {
  public :
    bool operator() (const Flavour_Pair & flpair1,
		     const Flavour_Pair & flpair2) const {
      if ((flpair1.first==flpair2.first &&
      	   flpair1.second==flpair2.second)) return false;
      if (ATOOLS::IsEqual(flpair1.first.HadMass()+flpair1.second.HadMass(),
			  flpair2.first.HadMass()+flpair2.second.HadMass())) {
	if (flpair1.first.Kfcode()>flpair2.first.Kfcode()) return false;
      }
      if (flpair1.first.HadMass()+flpair1.second.HadMass()<
	  flpair2.first.HadMass()+flpair2.second.HadMass()) return false;
      return true;
    }
  };

  typedef std::map<Flavour_Pair,double,
  		   Flavour_Pair_Sorting_Mass>         Double_Transition_List;

  typedef std::map<Flavour_Pair,
		   Double_Transition_List *>          Double_Transition_Map;
  
  class Double_Transitions {
  private:
    Double_Transition_Map m_transitions;
    double m_wtthres;

    void FillMap(Single_Transitions * singles);
    void Normalise();
  public:
    Double_Transitions(Single_Transitions * singles);
    ~Double_Transitions();
    
    Flavour_Pair             GetLightestTransition(const Flavour_Pair & flavs);
    Flavour_Pair             GetHeaviestTransition(const Flavour_Pair & flavs);
    double                   GetLightestMass(const Flavour_Pair & flavs);
    double                   GetHeaviestMass(const Flavour_Pair & flavs);
    Double_Transition_Map  * GetTransitions() { return &m_transitions; }
    Double_Transition_List * operator[](const Flavour_Pair & flavs);
    void                     Print(const bool& full=false);
  };
}
#endif
