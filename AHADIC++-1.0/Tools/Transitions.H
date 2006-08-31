#ifndef Transitions_H
#define Transitions_H

#include "Hadron_Multiplet.H"
#include "Hadron_Wave_Function.H"
#include "Constituents.H"
#include "Cluster.H"
#include "Flavour.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Flavour_Sorting_Mass {
  public :
    bool operator() (const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2) {
      if (fl1.Mass()>fl2.Mass()) return true;
      return false;
    }
  };

  typedef std::map<ATOOLS::Flavour,double,
		   Flavour_Sorting_Mass>         Single_Transition_List;
  typedef Single_Transition_List::iterator       Single_Transition_Siter;

  typedef std::map<Flavour_Pair,
		   Single_Transition_List *>     Single_Transition_Map;
  typedef Single_Transition_Map::iterator        Single_Transition_Miter;

  class Single_Transitions {
  private:
    Single_Transition_Map * p_transitions;
    double                  m_offset;

    ATOOLS::Flavour         m_lightest_constituent;
  public:
    Single_Transitions();
    ~Single_Transitions();
    bool                    MustDesintegrate(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    ATOOLS::Flavour         GetLightestTransition(const Flavour_Pair &);
    ATOOLS::Flavour         GetHeaviestTransition(const Flavour_Pair &);
    double                  GetLightestMass(const Flavour_Pair &);
    double                  GetHeaviestMass(const Flavour_Pair &);
    Single_Transition_Map * GetTransitions() { return p_transitions; }
    void                    PrintSingleTransitions();
  };


  class Flavour_Pair_Sorting_Mass {
  public :
    bool operator() (const Flavour_Pair & flpair1,
		     const Flavour_Pair & flpair2) {
      if ((flpair1.first.Mass()+flpair1.second.Mass())<
	  (flpair2.first.Mass()+flpair2.second.Mass())) return false;
      return true;
    }
  };

  typedef std::map<Flavour_Pair,double,
		   Flavour_Pair_Sorting_Mass>         Double_Transition_List;
  typedef Double_Transition_List::iterator            Double_Transition_Siter;
  
  typedef std::map<Flavour_Pair,
		   Double_Transition_List *>          Double_Transition_Map;
  typedef Double_Transition_Map::iterator             Double_Transition_Miter;
  
  class Double_Transitions {
  private:
    Double_Transition_Map * p_transitions;
    double                  m_offset;

  public:
    Double_Transitions();
    ~Double_Transitions();
    Flavour_Pair            GetLightestTransition(const Flavour_Pair &);
    Flavour_Pair            GetHeaviestTransition(const Flavour_Pair &);
    double                  GetLightestMass(const Flavour_Pair &);
    double                  GetHeaviestMass(const Flavour_Pair &);
    Double_Transition_Map * GetTransitions() { return p_transitions; }
    void                    PrintDoubleTransitions();
  };

}

#endif
