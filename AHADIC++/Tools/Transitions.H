#ifndef AHADIC_Tools_Transitions_H
#define AHADIC_Tools_Transitions_H

#include "AHADIC++/Tools/Hadron_Multiplet.H"
#include "AHADIC++/Tools/Hadron_Wave_Function.H"
#include "AHADIC++/Tools/Constituents.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Return_Value.H"
#include <map>


namespace AHADIC {
  class Flavour_Sorting_Mass {
  public :
    bool operator() (const ATOOLS::Flavour & fl1,const ATOOLS::Flavour & fl2) {
      if (fl1.HadMass()>fl2.HadMass()) return true;
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
    ATOOLS::Flavour         GetLightestConstituent() { return m_lightest_constituent; }
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
      if (flpair1.first==flpair2.first &&
	  flpair1.second==flpair2.second) return true;
      if ((flpair1.first.HadMass()+flpair1.second.HadMass())<=
	  (flpair2.first.HadMass()+flpair2.second.HadMass())) return false;
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


  /*!
    \file Transitions.H
    \brief Contains the classes AHADIC::Flavour_Sorting_Mass, AHADIC::Flavour_Pair_Sorting_Mass, 
    AHADIC::Single_Transitions, and AHADIC::Double_Transitions as well as some type definitions
    employing maps with the sorting criteria defined here.  
  */

  /*!
    \class Flavour_Sorting_Mass
    \brief Sorts flavours according to their mass (heavies upfront)
  */

  /*!
    \class Flavour_Pair_Sorting_Mass
    \brief Sorts pairs of flavours according to their added mass (heavies upfront)
  */

  /*!
    \class Single_Transitions
    \brief Basically a map of cluster constituents (a triplet-antitriplet flavour pair) and
    the corresponding list of all allowed hadrons, into which this cluster could transit.  The
    class Single_Transitions just provides methods to fill this map and some simple methods
    to act on it.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Single_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Single_Transition_Map * Single_Transitions::p_transitions
    \brief A Single_Transition_Map is defined to be a map of a AHADIC::Flavour_Pair and
    a Single_Transition_List, which in turn is a map of individual flavours, the hadrons,
    and the overlap their wave function has with a given flavour pair.    
  */
  /*!
    \var double Single_Transitions::m_offset
    \brief The value for the transition offset \f$\delta_1\f$, given by "TRANSITION_OFFSET".
  */
  /*!
    \var ATOOLS::Flavour Single_Transitions::m_lightest_constituent
    \brief The lightest active constituent, usually the \f$u\f$ or \f$d\f$ quark.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member functions for Single_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Single_Transitions::Single_Transitions()
    \brief In the constructor, mainly, the Single_Transition_Map is filled.  

    To that end, the method iterates over all hadrons in all multiplets.  For each hadron,
    all wave function components are iterated over.  Depending of whether this is a new
    flavour constellation or not, a new entry, i.e. a new Single_Transition_List, for this 
    flavour combination is initialised or not.  In both cases then the hadron and the 
    corresponding weight (wave function times general multiplet weight) are added to
    the Single_Transition_List, which in fact is a map of the hadron flavours and their 
    respective total weight.
  */
  /*!
    \fn Single_Transitions::~Single_Transitions()
    \brief Delete the Single_Transition_Map Single_Transitions::p_transitions. 
  */
  /*!
    \fn bool Single_Transitions::MustDesintegrate(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Chcks if a cluster made of diquarks can decay into two baryons or must
    decay into two mesons by breaking up te diquarks.

    If the cluster consists of two diquarks, and if the two lightest possible baryons
    are heavier than the cluster, then the two lightest mesons consisting of the quarks
    and antiquarks, that build the diquarks in the cluster, are selected and "true" is
    returned. 
  */
  /*!
    \fn ATOOLS::Flavour Single_Transitions::GetLightestConstituent() 
    \brief Returns the lightest constituent.
  */
  /*!
    \fn ATOOLS::Flavour Single_Transitions::GetLightestTransition(const Flavour_Pair &)
    \brief Returns the lightest hadron, into which a cluster with a certain
    flavour composition can transit.
  */
  /*!
    \fn ATOOLS::Flavour Single_Transitions::GetHeaviestTransition(const Flavour_Pair &)
    \brief Returns the heaviest hadron, into which a cluster with a certain
    flavour composition can transit.
  */
  /*!
    \fn double Single_Transitions::GetLightestMass(const Flavour_Pair &)
    \brief Returns the mass of the lightest hadron, into which a cluster
    with a certain flavour combination can transit.
  */
  /*!
    \fn double Single_Transitions::GetHeaviestMass(const Flavour_Pair &)
    \brief Returns the mass of the heaviest hadron, into which a cluster
    with a certain flavour combination can transit.
  */
  /*!
    \fn Single_Transition_Map * Single_Transitions::GetTransitions()
    \brief Returns a pointer to the full Single::Transitions::Single_Transition_Map.
  */
  /*!
    \fn void Single_Transitions::PrintSingleTransitions
    \brief Prints all allowed cluster->hadron transitions.
  */

  /*!
    \class Double_Transitions 
    \brief Basically a map of cluster constituents (a triplet-antitriplet flavour pair) and
    the corresponding list of all allowed hadron pairs, into which this cluster could decay.  The
    class Double_Transitions just provides methods to fill this map and some simple methods
    to act on it.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Double_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var Double_Transition_Map * Double_Transitions::p_transitions
    \brief A Double_Transition_Map is defined to be a map of a AHADIC::Flavour_Pair and
    a Double_Transition_List, which in turn is a map of individual pairs of flavours, the hadrons,
    and the overlap their combined wave function has with a given flavour pair (the cluster
    constituents), after a new flavour pair has been popped.    
  */
  /*!
    \var double Double_Transitions::m_offset
    \brief The value for the transition offset \f$\delta_2\f$, given by "DECAY_OFFSET".
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member functions for Double_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Double_Transitions::Double_Transitions()
    \brief Constructs the Double_Transition_Map.
  */
  /*!
    \fn Double_Transitions::~Double_Transitions()
    \brief Deletes the map of all double transitions, Double_Transitions::p_transitions.
  */
  /*!
    \fn Flavour_Pair Double_Transitions::GetLightestTransition(const Flavour_Pair &)
    \brief Returns the lightest hadron pair, into which a cluster with a certain flavour 
    combination can decay.
  */
  /*!
    \fn Flavour_Pair Double_Transitions::GetHeaviestTransition(const Flavour_Pair &)
    \brief Returns the heaviest hadron pair, into which a cluster with a certain flavour 
    combination can decay.
  */
  /*!
    \fn double Double_Transitions::GetLightestMass(const Flavour_Pair &)
    \brief Returns the mass of the lightest hadron pair, into which a cluster
    with a certain flavour combination can decay.
  */
  /*!
    \fn double Double_Transitions::GetHeaviestMass(const Flavour_Pair &)
    \brief Returns the mass of the heaviest hadron pair, into which a cluster
    with a certain flavour combination can decay.
  */
  /*!
    \fn Double_Transition_Map * Double_Transitions::GetTransitions()
    \brief Returns the list of all allowed cluster->hadron+hadron decays.
  */
  /*!
    \fn void Double_Transitions::PrintDoubleTransitions()
    \brief Prints all allowed cluster->hadron+hadron transitions (decays).
  */
}
#endif
