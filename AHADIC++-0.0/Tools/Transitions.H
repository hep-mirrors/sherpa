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

    Single_Transition_Map * GetSingleTransitions() { return p_transitions; }
    bool                    MustTransit(Cluster *,ATOOLS::Flavour &,
					const double=0.,const bool=false);
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
  public:
    All_Double_Transitions(All_Hadron_Multiplets *);
    ~All_Double_Transitions();

    bool                    IsoDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &);
    Double_Transition_Map * GetDoubleTransitions() { return p_transitions; }
    void                    PrintDoubleTransitions();
  };

  /*!
    \file 
    \brief Contains the classes AHADIC::All_Single_Transitions and  classes AHADIC::All_Double_Transitions
    as well as the little helpers AHADIC::Flavour_Sorting_Mass, AHADIC::stm, and AHADIC::dtm.
  */

  /*!
    \class All_Single_Transitions
    \brief This class is central to the management of transitions of clusters into single hadrons. It
    is employed by, e.g., AHADIC::Cluster_Transformer. A pointer to it can be obtained from the
    AHADIC::Hadronisation_Parameters through AHADIC::Hadronisation_Parameters::GetSingleTransitions().
    The Hadronisation_Parameters also initialise this object.

    So far, there are two modes implemented to select a hadron for a cluster-to-hadron transition:
    - a simple mass criterion, All_Single_Transitions::m_stmode=1. There the hadron is selected, for which
      \f$ |m_{\rm hadron}-m_{\rm cluster}|\f$ is minimal.
    - a criterion incorporating both mass and wave function of the hadron, m_stmode=2. There the hadron is
      selected for which the weight
      \f[
      w = \frac{\left|\Psi_{3\bar 3}^{\rm hadron}\right|^2}{\left|m_{\rm hadron}-m_{\rm cluster}^2\right|}
      \f]
      is maximal.
  */
  /*!
    \class All_Double_Transitions
    \brief This class is central to the management of transitions of clusters into hadron pairs. It
    is employed by, e.g., AHADIC::Hadron_Part.  A pointer to it can be obtained from the
    AHADIC::Hadronisation_Parameters through AHADIC::Hadronisation_Parameters::GetDoubleTransitions().
    The Hadronisation_Parameters also initialise this object.
  */
  /*!
    \class Flavour_Sorting_Mass
    \brief Used to build a AHADIC::Single_Transition_List, which is a set of (hadron-) flavours, sorted
    with decreasing mass through this class. Thus, the only method is an operator, which compares two flavours 
    according to their mass and returns either a zero or one.
  */
  /*!
    \struct stm
    \brief A tag which steers the mode of the AHADIC::All_Single_Transitions. So far included are
    two modes, namely a criterion based on masses only - simplemass(1) - and a more involved one, based on
    both masses and wave functions of the hadrons in question - massXwaves(2).
  */
  /*!
    \struct dtm
    \brief A tag which steers the mode of the AHADIC::All_Double_Transitions. So far included are
    a number of modes, namely phase space only(0), or wave function times phase space(1), phase space
    times flavour popping probability (10), or a combination of all three (11).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for All_Single_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var stm::code All_Single_Transitions::m_stmode 
    \brief The mode with which the cluster-to-hadron transitions are performed.
  */
  /*!
    \var All_Hadron_Multiplets * All_Single_Transitions::p_multiplets 
    \brief A pointer to the hadron multiplets, filled before by the AHADIC::Hadronisation_Parameters.
  */
  /*!
    \var Single_Transition_Map * All_Single_Transitions::p_transitions
    \brief A map of flavour pairs (triplet-antitriplet), i.e. potential consituents of the hadrons,
    and all possible hadrons whose wave function has a component equal to the flavour pair. It is the
    heart of the class and filled in the constructor.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for All_Single_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn All_Single_Transitions::All_Single_Transitions(All_Hadron_Multiplets *)
    \brief In the constructor, the Single_Transition_Map p_transitions is initialised and filled.

    To this end, all AHADIC::Hadron_Wave_Function, accessible through All_Hadron_Multiplets::GetWaveFunctions(),
    are iterated over. Whenever a Hadron_Wave_Function has a component (triplet-antitriplet), not included into
    the Single_Transition_Map p_transitions so far, a new entry is created and hadrons are filled in.
    By default, the method mass times wave function is selected for the single hadron transition of clusters.
  */
  /*!
    \fn All_Single_Transitions::~All_Single_Transitions()
    \brief Deletes the Single_Transition_Map p_transitions.
  */
  /*!
    \fn Single_Transition_Map * All_Single_Transitions::GetSingleTransitions()
    \brief Returns the Single_Transition_Map p_transitions.
  */
  /*!
    \fn bool All_Single_Transitions::MustTransit(Cluster *,ATOOLS::Flavour &,const double=0.)
    \brief Checks, if a cluster must transit into a hadron, which in case will override the argument.
    A cluster must transit if its mass is smaller than the heaviest hadron with a wave function component
    equal to the two constituents plus an offset.

    For a given cluster, the heaviest hadron with a wave function component equal to the two constituents
    is selected. If it is heavier than the cluster mass plus an offset, the cluster must decay into a hadron.
    This hadron is selected according to either of the methods steered by m_stmode, i.e. by one of the methods 
    SimpleMassCriterion(Cluster *,Single_Transition_List *,ATOOLS::Flavour &) or 
    MassTimesWavefunction(Cluster *,Single_Transition_List *,ATOOLS::Flavour &). At the moment, clusters
    consisting of two diquarks are not yet considered. 

    The method is called by AHADIC::Cluster_Transformer and AHADIC::Cluster_Decayer_Base.
  */
  /*!
    \fn bool All_Single_Transitions::NextLightest(Cluster *,ATOOLS::Flavour &)
    \brief For a given flavour constellation inside a cluster and a given hadron, the next lightest hadron
    with suitable flavour composition will be produced. If this succeeds, "true" is returned.

    To be more specific, if so far no hadron has been specified the first entry of the corresponding 
    Single_Transition_List, given by the cluster's two constituents, will override the argument and "true" will be 
    returned. If a hadron has been specified, the next one of the list will be looked for. If this is successful, 
    the flavour in the argument will be fixed accordingly and "true" will be returned. If in contrast, there is no 
    next one, i.e., if there is no hadron with less mass than the previous one, "false" will be returned. This method
    will lead to an abortion of the program, if no suitable Single_Transition_List for the two constituents is found. 

    The method is called by AHADIC::Cluster_Transformer and AHADIC::Hadron_Part.
  */
  /*!
    \fn void All_Single_Transitions::PrintSingleTransitions()
    \brief Prints all possible single transitions, including some check outputs.
  */
  /*!
    \fn void All_Single_Transitions::SimpleMassCriterion(Cluster *,Single_Transition_List *,ATOOLS::Flavour &)
    \brief Yields the hadron for the cluster-to-hadron transition, that has the smallest mass difference
    w.r.t. the cluster, i.e. for which \f$ |m_{\rm hadron}-m_{\rm cluster}|\f$ is minimal.
  */
  /*!
    \fn void All_Single_Transitions::MassTimesWavefunction(Cluster *,Single_Transition_List *,ATOOLS::Flavour &)
    \brief Yields the hadron for the cluster-to-hadron transition, that has the maximal transition weight,
    which is given by
    \f[
    w = \frac{\left|\Psi_{3\bar 3}^{\rm hadron}\right|^2}{\left|m_{\rm hadron}-m_{\rm cluster}^2\right|}
    \f]
  */


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for All_Double_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var stm::code All_Double_Transitions::m_dtmode 
    \brief The mode with which the cluster decays into hadron pairs are performed.
  */
  /*!
    \var All_Hadron_Multiplets * All_Double_Transitions::p_multiplets 
    \brief A pointer to the hadron multiplets, filled before by the AHADIC::Hadronisation_Parameters.
  */
  /*!
    \var Double_Transition_Map * All_Double_Transitions::p_transitions
    \brief A map of flavour pairs (triplet-antitriplet), i.e. potential consituents of clusters,
    and all possible hadron pairs into which they could decay. It is the heart of the class and filled 
    in the constructor. 
  */  
  /*!
    \var FlavCCMap * All_Double_Transitions::p_alloweds
    \brief A list of constituents that could be used as flavours "popped out", when the cluster splits. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for All_Double_Transitions
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn All_Double_Transitions::All_Double_Transitions(All_Hadron_Multiplets *)
    \brief In the constructor, the Double_Transition_Map p_transitions is initialised and filled.

    To this end, all allowed combinations of AHADIC::Hadron_Wave_Function are checked, accessible by
    All_Hadron_Multiplets::GetWaveFunctions(). Then the individual wave components of each wave function are 
    iterated over, to check for a legal decay of two cluster consitutent flavours into two hadrons (given 
    by the wave functions). Any such decay is legal, if the second flavour of the first wave component 
    equals the barred first flavour of the second wave component, i.e. if there is a common flavour that could
    be popped out of the vacuum in the cluster decay. If such a combination exists, the corresponding weights 
    (multiplet weight times wave component) from each Hadron_Wave_Function are combined to yield the constant 
    transition weight for the potential transition.
  */
  /*!
    \fn All_Double_Transitions::~All_Double_Transitions()
    \brief Deletes the Double_Transition_Map.
  */
  /*!
    \fn bool All_Double_Transitions::IsoDecay(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &)
    \brief Gives a decay mode for a cluster into two hadrons according to the phase space weight times
    the weight from the Double_Transition_Map All_Double_Transitions::p_transitions. If two suitable 
    flavours are found, "true" is returned.

    For each combination of two hadrons, the phase space weight is given by
    \f[
    w^2_{\rm PS} = [M_{\rm cluster}^2-(m_{{\rm had}_1}+m_{{\rm had}_2})^2]\cdot
    [M_{\rm cluster}^2-(m_{{\rm had}_1}-m_{{\rm had}_2})^2]\,.
    \f]
    It is multiplied with the corresponding \f$C\to HH\f$ transition weight, a combination of the fixed
    multiplet weight and the squares of the hadron wave function weights, i.e. the corresponding flavour wave 
    function components squared:
    \f[
    w_{\rm tot} = w_{\rm multi}^{(1)}\cdot w_{\rm multi}^{(2)} 
    \cdot \left|\Psi_{3\bar 3}^{(1)}\right|^2 \cdot \left|\Psi_{3\bar 3}^{(2)}\right|^2 \cdot w_{\rm PS}
    \f]
    The actual decay channel is then selected according to this combined weight.<br>
    <b>Note that this does not include any further suppression on the popping of the quark- or diquark pair
    necessary to produce the two hadrons</b>.

    The method is called by AHADIC::Hadron_Part.
  */
  /*!
    \fn Double_Transition_Map * All_Double_Transitions::GetDoubleTransitions()
    \brief Returns the pointer to the Double_Transition_Map.
  */
  /*!
    \fn void All_Double_Transitions::PrintDoubleTransitions()
    \brief Prints all possible channels for \f$C\to HH\f$ transitions, including some check outputs.
  */
}

#endif
