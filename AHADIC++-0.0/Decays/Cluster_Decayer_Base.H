#ifndef Cluster_Decayer_Base_H
#define Cluster_Decayer_Base_H

#include "Cluster_Part.H"
#include "Hadron_Part.H"
#include "Transitions.H"
#include "Cluster.H"


namespace AHADIC {
  class Cluster_Decayer_Base {
  private:
    All_Single_Transitions * p_stransitions;
    Cluster_Part           * p_cdecs;
    Hadron_Part            * p_chads;
 
    int             m_test;
    double          m_offset;
    ATOOLS::Flavour m_had1, m_had2;
    
    void TestOffSprings(Cluster *);
    void TreatHadDecay(Cluster *,Part_List *);
  public:
    Cluster_Decayer_Base(Cluster_Part *,Hadron_Part *);
    ~Cluster_Decayer_Base();
    bool Treat(Cluster *,Part_List *);
  };
  /*!
    \file Cluster_Decayer_Base.H
    \brief Contains the class AHADIC::Cluster_Decayer_Base. 
  */

  /*! 
    \class Cluster_Decayer_Base
    \brief This class is responsible for the organisation of individual cluster decays, for which it owns pointers
    to a AHADIC::Cluster_Part and a AHADIC::Hadron_Part object.

    The cluster decays are organised through a single method, namely Treat(Cluster *,Part_List *), which is
    called by AHADIC::Cluster_Decay_Handler::DecayIt(Cluster * cluster).
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Decayer_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var All_Single_Transitions * Cluster_Decayer_Base::p_stransitions
    \brief To test whether clusters are in fact hadrons.

    After a test decay with kinematics according to the strategy implemented in the AHADIC::Cluster_Part object,
    the emerging clusters are testeed whether they fall into the hadron regime. This is done by the method
    All_Single_Transitions::MustTransit(Cluster *,ATOOLS::Flavour &,double) in the method TestOffSprings(Cluster *).
  */
  /*!
    \var Cluster_Part * Cluster_Decayer_Base::p_cdecs
    \brief A pointer to the AHADIC::Cluster_Part, which produces test decays of a single cluster according to a
    certain strategy. 
  */
  /*!
    \var Hadron_Part * Cluster_Decayer_Base::p_chads
    \brief A pointer to the AHADIC::Hadron_Part, which is responsible for cluster decays involving hadrons.
  */
  /*!
    \var int Cluster_Decayer_Base::m_test
    \brief For each cluster, m_test distinguishes, how many hadrons are involved in the cluster decay. 

    This counter is set in TestOffSprings(Cluster *) and passed on to the relevant method in the Hadron_Part,
    Hadron_Part::RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &). In particular:
    - m_test = 0: no hadron
    - m_test = 1: the left cluster must be a hadron
    - m_test = 2: the right cluster must be a hadron
    - m_test = 3: both clusters must be hadrons
  */
  /*!
    \var double Cluster_Decayer_Base::m_offset
    \brief The offset parameter for the single transitions, initialised in the constructor by looking up the value
    in the AHADIC::Hadronisation_Parameters.
  */
  /*!
    \var ATOOLS::Flavour Cluster_Decayer_Base::m_had1
    \brief Hadron flavour for the left cluster emerging in the test decay.
  */
  /*!
    \var ATOOLS::Flavour Cluster_Decayer_Base::m_had2
    \brief Hadron flavour for the right cluster emerging in the test decay.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Decayer_Base
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Decayer_Base::Cluster_Decayer_Base(Cluster_Part *,Hadron_Part *)
    \brief Constructor, fills in the various variables of the class.
  */
  /*!
    \fn Cluster_Decayer_Base::~Cluster_Decayer_Base()
    \brief Destructor, deletes the Hadron_Part and the Cluster_Part, which have been constructed in the
    Cluster_Decay_Handler and passed on to the Cluster_Decayer_Base as argument of the constructor.
  */
  /*!
    \fn bool Cluster_Decayer_Base::Treat(Cluster *,Part_List *)
    \brief Organises the decay of a single cluster and returns "true" if a valid solution has been found.

    This is done in the following steps:
    - A test decay kinematics is produced by the Cluster_Part::TestDecay(Cluster *,Part_List *) method. It will
      return "true" as long as it can produce a decay kinematics plus a valid quark-antiquark (or diquark-antidiquark) 
      pair, otherwise it will return "false".
    - In the former case, i.e. if a valid cluster decay has been produced, the offsprings are tested individually, 
      whether they fall into the hadron regime. 
      -# If this is the case (m_test>0) the TreatHadDecay(Cluster *,Part_List *) method is being called, which will, 
         effectively, call Hadron_Part::RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &) to 
         add one or two hadrons to the particle list and delete the corresponding clusters. 
      -# Then, or if no hadron was produced, the method Treat returns "true".
    - If no valid cluster decay has been produced, a forced decay is imposed through a call to the method
      Hadron_Part::ForcedDecay(Cluster *,Part_List *). If this works out, "true" is returned by the Treat method.
    - Only if both strategies fail, "false" is returned.
  */
  /*!
    \fn void Cluster_Decayer_Base::TestOffSprings(Cluster *)
    \brief Checks the left and the right offspring in the test cluster decay and decides whether they fall into the
    hadron regime through calls to All_Single_Transitions::MustTransit(Cluster *,ATOOLS::Flavour &,double). 
    Correspondinlgy m_test is set in this method.
  */
  /*!
    \fn void Cluster_Decayer_Base::TreatHadDecay(Cluster *,Part_List *)
    \brief Deals with those cluster decays that involve one or two hadrons. 

    This method first calls Hadron_Part::RedoDecay(Cluster *,Part_List *,int,ATOOLS::Flavour &,ATOOLS::Flavour &)
    and then deletes those clusters that have been replaced by hadrons.
  */
}
#endif
