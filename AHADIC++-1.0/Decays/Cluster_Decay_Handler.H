#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Cluster_Part.H"
#include "Soft_Cluster_Handler.H"
#include "Blob_List.H"
#include "Cluster.H"
#include "Cluster_Decay_Analysis.H"
#include "Return_Value.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    bool                     m_cib;

    Soft_Cluster_Handler   * p_softclusters;
    Cluster_Part           * p_clus;

    Cluster_Decay_Analysis * p_analysis;
 
    ATOOLS::Blob *  DecayIt(Cluster *);
    void            AnalyseThis(ATOOLS::Blob *);
  public:
    Cluster_Decay_Handler(Soft_Cluster_Handler *,bool=false,bool=false);
    ~Cluster_Decay_Handler();

    ATOOLS::Return_Value::code DecayClusters(Cluster_List *,ATOOLS::Blob_List * = NULL);
  };

  /*!
    \file Cluster_Decay_Handler.H
    \brief Contains the class AHADIC::Cluster_Decay_Handler and the struct cdm.
  */
  /*! 
    \struct cdm
    \brief Steers the way the Cluster_Decay_Handler works, in particular the interplay of 
    AHADIC::Cluster_Part and AHADIC::Hadron_Part.

    For the AHADIC::Cluster_Part there are, so far, the following options
    - Simple_Q_over_M : 10
    - Four_Fermion    : 40

    and for the AHADIC::Hadron_Part, there are
    - Isotropic : 1
    - Retain    : 2.
    
    These numbers are added.
  */
  /*! 
    \class Cluster_Decay_Handler
    \brief Organises and steers the decays of a list of clusters into hadrons, where, intermediately, more clusters
    may be added to the list.

    The two main methods in this class are DecayClusters(Cluster_List *,ATOOLS::Blob *), which iterates over all
    clusters in the list, which has been filled before by the Cluster_Formation_Handler. For each single cluster,
    the method DecayIt(Cluster *) is called, which decays the cluster, and eventually, the emerging clusters
    until only hadrons are left. These hadrons are filled into the internal Part_List and, when the cluster decays have
    terminated, added to the fragmentation blob as outgoing particles. 
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member variables for Cluster_Decay_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \var cdm::code Cluster_Decay_Handler::m_cdm
    \brief The modus operandi of the cluster decays, steers which specific Cluster_Part and Hadron_Part
    implementation is initialised and passed as arguments to the Cluster_Decayer_Base. The default value
    is cdm::SimpleQoverM_Isotropic.
  */
  /*!
    \var Cluster_Decayer_Base * Cluster_Decay_Handler::p_decayer
    \brief The pointer to the Cluster_Decayer_Base which is initialised in the Cluster_Decay_Handler. 
  */
  /*!
    \var Cluster_Decay_Analysis * Cluster_Decay_Handler::p_analysis
    \brief A simple analysis of what happens in cluster decays. Its initialisation is steered by the arguments
    in the constructor.
  */
  /*!
    \var Cluster_Transformer * Cluster_Decay_Handler::p_transformer
    \brief To treat single clusters that cannot be dealt with by the Cluster_Decayer_Base. It is not yet functional
    in the form needed here.

    \todo Work on the way the Cluster_Transformer interacts with the Cluster_Decay_Handler.
  */
  /*!
    \var Part_List * Cluster_Decay_Handler::p_partlist
    \brief The internal particle list of hadron emerging in the decays of the clusters.
  */
  /*!
    \var ATOOLS::Blob * Cluster_Decay_Handler::p_blob
    \brief The fragmentation blob.
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Description of member methods for Cluster_Decay_Handler
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /*!
    \fn Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_Transformer *,bool=false)
    \brief Constructor, initialises the Cluster_Decayer_Base according to the strategy in m_cdm and the
    Cluster_Decay_Analysis object, if the flag is set to "true".
  */
  /*!
    \fn Cluster_Decay_Handler::~Cluster_Decay_Handler()
    \brief Destructor, deletes the p_decayer, the p_analysis, and the p_partlist. 
  */
  /*!
    \fn bool Cluster_Decay_Handler::DecayIt(Cluster *)
    \brief Manages the decays of a cluster and its left and right offsprings, if they are clusters. The method
    returns "true", if this was successful; if not, i.e. if the Cluster_Transformer had to deal with single clusters
    that could not decay according to standard strategies, it returns "false".

    The procedure in this method is:
    - If the Cluster_Decayer_Base managed to produce a valid cluster decay, irrespective of whether the decay was in 
      hadrons or clusters, through the method Cluster_Decayer_Base::Treat(Cluster *,Part_List *), the offsprings
      of the cluster, i.e. the left and right pointer are tested. If they exist, i.e. if they are clusters
      and not hadrons, DecayIt(Cluster *) is called again, with them as argument. Finally, "true" is returned.
    - If there was no valid decay produced by Cluster_Decayer_Base::Treat(Cluster *,Part_List *), the
      Cluster_Transformer::TreatSingleCluster(Cluster *,Part_List *) is called and "false" is returned.
      In this case, a forced \f$\pi^0\f$ emission is tried in order to make a hadron out of the cluster
      by absorbing the excess energy in this emission.
  */
  /*!
    \fn void Cluster_Decay_Handler::AnalyseThis(ATOOLS::Blob *)
    \brief Analyses the hadrons emerging from the cluster decays.
  */
  /*!
    \fn void Cluster_Decay_Handler::DecayClusters(Cluster_List *,ATOOLS::Blob *)
    \brief Decays all clusters into hadron and adds them as outgoing particles to the fragmentation blob.

    The key element in this is the call to DecayIt(Cluster *). If this was successful the iterator is
    just incremented, if not, the cluster in question must be deleted, since it has been replaced by one or
    two hadrons.
  */
}
#endif
